// This file has been modified from Ken Silverman's original release

//PND3D current octree format:
// ls: sol chi|
//+-----------+-------
//|>0   0   0 |pure air
//|>0   0   1 |mix of air&sol ; some surfs ; child on octv_t
//|>0   1   0 |pure solid     ; no surfs
//|>0   1   1 |pure solid     ; some surfs ; child on octv_t
//+-----------+-------
//| 0   0   0 |air voxel
//| 0   0   1 |(N/A)
//| 0   1   0 |interior voxel
//| 0   1   1 |surface voxel               ; child on surf_t (col/norm/etc.)
//+-----------+-------
//
//Targetted apps: VCD, SLAB6, PND3D, AOS, GAME_TD
//
//TODO:
// ! drawcone(),drawsph(),drawpol() should be affected by oct_fogdist in /cpu mode
// ! FIZEX:small objects that validly get stuck start daffy ducking (once gravity incs vel) - need verlet
// ! FIZEX:hovering objects like to get stuck.
// * oct_modnew(): should call recvoctfunc() multiple times for disconnected regions
// * FIZEX:use early out ohit 1v1.
// * optimization: oct_touch_oct(): change cr.x = 2.0 to 1.0, etc.. so MAT0ISIDENTITY can be 1
// * support mip-map for /ils=? mode
// * bad rounding for rasterization: think borders are visible around cubes :/
// * GAME_TD: write animation editor
// * support highly sheared (non-orthogonal) axes
// * get ARB ASM running on at least 1 old GPU that doesn't already support GLSL.
// * tiles instancing editor.
// * marching cubes.
// * edgeiswrap/edgeissol: support!
// * oct_rendslice (oct_t *loct, point3f *pp, point3f *pr, point3f *pd, float mixval, int imulcol);
// . optimization: oct_modnew(): use consecutive allocation instead of bitalloc() for newoct
// . optimization: oct_modnew(): support simultaneous chop off original oct as mode option
// . optimization: oct_mod(): special case for empty tree (may speed up loading)
// . oct_sol2bit(): make it respect loct->edgeiswrap/edgeissol
// . write oct_oct_sweep() - translational sweep collision of octree.
// . Alt+H slow with memcpy because writing directly to GPU texture?
// . slice editor
// . symmetry copy
// . double/halve dimension tool
// . hollow fill?
// . new PIXMETH idea: store x/y/z/chi at head of 2x2x2 surf list; ~25% more GRAM, but can do larger octree and maybe faster (incompatible w/instancing :/).
// . oct_rendclosecubes(): can't use cornmin[]; must render all cubes within manhattan dist; new code fails for high FOV :/
// . bstatus doesn't work on first click (likely winmain dinput acquire/gbstatus problem)
// x research voxel tearing (STRESSGRID.KC)
// x 2 voxels on diagonal disappear when pos exactly on integer coord (must be in *_dorect); solved with temp hack
// x optimize: drawcone brush

//08/21/2010: FPS comparison; tests use: vxl/untitled.vxl /1024x768, no effects (Ken's C2Q@2.66Ghz)
//                        MT:1  MT:4
//   voxlap6:DrawPaint     19     -  (cheat:limited scan dist)
//   voxlap6:DrawPoint     25     -  (cheat:limited scan dist)
//   voxlap6:Drawcast6D    37    64  (cheat:limited scan dist)
//   voxed (voxlap5)       40     -  (cheat:limited scan dist)
//   pnd3d                 65    96  <-bbox suspected?
//   pnd3d                 71   107  (cheats:using rcpps & z-buffer disabled)
//  (voxlap6:Drawcast4D)  (52) (135) (cheats:limited scan dist & 4DOF!)
//   game (voxlap5)        60     -  (cheat:mip-mapping)
//
//05/29/2011: Same settings as ^ but on KJSI7:
//                        MT:1  MT:4  MT:8
//   voxlap6:DrawPaint     21     -     -  (cheat:limited scan dist)
//   voxlap6:DrawPoint     32     -     -  (cheat:limited scan dist)
//   voxlap6:Drawcast6D    47    94    94  (cheat:limited scan dist)
//   voxed (voxlap5)       59     -     -
//   voxed (voxlap5)       64     -     -  (cheat:limited scan dist)
//   pnd3d                 35   100   150
//   pnd3d                 47   130   177  (cheats:using bbox)
//   pnd3d                 47   137   178  (cheats:using bbox, rcpps, no zbuf)
//  (voxlap6:Drawcast4D)  (74) (233) (330) (cheats:limited scan dist & 4DOF!)
//   game (voxlap5)       116     -     -  (cheat:mip-mapping)

//08/18/2011: After enabling early gcbit check for ls==0:
//                opnd3d_bak3:    pnd3d:    |09/04/2011:
//                 MT:1  MT:8   MT:1  MT:8  |MT8,ASM,PBO,GLSL,PIXMETH=1 (texmap!):
//   CANYON.VXL   23.63 135.7  25.08 141.9  | 222
//   coreef.vxl   23.68 123.9  24.63 127.5  | 228
//   lung.vxl     17.64  79.6  20.42 106.5  | 168
//   test.vxl     20.52  65.0  24.52 105.0  | 160
//   tomland.png  20.54 104.6  21.63 106.8  | 329
//   untitled.vxl 28.62 159.6  28.71 160.0  | 254

//02/10/2012: mod speed, untitled.vxl:
//                oldmod_0: newmod_2:
//   load time:      1325ms     931ms
//   #getsurf:    5,512,197 5,512,197
//   #isins_ls0: 52,828,714   982,008
//   #isins_ls1: 20,134,389   380,296
//   #isins_ls2:  4,824,145   106,584
//   #isins_ls3:  1,136,307    25,712
//   #isins_ls4:    262,982     5,808
//   #isins_ls5:     59,472     1,368
//   #isins_ls6:     12,033       296
//   #isins_ls7:      2,298       512
//   #isins_ls8:         64        64
//   #isins_ls9:          8         8

//03/02/2012: testing both mod & rend (big GPUSEBO=1 speedup today!)
//CPU pure:     | GPU,GPUSEBO=0 | GPU,GPUSEBO=1
//mod:   rend:  | mod:   rend:  | mod:   rend:
//143ms  134fps | 140ms  146fps | 147ms  499fps
//143ms  128fps | 548ms  142fps | 149ms  433fps
//              |               |
//..            | ..            | ..
//7-14ms 106fps | 762ms  153fps | 8-13ms 546fps
// :)     :/    |  :(     :/    |  :)     :)

//03/14/2012: test load speed&stats    uncomp: kzip: loadms: nod.num: sur.num:
//pnd3d ../voxlap/vxl/tomland.png      2088348          1189  2212518  6042719
//pnd3d ../aos/vxl/tomland.vxl         2101280           250   333153   936306
//pnd3d ../aos/vxl/grandcan.vxl        2497156           282   382338  1073619
//pnd3d ../aos/vxl/castlesandshit.vxl  2673600           283   392134  1087500
//pnd3d ../aos/vxl/world12.vxl         3253900           291   401378  1137110
//pnd3d ../aos/vxl/metropolis.vxl      6532084           372   556221  1670494
//pnd3d ../voxlap/vxl/untitled.vxl    12425384  3288588 1142  1862816  5512198
//pnd3d ../voxlap/vxl/untitled.vxl    12425384  3288588  712   650132  1852926 (edgeissol=61)
//pnd3d ../voxlap/vxl/untitled.kvo     5788180  2265217  501   650132  1852926 (edgeissol=61)
//pnd3d ../cavedemo/data/seed12t1.vxl 13501288          1600  2242535  5826251
//pnd3d ../voxlap/vxl/canyon.vxl      14321980          1642  2343131  6200467
//pnd3d ../cavedemo/data/seed7t3.vxl  20994016          2029  2720109  6936963
//pnd3d ../cavedemo/data/seed5t1.vxl  21624788          1938  2697394  7471136
//pnd3d ../cavedemo/data/seed9t1.vxl  22572000          2241  2914854  7228081
//pnd3d ../cavedemo/data/seed7t2.vxl  23240136          2285  2963912  7343796
//pnd3d ../cavedemo/data/seed7t1.vxl  23667988          2316  2980566  7395794
//pnd3d ../voxlap/vxl/coreef.vxl      25765080          2339  3175828  8157288
//pnd3d ../voxlap/vxl/test.vxl        28534696          2636  3426678  8557623
//pnd3d ../voxlap/vxl/lung.vxl        33680036          2898  3560093  9243838
//pnd3d ../voxlap/vxl/voxelstein.vxl  56591876          3488  4498673 13687906
//pnd3d ../voxlap/vxl/voxelstein.kvo  33458913
//
//rule of thumb: sur.num = nod.num*~{2.5-3}

//03/15/2012: idea for new PIXMETH?
//   //pix buffer:
//psurf:29 (points to start of 2x2x2 surf list)
//popcnt:3 (psurf offset)
//
//   //surf format:
//short x, y, z; char chi, filler;
//unsigned char b, g, r, a; signed char norm[3], tex;
//unsigned char b, g, r, a; signed char norm[3], tex;
//unsigned char b, g, r, a; signed char norm[3], tex;
//..
//
//pros:
//* 32bit transfer is nice
//* supports up to 16-bit lsid
//* generation of pix val in dorect may be negligibly faster
//
//cons:
//* requires ~25% more memory on GPU/surf list texture ://
//* surf list limited to 4GBy, but still 2x better than other pixmeth

//--------------------------------------------------
//   clearscreen(oct_fogcol);
//   zbuf: set all to 0x7f7f7f7f
//   {
//      if (first sprite) { gcbit:set cubearea to 1 }
//                   else { gcbit:set cubearea to 1 if zbuf>front }
//      *_dorect: if (gcbit[?]) gcbit[?] = 0; if (z<zbuf[?]) .. )
//   }
//   render other stuff here
//--------------------------------------------------

#include "octree_renderer.hpp"

#include "octree_modify.hpp"
#include "octree_physics.hpp"

#include "graphics.hpp"

#include "utilities/bits.hpp"
#include "utilities/colour.hpp"
#include "utilities/gl.hpp"
#include "utilities/maths.hpp"
#include "utilities/popcount.hpp"
#include "utilities/self_modifying_code.hpp"
#include "utilities/timestamp.hpp"
#include "utilities/thread_pool.hpp"
#include "utilities/window.hpp"

#include <algorithm>
#include <cstdarg>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <emmintrin.h>
#include <malloc.h>

#define MAXYDIM 2560

static int gdrawclose = 1;
static HINSTANCE ghinst;

static tiletype gdd;

static int gzbufoff, gddz;

static float ghx, ghy, ghz;
static int ighx, ighy, ighz;
static point3f girig, gidow, gifor, gipos;

float oct_fogdist = 16.0; //1e32;
unsigned int oct_fogcol = 0x605040;

//glBindTexture   id's (see also oct_t.octid, oct_t.tilid)
unsigned int gpixtexid;

//INT_PTR to support CPU as pointer & GPU as index for drawpol()
INT_PTR gskyid = 0;

unsigned int gfont6x8id;
unsigned int gtextbufid;

//glBufferTexture id's (see also oct_t.gbufid)
unsigned int gpixbufid;

static point3f fgxmul;
static point3f fgymul;
static point3f fgzmul;
static point3f fgnadd;
static point3f gxmul;
static point3f gymul;
static point3f gzmul;
static point3f gnadd;

static point3f gorig;
static point3i glorig;

static unsigned int *gcbit;
static unsigned int gcbitmax = 0;
static unsigned int gcbpl;
static unsigned int gcbitpl;
static unsigned int gcbitplo5;
static unsigned int gcbitplo3;
static unsigned int *gtbit;
static unsigned int gtbitmal = 0;

INT_PTR gcbitmal = 0;

char *gpixbuf = 0;
int gpixbufmal = 0;
int gpixxdim = 2048;
int  gpixydim = 2048;

void (*shaderfunc)(int,void*) = 0;

//-1=uninited, 0=drawtext6x8,etc.., 1=drawoct()
int gshadermode = -1;
//0=max, 1=60fps, 2=30fps, ..
int swapinterval = 1;

double gznear = 1.0/32768.0;
double gzfar = 256.0;

static int *zbuf = 0;
static int zbufmal = 0;
static INT_PTR zbufoff;

//temp for cpu shader only
static int gimixval;
static int gimulcol;

static int gtiloff[16];
alignas(16) static int dqtiloff[256][4];

//intermediate data between 1st pass & shader
static tiletype gxd;

typedef struct { INT_PTR f; int p, x, y, ltilesid; } tiles_t;
static tiles_t tiles[16]; //mips

#define cvttss2si(f) _mm_cvtt_ss2si(_mm_set_ss(f))

//--------------------------------------------------------------------------------------------------

HDC glhDC;
HGLRC glhRC;

unsigned int shadcur = 0;
//[0]=oct_*, [1]=drawcone, [2]=drawsky, [3]=drawtext, [4]=drawpol(need shader to brighten :/)
unsigned int shadprog[SHADNUM], shadvert[SHADNUM], shadfrag[SHADNUM];

void oct_freetex (INT_PTR ptr)
{
    if constexpr (oct_usegpu) {
        ((PFNGLDELETETEXTURES)glfp[glDeleteTextures])(1,(unsigned int *)&ptr);
    }
    else {
        free((void *)((tiles_t *)ptr)->f);
        free((void *)ptr);
    }
}

INT_PTR oct_loadtex (INT_PTR ptr, int xsiz, int ysiz, int flags)
{
    if constexpr (oct_usegpu) {
        unsigned int hand;
        ((PFNGLGENTEXTURES)glfp[glGenTextures])(1,(unsigned int *)&hand);
        kglalloctex(hand,(void *)ptr,xsiz,ysiz,1,flags);
        return((INT_PTR)hand);
    }
    else {
        tiles_t *ptt;
        int y, xtraline4bilin;

        ptt = (tiles_t *)malloc(sizeof(tiles_t));
        ptt->x = xsiz; ptt->y = ysiz;

        xtraline4bilin = (!(flags&KGL_NEAREST));
        ptt->p = (ptt->x+xtraline4bilin)*sizeof(int); ptt->f = (INT_PTR)malloc((ptt->y+xtraline4bilin)*ptt->p + xtraline4bilin*16/*for movups in bilinear hack*/);
        for(y=0;y<ptt->y;y++) memcpy((void *)(ptt->p*y + ptt->f),(void *)(y*xsiz*4 + ptr),xsiz*sizeof(int));

        if (!(flags&(KGL_CLAMP|KGL_CLAMP_TO_EDGE)))
        {
            for(y=0;y<ptt->y;y++) *(int *)(ptt->p*y + (ptt->x<<2) + ptt->f) = *(int *)(ptt->p*y + ptt->f);
            memcpy((void *)(ptt->p*ptt->y + ptt->f),(void *)ptt->f,(ptt->x+1)<<2);
        }
        else
        {
            for(y=0;y<ptt->y;y++) *(int *)(ptt->p*y + (ptt->x<<2) + ptt->f) = *(int *)(ptt->p*y + ((ptt->x-1)<<2) + ptt->f);
            memcpy((void *)(ptt->p*ptt->y + ptt->f),(void *)(ptt->p*(ptt->y-1) + ptt->f),(ptt->x+1)<<2);
        }

        for(ptt->ltilesid=1;(1<<ptt->ltilesid)<ptt->x;ptt->ltilesid<<=1);
        return((INT_PTR)ptt);
    }
}

char vshad_drawsky[] = R"(
    varying vec4 t;
    void main()
    {
        gl_Position = ftransform();
        t = gl_MultiTexCoord0;
    }
)";

char fshad_drawsky[] = R"(
    varying vec4 t;
    uniform samplerCube tex0;
    uniform vec4 fogcol;
    void main()
    {
        gl_FragColor = textureCube(tex0, t.xyz) * fogcol;
    }
)";

char vshadasm_drawsky[] = R"(!!ARBvp1.0
    PARAM ModelViewProj[4] = {state.matrix.mvp};
    TEMP t;
    DP4 t.x, ModelViewProj[0], vertex.position;
    DP4 t.y, ModelViewProj[1], vertex.position;
    DP4 t.z, ModelViewProj[2], vertex.position;
    DP4 t.w, ModelViewProj[3], vertex.position;
    MOV result.position, t;
    MOV result.texcoord[0], vertex.texcoord[0];
    END
)";

char fshadasm_drawsky[] = R"(!!ARBfp1.0
    TEMP t;
    TEX t, fragment.texcoord[0], texture[0], CUBE;
    MUL result.color, t, program.local[0];
    END
)";

enum class platform_type {
    cpu_cpp,
    cpu_asm,
    gpu_glsl,
    gpu_asm
};

template <platform_type platform>
static void drawsky(int xres, int yres);

template <>
void drawsky<platform_type::cpu_cpp>(int xres, int yres) {
    auto sky_draw_mt = [](int sy, void *) -> void {
        INT_PTR wptr, rptr2, picf;
        float f, fx0, fy0, fz0, fx, fy, fz, fx2, fy2, fz2, fxi, fyi, fzi, fhpic;
        int i, j, sx, sx0, sx1, sxe2, sxe, iu, iv, iu2, iv2, iui, ivi, picsizm1, icross[6], icrossn, fac, hpic, col, nfogcol;
        #define LINTERPSIZ 4

        tiles_t *skypic = (tiles_t *)gskyid;

        hpic = (skypic->x<<(16-1)); fhpic = (float)hpic;
        picsizm1 = ((skypic->x)<<16)-1;

        wptr = gdd.p*sy + gdd.f;
        nfogcol = ((oct_fogcol&0xfefefe)>>1);

        fx0 = (0-ghx)*girig.x + (sy-ghy)*gidow.x + ghz*gifor.x;
        fy0 = (0-ghx)*girig.y + (sy-ghy)*gidow.y + ghz*gifor.y;
        fz0 = (0-ghx)*girig.z + (sy-ghy)*gidow.z + ghz*gifor.z;

        icrossn = 0;
        f =  (fx0-fy0) / (girig.y - girig.x); if (fabs(fz0+girig.z*f) < fabs(fx0+girig.x*f)) { icross[icrossn] = cvttss2si(f)+1; if ((unsigned)icross[icrossn] < (unsigned)gdd.x) icrossn++; }
        f =  (fx0-fz0) / (girig.z - girig.x); if (fabs(fy0+girig.y*f) < fabs(fx0+girig.x*f)) { icross[icrossn] = cvttss2si(f)+1; if ((unsigned)icross[icrossn] < (unsigned)gdd.x) icrossn++; }
        f =  (fy0-fz0) / (girig.z - girig.y); if (fabs(fx0+girig.x*f) < fabs(fy0+girig.y*f)) { icross[icrossn] = cvttss2si(f)+1; if ((unsigned)icross[icrossn] < (unsigned)gdd.x) icrossn++; }
        f = -(fx0+fy0) / (girig.y + girig.x); if (fabs(fz0+girig.z*f) < fabs(fx0+girig.x*f)) { icross[icrossn] = cvttss2si(f)+1; if ((unsigned)icross[icrossn] < (unsigned)gdd.x) icrossn++; }
        f = -(fx0+fz0) / (girig.z + girig.x); if (fabs(fy0+girig.y*f) < fabs(fx0+girig.x*f)) { icross[icrossn] = cvttss2si(f)+1; if ((unsigned)icross[icrossn] < (unsigned)gdd.x) icrossn++; }
        f = -(fy0+fz0) / (girig.z + girig.y); if (fabs(fx0+girig.x*f) < fabs(fy0+girig.y*f)) { icross[icrossn] = cvttss2si(f)+1; if ((unsigned)icross[icrossn] < (unsigned)gdd.x) icrossn++; }
        for(i=1;i<icrossn;i++)
            for(j=0;j<i;j++)
                if (icross[i] < icross[j]) { sx = icross[i]; icross[i] = icross[j]; icross[j] = sx; }

        sx0 = 0; sx1 = gdd.x;

        while ((icrossn > 0) && (icross[icrossn-1] >= sx1)) icrossn--;
        sx = sx1-1;
        do
        {
            sxe = sx0;
            if ((icrossn > 0) && (icross[icrossn-1] > sx0)) { icrossn--; sxe = icross[icrossn]; }

            f = (float)(sx+sxe) * 0.5f;
            fx = girig.x*f + fx0; fx2 = (sx-1)*girig.x + fx0;
            fy = girig.y*f + fy0; fy2 = (sx-1)*girig.y + fy0;
            fz = girig.z*f + fz0; fz2 = (sx-1)*girig.z + fz0;
            if (fabs(fz) > std::max(fabs(fx),fabs(fy))) fac = (*(int *)&fz<0)*2+0;
            else if (fabs(fy) > fabs(fx))          fac = (*(int *)&fy>0)+4  ;
            else                                   fac = (*(int *)&fx<0)*2+1;
            picf = skypic->x*fac*skypic->p + skypic->f;

            switch(fac)
            {
                default:
                case 0: fx = fx2; fy = fy2; fz = fz2; fxi = girig.x; fyi = girig.y; fzi = girig.z; break;
                case 1: fx =-fz2; fy = fy2; fz = fx2; fxi =-girig.z; fyi = girig.y; fzi = girig.x; break;
                case 2: fx = fx2; fy =-fy2; fz = fz2; fxi = girig.x; fyi =-girig.y; fzi = girig.z; break;
                case 3: fx =-fz2; fy =-fy2; fz = fx2; fxi =-girig.z; fyi =-girig.y; fzi = girig.x; break;
                case 4: fx =-fx2; fy =-fz2; fz = fy2; fxi =-girig.x; fyi =-girig.z; fzi = girig.y; break;
                case 5: fx = fx2; fy =-fz2; fz = fy2; fxi = girig.x; fyi =-girig.z; fzi = girig.y; break;
            }

            f = fhpic/fz;
            iu = cvttss2si(fx*f) + hpic; iu = std::min(std::max(iu,0),picsizm1);
            iv = cvttss2si(fy*f) + hpic; iv = std::min(std::max(iv,0),picsizm1);
            while (sx-sxe >= (1<<LINTERPSIZ))
            {
                fx -= fxi*(1<<LINTERPSIZ); fy -= fyi*(1<<LINTERPSIZ); fz -= fzi*(1<<LINTERPSIZ);
                f = fhpic/fz;
                iu2 = cvttss2si(fx*f) + hpic; iu2 = std::min(std::max(iu2,0),picsizm1); iui = ((iu2-iu)>>LINTERPSIZ);
                iv2 = cvttss2si(fy*f) + hpic; iv2 = std::min(std::max(iv2,0),picsizm1); ivi = ((iv2-iv)>>LINTERPSIZ);
                for(sxe2=std::max(sx-(1<<LINTERPSIZ),sxe);sx>sxe2;sx--,iu+=iui,iv+=ivi)
                {
                    rptr2 = (iv>>16)*skypic->p + (iu>>16)*4 + picf;
                    col = *(int *)rptr2; //nearest
                    *(int *)((sx<<2)+wptr) = mulcol(col,nfogcol);
                }
            }
            for(;sx>=sxe;sx--,fx-=fxi,fy-=fyi,fz-=fzi)
            {
                f = fhpic/fz;
                iu = cvttss2si(fx*f) + hpic; iu = std::min(std::max(iu,0),picsizm1);
                iv = cvttss2si(fy*f) + hpic; iv = std::min(std::max(iv,0),picsizm1);
                rptr2 = (iv>>16)*skypic->p + (iu>>16)*4 + picf;
                col = *(int *)rptr2; //nearest
                *(int *)((sx<<2)+wptr) = mulcol(col,nfogcol);
            }
        } while (sx >= sx0);
    };
    htrun(sky_draw_mt,0,0,yres,true);
}

template <>
void drawsky<platform_type::cpu_asm>(int xres, int yres) {
    drawsky<platform_type::cpu_cpp>(xres, yres);
}

template <>
void drawsky<platform_type::gpu_glsl>(int xres, int yres) {
    // Set shader id.
    shadcur = 2;

    // Enable shader.
    ((PFNGLUSEPROGRAMPROC)glfp[glUseProgram])(shadprog[shadcur]);

    // Load fog colour.
    const float rgba[4] = {
        static_cast<float>((oct_fogcol >> 16) & 255) / 128.0f,
        static_cast<float>((oct_fogcol >> 8) & 255) / 128.0f,
        static_cast<float>((oct_fogcol >> 0) & 255) / 128.0f,
        static_cast<float>((oct_fogcol >> 24) & 255) / 128.0f
    };
    kglUniform4f(kglGetUniformLoc(shadprog[shadcur], "fogcol"),rgba[0], rgba[1], rgba[2], rgba[3]);

    // Configure projection and modelview matrices.
    ((PFNGLMATRIXMODE)glfp[glMatrixMode])(GL_PROJECTION);
    ((PFNGLLOADIDENTITY)glfp[glLoadIdentity])();
    ((PFNGLFRUSTUM)glfp[glFrustum])(-gznear, gznear, -(float)yres / (float)xres * gznear, (float)yres / (float)xres * gznear, gznear, gzfar);
    ((PFNGLMATRIXMODE)glfp[glMatrixMode])(GL_MODELVIEW);
    ((PFNGLLOADIDENTITY)glfp[glLoadIdentity])();

    // Enable depth testing.
    ((PFNGLENABLE)glfp[glEnable])(GL_DEPTH_TEST);

    kglActiveTexture(0);
    ((PFNGLBLENDFUNC)glfp[glBlendFunc])(GL_ONE, GL_NONE);
    ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_CUBE_MAP, gskyid);

    float g = ghz / ghx;
    // Highest distance without falling behind Z .. but why 221/256?
    float f = static_cast<float>(gzfar * (220.999992 / 256.0));

    ((PFNGLCOLOR4UB)glfp[glColor4ub])(255, 255, 255, 255);

    ((PFNGLBEGIN)glfp[glBegin])(GL_QUADS);
    ((PFNGLTEXCOORD3F)glfp[glTexCoord3f])(- girig.x - gidow.x + gifor.x * g, + girig.y + gidow.y - gifor.y * g, - girig.z - gidow.z + gifor.z * g);
    ((PFNGLVERTEX3F)glfp[glVertex3f])(-f, +f, -f);
    ((PFNGLTEXCOORD3F)glfp[glTexCoord3f])(+ girig.x - gidow.x + gifor.x * g, - girig.y + gidow.y - gifor.y * g, + girig.z - gidow.z + gifor.z * g);
    ((PFNGLVERTEX3F)glfp[glVertex3f])(+f, +f, -f);
    ((PFNGLTEXCOORD3F)glfp[glTexCoord3f])(+ girig.x + gidow.x + gifor.x * g, - girig.y - gidow.y - gifor.y * g, + girig.z + gidow.z + gifor.z * g);
    ((PFNGLVERTEX3F)glfp[glVertex3f])(+f, -f, -f);
    ((PFNGLTEXCOORD3F)glfp[glTexCoord3f])(- girig.x + gidow.x + gifor.x * g, + girig.y - gidow.y - gifor.y * g, - girig.z + gidow.z + gifor.z * g);
    ((PFNGLVERTEX3F)glfp[glVertex3f])(-f, -f, -f);
    ((PFNGLEND)glfp[glEnd])();

    shadcur = 0;
    ((PFNGLUSEPROGRAMPROC)glfp[glUseProgram])(shadprog[shadcur]);
}

template <>
void drawsky<platform_type::gpu_asm>(int xres, int yres) {
    shadcur = 2;
    ((PFNGLENABLE)glfp[glEnable])(GL_VERTEX_PROGRAM_ARB);
    ((PFNGLENABLE)glfp[glEnable])(GL_FRAGMENT_PROGRAM_ARB);
    ((PFNGLBINDPROGRAMARBPROC)glfp[glBindProgramARB])(GL_VERTEX_PROGRAM_ARB  ,shadvert[shadcur]);
    ((PFNGLBINDPROGRAMARBPROC)glfp[glBindProgramARB])(GL_FRAGMENT_PROGRAM_ARB,shadfrag[shadcur]);

    // Load fog colour.
    const float rgba[4] = {
        static_cast<float>((oct_fogcol >> 16) & 255) / 128.0f,
        static_cast<float>((oct_fogcol >> 8) & 255) / 128.0f,
        static_cast<float>((oct_fogcol >> 0) & 255) / 128.0f,
        static_cast<float>((oct_fogcol >> 24) & 255) / 128.0f
    };
    kglProgramLocalParam(0, rgba[0], rgba[1], rgba[2], rgba[3]);

    // Configure projection and modelview matrices.
    ((PFNGLMATRIXMODE)glfp[glMatrixMode])(GL_PROJECTION);
    ((PFNGLLOADIDENTITY)glfp[glLoadIdentity])();
    ((PFNGLFRUSTUM)glfp[glFrustum])(-gznear, gznear, -(float)yres / (float)xres * gznear, (float)yres / (float)xres * gznear, gznear, gzfar);
    ((PFNGLMATRIXMODE)glfp[glMatrixMode])(GL_MODELVIEW);
    ((PFNGLLOADIDENTITY)glfp[glLoadIdentity])();

    // Enable depth testing.
    ((PFNGLENABLE)glfp[glEnable])(GL_DEPTH_TEST);

    kglActiveTexture(0);
    ((PFNGLBLENDFUNC)glfp[glBlendFunc])(GL_ONE,GL_NONE);
    ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_CUBE_MAP,gskyid);

    float g = ghz/ghx;
    // highest distance without falling behind Z .. but why 221/256?
    float f = static_cast<float>(gzfar * (220.999992 / 256.0));

    ((PFNGLCOLOR4UB)glfp[glColor4ub])(255, 255, 255, 255);

    ((PFNGLBEGIN)glfp[glBegin])(GL_QUADS);
    ((PFNGLTEXCOORD3F)glfp[glTexCoord3f])(- girig.x - gidow.x + gifor.x * g, + girig.y + gidow.y - gifor.y * g, - girig.z - gidow.z + gifor.z * g);
    ((PFNGLVERTEX3F)glfp[glVertex3f])(-f, +f, -f);
    ((PFNGLTEXCOORD3F)glfp[glTexCoord3f])(+ girig.x - gidow.x + gifor.x * g, - girig.y + gidow.y - gifor.y * g, + girig.z - gidow.z + gifor.z * g);
    ((PFNGLVERTEX3F)glfp[glVertex3f])(+f, +f, -f);
    ((PFNGLTEXCOORD3F)glfp[glTexCoord3f])(+ girig.x + gidow.x + gifor.x * g, - girig.y - gidow.y - gifor.y * g, + girig.z + gidow.z + gifor.z * g);
    ((PFNGLVERTEX3F)glfp[glVertex3f])(+f, -f, -f);
    ((PFNGLTEXCOORD3F)glfp[glTexCoord3f])(- girig.x + gidow.x + gifor.x * g, + girig.y - gidow.y - gifor.y * g, - girig.z + gidow.z + gifor.z * g);
    ((PFNGLVERTEX3F)glfp[glVertex3f])(-f, -f, -f);
    ((PFNGLEND)glfp[glEnd])();
}

static void drawsky(int xres, int yres)
{
    if constexpr (oct_usegpu && oct_useglsl) {
        drawsky<platform_type::gpu_glsl>(xres, yres);
    }
    else if constexpr (oct_usegpu) {
        drawsky<platform_type::gpu_asm>(xres, yres);
    }
    else if constexpr (oct_useasm) {
        drawsky<platform_type::cpu_asm>(xres, yres);
    }
    else {
        drawsky<platform_type::cpu_cpp>(xres, yres);
    }
}

//--------------------------------------------------------------------------------------------------

//  0: 32bit: psurf:32                           :)    :)   :)   col_only

char vshad_drawoct[] = R"(
    varying vec2 tv;
    varying vec4 v;
    uniform vec4 vmulx, vmuly, vadd;
    void main()
    {
        gl_Position = ftransform();
        tv = gl_MultiTexCoord0.xy - vec2(1.0 / 8192.0, 0.0);
        v = gl_MultiTexCoord0.x * vmulx + gl_MultiTexCoord0.y * vmuly + vadd;
    };
)";

char fshad_drawoct[] = R"(
    varying vec2 tv;
    varying vec4 v;
    uniform sampler2D tex0, tex1, tex2, tex3;
    uniform vec4 pos, vmulx, vmuly, vadd, mulcol, fogcol;
    uniform float rxsid, rysid, mipoffs, maxmip, mixval, gsideshade[6];
    uniform int axsidm1, lxsidm1; // Shader fails if axsidm1 is renamed to xsidm1!
    uniform float depthmul, depthadd, fogdist;
    void main()
    {
        // Calculate i, the psurf stored in the first texture.
        ivec4 o1 = ivec4(texture2D(tex1,tv.xy) * 255.5);
        int i = (o1.r + o1.g*256 + o1.b*65536 + o1.a*16777216);
        if (i == 0) {
            discard;
        }

        //i = packUnorm4x8(texture2D(tex1, tv.xy));
        //if (i == 0u) {
        //    discard;
        //}

        // Calculate the location of the psurf in the second texture.
        vec2 v2 = vec2(float(i & axsidm1) * 2.0 * rxsid + 0.5 * rxsid, float(i >> lxsidm1) * rysid + 0.5 * rysid);

        // The surface holds the colour.
        vec4 c = texture2D(tex2, v2).bgra;
        c.a = 1.0;

        // The surface also holds the normal.
        vec4 p = texture2D(tex2, v2 + vec2(rxsid, 0.0));

        // Convert the surface normal into +-1.
        vec3 n = p.xyz - vec3(lessThan(vec4(0.5), p));

        // Apply the normal to light surfaces.
        c *= dot(n, vec3(1.0)) + 1.0;

        // Fix the alpha.
        c.a = 1.0;

        gl_FragDepth = depthmul + depthadd;
        gl_FragColor = mix(c * mulcol, fogcol, min(fogdist, 1.0));
        //gl_FragColor = c * mulcol; //NVidia shader compiler bug?
    }
)";

char fshad_drawoct_type2[] = R"(
    varying vec2 tv;
    varying vec4 v;
    uniform sampler2D tex0, tex1, tex2, tex3;
    uniform vec4 pos, vmulx, vmuly, vadd, mulcol, fogcol;
    uniform float rxsid, rysid, mipoffs, maxmip, mixval, gsideshade[6];
    uniform int axsidm1, lxsidm1; // Shader fails if axsidm1 is renamed to xsidm1!
    uniform float depthmul, depthadd, fogdist;
    void main()
    {
        vec4 c, n, p, w;
        ivec4 o0, o1, ip;
        vec2 v2;
        float f, g, h, gotf;
        int i;

        o0 = ivec4(texture2D(tex1,tv                     )*255.5);
        o1 = ivec4(texture2D(tex1,tv+vec2(1.0/4096.0,0.0))*255.5);
        i = (o1.r>>4) + (o1.g<<4) + (o1.b<<12) + (o1.a<<20);
        if (i == 0) {
            discard;
        }
        ip.x = o0.r + ((o0.g&15)<<8);
        ip.y = (o0.g>>4) + (o0.b<<4);
        ip.z = o0.a + ((o1.r&15)<<8);

        p = (vec4(lessThan(v,vec4(0.0)))+ip-pos)/v;
        n = vec4(lessThan(max(p.yxxw,p.zzyw),p));
        gotf = dot(p.xyz,n.xyz);
        w = v*gotf+pos-vec4(ip); w.xy = mix(w.xy,w.zz,n.xy);
        //   if (max(abs(w.x-.5),abs(w.y-.5)) > 0.5) discard;

        v2 = vec2(float(i&axsidm1)*2.0*rxsid + 0.5*rxsid,float(i>>lxsidm1)*rysid + 0.5*rysid);
        c = texture2D(tex2,v2                ).bgra;
        c.a = 1.0;
        p = texture2D(tex2,v2+vec2(rxsid,0.0));
        g = p.w*255.0;
        c *= dot(p.xyz-vec3(lessThan(vec4(.5),p)),vec3(1.0))+1.0;

        #if 0
            #ifdef GL_ARB_shader_texture_lod
                if (true) {
                    c = mix(c,texture2DLod(tex0,vec2(clamp(w.x,0.0,1.0)/32.0 +  mod(g,16.0)/16.0 + 1.0/64.0, clamp(w.y,0.0,1.0)/32.0 +floor(g/16.0)/16.0 + 1.0/64.0),min(log2(gotf/abs(dot(v.xyz,n.xyz)))+mipoffs,maxmip)),mixval);
                }
                else {
                    c *= texture2DLod(tex0,vec2(clamp(w.x,0.0,1.0)/32.0 +  mod(g,16.0)/16.0 + 1.0/64.0, clamp(w.y,0.0,1.0)/32.0 +floor(g/16.0)/16.0 + 1.0/64.0),min(log2(gotf/abs(dot(v.xyz,n.xyz)))+mipoffs,maxmip));
                }
            #else
                c = mix(c,texture2D   (tex0,vec2(clamp(w.x,0.0,1.0)/32.0 +  mod(g,16.0)/16.0 + 1.0/64.0, clamp(w.y,0.0,1.0)/32.0 +floor(g/16.0)/16.0 + 1.0/64.0)),mixval);
            #endif
        #else
            //Test code for (oct_usefilter == 3); not needed for newer GLSL..
            float lod = clamp(log2(gotf/abs(dot(v.xyz,n.xyz)))+mipoffs,0.0,maxmip);
            vec2 tx = vec2(clamp(w.x,0.0,1.0)/64.0 +  mod(g,16.0)/32.0 + 65.0/128.0, clamp(w.y,0.0,1.0)/32.0 +floor(g/16.0)/16.0 +  1.0/ 64.0) * pow(0.5,floor(lod));
            w = mix(texture2D(tex0,tx),texture2D(tex0,tx*.5),mod(lod,1.0));
            if (true)
                c = mix(c,w,mixval);
            else
                c *= w;
            c = mix(c,fogcol,mixval);
        #endif

        c *= mulcol;
        c.rgb -= vec3(gsideshade[int(dot(vec3(lessThan(v.xyz,vec3(0.0)))+vec3(0.0,2.0,4.0),n.xyz))]);
        c = mix(c,fogcol,min(gotf*fogdist,1.0));
        gl_FragDepth = depthmul/gotf + depthadd;
        gl_FragColor = c;
    }
)";

char vshadasm_drawoct[] = R"(!!ARBvp1.0
    PARAM ModelViewProj[4] = {state.matrix.mvp};
    TEMP t;
    DP4 t.x, ModelViewProj[0], vertex.position;
    DP4 t.y, ModelViewProj[1], vertex.position;
    DP4 t.z, ModelViewProj[2], vertex.position;
    DP4 t.w, ModelViewProj[3], vertex.position;
    MOV result.position, t;
    MOV result.color, vertex.color;
    MOV result.texcoord[0], vertex.texcoord[0];
    END
)";

char fshadasm_drawoct[] = R"(!!ARBfp1.0
    ATTRIB tv = fragment.texcoord[0];
    PARAM sidmul0 = program.local[0]; #{ 1/xsid, 2/xsid,0,0}
    PARAM sidmul1 = program.local[1]; #{ 2     , 1/ysid,0,0}
    PARAM sidadd1 = program.local[2]; #{.5/xsid,.5/ysid,0,0}
    PARAM depth   = program.local[3]; #{depthmul,depthadd,0,0}
    PARAM vmulx   = program.local[4]; #{vmulx.x,vmulx.y,vmulx.z,0}
    PARAM vmuly   = program.local[5]; #{vmuly.x,vmuly.y,vmuly.z,0}
    PARAM vadd    = program.local[6]; #{vadd.x ,vadd.y ,vadd.z ,0}
    PARAM pos     = program.local[7]; #{pos.x  ,pos.y  ,pos.z  ,0}
    PARAM mixnfog = program.local[8]; #{0,1/oct_fogdist,mixval}
    PARAM mulcol  = program.local[9]; #mulcol
    PARAM fogcol  = program.local[10]; #fogcol
    PARAM gsidsh0 = program.local[11]; #{gsideshade[0],gsideshade[1],gsideshade[2],0}
    PARAM gsidsh1 = program.local[12]; #{gsideshade[3],gsideshade[4],gsideshade[5],0}
    TEMP t, u, i, c, p, q;

    TEX u, tv, texture[1], 2D;                #u = floor(texture2D(tex1,t.xy)*255.5);
    MUL u, u, {255.5,255.5,255.5,255.5};
    FLR u, u;
    DP4 i, u, {1.0,256.0,65536.0,16777216.0}; #i = t.r + t.g*256.0 + t.b*65536.0 + t.a*16777216.0;
    SGE t, -i, {0.0};                         #if (i == 0.0) discard;
    KIL -t.x;
    MUL t.y, i.x, sidmul0;                    #t = vec2(frac(f*rxsid)*2.0 + 0.5*rxsid,floor(f*rxsid*2)*rysid + 0.5*rysid);
    DP4 i, u, {1.0,256.0,0.0,0.0};            #Trick required for 32-bit precision: limits i to 16-bit precision to avoid float truncation
    MUL t.x, i.x, sidmul0.x;
    FLR t.y, t.y;
    MAD t, t, sidmul1, sidadd1;
    ADD q, t, sidmul0.xwww;
    TEX c, t, texture[2], 2D;                 #c = texture2D(tex2,t).bgra;
    TEX p, q, texture[2], 2D;                 #p = texture2D(tex2,t+vec2(rxsid,0.0));
    SLT t, {.5,.5,.5,0.0}, p;                 #c *= dot(p.xyz-vec3(lessThan(vec4(.5),p)),vec3(1.0))+1.0;
    SUB t, p, t;
    DPH t, t, {1.0,1.0,1.0,1.0};
    MUL c, t, c.bgra;
    MOV c.a, {1.0};
    #MOV result.depth, depth.y;                #gl_FragDepth = depthadd;
    MOV result.depth.z, {0.0,0.0,0.0,0.0};    #gl_FragDepth = depthadd;
    MUL result.color, c, mulcol;              #gl_FragColor = c*mulcol;
    END
)";



//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------


static void oct_debugprintftree (oct_t *loct)
{
    int i, ls, ind, stkind[OCT_MAXLS], stknum[OCT_MAXLS];

    ls = loct->lsid-1; stkind[ls] = loct->head; stknum[ls] = 1;
    while (1)
    {
        ind = stkind[ls]; stkind[ls]++; stknum[ls]--; //2sibly
        if (ls >= 0) //2child
        {
            octv_t *ptr = (octv_t *)loct->nod.buf;
            printf("nod[%5d]: ls=%d, chi:%02x, sol:%02x, ind:%d\n",ind,ls,ptr[ind].chi,ptr[ind].sol,ptr[ind].ind);
            i = popcount8(ptr[ind].chi);
            ls--; stkind[ls] = ptr[ind].ind; stknum[ls] = i;
        } else { printf("sur[%5d]\n",ind); } //surf
        while (stknum[ls] <= 0) { ls++; if (ls >= loct->lsid) return; } //2parent
    }
}



#if 0
#define RUBBLE_NUMREG 4
static unsigned short rubble_boxsum[RUBBLE_NUMREG][32+1][32+1][32+1];
static int rubble_inited = 0;
static float rubble_getsc (point3i *pt, int numpt, int x, int y, int z, int *bestj)
{
    float f, bestf, dx, dy, dz;
    int i, j;

    bestf = 1e32f; (*bestj) = 0;
    for(j=0;j<RUBBLE_NUMREG;j++)
    {
        f = 0.f;
        for(i=0;i<numpt;i++)
        {
            dx = (float)(pt[j*numpt+i].x-x)+.5;
            dy = (float)(pt[j*numpt+i].y-y)+.5;
            dz = (float)(pt[j*numpt+i].z-z)+.5;
            f += 1.f/(dx*dx + dy*dy + dz*dz);
        }
        if (f < bestf) { bestf = f; (*bestj) = j; }
    }
    return(f);
}
static void rubble_genboxsums() //Voxlap rubble.. generates 4 boxsums
{
    #define RUBBLE_NUMPT 7
    char bmp[32][32][32];
    point3i pt[RUBBLE_NUMREG*RUBBLE_NUMPT];
    float f, maxf;
    int i, j, k, x, y, z, xi;

        //for each region, generate n points within 11 of the center
    srand(1);
    for(j=0;j<RUBBLE_NUMREG;j++)
        for(i=0;i<RUBBLE_NUMPT;i++)
        {
            k = j*RUBBLE_NUMPT+i;
            do
            {
                pt[k].x = (rand()&31)-16;
                pt[k].y = (rand()&31)-16;
                pt[k].z = (rand()&31)-16;
            } while (pt[k].x*pt[k].x + pt[k].y*pt[k].y + pt[k].z*pt[k].z >= 11*11);
            pt[k].x += 16; pt[k].y += 16; pt[k].z += 16;
        }

        //find highest score along edges of cube
    maxf = 0.f;
    for(z=32-1;z>=0;z--)
        for(y=32-1;y>=0;y--)
        {
            if ((y == 0) || (z == 0) || (y == 32-1) || (z == 32-1)) xi = 1; else xi = 32-1;
            for(x=32-1;x>=0;x-=xi) maxf = std::max(maxf,rubble_getsc(pt,RUBBLE_NUMPT,x,y,z,&j));
        }

        //write region buffer; 0 for blank, 1-NUMREG for object
    memset(bmp,255,sizeof(bmp));
    for(z=32-1;z>=0;z--)
        for(y=32-1;y>=0;y--)
            for(x=32-1;x>=0;x--)
                { f = rubble_getsc(pt,RUBBLE_NUMPT,x,y,z,&j); if (f > maxf) bmp[z][y][x] = j; }

    for(k=0;k<RUBBLE_NUMREG;k++) brush_bmp_calcsum(&rubble_boxsum[k][0][0][0],&bmp[0][0][0],k,32,32,32);
}
#endif

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------

//WARNING:oct_drawoct() asm depends on this exact structure!
typedef struct
{
    float x0, y0, z0, dumzero;
    float xv, yv, zv, cosang2;
} isintcs_t;

//cx0,cy0,cz0: cone focus
//cxv,cyv,czv: cone center vector direction, normalized
//       crad: cone angle (center vector to edge)
//         sr: sphere radius
static void isint_cone_sph_init (isintcs_t *iics, float cx0, float cy0, float cz0, float cxv, float cyv, float czv, float crad, float sr)
{
    float f = sr/sin(crad);
    iics->x0 = cx0-f*cxv;
    iics->y0 = cy0-f*cyv;
    iics->z0 = cz0-f*czv;
    iics->dumzero = 0.f;
    iics->xv = cxv;
    iics->yv = cyv;
    iics->zv = czv;
    f = cos(crad); iics->cosang2 = f*f;
}
    //sx,sy,sz: sphere center
static int isint_cone_sph (isintcs_t *iics, float sx, float sy, float sz)
{
    float f;
    sx -= iics->x0; sy -= iics->y0; sz -= iics->z0;
    f = sx*iics->xv + sy*iics->yv + sz*iics->zv; if (f <= 0.f) return(0);
    return((sx*sx + sy*sy + sz*sz)*iics->cosang2 < f*f);
}
//--------------------------------------------------------------------------------------------------

typedef struct { int x0, y0, x1, y1; } rect_t;
typedef struct { float z, z2, x, y; } lpoint4d;
typedef struct { short x, y, z, dum; } stk2_t;
alignas(16) static rect_t grect[(4+4)*(4+4)+1024];
alignas(16) static lpoint4d gposadd[8][OCT_MAXLS];
alignas(8) static stk2_t gklslut[8][OCT_MAXLS], glandlut[OCT_MAXLS], gldirlut[OCT_MAXLS], glcenadd, glcensub, gkls0[8];
static float cornmin[OCT_MAXLS], cornmax[OCT_MAXLS], scanmax[OCT_MAXLS];
alignas(16) static isintcs_t giics[OCT_MAXLS];
static int ordlut[8][256]; //ordered list, 4 bits per cell, compacted to skip unused cells

//  4------5       z
//  |\     |\       \
//  | 0----+-1       *----x
//  | |    | |       |
//  6-+----7 |       |
//   \|     \|       y
//    2------3
static const char vanlut[27][3] =
{
    1,2,3, 3,2,1, 3,0,1, // 0 -  2  // 4--5  z
    2,1,3, 3,2,2, 3,0,2, // 3 -  5  // 0--1   \--x
    2,1,0, 2,3,0, 0,3,2, // 6 -  8  // 6--7   |
    1,3,2, 2,2,1, 2,0,1, // 9 - 11  // 2--3   y
    1,3,1, 0,0,0, 0,2,0, //12 - 14
    3,1,0, 0,0,3, 0,2,3, //15 - 17
    0,3,2, 0,1,2, 2,1,0, //18 - 20
    0,3,1, 1,0,0, 1,2,0, //21 - 23
    3,0,1, 1,0,3, 1,2,3, //24 - 26
};

static const char vanflip[27] = {7,6,6, 5,2,6, 5,5,4, 3,4,6, 2,0,0, 5,0,0, 3,3,2, 3,0,0, 1,0,0};
int oct_sideshade[6] = {0x000000,0x030303,0x060606,0x090909,0x0c0c0c,0x0f0f0f};
alignas(16) static int dqgsideshade[6][4];
alignas(16) static float gvanx[4], gvany[4], gvanxmh[4], gvanymh[4], cornvadd[6][4], cornvadd2[27][16] = {0};
alignas(16) static int vause[27][4] = {0};

// Render nearby cubes with classic brute force algo..
static void drawcube_close (oct_t *loct, int xx, int yy, int zz, int ipsurf, int ind, int facemask)
{
    static_cast<void>(loct);
    //                              left     right     up      down    front     back
    static const char fac[6][4] = {4,0,2,6, 1,5,7,3, 4,5,1,0, 2,3,7,6, 0,1,3,2, 5,4,6,7};
    point3f vt[8], vt2[5];
    float f, t, fx, fxinc, sx[5], sy[5];
    INT_PTR sptr;
    alignas(16) char val[16];
    unsigned int *gcptr;
    int i, j, k, kk, i0, i1, sx0, sy0, sx1, sy1, sxe, my0, my1;
    int n, x, y, *iptr, isy[5], xc[2][MAXYDIM], dir[3];
    #define SCISDIST2 1e-6f

        //Rotate&Translate
    vt[0].x = xx*gxmul.x + yy*gymul.x + zz*gzmul.x + gnadd.x;
    vt[0].y = xx*gxmul.y + yy*gymul.y + zz*gzmul.y + gnadd.y;
    vt[0].z = xx*gxmul.z + yy*gymul.z + zz*gzmul.z + gnadd.z;
                             vt[1].x = vt[0  ].x+gxmul.x; vt[1].y = vt[0  ].y+gxmul.y; vt[1].z = vt[  0].z+gxmul.z;
    for(i=2;i<4;i++) { vt[i].x = vt[i-2].x+gymul.x; vt[i].y = vt[i-2].y+gymul.y; vt[i].z = vt[i-2].z+gymul.z; }
    for(   ;i<8;i++) { vt[i].x = vt[i-4].x+gzmul.x; vt[i].y = vt[i-4].y+gzmul.y; vt[i].z = vt[i-4].z+gzmul.z; }

    dir[0] = xx-glorig.x;
    dir[1] = yy-glorig.y;
    dir[2] = zz-glorig.z;

    *(int *)&val[0] = ipsurf;

    for(k=3-1;k>=0;k--)
    {
            //Back-face cull
        if (!dir[k]) continue;
        kk = (dir[k]<0)+k*2;
        if (!(facemask&(1<<kk))) continue;

            //Clip near plane
        n = 0;
        for(i=4-1,j=0;j<4;i=j,j++)
        {
            i0 = fac[kk][i]; i1 = fac[kk][j];
            if (vt[i0].z >= SCISDIST2) { vt2[n] = vt[i0]; n++; }
            if ((vt[i0].z >= SCISDIST2) != (vt[i1].z >= SCISDIST2))
            {
                t = (SCISDIST2-vt[i0].z)/(vt[i1].z-vt[i0].z);
                vt2[n].x = (vt[i1].x-vt[i0].x)*t + vt[i0].x;
                vt2[n].y = (vt[i1].y-vt[i0].y)*t + vt[i0].y;
                vt2[n].z = SCISDIST2; n++;
            }
        } if (n < 3) continue;

        //Project
        my0 = 0x7fffffff; my1 = 0x80000000;
        for(i=0;i<n;i++)
        {
            f = ghz/vt2[i].z;
            sx[i] = vt2[i].x*f + ghx;
            sy[i] = vt2[i].y*f + ghy;
            isy[i] = std::min(std::max((int)sy[i],grect[ind].y0),grect[ind].y1);
            if (isy[i] < my0) my0 = isy[i];
            if (isy[i] > my1) my1 = isy[i];
            sy[i] -= .5;
        }

        //Rasterize
        for(i=n-1,j=0;j<n;i=j,j++)
        {
            if (isy[i] == isy[j]) continue;
            if (isy[i] < isy[j]) { sy0 = isy[i]; sy1 = isy[j]; iptr = &xc[1][0]; }
                                 else { sy0 = isy[j]; sy1 = isy[i]; iptr = &xc[0][0]; }
            fxinc = (sx[j]-sx[i])/(sy[j]-sy[i]); fx = ((float)sy0-sy[i])*fxinc + sx[i];
            for(y=sy0;y<sy1;y++,fx+=fxinc) iptr[y] = (int)fx;
        }

        //Draw hlines
        sptr = my0*gxd.p + gxd.f;
        gcptr = (unsigned int *)(my0*gcbitplo3 + (INT_PTR)gcbit);
        for(y=my0;y<my1;y++,sptr+=gxd.p,gcptr+=gcbitplo5)
            for(sx1=std::max(xc[0][y],grect[ind].x0),sxe=std::min(xc[1][y],grect[ind].x1);sx1<sxe;)
            {
                sx0 = dntil1(gcptr,sx1,sxe); if (sx0 >= sxe) break;
                sx1 = dntil0(gcptr,sx0,sxe); if (sx1 > sxe) sx1 = sxe;
                setzrange0(gcptr,sx0,sx1);
                if constexpr (false) {
                    for(x=sx0;x<sx1;x++) {
                        //WARNING:memcpy VERY SLOW to GPU memory?
                        memcpy((void *)(sptr+x*PIXBUFBYPP),val,PIXBUFBYPP);
                    }
                } else if constexpr (PIXBUFBYPP == 8) {
                    #if !defined(_WIN32)
                        for(x=sx0;x<sx1;x++) {
                            //WARNING:memcpy VERY SLOW to GPU memory?
                            memcpy((void *)(sptr+x*PIXBUFBYPP),val,PIXBUFBYPP);
                        }
                    #else
                        _asm
                        {
                            movaps xmm0, val
                            mov ecx, sptr
                            mov eax, sx0
                            mov edx, sx1
                            lea ecx, [ecx+edx*8]
                            sub eax, edx
                            jge short endit
                        begit:
                            movlps [ecx+eax*8], xmm0
                            add eax, 1
                            jl short begit
                        endit:
                        }
                    #endif
                }
                else {
                    sx0 *= PIXBUFBYPP;
                    sx1 *= PIXBUFBYPP;
                    for(x=sx0;x<sx1;x+=4) {
                        *(int *)(sptr+x) = *(int *)&val[x&(PIXBUFBYPP-1)];
                    }
                }
            }
    }
}

static void oct_rendclosecubes (oct_t *loct, point3i *cent)
{
    typedef struct { octv_t *ptr; int x, y, z, j, ord; } stk_t;
    stk_t stk[OCT_MAXLS];
    octv_t *ptr;
    int i, j, k, ls, s, x, y, z, nx, ny, nz, ord, vis, surftyp;
    #define BITDIST 6 //Must be <= 6
    static unsigned int bitvis[16*16], *pbitvis;

    surftyp = oct_getvox(loct,cent->x,cent->y,cent->z,0);
    if (surftyp != 1)
    {
        oct_sol2bit(loct,bitvis,cent->x-BITDIST-1,cent->y-BITDIST-1,cent->z-BITDIST-1,16,16,16,1);
    }

    ls = loct->lsid-1; s = (1<<ls); ptr = &((octv_t *)loct->nod.buf)[loct->head]; x = 0; y = 0; z = 0; j = 8-1;
    ord = (z+s > cent->z)*4 + (y+s > cent->y)*2 + (x+s > cent->x);
    while (1)
    {
        k = j^ord;

        i = (1<<k); if (!(ptr->chi&i)) goto tosibly;

        nx = (( k    &1)<<ls)+x;
        ny = (((k>>1)&1)<<ls)+y;
        nz = (((k>>2)&1)<<ls)+z;

             //skip if > DIST manhattans
        if ((nx > cent->x+BITDIST) || (nx+s <= cent->x-BITDIST)) goto tosibly;
        if ((ny > cent->y+BITDIST) || (ny+s <= cent->y-BITDIST)) goto tosibly;
        if ((nz > cent->z+BITDIST) || (nz+s <= cent->z-BITDIST)) goto tosibly;

            //skip if not in frustum pyramid
        if (!isint_cone_sph(&giics[ls],(float)nx,(float)ny,(float)nz)) goto tosibly;

            //skip if not in front //FIXFIXFIXFIX:comment this out and fix do_rect to not draw manhattan!
        //FIXFIXFIXFIX
        if (nx*gxmul.z + ny*gymul.z + nz*gzmul.z + gnadd.z > cornmin[ls]) goto tosibly;

        if (ls <= 0)
        {
            if (surftyp == 1) vis = 63;
            else
            {
                vis = 0; k = (1<< (nx-(cent->x-BITDIST-1)) );
                pbitvis = &bitvis[(ny-(cent->y-BITDIST-1)) +
                                        (nz-(cent->z-BITDIST-1))*16];
                if (!(pbitvis[  0]&(k>>1))) vis |= (1<<0);
                if (!(pbitvis[  0]&(k<<1))) vis |= (1<<1);
                if (!(pbitvis[ -1]&(k   ))) vis |= (1<<2);
                if (!(pbitvis[ +1]&(k   ))) vis |= (1<<3);
                if (!(pbitvis[-16]&(k   ))) vis |= (1<<4);
                if (!(pbitvis[+16]&(k   ))) vis |= (1<<5);
                if (surftyp) vis ^= 63;
            }
            drawcube_close(loct,nx,ny,nz,popcount8(ptr->chi&(i-1)) + ptr->ind,0,vis);
            goto tosibly;
        }

        stk[ls].ptr = ptr; stk[ls].x = x; stk[ls].y = y; stk[ls].z = z; stk[ls].j = j; stk[ls].ord = ord; ls--; s >>= 1; //2child
        ptr = &((octv_t *)loct->nod.buf)[popcount8(ptr->chi&(i-1)) + ptr->ind]; x = nx; y = ny; z = nz; j = 8-1;
        ord = (z+s > cent->z)*4 + (y+s > cent->y)*2 + (x+s > cent->x);
        continue;

tosibly:;
        j--; if (j >= 0) continue;
        do { ls++; s <<= 1; if (ls >= loct->lsid) return; j = stk[ls].j-1; } while (j < 0); //2parent
        ptr = stk[ls].ptr; x = stk[ls].x; y = stk[ls].y; z = stk[ls].z; ord = stk[ls].ord;
    }
}

// TO REMOVE
// This is only used by a some ASM code below.
#define B2(n)    n ,    n+1 ,    n+1 ,    n+2
#define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
constexpr static const int popcount8_array[256] = {
    B6(0), B6(1), B6(1), B6(2)
};
#undef B6
#undef B4
#undef B2
// END REMOVE

static void oct_ftob2_dorect (int ind, void *_)
{
    alignas(16) lpoint4d stk[OCT_MAXLS];
    alignas(16) float psmin[4], cornadd[OCT_MAXLS][8];
    alignas(16) float vxi[8], vyi[8];
    int i, x, y, z, x0, y0, z0, x1, y1, z1, ls, stkind[OCT_MAXLS];
    unsigned int stkord[OCT_MAXLS];
    oct_t *loct = (oct_t *)_;

    // variables for asm code
    alignas(16) static const int dq0fff[4] = {-1,-1,-1, 0};
    alignas(16) static short dqpitch[8];
    alignas(8) static const short cubvmul[4] = {1*16,3*16,9*16,0};
    static int gxdf = 0, gxdp = 0;
    static const float qzero = 0.f;

    if constexpr (oct_useasm) {
        if (ind < 0) {
            static int ognodbuf = 0, ogcbit = 0, ongzbufoff = 0, oglsid = 0, init = 1;
            if (ognodbuf != (INT_PTR)loct->nod.buf) //FIX:incompatible with MT rendering of separate models
            {
                ognodbuf = (INT_PTR)loct->nod.buf;
                i = (INT_PTR)&((octv_t *)loct->nod.buf)->chi;
                SELFMODVAL(init,selfmod_goctchi0-4,i);
                SELFMODVAL(init,selfmod_goctchi1-4,i);
                SELFMODVAL(init,selfmod_goctchi2-4,i);
                i = (INT_PTR)&((octv_t *)loct->nod.buf)->ind;
                SELFMODVAL(init,selfmod_goctind0-4,i);
                SELFMODVAL(init,selfmod_goctind1-4,i);
            }
            if (ogcbit != (INT_PTR)gcbit)
            {
                ogcbit = (INT_PTR)gcbit;
                SELFMODVAL(init,selfmod_gcbit0-4,ogcbit);
                SELFMODVAL(init,selfmod_gcbit1c-4,ogcbit);
                SELFMODVAL(init,selfmod_gcbit2c-4,ogcbit);
                SELFMODVAL(init,selfmod_gcbit3c-4,ogcbit);
            }
            if (oglsid != loct->lsid)
            {
                oglsid = loct->lsid;
                SELFMODVAL(init,selfmod_maxls4-1,(char)(loct->lsid*4));
            }

            dqpitch[0] = 1; dqpitch[1] = gcbitpl; dqpitch[2] = PIXBUFBYPP; dqpitch[3] = gxd.p; dqpitch[4] = 0; dqpitch[5] = 0; dqpitch[6] = 0; dqpitch[7] = 0;
            gxdf = gxd.f;
            gxdp = gxd.p;

            init = 0;
            return;
        }
    }

    psmin[0] = (float)grect[ind].x0; psmin[1] = (float)grect[ind].y0; psmin[2] = (float)-grect[ind].x1; psmin[3] = (float)-grect[ind].y1;

    x = grect[ind].x0+grect[ind].x1 - ighx*2; y = grect[ind].y0+grect[ind].y1 - ighy*2; z = ighz*2;
    x0 = (gxmul.x*z >= gxmul.z*x); y0 = (gymul.x*z >= gymul.z*x); z0 = (gzmul.x*z >= gzmul.z*x);
    cornadd[0][2] = x0*fgxmul.x + y0*fgymul.x + z0*fgzmul.x; cornadd[0][0] = fgxmul.x + fgymul.x + fgzmul.x - cornadd[0][2];
    cornadd[0][6] = x0*fgxmul.z + y0*fgymul.z + z0*fgzmul.z; cornadd[0][4] = fgxmul.z + fgymul.z + fgzmul.z - cornadd[0][6];
    x1 = (gxmul.y*z >= gxmul.z*y); y1 = (gymul.y*z >= gymul.z*y); z1 = (gzmul.y*z >= gzmul.z*y);
    cornadd[0][3] = x1*fgxmul.y + y1*fgymul.y + z1*fgzmul.y; cornadd[0][1] = fgxmul.y + fgymul.y + fgzmul.y - cornadd[0][3];
    cornadd[0][7] = x1*fgxmul.z + y1*fgymul.z + z1*fgzmul.z; cornadd[0][5] = fgxmul.z + fgymul.z + fgzmul.z - cornadd[0][7];
    cornadd[0][6] = -cornadd[0][6]; cornadd[0][7] = -cornadd[0][7]; //Negation is asm trick to eliminate minps
    for(ls=1;ls<loct->lsid;ls++) for(i=0;i<8;i++) cornadd[ls][i] = cornadd[ls-1][i]*2;

    ls = loct->lsid-1;
    stk[ls].x = fgnadd.x; stk[ls].y = fgnadd.y; stk[ls].z = fgnadd.z;

    //--------------------------------------------------------------------------------

    if constexpr (oct_useasm) {
        stk[ls].z2 = -fgnadd.z;
        stkind[loct->lsid-1] = loct->head;

#if defined(_WIN32)
        __asm
        {
            //eax:ls*4                          mm0:temp:short life/color  xmm0:[   -sy1    -sx1     sy0     sx0]
            //ebx:temp:gcbitptr                 mm1:-                      xmm1:[sy0-sy1 sx0-sx1 sy0-sy1 sx0-sx1]/temp
            //ecx:temp:many places;cl as shift  mm2:[- z y x]              xmm2:[      y       x      z2       z] (3D screen space coord)
            //edx:temp:many places              mm3:temp                   xmm3:[      0       0    zptr gcbtptr] (pmaddwd results)
            //esi:temp:many places              mm4:temp:holds gddz/gxd    xmm4:-
            //edi:ord (stkord[ls] cache)        mm5:temp                   xmm5:-
            //ebp:(stack offset)                mm6:-                      xmm6:-
            //esp:dec y/jnz                     mm7:[0 0 0 0]              xmm7:(temp)

            push ebx
            push esi
            push edi
            mov i, esp

            mov eax, ls
            shl eax, 2

            pxor mm2, mm2
            pxor mm7, mm7

            mov esi, stkind[eax]
            jmp short in2it_asm

        align 16
        tochild_asm:
            mov stkord[eax], edi          ;stkord[ls] = ord;
            sub eax, 4                    ;ls--;

            movaps stk[eax*4], xmm2       ;stk[ls].x = pt.x; stk[ls].y = pt.y; stk[ls].z = pt.z;

            and edi, 7

            mov esi, stkind[eax+4]        ;stkind[ls] = popcount8(loct->nod.buf[stkind[ls+1]].chi&pow2m1[k]) + loct->nod.buf[stkind[ls+1]].ptr; //2child
            movzx ecx, byte ptr [esi*8+0x88888888] _asm selfmod_goctchi2: ;octv_t.chi
            and ecx, pow2m1[edi*4]
            mov esi, [esi*8+0x88888888] _asm selfmod_goctind0: ;octv_t.ind
            add esi, popcount8_array[ecx*4]
            mov stkind[eax], esi

            pand mm2, glandlut[eax*2+8]   ;stk2.x &= glandlut[ls+1].x; stk2.y &= glandlut[ls+1].y; stk2.z &= glandlut[ls+1].z;
            shl edi, LOCT_MAXLS+3         ;stk2.x |= gklslut[k][ls+1].x; stk2.y |= gklslut[k][ls+1].y; stk2.z |= gklslut[k][ls+1].z;
            por mm2, gklslut[edi+eax*2+8]

        in2it_asm:
            movq mm0, mm2                 ;ord = ordlut[(stk2.x > gldirlut[ls].x) + (stk2.y > gldirlut[ls].y)*2 + (stk2.z > gldirlut[ls].z)*4][loct->nod.buf[stkind[ls]].chi];
            pcmpgtw mm0, gldirlut[eax*2]
            packsswb mm0, mm7
            pmovmskb edx, mm0
            movzx esi, byte ptr [esi*8+0x88888888] _asm selfmod_goctchi0: ;octv_t.chi
            shl edx, 10
            mov edi, ordlut[edx+esi*4]

        top_asm:
            mov edx, edi                  ;k = (ord&7);
            and edx, 7

            ;Start divide calculation early..
            movaps xmm2, stk[eax*4]
            mov ebx, edx                  ;backup for code at covskip_asm:
            shl edx, LOCT_MAXLS+4         ;ptr = &gposadd[k][ls]; pt.* = stk[ls].* + ptr->*;
            addps xmm2, gposadd[edx+eax*4]

            pshufd xmm0, xmm2, 0xee       ;xmm0: [pt.y  pt.x  pt.y pt.x]
            pshufd xmm1, xmm2, 0x50       ;xmm1: [pt.z2 pt.z2 pt.z pt.z]
            addps xmm0, cornadd[eax*8]    ;xmm0: [pt.y +fptr[3] pt.x +fptr[2] pt.y+fptr[1] pt.x+fptr[0]]
            addps xmm1, cornadd[eax*8+16] ;xmm1: [pt.z2+fptr[7] pt.z2+fptr[6] pt.z+fptr[5] pt.z+fptr[4]]

            #if 1
                    divps xmm0, xmm1
            #elif 0
                    rcpps xmm7, xmm1 ;newton-raphson formula (doesn't seem any faster than divps :/)
                    mulps xmm1, xmm7
                    mulps xmm1, xmm7
                    addps xmm7, xmm7
                    subps xmm7, xmm1
                    mulps xmm0, xmm7
            #else
                    rcpps xmm1, xmm1 ;noticable artifacts (random missing pixels)
                    mulps xmm0, xmm1
            #endif

            ;Lo probability out..
            ucomiss xmm2, cornmin[eax]    ;if (pt.z <= cornmin[ls]) { //WARNING:can't use int cmp: cornmin[] may be +/-
            ja short cornskip
            ucomiss xmm2, cornmax[eax]    ;if (pt.z <= cornmax[ls]) goto tosibly_asm; else goto tochild_asm; //WARNING:can't use int cmp: cornmax[] may be +/-
            jbe short tosibly_asm

            ;FIXFIXFIX:convert cone algo to use xmm2/screen space 3d coord!
            ;if (!isint_cone_sph(&giics[ls],(stk2.x&glandlut[ls].x)|gklslut[k][ls].x,
            ;                               (stk2.y&glandlut[ls].y)|gklslut[k][ls].y,
            ;                               (stk2.z&glandlut[ls].z)|gklslut[k][ls].z)) goto tosibly_asm;
            movq mm0, mm2
            pand mm0, glandlut[eax*2]     ;stk2.x &= glandlut[ls].x; stk2.y &= glandlut[ls].y; stk2.z &= glandlut[ls].z;
            shr edx, 1                    ;edx = (ord&7)<<(LOCT_MAXLS+3)
            por mm0, gklslut[edx+eax*2]   ;stk2.x |= gklslut[k][ls].x; stk2.y |= gklslut[k][ls].y; stk2.z |= gklslut[k][ls].z;
                                                    ;sx = (float)stk2.x; sy = (float)stk2.y; sz = (float)stk2.z;
            movq2dq xmm4, mm0             ;xmm4:[0  0| 0  0| 0 sz| sy sx]
            pshufd xmm4, xmm4, 0xd8       ;xmm4:[0  0| 0 sz| 0  0| sy sx]
            pshuflw xmm4, xmm4, 0xd8      ;xmm4:[0  0| 0 sz| 0 sy|  0 sx]
            cvtdq2ps xmm4, xmm4
            subps xmm4, giics[eax*8]      ;sx -= iics->x0; sy -= iics->y0; sz -= iics->z0;
            movaps xmm5, xmm4
            mulps xmm4, giics[eax*8+16]   ;sx *= iics->xv; sy *= iics->yv; sz *= iics->zv;
            mulps xmm5, xmm5              ;sx *= sx;       sy *= sy;       sz *= sz;
            movhlps xmm6, xmm4
            movhlps xmm7, xmm5
            addps xmm4, xmm6
            addps xmm5, xmm7
            pshuflw xmm6, xmm4, 0xe
            pshuflw xmm7, xmm5, 0xe
            addss xmm4, xmm6              ;xmm4[0] = sx*iics->xv + sy*iics->yv + sz*iics->zv
            addss xmm5, xmm7              ;xmm5[0] = sx*sx + sy*sy + sz*sz
            ucomiss xmm4, qzero           ;if (xmm4[0] <= 0.f) goto sibly;
            jbe short tosibly_asm
            mulss xmm4, xmm4              ;if (xmm4[0]*xmm4[0] <= xmm5[0]*iics->cosang2) goto tosibly_asm;
            mulss xmm5, giics[eax*8+28]
            ucomiss xmm4, xmm5
            jbe short tosibly_asm

            jmp short tochild_asm

        align 16
        cornskip:
            ucomiss xmm2, scanmax[eax]    ;if (pt.z >= scanmax[ls]) goto tosibly_asm;
            jae short tosibly_asm

            maxps xmm0, psmin             ;xmm0: [std::max(-sy1,-ymax) std::max(-sx1,-xmax) std::max(sy0,ymin) std::max(sx0,xmin)]
            cvtps2dq xmm0, xmm0           ;xmm0: [-sy1 -sx1  sy0  sx0]
            pshufd xmm1, xmm0, 0x4e       ;xmm1: [ sy0  sx0 -sy1 -sx1]
            paddd xmm1, xmm0              ;xmm1: [sy0-sy1 sx0-sx1 sy0-sy1 sx0-sx1]
            movmskps edx, xmm1
            cmp edx, 15
            jne short tosibly_asm
            pxor xmm1, dq0fff             ;xmm1: [sy0-sy1 sx1-sx0-1 sy1-sy0-1 sx1-sx0-1]

            ;zptr = sy0*gdd.p + (sx0<<2); i = sy0*gcbitpl + sx0;
            pshuflw xmm3, xmm0, 0x88      ;   xmm3:[? ? ? ?   sy0 sx0     sy0 sx0]
            pmaddwd xmm3, dqpitch         ;dqpitch:[0 0 0 0 gdd.p   4 gcbitpl   1]
            movd esp, xmm3

            ;NOTE:always checking gcbit is faster than only doing it for ls>0
            movd ecx, xmm0 ;sx0        ;v = (sx0&7)+sx1;
            and ecx, 7
            movd esi, xmm1 ;sx1
            add esi, ecx
            cmp esi, 31                ;if (v >= 32) goto tochild_asm; //always recurse for large regions (don't bother to check cover map)
            jae short covskip_asm
            shr esp, 3                 ;uptr = (unsigned int *)(((int)gcbit) + (i>>3));
            pextrw edx, xmm1, 2 ;sy1   ;v = npow2[sx0&7]&pow2m1[v];
            mov esi, pow2m1[esi*4+4]
            and esi, npow2[ecx*4]      ;for(;sy1>0;sy1--,uptr+=gcbitplo5) if (uptr[0]&v) goto tochild_asm;
        covbegy:
            test esi, [esp+0x88888888] _asm selfmod_gcbit0:
            jne short covskip_asm
            add esp, gcbitplo3
            sub edx, 1
            jge short covbegy
            jmp short tosibly_asm          ;goto tosibly_asm;

        align 16                         ;/*}*/
        covskip_asm:
            test eax, eax
            jnz short tochild_asm

            mov edx, stkind[eax]          ;mm3: ipsurf = popcount8(loct->nod.buf[stkind[ls]].chi&pow2m1[k]) + loct->nod.buf[stkind[ls]].ptr;
            movzx esi, byte ptr [edx*8+0x88888888] _asm selfmod_goctchi1: ;octv_t.chi
            and esi, pow2m1[ebx*4]
            mov edx, [edx*8+0x88888888] _asm selfmod_goctind1: ;octv_t.ind
            add edx, popcount8_array[esi*4]
            movd mm3, edx

        ;---------------------------------------------------------------------------------------------------

            movq mm0, gkls0[ebx*8]        ;mm0:[- z=stk2.z+gkls0[k].z y=stk2.y+gkls0[k].y x=stk2.x+gkls0[k].x]
            paddw mm0, mm2

            ;if ((labs(stk2.x+gkls0[k].x-glorig.x) <= 6) &&
            ;    (labs(stk2.y+gkls0[k].y-glorig.y) <= 6) &&
            ;    (labs(stk2.z+gkls0[k].z-glorig.z) <= 6)) goto tosibly_asm;
            ;movq mm5, mm0
            ;psubw mm5, ?                 ;[0 -32768+ 6-glorig.z -32768+ 6-glorig.y -32768+ 6-glorig.x]
            ;pcmpgtw mm5, ?               ;[0 -32768+12-glorig.z -32768+12-glorig.y -32768+12-glorig.x]
            ;pmovmskb ecx, mm5
            ;test ecx, 0x3f
            ;jnz short tochild_asm

            movq mm5, mm0
            paddsw mm5, glcenadd
            psubusw mm5, glcensub         ;mm5:[0 sgn(z)+1][sgn(y)+1 sgn(x)+1]
            pmaddwd mm5, cubvmul          ;   *[0     9*16][    3*16     1*16]
            movd ecx, mm5
            pextrw esp, mm5, 2
            add ecx, esp

            movd ebx, xmm3                ;ebx: gcbitptr

            pshuflw xmm3, xmm3, 0xe       ;eax: sptr = sy0*gxd.p + (sx0<<1) + gxd.f;
            movd eax, xmm3
            add eax, gxdf

        ;===================================================================================================

            //Cube rasterization. Register status at this point:
            //eax:gddf                    mm0:xyzsurf data  xmm0:[   -sy1    -sx1     sy0     sx0]
            //ebx:temp:gcbitptr           mm1:-             xmm1:[      -       -   sy1-1   sx1-1]
            //ecx:temp:ind27*16           mm2:[- z y x]     xmm2:[      y       x      z2       z] (3D screen space coord)
            //edx:-                       mm3:ipsurf        xmm3:-
            //esi:-                       mm4:-             xmm4:-
            //edi:ord (stkord[ls] cache)  mm5:-             xmm5:-
            //ebp:(stack offset)          mm6:-             xmm6:-
            //esp:-                       mm7:[0 0 0 0]     xmm7:-

            pshufd xmm3, xmm2, 0x00       ;xmm3: [pt.z pt.z pt.z pt.z]
            movaps xmm5, xmm3             ;xmm5: [pt.z pt.z pt.z pt.z]
            mulps xmm3, gvany[0]
            mulps xmm5, gvanx[0]
            pshufd xmm4, xmm2, 0xff       ;xmm4: [pt.y pt.y pt.y pt.y]
            subps xmm3, xmm4
            pshufd xmm4, xmm2, 0xaa       ;xmm4: [pt.x pt.x pt.x pt.x]
            subps xmm4, xmm5
            pand xmm3, vause[ecx]         ;xmm3: [   0   k2   k1   k0]
            pand xmm4, vause[ecx]         ;xmm4: [   0   k6   k5   k4]

            pshufd xmm0, xmm0, 0x6a       ;xmm0:[ y0   -x1   -x1   -x1]
            pxor xmm0, dq0fff             ;xmm0:[ y0  x1-1  x1-1  x1-1]
            cvtdq2ps xmm5, xmm0           ;xmm5:[fy0 fx1-1 fx1-1 fx1-1]
            pshufd xmm6, xmm5, 0xff       ;xmm6:[ -    fy0   fy0   fy0]
            subps xmm5, gvanxmh[0]        ;xmm5:[ -     ka    k9    k8]
            subps xmm6, gvanymh[0]        ;xmm6:[ -     ke    kd    kc]

            movaps xmm0, cornvadd2[ecx*4+32]
            movaps xmm7, cornvadd2[ecx*4+48]
            subps xmm0, xmm3
            subps xmm7, xmm4
            addps xmm3, cornvadd2[ecx*4]
            addps xmm4, cornvadd2[ecx*4+16]
            movaps vxi[ 0], xmm3
            movaps vxi[16], xmm0
            movaps vyi[ 0], xmm4
            movaps vyi[16], xmm7

            mulps xmm0, xmm5
            mulps xmm5, xmm3
            mulps xmm4, xmm6
            mulps xmm6, xmm7
            addps xmm5, xmm4              ;xmm5 = xmm5*vxi[ 0] + xmm6*vyi[ 0] ;<- oval[ 0]
            addps xmm6, xmm0              ;xmm6 = xmm5*vxi[16] + xmm6*vyi[16] ;<- oval[16]
            pand xmm5, vause[ecx]
            pand xmm6, vause[ecx]

            ;addps xmm5, dq?               ;oval[ 0]
            ;addps xmm6, dq?               ;oval[16]

        ;===================================================================================================

            ; val: oval: | [3] [2] [1] [0]
            ;------------+----------------
            ; xmm0  xmm5 |  0  v2  v1  v0
            ; xmm3  xmm6 |  0  v6  v5  v4

        ;---------------------------------------------------------------------------------------------------

            ;eax:gxdptr, ebx:gcbitptr, ecx:temp, edx:temp, esi:xcnt, edi:(ord), ebp:(stack offset), esp:ycnt
            pextrw esp, xmm1, 2 ;sy1-1    ;sy = sy1-1;
        begpixyc:
            movd esi, xmm1 ;sx1-1      ;sx = sx1-1;
            movaps xmm0, xmm5
            movaps xmm3, xmm6
            addps xmm5, vyi[ 0]
            addps xmm6, vyi[16]
        begpixxc:

        #if 0
            movmskps ecx, xmm0   ;3,1
            movmskps edx, xmm3   ;3,1
            subps xmm0, vxi[ 0]
            subps xmm3, vxi[16]
            xor ecx, edx
            jnz short skippixc
        #else
            movaps xmm7, xmm0    ;1,0.33
            pxor xmm7, xmm3      ;1,0.33
            movmskps edx, xmm7   ;3,1
            subps xmm0, vxi[ 0]
            subps xmm3, vxi[16]
            test edx, edx
            jnz short skippixc
        #endif

        lea ecx, [ebx+esi]      ;v = i+sx;
        mov edx, 1              ;if (gcbit[v>>5]&(1<<v)) {
        shl edx, cl
        shr ecx, 5
        test edx, [ecx*4+0x88888888] _asm selfmod_gcbit1c:
        jz short skippixc
        xor edx, [ecx*4+0x88888888] _asm selfmod_gcbit2c:            ;gcbit[v>>5] ^= (1<<v);
        mov [ecx*4+0x88888888], edx _asm selfmod_gcbit3c:
        movd [eax+esi*4], mm3

        skippixc:
            sub esi, 1              ;sx--;
            jge short begpixxc      ;} while (sx >= 0);
            add eax, gxdp              ;sptr += gxdp;
            add ebx, gcbitpl           ;i += gcbitpl;
            sub esp, 1                 ;sy--;
            jge short begpixyc         ;while (sy >= 0);

            xor eax, eax                  ;clear eax back to 0 for use as ls*4

        ;---------------------------------------------------------------------------------------------------

        tosibly_asm:
            shr edi, 4                    ;ord >>= 4; if (ord) goto top_asm;
            jnz short top_asm
            add eax, 4                    ;ls++; //2parent
            mov edi, stkord[eax]          ;ord = stkord[ls];
            cmp eax, 0x78 _asm selfmod_maxls4:
            jb short tosibly_asm              ;} while (ls < loct->lsid);

            mov esp, i
            pop edi
            pop esi
            pop ebx

            emms
        }
#else
        fprintf(stderr, "Unimplemented\n");
        std::abort();
#endif
    }
    //--------------------------------------------------------------------------------------------
    else {
        alignas(16) lpoint4d pt;
        alignas(16) rect_t smax;
        alignas(8)  stk2_t stk2;
        alignas(16) float oval[8] = {0}, val[8] = {0};
        float fx, fy, k0, k1, k2, k4, k5, k6, k8, k9, ka, kc, kd, ke;
        lpoint4d *pptr;
        float *fptr;
        unsigned int *uptr, ord;
        int k, v, sx, ind27;
        //alignas(8) stk2_t col2;
        int sptr, ipsurf;

        stk2.x = 0;
        stk2.y = 0;
        stk2.z = 0;
        stkind[ls] = loct->head;
        goto in2it_cpp;

    tochild_cpp:;
        stkord[ls] = ord; ls--;
        stk[ls].x = pt.x; stk[ls].y = pt.y; stk[ls].z = pt.z;
        stkind[ls] = popcount8(((octv_t *)loct->nod.buf)[stkind[ls+1]].chi&pow2m1[k]) + ((octv_t *)loct->nod.buf)[stkind[ls+1]].ind; //2child
        stk2.x = (stk2.x&glandlut[ls+1].x)|gklslut[k][ls+1].x;
        stk2.y = (stk2.y&glandlut[ls+1].y)|gklslut[k][ls+1].y;
        stk2.z = (stk2.z&glandlut[ls+1].z)|gklslut[k][ls+1].z;
    in2it_cpp:;
        ord = ordlut[(stk2.x > gldirlut[ls].x) + (stk2.y > gldirlut[ls].y)*2 + (stk2.z > gldirlut[ls].z)*4][((octv_t *)loct->nod.buf)[stkind[ls]].chi];
    top_cpp:;
        k = (ord&7);

        pptr = &gposadd[k][ls];
        pt.x = stk[ls].x + pptr->x;
        pt.y = stk[ls].y + pptr->y;
        pt.z = stk[ls].z + pptr->z;

        if (pt.z <= cornmin[ls]) //WARNING:can't use int cmp: cornmin[] may be +/-
        {
            if (pt.z <= cornmax[ls]) goto tosibly_cpp; //WARNING:can't use int cmp: cornmax[] may be +/-
            if (!isint_cone_sph(&giics[ls],(stk2.x&glandlut[ls].x)|gklslut[k][ls].x,(stk2.y&glandlut[ls].y)|gklslut[k][ls].y,(stk2.z&glandlut[ls].z)|gklslut[k][ls].z)) goto tosibly_cpp;
            goto tochild_cpp;
        }
        if (pt.z >= scanmax[ls]) goto tosibly_cpp;

        fptr = &cornadd[ls][0];
        smax.x0 = std::max((int)((pt.x+fptr[0]) / (pt.z+fptr[4])),grect[ind].x0);
        smax.y0 = std::max((int)((pt.y+fptr[1]) / (pt.z+fptr[5])),grect[ind].y0);
        smax.x1 = std::min((int)((pt.x+fptr[2]) / (pt.z-fptr[6])),grect[ind].x1);
        smax.y1 = std::min((int)((pt.y+fptr[3]) / (pt.z-fptr[7])),grect[ind].y1);
        smax.x1 -= smax.x0; smax.y1 -= smax.y0; if ((smax.x1 <= 0) || (smax.y1 <= 0)) goto tosibly_cpp;
        i = smax.y0*gcbitpl+smax.x0;

        v = (smax.x0&7)+smax.x1; if (v >= 32) goto covskip_cpp; //always recurse for large regions (don't bother to check cover map)
        uptr = (unsigned int *)(((INT_PTR)gcbit) + (i>>3)); v = npow2[smax.x0&7]&pow2m1[v];
        for(sx=smax.y1;sx>0;sx--,uptr+=gcbitplo5) if (uptr[0]&v) goto covskip_cpp;
        goto tosibly_cpp;
    covskip_cpp:;
        if (ls) goto tochild_cpp;

        sptr = smax.y0*gxd.p + smax.x0*PIXBUFBYPP + gxd.f;
        ipsurf = popcount8(((octv_t *)loct->nod.buf)[stkind[0]].chi&pow2m1[k]) + ((octv_t *)loct->nod.buf)[stkind[0]].ind;
        {
    #if 1    //FIXFIXFIXFIX:should be 0!
            drawcube_close(loct,stk2.x+gkls0[k].x,stk2.y+gkls0[k].y,stk2.z+gkls0[k].z,ipsurf,ind,63); goto tosibly_cpp; //fun with brute force
    #endif

            fx = (float)(smax.x0+smax.x1-1);
            fy = (float)(smax.y0          );

            //if ((labs(stk2.x+gkls0[k].x-glorig.x) <= 6) &&
            //    (labs(stk2.y+gkls0[k].y-glorig.y) <= 6) &&
            //    (labs(stk2.z+gkls0[k].z-glorig.z) <= 6)) goto tosibly_cpp;

            ind27 = std::min(std::max(stk2.x+gkls0[k].x-glorig.x+1,0),2)*1 + std::min(std::max(stk2.y+gkls0[k].y-glorig.y+1,0),2)*3 + std::min(std::max(stk2.z+gkls0[k].z-glorig.z+1,0),2)*9;

            k0 = pt.z*gvany[0] - pt.y; *(int *)&k0 &= vause[ind27][0]; k4 = pt.x - pt.z*gvanx[0]; *(int *)&k4 &= vause[ind27][0];
            k1 = pt.z*gvany[1] - pt.y; *(int *)&k1 &= vause[ind27][1]; k5 = pt.x - pt.z*gvanx[1]; *(int *)&k5 &= vause[ind27][1];
            k2 = pt.z*gvany[2] - pt.y; *(int *)&k2 &= vause[ind27][2]; k6 = pt.x - pt.z*gvanx[2]; *(int *)&k6 &= vause[ind27][2];
            k8 = fx-gvanxmh[0]; kc = fy-gvanymh[0];
            k9 = fx-gvanxmh[1]; kd = fy-gvanymh[1];
            ka = fx-gvanxmh[2]; ke = fy-gvanymh[2];
            vxi[0] = cornvadd2[ind27][ 0] + k0; vyi[0] = cornvadd2[ind27][ 4] + k4; oval[0] = k8*vxi[0] + kc*vyi[0] + 1e-1; //FIX:+1e-1's not yet implemented in ASM!
            vxi[1] = cornvadd2[ind27][ 1] + k1; vyi[1] = cornvadd2[ind27][ 5] + k5; oval[1] = k9*vxi[1] + kd*vyi[1] + 1e-1;
            vxi[2] = cornvadd2[ind27][ 2] + k2; vyi[2] = cornvadd2[ind27][ 6] + k6; oval[2] = ka*vxi[2] + ke*vyi[2] + 1e-1;
            vxi[4] = cornvadd2[ind27][ 8] - k0; vyi[4] = cornvadd2[ind27][12] - k4; oval[4] = k8*vxi[4] + kc*vyi[4] + 1e-1;
            vxi[5] = cornvadd2[ind27][ 9] - k1; vyi[5] = cornvadd2[ind27][13] - k5; oval[5] = k9*vxi[5] + kd*vyi[5] + 1e-1;
            vxi[6] = cornvadd2[ind27][10] - k2; vyi[6] = cornvadd2[ind27][14] - k6; oval[6] = ka*vxi[6] + ke*vyi[6] + 1e-1;
            do
            {
                val[0]  = oval[0]; val[1]  = oval[1]; val[2]  = oval[2]; val[4]  = oval[4]; val[5]  = oval[5]; val[6]  = oval[6];
                oval[0] += vyi[0]; oval[1] += vyi[1]; oval[2] += vyi[2]; oval[4] += vyi[4]; oval[5] += vyi[5]; oval[6] += vyi[6];
                sx = smax.x1;
                do
                {
                    if (((*(int *)&val[0]^*(int *)&val[4]) >= 0) && ((*(int *)&val[1]^*(int *)&val[5]) >= 0) && ((*(int *)&val[2]^*(int *)&val[6]) >= 0))
                    {
                        v = i+sx-1;
                        if (gcbit[v>>5]&(1<<v))
                        {
                            gcbit[v>>5] ^= (1<<v);
                            *(int *)((sx-1)*4 + sptr) = ipsurf;
                        }
                    }
                    val[0] -= vxi[0]; val[1] -= vxi[1]; val[2] -= vxi[2]; val[4] -= vxi[4]; val[5] -= vxi[5]; val[6] -= vxi[6];
                    sx--;
                } while (sx > 0);
                sptr += gxd.p;
                i += gcbitpl; smax.y1--;
            } while (smax.y1 > 0);
        }

    tosibly_cpp:;
        do
        {
            ord >>= 4; if (ord) goto top_cpp;
            ls++; if (ls >= loct->lsid) return; //2parent
            ord = stkord[ls];
        } while (1);
    }
    //-------------------------------------------------------------------------------------------
}

static void oct_setcam (int xres, int yres, tiletype *dd, int zbufoff, point3f *p, point3f *r, point3f *d, point3f *f, float hx, float hy, float hz)
{
    static int inited = 0;

    if (!inited)
    {
        int i, ls;
        inited = 1;
        for(ls=0;ls<OCT_MAXLS;ls++)
        {
            glandlut[ls].x = ((-2)<<(ls-0));
            glandlut[ls].y = ((-2)<<(ls-0));
            glandlut[ls].z = ((-2)<<(ls-0));
            glandlut[ls].dum = 0;
        }
        for(ls=0;ls<OCT_MAXLS;ls++)
            for(i=0;i<8;i++)
            {
                gklslut[i][ls].x = (((i   )&1)<<ls);
                gklslut[i][ls].y = (((i>>1)&1)<<ls);
                gklslut[i][ls].z = (((i>>2)&1)<<ls);
                gklslut[i][ls].dum = 0;
            }
        for(i=0;i<8;i++) gkls0[i] = gklslut[i][0];

        int ord, msk, v, n;
        for(ord=0;ord<8;ord++)
            for(msk=0;msk<256;msk++)
            {
                v = 0; n = 0;
                for(i=0;i<8;i++) if (msk&(1<<(ord^i^7))) { v += ((ord^i^15)<<n); n += 4; }
                ordlut[ord][msk] = v;
            }
        if constexpr (!oct_usegpu) {
            gxd.f = 0; gxd.p = 0; gxd.x = 0; gxd.y = 0;
        }
    }
    if constexpr (!oct_usegpu) {
        if (dd->x*dd->y > gxd.x*gxd.y)
        {
            if (gxd.f) free((void *)gxd.f);
            gxd.f = (INT_PTR)malloc(dd->x*dd->y*PIXBUFBYPP);
            if (!gxd.f) {
                fprintf(stderr, "oct_setcam():malloc failed");
                std::abort();
            }
        }
        gxd.p = dd->x*PIXBUFBYPP; gxd.x = dd->x; gxd.y = dd->y;
    }

    if constexpr (!oct_usegpu) {
        gdd = (*dd);
        gzbufoff = zbufoff;
        gddz = gdd.f+zbufoff;
    }
    else if constexpr (!oct_usegpubo)
    {
        if (gpixbufmal < gpixxdim*gpixydim*PIXBUFBYPP)
        {
            gpixbufmal = gpixxdim*gpixydim*PIXBUFBYPP;
            gpixbuf = (char *)realloc(gpixbuf,gpixbufmal);
        }
        gxd.f = (INT_PTR)gpixbuf; gxd.p = xres*PIXBUFBYPP; gxd.x = xres; gxd.y = yres;
    }

    gipos = (*p); girig = (*r); gidow = (*d); gifor = (*f);
    ghx = hx; ghy = hy; ghz = hz;
    ighx = (int)hx; ighy = (int)hy; ighz = (int)hz;
}

static void oct_setcam (int xres, int yres, tiletype *dd, int zbufoff, point3d *p, point3d *r, point3d *d, point3d *f, double hx, double hy, double hz)
{
    point3f fp, fr, fd, ff;
    fp.x = p->x; fp.y = p->y; fp.z = p->z;
    fr.x = r->x; fr.y = r->y; fr.z = r->z;
    fd.x = d->x; fd.y = d->y; fd.z = d->z;
    ff.x = f->x; ff.y = f->y; ff.z = f->z;
    oct_setcam(xres, yres, dd,zbufoff,&fp,&fr,&fd,&ff,(float)hx,(float)hy,(float)hz);
}

//--------------------------------------------------------------------------------------------------
static int glmost[MAXYDIM];
static int grmost[MAXYDIM];
static int gxmin;
static int gymin;
static int gxmax;
static int gymax;
static void rastquad (int xres, int yres, point3f vt[8], const int ind[4], int *lmost, int *rmost, int *xmin, int *ymin, int *xmax, int *ymax)
{
    point3f vt2[8];
    float f, fi, fsx[8], fsy[8];
    int i, j, sy, imost[2], scry[8], n2;
    #define SCISDIST 5e-4f //NOTE:1e-6f fails on some /ils>0

    (*xmin) = 0x7fffffff; (*ymin) = 0x7fffffff;
    (*xmax) = 0x80000000; (*ymax) = 0x80000000;
    for(i=4-1,j=0,n2=0;j<4;i=j,j++)
    {
        if (vt[ind[i]].z >= SCISDIST) { vt2[n2] = vt[ind[i]]; n2++; }
        if ((vt[ind[i]].z >= SCISDIST) != (vt[ind[j]].z >= SCISDIST))
        {
            f = (SCISDIST-vt[ind[j]].z)/(vt[ind[i]].z-vt[ind[j]].z);
            vt2[n2].x = (vt[ind[i]].x-vt[ind[j]].x)*f + vt[ind[j]].x;
            vt2[n2].y = (vt[ind[i]].y-vt[ind[j]].y)*f + vt[ind[j]].y;
            vt2[n2].z = SCISDIST; n2++;
        }
    } if (n2 < 3) return;

    for(i=n2-1;i>=0;i--)
    {
        f = ghz/vt2[i].z;
        fsx[i] = vt2[i].x*f + ghx;
        fsy[i] = vt2[i].y*f + ghy;

        j = std::min(std::max(-cvttss2si(-fsx[i]),0),xres);
        if (j < (*xmin)) { (*xmin) = j; }
        if (j > (*xmax)) { (*xmax) = j; }

        scry[i] = std::min(std::max(-cvttss2si(-fsy[i]),0),yres);
        if (scry[i] < (*ymin)) { (*ymin) = scry[i]; imost[0] = i; }
        if (scry[i] > (*ymax)) { (*ymax) = scry[i]; imost[1] = i; }
    }
    i = imost[1]; sy = (*ymax);
    do
    {
        j = i+1; if (j >= n2) j = 0;
        if (sy > scry[j])
        {
            fi = (fsx[j]-fsx[i])/(fsy[j]-fsy[i]); f = (sy-fsy[j])*fi + fsx[j] + 0.5; //make CPU render more than necessary; don't trust GPU rounding
            do { f -= fi; sy--; lmost[sy] = cvttss2si(f); } while (sy > scry[j]);
        }
        i = j;
    } while (i != imost[0]);
    do
    {
        j = i+1; if (j >= n2) j = 0;
        if (sy < scry[j])
        {
            fi = (fsx[j]-fsx[i])/(fsy[j]-fsy[i]); f = (sy-fsy[i])*fi + fsx[i] + 1.5; //make CPU render more than necessary; don't trust GPU rounding
            do { rmost[sy] = cvttss2si(f); sy++; f += fi; } while (sy < scry[j]);
        }
        i = j;
    } while (i != imost[1]);
}

//set gcbit if rasterized quad depth (front faces of current sprite) is behind value in zbuf
static void setgcbitifz (int xres, int yres, oct_t *loct, int face, int *lmost, int *rmost, int ymin, int ymax)
{
    alignas(16) constexpr static const float dq0123[4] = {0.f,1.f,2.f,3.f};
    alignas(16) constexpr static const float dq1111[4] = {1.f,1.f,1.f,1.f};
    alignas(16) constexpr static const float dq4444[4] = {4.f,4.f,4.f,4.f};
    float f, fx, fy, fz, topt;
    int sy, sx0, sx1;

    //intx = vx*t = gxmul.x*a + gymul.x*b + gzmul.x*c + gnadd.x
    //inty = vy*t = gxmul.y*a + gymul.y*b + gzmul.y*c + gnadd.y
    //intz = vz*t = gxmul.z*a + gymul.z*b + gzmul.z*c + gnadd.z
    //(-vx)*t + gymul.x*b + gzmul.x*c = gxmul.x-gnadd.x
    //(-vy)*t + gymul.y*b + gzmul.y*c = gxmul.y-gnadd.y
    //(-vz)*t + gymul.z*b + gzmul.z*c = gxmul.z-gnadd.z
    topt = 0.0; if (face&1) f = 0.f; else f = (float)loct->sid;
    switch(face>>1)
    {
        case 0: fx = gzmul.y*gymul.z - gymul.y*gzmul.z; topt  = (gxmul.x*f + gnadd.x)*fx;
                  fy = gzmul.z*gymul.x - gymul.z*gzmul.x; topt += (gxmul.y*f + gnadd.y)*fy;
                  fz = gzmul.x*gymul.y - gymul.x*gzmul.y; topt += (gxmul.z*f + gnadd.z)*fz; break;
        case 1: fx = gxmul.y*gzmul.z - gzmul.y*gxmul.z; topt  = (gymul.x*f + gnadd.x)*fx;
                  fy = gxmul.z*gzmul.x - gzmul.z*gxmul.x; topt += (gymul.y*f + gnadd.y)*fy;
                  fz = gxmul.x*gzmul.y - gzmul.x*gxmul.y; topt += (gymul.z*f + gnadd.z)*fz; break;
        case 2: fx = gymul.y*gxmul.z - gxmul.y*gymul.z; topt  = (gzmul.x*f + gnadd.x)*fx;
                  fy = gymul.z*gxmul.x - gxmul.z*gymul.x; topt += (gzmul.y*f + gnadd.y)*fy;
                  fz = gymul.x*gxmul.y - gxmul.x*gymul.y; topt += (gzmul.z*f + gnadd.z)*fz; break;
    }
    f = (float)loct->sid/(topt*512.0); fx *= f; fy *= f; fz *= f;
    //z_buffer_depth = 1.0 / ((sx-ghx)*fx + (sy-ghy)*fy + ghz*fz)  (true at this line of code)

    for(sy=ymin;sy<ymax;sy++)
    {
        sx0 = std::max(lmost[sy],   0);
        sx1 = std::min(rmost[sy],xres); if (sx0 >= sx1) continue;
#if 1    //safe
        setzrange1(&gcbit[sy*gcbitplo5],sx0,sx1); //just let it overdraw - no savings on rendering 2nd+ sprite
#elif 0
        //WARNING: this block may cause a rare crash in the CPU shader, trying to load an invalid surf due to gcbit being already 0! Z-fighting?
        //Pure C version of below
        f = (sx0-ghx)*fx + (sy-ghy)*fy + ghz*fz;
        fptr = (float *)(gdd.p*sy + gddz);
        gcptr = &gcbit[sy*gcbitplo5];
        for(sx=sx0;sx<sx1;sx++,f+=fx)
        {
            g = fptr[sx]*f;
            if (*(int *)&g > 0x3f800000) //if (g > 1.f)
            {
                gcptr[sx>>5] |= (1<<sx);
                //*(int *)(gdd.f + gdd.p*sy + (sx<<2)) = (((int)f)&255)*0x10101; //Debug only!
            }
        }
#else
        //WARNING: this block may cause a rare crash in the CPU shader, trying to load an invalid surf due to gcbit being already 0! Z-fighting?
        //SSE2 version of above
        f = (sx0-ghx)*fx + (sy-ghy)*fy + ghz*fz;
        fptr = (float *)(gdd.p*sy + gddz);
        gcptr = &gcbit[sy*gcbitplo5];
        while (sx0&7) {        g = fptr[sx0]*(             f); if (*(int *)&g > 0x3f800000) gcptr[sx0>>5] |= (1<<sx0); sx0++; f += fx; if (sx0 >= sx1) goto cont; }
        while (sx1&7) { sx1--; g = fptr[sx1]*((sx1-sx0)*fx+f); if (*(int *)&g > 0x3f800000) gcptr[sx1>>5] |= (1<<sx1);                 if (sx0 >= sx1) goto cont; }
        fptr += sx0; gcptr = (unsigned int *)((sx0>>3) + (INT_PTR)gcptr); sx1 = ((sx1-sx0)>>3);
        __asm
        {
            push esi
            push edi

            movss xmm0, f         ;xmm0: [     0      0      0      f]
            movss xmm1, fx        ;xmm1: [     0      0      0     fx]
            pshufd xmm0, xmm0, 0  ;xmm0: [     f      f      f      f]
            pshufd xmm1, xmm1, 0  ;xmm1: [    fx     fx     fx     fx]
            movaps xmm2, xmm1     ;xmm2: [    fx     fx     fx     fx]
            mulps xmm1, dq4444    ;xmm1: [  fx*4   fx*4   fx*4   fx*4]
            mulps xmm2, dq0123    ;xmm2: [  fx*3   fx*2   fx*1   fx*0]
            addps xmm0, xmm2      ;xmm0: [f+fx*3 f+fx*2 f+fx*1 f+fx*0]
            movaps xmm7, dq1111   ;xmm7: [     1      1      1      1]
            mov esi, fptr
            mov edi, gcptr
            mov ecx, sx1
            add edi, ecx
            neg ecx
         beg:
            movaps xmm2, [esi]
            movaps xmm3, [esi+16]
            mulps xmm2, xmm0
            addps xmm0, xmm1
            mulps xmm3, xmm0
            addps xmm0, xmm1
            pcmpgtd xmm2, xmm7
            pcmpgtd xmm3, xmm7
            movmskps eax, xmm2
            movmskps edx, xmm3
            shl edx, 4
            add eax, edx
            mov byte ptr [edi+ecx], al
            add esi, 32
            add ecx, 1
            jnz short beg

            pop edi
            pop esi
        }
    cont:;
#endif
    }
}

static void setgcbit (int xres, int yres, int *lmost, int *rmost, int ymin, int ymax)
{
    int sy, sx0, sx1;
    for(sy=ymin;sy<ymax;sy++)
    {
        sx0 = std::max(lmost[sy],   0);
        sx1 = std::min(rmost[sy],xres); if (sx0 >= sx1) continue;
        setzrange1(&gcbit[sy*gcbitplo5],sx0,sx1);
    }
}

static void clearbackpixes (int sy, void* loct)
{
    int sx, sxs, sx0, sx1;
    unsigned int *gcptr;
    INT_PTR rptr;

    if ((sy < gymin) || (sy >= gymax)) return;
    sx0 = std::max(glmost[sy],0);
    sx1 = std::min(grmost[sy],((oct_t*)loct)->xres);
    if (sx0 >= sx1) return;

    rptr = sy*gxd.p + gxd.f;
    gcptr = &gcbit[sy*gcbitplo5];

    sxs = sx1;
    while (1)
    {
        sx  = uptil1(gcptr,sxs-1); if (sx <= sx0) return;
        sxs = uptil0(gcptr,sx -1); sxs = std::max(sxs,sx0);
        memset4((void *)(rptr+sxs*PIXBUFBYPP),0x00000000,(sx-sxs)*PIXBUFBYPP);
    }
}

static point3f grxmul, grymul, grzmul;
static float gzscale, gmipoffs;
#define LTEXPREC 7
#define TEXPREC (1<<LTEXPREC)

// Unused variables but required for compiling
alignas(16) static float dqgrxmul[4], dqgorig[4], dqfwmul[4];
alignas(16) static int dqiwmsk[4];
alignas(16) static short dqfogcol[8], dqimulcol[8], dqtilepitch4[16][8], dqmixval[8];
static int tilesf[16];

//pixmeth 0..
void cpu_shader_solcol (int sy, void *_)
{
    int sx, sxs, sx0, sx1, *cptr, *rptr;
    surf_t *psurf;
    unsigned int *gcptr;
    oct_t *loct;

    if ((sy < gymin) || (sy >= gymax)) return;
    sx0 = std::max(glmost[sy],0); sx1 = std::min(grmost[sy],gdd.x); if (sx0 >= sx1) return;

    rptr = (int *)(sy*gxd.p + gxd.f);
    cptr = (int *)(sy*gdd.p + gdd.f);
    gcptr = &gcbit[sy*gcbitplo5];
    loct = (oct_t *)_;

    sxs = sx1;
    while (1)
    {
        sx  = uptil0(gcptr,sxs-1); if (sx <= sx0) return;
        sxs = uptil1(gcptr,sx -1); sxs = std::max(sxs,sx0);
        for(sx--;sx>=sxs;sx--)
        {
            psurf = &((surf_t *)loct->sur.buf)[rptr[sx]];
            cptr[sx] = mulsc(*(int *)&psurf->b,psurf->norm[0]+psurf->norm[1]+psurf->norm[2]+256);
        }
    }
}
void cpu_shader_znotex (int sy,void *_) { cpu_shader_solcol(sy,_); }
void cpu_shader_texmap (int sy,void *_) { cpu_shader_solcol(sy,_); }
void cpu_shader_texmap_mc (int,void *){}

//--------------------------------------------------------------------------------------------------

void oct_drawoct (oct_t *loct, point3f *pp, point3f *pr, point3f *pd, point3f *pf, float mixval, int imulcol)
{
    point3f vt[8];
    float f, g, h, fx, fy, fz, mat[9];
    int i, j, k, x, y, z, ls, split[2][4+4+1+512], splitn[2], rectn, vaneg[24], dorast;
    static int lmost[MAXYDIM], rmost[MAXYDIM], xmin, ymin, xmax, ymax;
    static const int ind[6][4] = {1,5,7,3, 4,0,2,6, 2,3,7,6, 4,5,1,0, 5,4,6,7, 0,1,3,2};

    if constexpr (oct_usegpu && oct_usegpubo) {
        if (loct->gsurf) {
            loct->gsurf = 0;
            ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_2D,loct->octid);
            bo_end(loct->bufid,0,0,loct->gxsid,loct->gysid,GL_RGBA,GL_UNSIGNED_BYTE,0);
        }
    }

    if (!((octv_t *)loct->nod.buf)[loct->head].chi) return;

    invert3x3(pr,pd,pf,mat);
    fx = gipos.x-pp->x; fy = gipos.y-pp->y; fz = gipos.z-pp->z;
    gorig.x = fx*mat[0] + fy*mat[1] + fz*mat[2];
    gorig.y = fx*mat[3] + fy*mat[4] + fz*mat[5];
    gorig.z = fx*mat[6] + fy*mat[7] + fz*mat[8];

        //FIXFIXFIX:ugly hack to avoid visual artifact :/
    #define HACKTHRESH 1e-4
    f = gorig.x-floor(gorig.x+.5); if (fabs(f) < HACKTHRESH) { if (f >= 0) gorig.x += HACKTHRESH-f; else gorig.x -= HACKTHRESH+f; }
    f = gorig.y-floor(gorig.y+.5); if (fabs(f) < HACKTHRESH) { if (f >= 0) gorig.y += HACKTHRESH-f; else gorig.y -= HACKTHRESH+f; }
    f = gorig.z-floor(gorig.z+.5); if (fabs(f) < HACKTHRESH) { if (f >= 0) gorig.z += HACKTHRESH-f; else gorig.z -= HACKTHRESH+f; }

    gxmul.x = girig.x*pr->x + girig.y*pr->y + girig.z*pr->z;
    gxmul.y = gidow.x*pr->x + gidow.y*pr->y + gidow.z*pr->z;
    gxmul.z = gifor.x*pr->x + gifor.y*pr->y + gifor.z*pr->z;
    gymul.x = girig.x*pd->x + girig.y*pd->y + girig.z*pd->z;
    gymul.y = gidow.x*pd->x + gidow.y*pd->y + gidow.z*pd->z;
    gymul.z = gifor.x*pd->x + gifor.y*pd->y + gifor.z*pd->z;
    gzmul.x = girig.x*pf->x + girig.y*pf->y + girig.z*pf->z;
    gzmul.y = gidow.x*pf->x + gidow.y*pf->y + gidow.z*pf->z;
    gzmul.z = gifor.x*pf->x + gifor.y*pf->y + gifor.z*pf->z;
    gnadd.x = -(gorig.x*gxmul.x + gorig.y*gymul.x + gorig.z*gzmul.x);
    gnadd.y = -(gorig.x*gxmul.y + gorig.y*gymul.y + gorig.z*gzmul.y);
    gnadd.z = -(gorig.x*gxmul.z + gorig.y*gymul.z + gorig.z*gzmul.z);

    f = 1.0/(float)loct->sid;
    grxmul.x = (girig.x*mat[0] + girig.y*mat[1] + girig.z*mat[2])*f;
    grymul.x = (gidow.x*mat[0] + gidow.y*mat[1] + gidow.z*mat[2])*f;
    grzmul.x = (gifor.x*mat[0] + gifor.y*mat[1] + gifor.z*mat[2])*f;
    grxmul.y = (girig.x*mat[3] + girig.y*mat[4] + girig.z*mat[5])*f;
    grymul.y = (gidow.x*mat[3] + gidow.y*mat[4] + gidow.z*mat[5])*f;
    grzmul.y = (gifor.x*mat[3] + gifor.y*mat[4] + gifor.z*mat[5])*f;
    grxmul.z = (girig.x*mat[6] + girig.y*mat[7] + girig.z*mat[8])*f;
    grymul.z = (gidow.x*mat[6] + gidow.y*mat[7] + gidow.z*mat[8])*f;
    grzmul.z = (gifor.x*mat[6] + gifor.y*mat[7] + gifor.z*mat[8])*f;

    f = 1.0f; //> 1 here shrinks cubes: weird effect, but slower fps

    // Enable this line for crazy cube shrinking effect :P
    if constexpr (false) {
        double d = time_now_seconds();
        f = 1.0 / (sin(d*4.0) * 0.5 + 0.5);
    }

    fx = ghx*f; fy = ghy*f; fz = ghz*f;
    fgxmul.x = gxmul.x*fz + gxmul.z*fx; fgxmul.y = gxmul.y*fz + gxmul.z*fy; fgxmul.z = gxmul.z*f;
    fgymul.x = gymul.x*fz + gymul.z*fx; fgymul.y = gymul.y*fz + gymul.z*fy; fgymul.z = gymul.z*f;
    fgzmul.x = gzmul.x*fz + gzmul.z*fx; fgzmul.y = gzmul.y*fz + gzmul.z*fy; fgzmul.z = gzmul.z*f;
    fgnadd.x = gnadd.x*fz + gnadd.z*fx; fgnadd.y = gnadd.y*fz + gnadd.z*fy; fgnadd.z = gnadd.z*f;

                       gposadd[  0][0].x =                        0; gposadd[  0][0].y =                        0; gposadd[  0][0].z =                        0;
                       gposadd[  1][0].x =                 fgxmul.x; gposadd[  1][0].y =                 fgxmul.y; gposadd[  1][0].z =                 fgxmul.z;
    for(j=0;j<2;j++) { gposadd[j+2][0].x = gposadd[j][0].x+fgymul.x; gposadd[j+2][0].y = gposadd[j][0].y+fgymul.y; gposadd[j+2][0].z = gposadd[j][0].z+fgymul.z; }
    for(j=0;j<4;j++) { gposadd[j+4][0].x = gposadd[j][0].x+fgzmul.x; gposadd[j+4][0].y = gposadd[j][0].y+fgzmul.y; gposadd[j+4][0].z = gposadd[j][0].z+fgzmul.z; }
    for(j=0;j<8;j++) gposadd[j][0].z2 = -gposadd[j][0].z;
    for(j=0;j<8;j++)
    {
#if !defined(_WIN32)
        for(i=1;i<loct->lsid;i++)
        {
            gposadd[j][i].z = gposadd[j][i-1].z *2;
            gposadd[j][i].z2= gposadd[j][i-1].z2*2;
            gposadd[j][i].x = gposadd[j][i-1].x *2;
            gposadd[j][i].y = gposadd[j][i-1].y *2;
        }
#else
        i = (int)&gposadd[j][           0];
        k = (int)&gposadd[j][loct->lsid-1];
        __asm
        {
            mov eax, i
            movaps xmm0, [eax]
        beg:
            add eax, 16
            addps xmm0, xmm0
            movaps [eax], xmm0
            cmp eax, k
            jb short beg
        }
#endif
    }

    glorig.x = (int)floor(gorig.x);
    glorig.y = (int)floor(gorig.y);
    glorig.z = (int)floor(gorig.z);
    for(ls=loct->lsid-1;ls>=0;ls--)
    {
        i = (1<<ls);
        gldirlut[ls].x = std::min(std::max(glorig.x - i,-32768),32767);
        gldirlut[ls].y = std::min(std::max(glorig.y - i,-32768),32767);
        gldirlut[ls].z = std::min(std::max(glorig.z - i,-32768),32767);
        gldirlut[ls].dum = 32767;
    }

    glcenadd.x = -std::min(std::max(glorig.x,-1),loct->sid)+32766;   //-( -1)+32766 = 32767
    glcenadd.y = -std::min(std::max(glorig.y,-1),loct->sid)+32766;   //-(256)+32766 = 32510
    glcenadd.z = -std::min(std::max(glorig.z,-1),loct->sid)+32766;
    glcenadd.dum = 32767;
    glcensub.x = 32765;
    glcensub.y = 32765;
    glcensub.z = 32765;
    glcensub.dum = -1;

    f = 3.0*sqrt(pr->z*pr->z + pd->z*pd->z + pf->z*pf->z); g = f;
    fx = std::min(fgxmul.z,0.0f) + std::min(fgymul.z,0.0f) + std::min(fgzmul.z,0.0f);
    fy = std::max(fgxmul.z,0.0f) + std::max(fgymul.z,0.0f) + std::max(fgzmul.z,0.0f);
    h = oct_fogdist;
    for(i=0;i<loct->lsid;i++,fx+=fx,fy+=fy)
    {
        cornmin[i] = f - fx;
        cornmax[i] = g - fy;
        scanmax[i] = h - fx; //0x7f7fffff;
    }
    cornmax[0] = 1e32f; //Prevent goto tochild when (!ls)

    if constexpr (oct_usegpu)
    {
        //FIX:worth it to read zbuf from GPU?
        //Clear only area that cube covers
                           vt[  0].x =         gnadd.x;           vt[  0].y =         gnadd.y;           vt[  0].z =         gnadd.z;
                           vt[  1].x = gnadd.x+gxmul.x*loct->sid; vt[  1].y = gnadd.y+gxmul.y*loct->sid; vt[  1].z = gnadd.z+gxmul.z*loct->sid;
        for(i=0;i<2;i++) { vt[i+2].x = vt[i].x+gymul.x*loct->sid; vt[i+2].y = vt[i].y+gymul.y*loct->sid; vt[i+2].z = vt[i].z+gymul.z*loct->sid; }
        for(i=0;i<4;i++) { vt[i+4].x = vt[i].x+gzmul.x*loct->sid; vt[i+4].y = vt[i].y+gzmul.y*loct->sid; vt[i+4].z = vt[i].z+gzmul.z*loct->sid; }

        dorast = 0;
        for(i=0;i<3;i++)
        {
            switch(i)
            {
                case 0: {
                    fx = gxmul.y * gymul.z - gxmul.z * gymul.y;
                    fy = gxmul.z * gymul.x - gxmul.x * gymul.z;
                    fz = gxmul.x * gymul.y - gxmul.y * gymul.x;
                } break;
                case 1: {
                    fx = gymul.y * gzmul.z - gymul.z * gzmul.y;
                    fy = gymul.z * gzmul.x - gymul.x * gzmul.z;
                    fz = gymul.x * gzmul.y - gymul.y * gzmul.x;
                } break;
                case 2: {
                    fx = gzmul.y * gxmul.z - gzmul.z * gxmul.y;
                    fy = gzmul.z * gxmul.x - gzmul.x * gxmul.z;
                    fz = gzmul.x * gxmul.y - gzmul.y * gxmul.x;
                } break;
            }
            f = std::min(-(vt[0].x*fx + vt[0].y*fy + vt[0].z*fz),(vt[7].x*fx + vt[7].y*fy + vt[7].z*fz));
            if (f >= 0.0) continue;
            if (f*f <= (fx*fx + fy*fy + fz*fz)*SCISDIST*SCISDIST*16*16) continue; //SCISDIST*16:assumes not zoomed out a ridiculous amount
            dorast = 1; break;
        }
        if (dorast)
        {
            gymin = 0x7fffffff; gymax = 0x80000000;
            for(j=3-1;j>=0;j--)
            {
                i = ((int *)&glorig)[j]; if ((unsigned)i < (unsigned)loct->sid) continue;
                i = (((unsigned)i)>>31)+j*2;
                if (gymin >= gymax)
                {
                    rastquad(loct->xres, loct->yres, vt,ind[i],glmost,grmost,&gxmin,&gymin,&gxmax,&gymax);
                }
                else
                {
                    rastquad(loct->xres, loct->yres, vt,ind[i],lmost,rmost,&xmin,&ymin,&xmax,&ymax);
                    if (ymin >= ymax) continue;

                        //take union of rasts
                    for(y=std::max(ymin,gymin),k=std::min(ymax,gymax);y<k;y++)
                    {
                        glmost[y] = std::min(glmost[y],lmost[y]);
                        grmost[y] = std::max(grmost[y],rmost[y]);
                    }
                    while (ymin < gymin) { gymin--; glmost[gymin] = lmost[gymin]; grmost[gymin] = rmost[gymin];          }
                    while (ymax > gymax) {          glmost[gymax] = lmost[gymax]; grmost[gymax] = rmost[gymax]; gymax++; }
                    gxmin = std::min(gxmin,xmin); gxmax = std::max(gxmax,xmax);
                }
            }
            if (gymin >= gymax) { gymin = 0; gymax = 0; }
        }
        else
        {
            gxmin = 0; gymin = 0; gxmax = loct->xres; gymax = loct->yres;
            memset4(glmost,         0,loct->yres*4);
            memset4(grmost,loct->xres,loct->yres*4);
        }
        setgcbit(loct->xres, loct->yres, glmost,grmost,gymin,gymax);
    }
    else
    {
        memset16((void *)gcbit,-1,gcbpl*loct->yres); //Clear cover buffer

        gxmin = 0; gymin = 0; gxmax = loct->xres; gymax = loct->yres;
        memset4(glmost,         0,loct->yres*4);
        memset4(grmost,loct->xres,loct->yres*4);
    }

    g = 65536.f*65536.f;
    if (fabs(gxmul.z)*g < ghz) { f = g; if (gxmul.z < 0) f = -g; } else f = ghz/gxmul.z;
    gvanx[0] = gxmul.x*f + ghx;
    gvany[0] = gxmul.y*f + ghy;
    if (fabs(gymul.z)*g < ghz) { f = g; if (gymul.z < 0) f = -g; } else f = ghz/gymul.z;
    gvanx[1] = gymul.x*f + ghx;
    gvany[1] = gymul.y*f + ghy;
    if (fabs(gzmul.z)*g < ghz) { f = g; if (gzmul.z < 0) f = -g; } else f = ghz/gzmul.z;
    gvanx[2] = gzmul.x*f + ghx;
    gvany[2] = gzmul.y*f + ghy;
    for(i=0;i<3;i++) { gvanxmh[i] = gvanx[i]-.5; gvanymh[i] = gvany[i]-.5; }

    for(i=0;i<4;i++)
    {
        fx = 0.f; fy = 0.f; fz = 0.f;
        if ((i  )&1) { fx += gxmul.x; fy += gxmul.y; fz += gxmul.z; }
        if ((i+1)&2) { fx += gymul.x; fy += gymul.y; fz += gymul.z; }
        if ((i  )&2) { fx += gzmul.x; fy += gzmul.y; fz += gzmul.z; }
        fx = fx*ghz + fz*ghx;
        fy = fy*ghz + fz*ghy;
        cornvadd[0][i] = gvanx[0]*fz - fx;
        cornvadd[1][i] = gvanx[1]*fz - fx;
        cornvadd[2][i] = gvanx[2]*fz - fx;
        cornvadd[3][i] = gvany[0]*fz - fy;
        cornvadd[4][i] = gvany[1]*fz - fy;
        cornvadd[5][i] = gvany[2]*fz - fy;
    }

    memset(vaneg,0,sizeof(vaneg));
    for(i=0;i<27;i++)
    {
        vause[i][0] = ((i != 14) && (i != 12)); //left/right
        vause[i][1] = ((i != 16) && (i != 10)); //up/down
        vause[i][2] = ((i != 22) && (i !=  4)); //front/back
        vause[i][3] = 0;

        cornvadd2[i][ 0] = cornvadd[3][vanlut[i][2]]*+vause[i][0];
        cornvadd2[i][ 1] = cornvadd[4][vanlut[i][0]]*+vause[i][1];
        cornvadd2[i][ 2] = cornvadd[5][vanlut[i][1]]*+vause[i][2];
        cornvadd2[i][ 3] = 0.f;
        cornvadd2[i][ 4] = cornvadd[0][vanlut[i][2]]*-vause[i][0];
        cornvadd2[i][ 5] = cornvadd[1][vanlut[i][0]]*-vause[i][1];
        cornvadd2[i][ 6] = cornvadd[2][vanlut[i][1]]*-vause[i][2];
        cornvadd2[i][ 7] = 0.f;
        cornvadd2[i][ 8] = cornvadd[3][vanlut[i][0]]*-vause[i][0];
        cornvadd2[i][ 9] = cornvadd[4][vanlut[i][1]]*-vause[i][1];
        cornvadd2[i][10] = cornvadd[5][vanlut[i][2]]*-vause[i][2];
        cornvadd2[i][11] = 0.f;
        cornvadd2[i][12] = cornvadd[0][vanlut[i][0]]*+vause[i][0];
        cornvadd2[i][13] = cornvadd[1][vanlut[i][1]]*+vause[i][1];
        cornvadd2[i][14] = cornvadd[2][vanlut[i][2]]*+vause[i][2];
        cornvadd2[i][15] = 0.f;

        vause[i][0] = -vause[i][0];
        vause[i][1] = -vause[i][1];
        vause[i][2] = -vause[i][2];
    }

    //--------------------------------------------------------------------------------------------------

    //Collect horz&vert splits
    split[0][0] = gxmin; splitn[0] = 1;
    split[1][0] = gymin; splitn[1] = 1;

    if (gxmul.z != 0.f)
    {
        i = (int)(gxmul.x/gxmul.z*ghz + ghx); if ((i > gxmin) && (i < gxmax)) { split[0][splitn[0]] = i; splitn[0]++; }
        i = (int)(gxmul.y/gxmul.z*ghz + ghy); if ((i > gymin) && (i < gymax)) { split[1][splitn[1]] = i; splitn[1]++; }
    }
    if (gymul.z != 0.f)
    {
        i = (int)(gymul.x/gymul.z*ghz + ghx); if ((i > gxmin) && (i < gxmax)) { split[0][splitn[0]] = i; splitn[0]++; }
        i = (int)(gymul.y/gymul.z*ghz + ghy); if ((i > gymin) && (i < gymax)) { split[1][splitn[1]] = i; splitn[1]++; }
    }
    if (gzmul.z != 0.f)
    {
        i = (int)(gzmul.x/gzmul.z*ghz + ghx); if ((i > gxmin) && (i < gxmax)) { split[0][splitn[0]] = i; splitn[0]++; }
        i = (int)(gzmul.y/gzmul.z*ghz + ghy); if ((i > gymin) && (i < gymax)) { split[1][splitn[1]] = i; splitn[1]++; }
    }

    //------------------------------------------------
    if constexpr (true) // (oct_numcpu >= 2)
    {
        #define XSPLITS 0
        #define YSPLITS 16
        for(i=1;i<XSPLITS;i++) { j = (i*loct->xres)/XSPLITS; if ((j > gxmin) && (j < gxmax)) { split[0][splitn[0]] = j; splitn[0]++; } }
        for(i=1;i<YSPLITS;i++) { j = (i*loct->yres)/YSPLITS; if ((j > gymin) && (j < gymax)) { split[1][splitn[1]] = j; splitn[1]++; } }
    }
    //------------------------------------------------

    for(k=2-1;k>=0;k--)
    {
        //Sort x's (split[0][?]) or y's (split[1][?])
        for(i=2;i<splitn[k];i++) {
            for(j=1;j<i;j++) {
                if (split[k][i] < split[k][j]) {
                    z = split[k][i];
                    split[k][i] = split[k][j];
                    split[k][j] = z;
                }
            }
        }

        //Eliminate duplicates
        for(i=1,j=1;i<splitn[k];i++) {
            if (split[k][i] != split[k][j-1]) {
                split[k][j] = split[k][i];
                j++;
            }
        }
        splitn[k] = j;
    }
    split[0][splitn[0]] = gxmax;
    split[1][splitn[1]] = gymax;

//--------------------------------------------------------------------------------------------------

    rectn = 1;
    grect[0].x0 = 0; grect[0].y0 = 0; grect[0].x1 = gxmax; grect[0].y1 = gymax;
    for(x=0;x<splitn[0];x++) {
        //NOTE: Keep y as inner loop - reduces chances of multithread contention (and glitches with xor gcbit in inner loop)
        for(y=0;y<splitn[1];y++) {
            grect[rectn].x0 = split[0][x]; grect[rectn].x1 = split[0][x+1];
            grect[rectn].y0 = split[1][y]; grect[rectn].y1 = split[1][y+1];
            rectn++;
        }
    }

    //f = cone radius angle
    fx = std::max(ghx-0.f,(float)loct->xres-ghx);
    fy = std::max(ghy-0.f,(float)loct->yres-ghy);
    f = atan(sqrt(fx*fx + fy*fy)/ghz);
    g = 1.0/sqrt(gxmul.z*gxmul.z + gymul.z*gymul.z + gzmul.z*gzmul.z); fx = gxmul.z*g; fy = gymul.z*g; fz = gzmul.z*g;
    g = sqrt(3.0)*0.5; //sphere radius
    h = 0.5;
    for(i=0;i<loct->lsid;i++,g+=g,h+=h) isint_cone_sph_init(&giics[i],gorig.x-h,gorig.y-h,gorig.z-h,fx,fy,fz,f,g);

    if constexpr (oct_usegpu)
    {
        x = loct->xres*(PIXBUFBYPP>>2);
        if constexpr (oct_usegpubo) {
            gxd.f = (INT_PTR)bo_begin(gpixbufid,x*loct->yres*PIXBUFBYPP);
            gxd.p = x*4;
            gxd.x = loct->xres;
            gxd.y = loct->yres;
        }
    }

    if (gdrawclose) oct_rendclosecubes(loct,&glorig);
    if (gdrawclose < 2)
    {
        if constexpr (oct_useasm) {
            oct_ftob2_dorect(-1,loct);
        }
        htrun(oct_ftob2_dorect,loct,1,rectn,true);
    }

    if constexpr (oct_usegpu)
    {
        htrun(clearbackpixes,loct,gymin,gymax,true);

        //memcpy xres*yres*PIXBUFBYPP bytes from CPU->GPU
        ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_2D,gpixtexid);
        if constexpr (oct_usegpubo) {
            bo_end(gpixbufid,0,gymin,x,gymax-gymin,GL_RGBA,GL_UNSIGNED_BYTE,gymin*gxd.p);
        }
        else {
            ((PFNGLTEXSUBIMAGE2D)glfp[glTexSubImage2D])(GL_TEXTURE_2D,0,0,gymin,x,gymax-gymin,GL_RGBA,GL_UNSIGNED_BYTE,(void *)(gymin*gxd.p + gxd.f));
        }

        //send vals to shader
        //v.x = sx*grxmul.x + sy*grymul.x + ((.5-ghx)*grxmul.x + (.5-ghy)*grymul.x + ghz*grzmul.x);
        //v.y = sx*grxmul.y + sy*grymul.y + ((.5-ghx)*grxmul.y + (.5-ghy)*grymul.y + ghz*grzmul.y);
        //v.z = sx*grxmul.z + sy*grymul.z + ((.5-ghx)*grxmul.z + (.5-ghy)*grymul.z + ghz*grzmul.z);
        //v = t.x*vmulx + t.y*vmuly + vadd; <- in shader
        if constexpr (oct_useglsl) {
            ((PFNGLUSEPROGRAMPROC)glfp[glUseProgram])(shadprog[shadcur]);

            kglUniform1f(kglGetUniformLoc(shadprog[shadcur], "depthmul"),-gznear*gzfar/(gzfar-gznear)*(float)loct->sid);
            kglUniform1f(kglGetUniformLoc(shadprog[shadcur], "depthadd"),gzfar/(gzfar-gznear));

            kglUniform4f(kglGetUniformLoc(shadprog[shadcur], "pos"),gorig.x,gorig.y,gorig.z,0.f);
            f = 1.0/ghz;
            kglUniform4f(kglGetUniformLoc(shadprog[shadcur], "vmulx"),grxmul.x*gpixxdim*f,
                                                                grxmul.y*gpixxdim*f,
                                                                grxmul.z*gpixxdim*f,0.f);
            kglUniform4f(kglGetUniformLoc(shadprog[shadcur], "vmuly"),grymul.x*gpixydim*f,
                                                                grymul.y*gpixydim*f,
                                                                grymul.z*gpixydim*f,0.f);
            kglUniform4f(kglGetUniformLoc(shadprog[shadcur], "vadd"),(-ghx*grxmul.x - ghy*grymul.x + ghz*grzmul.x)*f,
                                                              (-ghx*grxmul.y - ghy*grymul.y + ghz*grzmul.y)*f,
                                                              (-ghx*grxmul.z - ghy*grymul.z + ghz*grzmul.z)*f,0.f);
            kglUniform1f(kglGetUniformLoc(shadprog[shadcur], "mipoffs"),-ghz/ghx + (float)tiles[0].ltilesid - loct->lsid*2 - log(pr->z*pr->z + pd->z*pd->z + pf->z*pf->z)/log(2.f) - 9.f); //-=sharp, +=blurry
            kglUniform1f(kglGetUniformLoc(shadprog[shadcur], "maxmip"),(float)tiles[0].ltilesid);
            kglUniform1f(kglGetUniformLoc(shadprog[shadcur], "mixval"),mixval);
            kglUniform4f(kglGetUniformLoc(shadprog[shadcur], "mulcol"),(float)((imulcol>>16)&255)/64.0,(float)((imulcol>>8)&255)/64.0,(float)(imulcol&255)/64.0,(float)(((unsigned)imulcol)>>24)/256.0);
            kglUniform4f(kglGetUniformLoc(shadprog[shadcur], "fogcol"),(float)((oct_fogcol>>16)&255)/255.0,(float)((oct_fogcol>>8)&255)/255.0,(float)(oct_fogcol&255)/255.0,(float)(((unsigned)oct_fogcol)>>24)/255.0);
            kglUniform1f(kglGetUniformLoc(shadprog[shadcur], "fogdist"),1.0/(oct_fogdist*(double)loct->sid));

            kglUniform1f(kglGetUniformLoc(shadprog[shadcur], "rxsid"),1.0/(float)loct->gxsid);
            kglUniform1f(kglGetUniformLoc(shadprog[shadcur], "rysid"),1.0/(float)loct->gysid);
            kglUniform1i(kglGetUniformLoc(shadprog[shadcur], "axsidm1"),loct->gxsid-1);
            kglUniform1i(kglGetUniformLoc(shadprog[shadcur], "lxsidm1"),loct->glxsid-1);

            for(i=0;i<6;i++) kglUniform1f(kglGetUniformLoc(shadprog[shadcur], "gsideshade")+i,(float)(oct_sideshade[i]&255)/255.f);
        }
        else
        {
            //kglProgramLocalParam(15?,drawconedat.h.x/drawconedat.h.z,0,0,0); //zoom
            shadcur = 0;
            ((PFNGLENABLE)glfp[glEnable])(GL_VERTEX_PROGRAM_ARB);
            ((PFNGLENABLE)glfp[glEnable])(GL_FRAGMENT_PROGRAM_ARB);
            ((PFNGLBINDPROGRAMARBPROC)glfp[glBindProgramARB])(GL_VERTEX_PROGRAM_ARB  ,shadvert[shadcur]);
            ((PFNGLBINDPROGRAMARBPROC)glfp[glBindProgramARB])(GL_FRAGMENT_PROGRAM_ARB,shadfrag[shadcur]);

            kglProgramLocalParam( 0,1.0/((float)loct->gxsid),2.0/((float)loct->gxsid),0.0,0.0); //sidmul0
            kglProgramLocalParam( 1,2.0,1.0/((float)loct->gysid),0.0,0.0); //sidmul1
            kglProgramLocalParam( 2,0.5*(float)loct->gxsid,0.5*(float)loct->gysid,0.0,0.0); //sidadd1
            kglProgramLocalParam( 3,-gznear*gzfar/(gzfar-gznear)*(float)loct->sid,gzfar/(gzfar-gznear),0.0,0.0); //depth

            f = 1.0/ghz;
            kglProgramLocalParam( 4,grxmul.x*gpixxdim*f,grxmul.y*gpixxdim*f,grxmul.z*gpixxdim*f,0.f); //vmulx
            kglProgramLocalParam( 5,grymul.x*gpixydim*f,grymul.y*gpixydim*f,grymul.z*gpixydim*f,0.f); //vmuly
            kglProgramLocalParam( 6,(-ghx*grxmul.x - ghy*grymul.x + ghz*grzmul.x)*f, (-ghx*grxmul.y - ghy*grymul.y + ghz*grzmul.y)*f, (-ghx*grxmul.z - ghy*grymul.z + ghz*grzmul.z)*f,0.f); //vadd
            kglProgramLocalParam( 7,gorig.x,gorig.y,gorig.z,0.f); //pos

            kglProgramLocalParam( 8,(float)1,1.0/(oct_fogdist*(double)loct->sid),mixval,0.0); //mixnfog
            kglProgramLocalParam( 9,(float)((imulcol>>16)&255)/64.0,(float)((imulcol>>8)&255)/64.0,(float)(imulcol&255)/64.0,(float)(((unsigned)imulcol)>>24)/256.0);    //mulcol
            kglProgramLocalParam(10,(float)((oct_fogcol>>16)&255)/255.0,(float)((oct_fogcol>>8)&255)/255.0,(float)(oct_fogcol&255)/255.0,(float)(((unsigned)oct_fogcol)>>24)/255.0); //fogcol

            kglProgramLocalParam(11,(float)(oct_sideshade[0]&255)/255.f,(float)(oct_sideshade[2]&255)/255.f,(float)(oct_sideshade[4]&255)/255.f,0.0); //gsideshade024
            kglProgramLocalParam(12,(float)(oct_sideshade[1]&255)/255.f,(float)(oct_sideshade[3]&255)/255.f,(float)(oct_sideshade[5]&255)/255.f,0.0); //gsideshade135

            kglProgramLocalParam(13,-ghz/ghx + (float)tiles[0].ltilesid - loct->lsid*2 - log(pr->z*pr->z + pd->z*pd->z + pf->z*pf->z)/log(2.f) - 9.f, (float)tiles[0].ltilesid,0.0,0.0); //mipdat (mipoffs,maxmip,0,0)
        }
        kglActiveTexture(0); ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_2D,loct->tilid);
        kglActiveTexture(1); ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_2D,gpixtexid   );
        kglActiveTexture(2); ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_2D,loct->octid);

        //render fullscreen quad using shader
        fx = (float)loct->xres/gpixxdim;
        fy = (float)loct->yres/gpixydim;
        g = 256.f;
        f = (float)loct->yres/(float)loct->xres*g;
        if (dorast)
        {
            for(j=3-1;j>=0;j--)
            {
                point3f vt2[8];
                float ff, gg, hh, fsx, fsy;
                int n2, m;
                #define SCISDIST 5e-4f

                i = ((int *)&glorig)[j]; if ((unsigned)i < (unsigned)loct->sid) continue;
                i = (((unsigned)i)>>31)+j*2;

                for(k=4-1,m=0,n2=0;m<4;k=m,m++)
                {
                    if (vt[ind[i][k]].z >= SCISDIST) { vt2[n2] = vt[ind[i][k]]; n2++; }
                    if ((vt[ind[i][k]].z >= SCISDIST) != (vt[ind[i][m]].z >= SCISDIST))
                    {
                        ff = (SCISDIST-vt[ind[i][m]].z)/(vt[ind[i][k]].z-vt[ind[i][m]].z);
                        vt2[n2].x = (vt[ind[i][k]].x-vt[ind[i][m]].x)*ff + vt[ind[i][m]].x;
                        vt2[n2].y = (vt[ind[i][k]].y-vt[ind[i][m]].y)*ff + vt[ind[i][m]].y;
                        vt2[n2].z = SCISDIST; n2++;
                    }
                }
                if (n2 < 3) continue;

                gg = ghz/ghx;
                hh = (float)loct->xres/(float)loct->yres;
                ((PFNGLBEGIN)glfp[glBegin])(GL_TRIANGLE_FAN);
                for(i=0;i<n2;i++)
                {
                    ff = gg/vt2[i].z;
                    fsx = vt2[i].x*ff;
                    fsy = vt2[i].y*ff*hh;
                    ((PFNGLTEXCOORD2F)glfp[glTexCoord2f])((fsx*.5+.5)*fx,(fsy*.5+.5)*fy); ((PFNGLVERTEX3F)glfp[glVertex3f])(g*fsx,-f*fsy,-g);
                }
                ((PFNGLEND)glfp[glEnd])();
            }
        }
        else
        {
            ((PFNGLBEGIN)glfp[glBegin])(GL_QUADS);
            ((PFNGLTEXCOORD2F)glfp[glTexCoord2f])( 0, 0); ((PFNGLVERTEX3F)glfp[glVertex3f])(-g,+f,-g);
            ((PFNGLTEXCOORD2F)glfp[glTexCoord2f])(fx, 0); ((PFNGLVERTEX3F)glfp[glVertex3f])(+g,+f,-g);
            ((PFNGLTEXCOORD2F)glfp[glTexCoord2f])(fx,fy); ((PFNGLVERTEX3F)glfp[glVertex3f])(+g,-f,-g);
            ((PFNGLTEXCOORD2F)glfp[glTexCoord2f])( 0,fy); ((PFNGLVERTEX3F)glfp[glVertex3f])(-g,-f,-g);
            ((PFNGLEND)glfp[glEnd])();
        }
    }
    else {
        gimixval = std::min(std::max((int)(mixval*32768.0),0),32767); gimulcol = imulcol;
        gzscale = ghz/loct->sid;
        gmipoffs = (float)klog2up7((tiles[0].x>>5)/(sqrt(pr->x*pr->x + pd->y*pd->y + pf->z*pf->z)*loct->sid));
        static int otiles0x = 0;
        if (tiles[0].x != otiles0x)
        {
            otiles0x = tiles[0].x;
            for(i=0;i<16;i++) gtiloff[i] = i*((tiles[0].x>>4)<<LTEXPREC) + ((tiles[0].x>>6)<<LTEXPREC);
            for(i=0;i<256;i++) { dqtiloff[i][0] = gtiloff[i&15]; dqtiloff[i][1] = gtiloff[i>>4]; }
        }
        for(i=0;i<6;i++) { dqgsideshade[i][0] = (oct_sideshade[i]&0xff) + ((oct_sideshade[i]&0xff00)<<8); dqgsideshade[i][1] = ((oct_sideshade[i]&0xff0000)>>16); }

        if ((shaderfunc == cpu_shader_znotex) || (shaderfunc == cpu_shader_texmap) || (shaderfunc == cpu_shader_texmap_mc))
        {
            dqgrxmul[0] = grxmul.x; dqgrxmul[1] = grxmul.y; dqgrxmul[2] = grxmul.z; dqgrxmul[3] = 0.f;
            dqgorig[0] = gorig.x; dqgorig[1] = gorig.y; dqgorig[2] = gorig.z; dqgorig[3] = 0.f;
            dqfogcol[0] = (oct_fogcol&255); dqfogcol[1] = ((oct_fogcol>>8)&255); dqfogcol[2] = ((oct_fogcol>>16)&255);
            dqimulcol[0] = (gimulcol&255); dqimulcol[1] = ((gimulcol>>8)&255); dqimulcol[2] = ((gimulcol>>16)&255);
            if ((shaderfunc == cpu_shader_texmap) || (shaderfunc == cpu_shader_texmap_mc))
            {
                float fwmul;
                int iwmsk, tilx5, ltilmsk;

                ltilmsk = std::max(((tiles[0].ltilesid-1)<<LTEXPREC) - 1,0); //ltilesid-1 gives 2x2
                tilx5 = (tiles[0].x>>5); iwmsk = (tilx5<<LTEXPREC)-1; fwmul = (float)(iwmsk+1);

                dqfwmul[0] = fwmul; dqfwmul[1] = fwmul; dqfwmul[2] = 0.f; dqfwmul[3] = 0.f;
                dqiwmsk[0] = iwmsk; dqiwmsk[1] = iwmsk; dqiwmsk[2] = 0; dqiwmsk[3] = 0;

                for(i=0;i<16;i++) { dqtilepitch4[i][0] = 4; dqtilepitch4[i][1] = tiles[i].p; tilesf[i] = tiles[i].f; } //FIXFIXFIXFIX:reduce count/optimize w/cache!
                dqmixval[0] = gimixval; dqmixval[1] = gimixval; dqmixval[2] = gimixval; dqmixval[3] = gimixval;
            }
        }
        htrun(shaderfunc,loct,gymin,gymax,true);
    }

#if 0
    //Debug only: show split lines
    static int bozocnt = 0; bozocnt++;
    for(i=0;i<splitn[1];i++) { for(x=(bozocnt&7);x<gdd.x;x+=8) *(int *)(gdd.p*split[1][i] + (x<<2) + gdd.f) = 0xc0a080; }
    for(i=0;i<splitn[0];i++) { for(y=(bozocnt&7);y<gdd.y;y+=8) *(int *)(gdd.p*y + (split[0][i]<<2) + gdd.f) = 0xc0a080; }
#endif
}

void oct_drawoct (oct_t *loct, point3d *dp, point3d *dr, point3d *dd, point3d *df, float mixval, int imulcol)
{
    point3f fp, fr, fd, ff;
    fp.x = dp->x; fp.y = dp->y; fp.z = dp->z;
    fr.x = dr->x; fr.y = dr->y; fr.z = dr->z;
    fd.x = dd->x; fd.y = dd->y; fd.z = dd->z;
    ff.x = df->x; ff.y = df->y; ff.z = df->z;
    oct_drawoct(loct,&fp,&fr,&fd,&ff,mixval,imulcol);
}

//--------------------------------------------------------------------------------------------------

static void oct_surf_normalize (surf_t *psurf)
{
    unsigned char *uptr;
    int t, bb, gg, rr, cnt;
    float dx, dy, dz;

    bb = psurf->b; gg = psurf->g; rr = psurf->r; cnt = 256; t = psurf->tex;
    do
    {
        //if (oct_usefilter != 3) uptr = (unsigned char *)((t&255)*4 + tiles[tiles[0].ltilesid+1].f);
        //                   else
        uptr = (unsigned char *)((t&15)*4 + ((t>>4)&15)*tiles[tiles[0].ltilesid+1].p + tiles[tiles[0].ltilesid+1].f);
        //surf->b(old) = cptr[0] * surf->b(new) * 2.0/256.0;
        psurf->b = std::min(std::max(bb*(256/2) / std::max((int)uptr[0],1),0),255);
        psurf->g = std::min(std::max(gg*(256/2) / std::max((int)uptr[1],1),0),255);
        psurf->r = std::min(std::max(rr*(256/2) / std::max((int)uptr[2],1),0),255);
        dx = ((float)uptr[0]) * (float)psurf->b * (2.0/256.0) - bb;
        dy = ((float)uptr[1]) * (float)psurf->g * (2.0/256.0) - gg;
        dz = ((float)uptr[2]) * (float)psurf->r * (2.0/256.0) - rr;
        if (dx*dx + dy*dy + dz*dz < 16*16) break;
        cnt--; if (cnt < 0) break;

        t = (t-1)&255;
    } while (1);
    psurf->tex = t;
}

void oct_normalizecols (oct_t *loct)
{
    typedef struct { octv_t *ptr; int j; } stk_t;
    stk_t stk[OCT_MAXLS];
    surf_t *psurf;
    octv_t *ptr;
    int i, j, ls, ind;

    if constexpr (oct_usegpu) {
        if constexpr (oct_usegpubo) {
            if (!loct->gsurf) loct->gsurf = (surf_t *)bo_begin(loct->bufid,0);
        }
        else {
            ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_2D,loct->octid);
        }
    }

    ls = loct->lsid-1; ptr = &((octv_t *)loct->nod.buf)[loct->head]; j = 8-1;
    while (1)
    {
        i = (1<<j); if (!(ptr->chi&i)) goto tosibly;

        if (ls <= 0)
        {
            ind = popcount8(ptr->chi&(i-1)) + ptr->ind; psurf = &((surf_t *)loct->sur.buf)[ind];
            oct_surf_normalize(psurf);
            if constexpr (oct_usegpu) {
                if constexpr (oct_usegpubo) {
                    memcpy(&loct->gsurf[ind],psurf,loct->sur.siz);
                }
                else {
                    ((PFNGLTEXSUBIMAGE2D)glfp[glTexSubImage2D])(GL_TEXTURE_2D,0,(ind&((loct->gxsid>>1)-1))<<1,ind>>(loct->glysid-1),2,1,GL_RGBA,GL_UNSIGNED_BYTE,(void *)&((surf_t *)loct->sur.buf)[ind]);
                }
            }
            goto tosibly;
        }

        stk[ls].ptr = ptr; stk[ls].j = j; ls--; //2child
        ptr = &((octv_t *)loct->nod.buf)[popcount8(ptr->chi&(i-1)) + ptr->ind]; j = 8-1;
        continue;

tosibly:
        j--; if (j >= 0) continue;
        do { ls++; if (ls >= loct->lsid) goto break2; j = stk[ls].j-1; } while (j < 0); //2parent
        ptr = stk[ls].ptr;
    }
break2:;
}

//Limitations: zs must be multiple of 32
//Based on pnd3dold.c:octgenmip()
typedef struct { int ys, zs, isor; unsigned int *rbit, *wbit; } genmip_voxbits_t;
static short unshuflut[256] = //unshuflut[i] = (i&1) + ((i&2)<<7) + ((i&4)>>1) + ((i&8)<<6) + ((i&16)>>2) + ((i&32)<<5) + ((i&64)>>3) + ((i&128)<<4); //see NEW/PACKBIT.BAS
{
    0x0000,0x0001,0x0100,0x0101,0x0002,0x0003,0x0102,0x0103,0x0200,0x0201,0x0300,0x0301,0x0202,0x0203,0x0302,0x0303,
    0x0004,0x0005,0x0104,0x0105,0x0006,0x0007,0x0106,0x0107,0x0204,0x0205,0x0304,0x0305,0x0206,0x0207,0x0306,0x0307,
    0x0400,0x0401,0x0500,0x0501,0x0402,0x0403,0x0502,0x0503,0x0600,0x0601,0x0700,0x0701,0x0602,0x0603,0x0702,0x0703,
    0x0404,0x0405,0x0504,0x0505,0x0406,0x0407,0x0506,0x0507,0x0604,0x0605,0x0704,0x0705,0x0606,0x0607,0x0706,0x0707,
    0x0008,0x0009,0x0108,0x0109,0x000a,0x000b,0x010a,0x010b,0x0208,0x0209,0x0308,0x0309,0x020a,0x020b,0x030a,0x030b,
    0x000c,0x000d,0x010c,0x010d,0x000e,0x000f,0x010e,0x010f,0x020c,0x020d,0x030c,0x030d,0x020e,0x020f,0x030e,0x030f,
    0x0408,0x0409,0x0508,0x0509,0x040a,0x040b,0x050a,0x050b,0x0608,0x0609,0x0708,0x0709,0x060a,0x060b,0x070a,0x070b,
    0x040c,0x040d,0x050c,0x050d,0x040e,0x040f,0x050e,0x050f,0x060c,0x060d,0x070c,0x070d,0x060e,0x060f,0x070e,0x070f,
    0x0800,0x0801,0x0900,0x0901,0x0802,0x0803,0x0902,0x0903,0x0a00,0x0a01,0x0b00,0x0b01,0x0a02,0x0a03,0x0b02,0x0b03,
    0x0804,0x0805,0x0904,0x0905,0x0806,0x0807,0x0906,0x0907,0x0a04,0x0a05,0x0b04,0x0b05,0x0a06,0x0a07,0x0b06,0x0b07,
    0x0c00,0x0c01,0x0d00,0x0d01,0x0c02,0x0c03,0x0d02,0x0d03,0x0e00,0x0e01,0x0f00,0x0f01,0x0e02,0x0e03,0x0f02,0x0f03,
    0x0c04,0x0c05,0x0d04,0x0d05,0x0c06,0x0c07,0x0d06,0x0d07,0x0e04,0x0e05,0x0f04,0x0f05,0x0e06,0x0e07,0x0f06,0x0f07,
    0x0808,0x0809,0x0908,0x0909,0x080a,0x080b,0x090a,0x090b,0x0a08,0x0a09,0x0b08,0x0b09,0x0a0a,0x0a0b,0x0b0a,0x0b0b,
    0x080c,0x080d,0x090c,0x090d,0x080e,0x080f,0x090e,0x090f,0x0a0c,0x0a0d,0x0b0c,0x0b0d,0x0a0e,0x0a0f,0x0b0e,0x0b0f,
    0x0c08,0x0c09,0x0d08,0x0d09,0x0c0a,0x0c0b,0x0d0a,0x0d0b,0x0e08,0x0e09,0x0f08,0x0f09,0x0e0a,0x0e0b,0x0f0a,0x0f0b,
    0x0c0c,0x0c0d,0x0d0c,0x0d0d,0x0c0e,0x0c0f,0x0d0e,0x0d0f,0x0e0c,0x0e0d,0x0f0c,0x0f0d,0x0e0e,0x0e0f,0x0f0e,0x0f0f,
};
static void genmip_voxbits (int x, void *vptr)
{
    genmip_voxbits_t *gm = (genmip_voxbits_t *)vptr;
    int i, j, i0, i1, y, z, rxpit, rypit, wxpit, wypit;

    wxpit = gm->zs;        rxpit = (wxpit<<1);
    wypit = gm->ys*gm->zs; rypit = (wypit<<2);
    if (gm->zs >= 16)
    {
        if (gm->isor)
        {
            for(y=0;y<gm->ys;y++)
            {
                i0 = (x<<1)*rxpit + (y<<1)*rypit;
                i1 =  x    *wxpit +  y    *wypit; i1 >>= 4;
                for(z=0;z<gm->zs;z+=16,i0+=32,i1++)
                {
                    i = gm->rbit[(i0      )>>5] | gm->rbit[(i0+rxpit      )>>5] |
                         gm->rbit[(i0+rypit)>>5] | gm->rbit[(i0+rxpit+rypit)>>5];
                    j = (((i>>1)|i)&0x55555555); j += (j>>15); j = (unshuflut[(j>>8)&255]<<4) + unshuflut[j&255]; //LUT: fastest on Ken's C2Q
                    ((short *)gm->wbit)[i1] = (short)j;
                }
            }
        }
        else
        {
            for(y=0;y<gm->ys;y++)
            {
                i0 = (x<<1)*rxpit + (y<<1)*rypit;
                i1 =  x    *wxpit +  y    *wypit; i1 >>= 4;
                for(z=0;z<gm->zs;z+=16,i0+=32,i1++)
                {
                    i = gm->rbit[(i0      )>>5] & gm->rbit[(i0+rxpit      )>>5] &
                         gm->rbit[(i0+rypit)>>5] & gm->rbit[(i0+rxpit+rypit)>>5];
                    j = (((i>>1)&i)&0x55555555); j += (j>>15); j = (unshuflut[(j>>8)&255]<<4) + unshuflut[j&255]; //LUT: fastest on Ken's C2Q
                    ((short *)gm->wbit)[i1] = (short)j;
                }
            }
        }
    }
    else
    {
        if (gm->isor)
        {
            for(y=0;y<gm->ys;y++)
            {
                i0 = (x<<1)*rxpit + (y<<1)*rypit;
                i1 =  x    *wxpit +  y    *wypit;
                for(z=0;z<gm->zs;z++,i0+=2,i1++)
                {
                    i = 0;
                    j = i0            ; if (gm->rbit[j>>5]&(3<<j)) i = 1;
                    j = i0+rxpit      ; if (gm->rbit[j>>5]&(3<<j)) i = 1;
                    j = i0+rypit      ; if (gm->rbit[j>>5]&(3<<j)) i = 1;
                    j = i0+rxpit+rypit; if (gm->rbit[j>>5]&(3<<j)) i = 1;
                    if (!i) gm->wbit[i1>>5] &=~(1<<i1);
                        else gm->wbit[i1>>5] |= (1<<i1);
                }
            }
        }
        else
        {
            for(y=0;y<gm->ys;y++)
            {
                i0 = (x<<1)*rxpit + (y<<1)*rypit;
                i1 =  x    *wxpit +  y    *wypit;
                for(z=0;z<gm->zs;z++,i0+=2,i1++)
                {
                    i = 1;
                    j = i0            ; if (((gm->rbit[j>>5]>>j)&3) != 3) i = 0;
                    j = i0+rxpit      ; if (((gm->rbit[j>>5]>>j)&3) != 3) i = 0;
                    j = i0+rypit      ; if (((gm->rbit[j>>5]>>j)&3) != 3) i = 0;
                    j = i0+rxpit+rypit; if (((gm->rbit[j>>5]>>j)&3) != 3) i = 0;
                    if (!i) gm->wbit[i1>>5] &=~(1<<i1);
                        else gm->wbit[i1>>5] |= (1<<i1);
                }
            }
        }
    }
}

//NOTE: font is stored vertically first! (like .ART files)
const uint64_t font6x8[] = //256 DOS chars, from: DOSAPP.FON (tab blank)
{
    0x3E00000000000000,0x6F6B3E003E455145,0x1C3E7C3E1C003E6B,0x3000183C7E3C1800,
    0x7E5C180030367F36,0x000018180000185C,0x0000FFFFE7E7FFFF,0xDBDBC3FF00000000,
    0x0E364A483000FFC3,0x6000062979290600,0x0A7E600004023F70,0x2A1C361C2A003F35,
    0x0800081C3E7F0000,0x7F361400007F3E1C,0x005F005F00001436,0x22007F017F090600,
    0x606060002259554D,0x14B6FFB614000060,0x100004067F060400,0x3E08080010307F30,
    0x08083E1C0800081C,0x0800404040407800,0x3F3C3000083E083E,0x030F3F0F0300303C,
    0x0000000000000000,0x0003070000065F06,0x247E247E24000307,0x630000126A2B2400,
    0x5649360063640813,0x0000030700005020,0x00000000413E0000,0x1C3E080000003E41,
    0x08083E080800083E,0x0800000060E00000,0x6060000008080808,0x0204081020000000,
    0x00003E4549513E00,0x4951620000407F42,0x3649494922004649,0x2F00107F12141800,
    0x494A3C0031494949,0x0305097101003049,0x0600364949493600,0x6C6C00001E294949,
    0x00006CEC00000000,0x2400004122140800,0x2241000024242424,0x0609590102000814,
    0x7E001E555D413E00,0x49497F007E111111,0x224141413E003649,0x7F003E4141417F00,
    0x09097F0041494949,0x7A4949413E000109,0x00007F0808087F00,0x4040300000417F41,
    0x412214087F003F40,0x7F00404040407F00,0x04027F007F020402,0x3E4141413E007F08,
    0x3E00060909097F00,0x09097F005E215141,0x3249494926006619,0x3F0001017F010100,
    0x40201F003F404040,0x3F403C403F001F20,0x0700631408146300,0x4549710007087008,
    0x0041417F00000043,0x0000201008040200,0x01020400007F4141,0x8080808080800402,
    0x2000000007030000,0x44447F0078545454,0x2844444438003844,0x38007F4444443800,
    0x097E080008545454,0x7CA4A4A418000009,0x0000007804047F00,0x8480400000407D00,
    0x004428107F00007D,0x7C0000407F000000,0x04047C0078041804,0x3844444438000078,
    0x380038444444FC00,0x44784400FC444444,0x2054545408000804,0x3C000024443E0400,
    0x40201C00007C2040,0x3C6030603C001C20,0x9C00006C10106C00,0x54546400003C60A0,
    0x0041413E0800004C,0x0000000077000000,0x02010200083E4141,0x3C2623263C000001,
    0x3D001221E1A11E00,0x54543800007D2040,0x7855555520000955,0x2000785554552000,
    0x5557200078545555,0x1422E2A21C007857,0x3800085555553800,0x5555380008555455,
    0x00417C0100000854,0x0000004279020000,0x2429700000407C01,0x782F252F78007029,
    0x3400455554547C00,0x7F097E0058547C54,0x0039454538004949,0x3900003944453800,
    0x21413C0000384445,0x007C20413D00007D,0x3D00003D60A19C00,0x40413C00003D4242,
    0x002466241800003D,0x29006249493E4800,0x16097F00292A7C2A,0x02097E8840001078,
    0x0000785555542000,0x4544380000417D00,0x007D21403C000039,0x7A0000710A097A00,
    0x5555080000792211,0x004E51514E005E55,0x3C0020404D483000,0x0404040404040404,
    0x506A4C0817001C04,0x0000782A34081700,0x0014080000307D30,0x0814000814001408,
    0x55AA114411441144,0xEEBBEEBB55AA55AA,0x0000FF000000EEBB,0x0A0A0000FF080808,
    0xFF00FF080000FF0A,0x0000F808F8080000,0xFB0A0000FE0A0A0A,0xFF00FF000000FF00,
    0x0000FE02FA0A0000,0x0F0800000F080B0A,0x0F0A0A0A00000F08,0x0000F80808080000,
    0x080808080F000000,0xF808080808080F08,0x0808FF0000000808,0x0808080808080808,
    0xFF0000000808FF08,0x0808FF00FF000A0A,0xFE000A0A0B080F00,0x0B080B0A0A0AFA02,
    0x0A0AFA02FA0A0A0A,0x0A0A0A0AFB00FF00,0xFB00FB0A0A0A0A0A,0x0A0A0B0A0A0A0A0A,
    0x0A0A08080F080F08,0xF808F8080A0AFA0A,0x08080F080F000808,0x00000A0A0F000000,
    0xF808F8000A0AFE00,0x0808FF00FF080808,0x08080A0AFB0A0A0A,0xF800000000000F08,
    0xFFFFFFFFFFFF0808,0xFFFFF0F0F0F0F0F0,0xFF000000000000FF,0x0F0F0F0F0F0FFFFF,
    0xFE00241824241800,0x01017F0000344A4A,0x027E027E02000003,0x1800006349556300,
    0x2020FC00041C2424,0x000478040800001C,0x3E00085577550800,0x02724C00003E4949,
    0x0030595522004C72,0x1800182418241800,0x2A2A1C0018247E24,0x003C02023C00002A,
    0x0000002A2A2A2A00,0x4A4A510000242E24,0x00514A4A44000044,0x20000402FC000000,
    0x2A08080000003F40,0x0012241224000808,0x0000000609090600,0x0008000000001818,
    0x02023E4030000000,0x0900000E010E0100,0x3C3C3C0000000A0D,0x000000000000003C,
};

void oct_drawpix (int xres, int yres, int x, int y, int col)
{
    if constexpr (oct_usegpu) {
        // Set shader mode for 2d
        if (!gshadermode) return;
        gshadermode = 0;

        if constexpr (oct_useglsl)
        {
            ((PFNGLUSEPROGRAMPROC)glfp[glUseProgram])(0);
        }
        else
        {
            ((PFNGLDISABLE)glfp[glDisable])(GL_VERTEX_PROGRAM_ARB);
            ((PFNGLDISABLE)glfp[glDisable])(GL_FRAGMENT_PROGRAM_ARB);
        }
        kglActiveTexture(0);

        ((PFNGLMATRIXMODE)glfp[glMatrixMode])(GL_PROJECTION); ((PFNGLLOADIDENTITY)glfp[glLoadIdentity])(); ((PFNGLORTHO)glfp[glOrtho])(0,xres,yres,0,-1,1);
        ((PFNGLMATRIXMODE)glfp[glMatrixMode])(GL_MODELVIEW); ((PFNGLLOADIDENTITY)glfp[glLoadIdentity])();

        ((PFNGLDISABLE)glfp[glDisable])(GL_DEPTH_TEST);
        ((PFNGLBLENDFUNC)glfp[glBlendFunc])(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        ((PFNGLENABLE)glfp[glEnable])(GL_TEXTURE_2D);

        // Render
        ((PFNGLDISABLE)glfp[glDisable])(GL_TEXTURE_2D);
        ((PFNGLCOLOR3UB)glfp[glColor3ub])((col>>16)&255,(col>>8)&255,col&255);
        ((PFNGLBEGIN)glfp[glBegin])(GL_POINTS); ((PFNGLTEXCOORD2F)glfp[glTexCoord2f])(0.5f,219.f/256.f);/*middle of solid white character*/
        ((PFNGLVERTEX2F)glfp[glVertex2f])((float)x,(float)y);
        ((PFNGLEND)glfp[glEnd])();
    }
    else {
        if (((unsigned)x >= (unsigned)gdd.x) || ((unsigned)y >= (unsigned)gdd.y)) {
            return;
        }
        *(int *)(gdd.p*y + (x<<2) + gdd.f) = col;
    }
}

void oct_drawline (int xres, int yres, float x0, float y0, float x1, float y1, int col)
{
    if constexpr (oct_usegpu)
    {
        // Set shader mode for 2d
        if (!gshadermode) return;
        gshadermode = 0;

        if constexpr (oct_useglsl)
        {
            ((PFNGLUSEPROGRAMPROC)glfp[glUseProgram])(0);
        }
        else
        {
            ((PFNGLDISABLE)glfp[glDisable])(GL_VERTEX_PROGRAM_ARB);
            ((PFNGLDISABLE)glfp[glDisable])(GL_FRAGMENT_PROGRAM_ARB);
        }
        kglActiveTexture(0);

        ((PFNGLMATRIXMODE)glfp[glMatrixMode])(GL_PROJECTION); ((PFNGLLOADIDENTITY)glfp[glLoadIdentity])(); ((PFNGLORTHO)glfp[glOrtho])(0,xres,yres,0,-1,1);
        ((PFNGLMATRIXMODE)glfp[glMatrixMode])(GL_MODELVIEW); ((PFNGLLOADIDENTITY)glfp[glLoadIdentity])();

        ((PFNGLDISABLE)glfp[glDisable])(GL_DEPTH_TEST);
        ((PFNGLBLENDFUNC)glfp[glBlendFunc])(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        ((PFNGLENABLE)glfp[glEnable])(GL_TEXTURE_2D);

        // Render
        ((PFNGLDISABLE)glfp[glDisable])(GL_TEXTURE_2D);
        ((PFNGLCOLOR4UB)glfp[glColor4ub])((col>>16)&255,(col>>8)&255,col&255,(col>>24)&255);
        ((PFNGLBEGIN)glfp[glBegin])(GL_LINES); ((PFNGLTEXCOORD2F)glfp[glTexCoord2f])(0.5f,219.f/256.f);/*middle of solid white character*/
        ((PFNGLVERTEX2F)glfp[glVertex2f])(x0,y0);
        ((PFNGLVERTEX2F)glfp[glVertex2f])(x1,y1);
        ((PFNGLEND)glfp[glEnd])();
    }
    else
    {
        float f;
        int i, ipx[2], ipy[2];

        if (x0 < 0.f) {
            if (x1 < 0.f) {
                return;
            }
            y0 = (0.f - x0) * (y1 - y0) / (x1 - x0) + y0;
            x0 = 0.f;
        }
        else if (x0 > (float)gdd.x) {
            if (x1 > (float)gdd.x) {
                return;
            }
            y0 = ((float)gdd.x - x0) * (y1 - y0) / (x1 - x0) + y0;
            x0 = (float)gdd.x;
        }
        if (y0 < 0.f) {
            if (y1 < 0.f)
                return;
            x0 = (0.f - y0) * (x1 - x0) / (y1 - y0) + x0;
            y0 = 0.f;
        }
        else if (y0 > (float)gdd.y) {
            if (y1 > (float)gdd.y) {
                return;
            }
            x0 = ((float)gdd.y - y0) * (x1 - x0) / (y1 - y0) + x0;
            y0 = (float)gdd.y;
        }
        if (x1 < 0.f) {
            y1 = (0.f - x1) * (y1 - y0) / (x1 - x0) + y1;
            x1 = 0.f;
        }
        else if (x1 > (float)gdd.x) {
            y1 = ((float)gdd.x - x1) * (y1 - y0) / (x1 - x0) + y1;
            x1 = (float)gdd.x;
        }
        if (y1 < 0.f) {
            x1 = (0.f - y1) * (x1 - x0) / (y1 - y0) + x1;
            y1 = 0.f;
        }
        else if (y1 > (float)gdd.y) {
            x1 = ((float)gdd.y - y1) * (x1 - x0) / (y1 - y0) + x1;
            y1 = (float)gdd.y;
        }

        x1 -= x0;
        y1 -= y0;
        i = cvttss2si(ceil(std::max(fabs(x1), fabs(y1))));
        if (!(i & 0x7fffffff)) {
            return;
        }
        f = 65536.0 / (float)i;
        ipx[0] = (int)(x0 * 65536.0);
        ipx[1] = (int)(x1 * f);
        ipy[0] = (int)(y0 * 65536.0);
        ipy[1] = (int)(y1 * f);
        for (; i > 0; i--) {
            oct_drawpix(xres,yres,ipx[0] >> 16, ipy[0] >> 16, col);
            ipx[0] += ipx[1];
            ipy[0] += ipy[1];
        }
    }
}

char vshad_drawtext[] = R"(
    varying vec4 t;
    void main()
    {
       gl_Position = ftransform();
       t = gl_MultiTexCoord0;
    }
)";

char fshad_drawtext[] = R"(
    varying vec4 t;
    uniform sampler1D tex0;
    uniform sampler2D tex1;
    uniform vec2 rtxtmul0, rtxtmul1, rtxtmul2, txtmod;
    uniform vec4 fcol, bcol;
    void main()
    {
       float ch = texture1D(tex0, dot(floor(t.xy * rtxtmul0), rtxtmul1)).x;
       ch *= (255.001/256.0);
       vec4 inchar = texture2D(tex1, mod(t.xy, txtmod) * rtxtmul2 + vec2(0.0,ch));
       gl_FragColor = mix(bcol, fcol, inchar.r);
    }
)";

char vshadasm_drawtext[] = R"(!!ARBvp1.0
    PARAM ModelViewProj[4] = {state.matrix.mvp};
    TEMP t;
    DP4 t.x, ModelViewProj[0], vertex.position;
    DP4 t.y, ModelViewProj[1], vertex.position;
    DP4 t.z, ModelViewProj[2], vertex.position;
    DP4 t.w, ModelViewProj[3], vertex.position;
    MOV result.position, t;
    MOV result.texcoord[0], vertex.texcoord[0];
    END
)";

char fshadasm_drawtext[] = R"(!!ARBfp1.0
    ATTRIB t = fragment.texcoord[0];
    PARAM rtxtmul0 = program.local[0];
    PARAM rtxtmul1 = program.local[1];
    PARAM rtxtmul2 = program.local[2];
    PARAM rtxtmod  = program.local[3];
    PARAM fcol     = program.local[4];
    PARAM bcol     = program.local[5];
    TEMP a, b;
    MUL a, t, rtxtmul0;
    FLR a, a;
    DP3 a, a, rtxtmul1;
    TEX a, a, texture[0], 1D;
    MUL a, a, {0.0,0.996,0.0,0.0}; #~255.001/256
    MUL b, t, rtxtmod;
    FRC b, b;
    MAD b, b, rtxtmul2, a;
    TEX b, b, texture[1], 2D;
    LRP result.color, b.r, fcol, bcol;
    END
)";

//FIX:dynamic!

//* Caller must update (0..dx-1, 0..dy-1) of GPU tex#0 before calling oct_drawtext6x8()
//* Full alpha is supported on GPU.

void oct_drawtext6x8 (int xres, int yres, int x0, int y0, int fcol, int bcol, const char *fmt, ...)
{
    va_list arglist;
    unsigned char st[1024];

    if (!fmt) return;
    va_start(arglist,fmt);
    if (_vsnprintf((char *)&st,sizeof(st)-1,fmt,arglist)) st[sizeof(st)-1] = 0;
    va_end(arglist);

    if constexpr (oct_usegpu) {
        #define FXSIZ 6
        #define FYSIZ 8
        unsigned char *cptr, txtbuf[TXTBUFSIZ];
        int i, x, y, xmax, ymax;

            //Calculate rectangle area
        x = 0; y = 0; xmax = -1; ymax = -1;
        for(cptr=st;cptr[0];x++)
        {
            i = *cptr++;
            if (i == 32) continue;
            if (i == 9) { x += 2; continue; }
            if (i == '\n') { x = -1; y++; continue; }
            if (x > xmax) xmax = x;
        }
        ymax = y; if ((xmax|ymax) < 0) return;
        xmax++; ymax++; if (xmax*ymax >= TXTBUFSIZ) return; //:/

            //Render string
        x = 0; y = 0; memset(txtbuf,32,xmax*ymax);
        for(cptr=st;cptr[0];x++)
        {
            i = *cptr++;
            if (i == 32) continue;
            if (i == 9) { x += 2; continue; }
            if (i == '\n') { x = -1; y++; continue; }
            txtbuf[y*xmax+x] = (unsigned char)i;
        }

        if (shadcur != 3)
        {
            shadcur = 3;
            ((PFNGLMATRIXMODE)glfp[glMatrixMode])(GL_PROJECTION); ((PFNGLLOADIDENTITY)glfp[glLoadIdentity])(); ((PFNGLORTHO)glfp[glOrtho])(0,xres,yres,0,-1,1);
            ((PFNGLBLENDFUNC)glfp[glBlendFunc])(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

            kglActiveTexture(1); ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_2D,gfont6x8id);
            kglActiveTexture(0); ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_1D,gtextbufid);

            if constexpr (oct_useglsl)
            {
                ((PFNGLUSEPROGRAMPROC)glfp[glUseProgram])(shadprog[shadcur]);
            }
            else
            {
                ((PFNGLENABLE)glfp[glEnable])(GL_VERTEX_PROGRAM_ARB);
                ((PFNGLENABLE)glfp[glEnable])(GL_FRAGMENT_PROGRAM_ARB);
                ((PFNGLBINDPROGRAMARBPROC)glfp[glBindProgramARB])(GL_VERTEX_PROGRAM_ARB  ,shadvert[shadcur]);
                ((PFNGLBINDPROGRAMARBPROC)glfp[glBindProgramARB])(GL_FRAGMENT_PROGRAM_ARB,shadfrag[shadcur]);
            }
        }

        ((PFNGLTEXSUBIMAGE1D)glfp[glTexSubImage1D])(GL_TEXTURE_1D,0,0,xmax*ymax,GL_LUMINANCE,GL_UNSIGNED_BYTE,(void *)txtbuf);

        if constexpr (oct_useglsl)
        {
            kglUniform2f(kglGetUniformLoc(shadprog[shadcur], "rtxtmul0"),1.f/(float)FXSIZ,1.f/(float)FYSIZ);
            kglUniform2f(kglGetUniformLoc(shadprog[shadcur], "rtxtmul1"),1.f/(float)TXTBUFSIZ,(float)xmax/(float)TXTBUFSIZ);
            kglUniform2f(kglGetUniformLoc(shadprog[shadcur], "rtxtmul2"),1.f/(float)8,1.f/(float)(FYSIZ*256));
            kglUniform2f(kglGetUniformLoc(shadprog[shadcur], "txtmod"),FXSIZ,FYSIZ);
            kglUniform4f(kglGetUniformLoc(shadprog[shadcur], "fcol"),(float)((fcol>>16)&255)*(1.f/255.f),(float)((fcol>> 8)&255)*(1.f/255.f),
                                                              (float)((fcol>> 0)&255)*(1.f/255.f),(float)((fcol>>24)&255)*(1.f/255.f));
            kglUniform4f(kglGetUniformLoc(shadprog[shadcur], "bcol"),(float)((bcol>>16)&255)*(1.f/255.f),(float)((bcol>> 8)&255)*(1.f/255.f),
                                                              (float)((bcol>> 0)&255)*(1.f/255.f),(float)((bcol>>24)&255)*(1.f/255.f));
        }
        else
        {
            kglProgramLocalParam(0,1.f/(float)FXSIZ,1.f/(float)FYSIZ,0,0);                                   //rtxtmul0
            kglProgramLocalParam(1,1.f/(float)TXTBUFSIZ,(float)xmax/(float)TXTBUFSIZ,0,0);                   //rtxtmul1
            kglProgramLocalParam(2,FXSIZ/(float)8,1.0/256.0,0,0);                                            //rtxtmul2*txtmod
            kglProgramLocalParam(3,1.f/(float)FXSIZ,1.0/(float)FYSIZ,0,0);                                   //rtxtmod
            kglProgramLocalParam(4,(float)((fcol>>16)&255)*(1.f/255.f),(float)((fcol>> 8)&255)*(1.f/255.f),  //fcol
                                          (float)((fcol>> 0)&255)*(1.f/255.f),(float)((fcol>>24)&255)*(1.f/255.f));
            kglProgramLocalParam(5,(float)((bcol>>16)&255)*(1.f/255.f),(float)((bcol>> 8)&255)*(1.f/255.f),  //bcol
                                          (float)((bcol>> 0)&255)*(1.f/255.f),(float)((bcol>>24)&255)*(1.f/255.f));
        }

        xmax *= FXSIZ; ymax *= FYSIZ;
        ((PFNGLBEGIN)glfp[glBegin])(GL_QUADS);
        ((PFNGLTEXCOORD2F)glfp[glTexCoord2f])(        0.f,        0.f); ((PFNGLVERTEX2I)glfp[glVertex2i])(x0     ,y0     );
        ((PFNGLTEXCOORD2F)glfp[glTexCoord2f])((float)xmax,        0.f); ((PFNGLVERTEX2I)glfp[glVertex2i])(x0+xmax,y0     );
        ((PFNGLTEXCOORD2F)glfp[glTexCoord2f])((float)xmax,(float)ymax); ((PFNGLVERTEX2I)glfp[glVertex2i])(x0+xmax,y0+ymax);
        ((PFNGLTEXCOORD2F)glfp[glTexCoord2f])(        0.f,(float)ymax); ((PFNGLVERTEX2I)glfp[glVertex2i])(x0     ,y0+ymax);
        ((PFNGLEND)glfp[glEnd])();

        shadcur = 0;
        if (oct_useglsl) ((PFNGLUSEPROGRAMPROC)glfp[glUseProgram])(shadprog[shadcur]);
    }
    else {
        unsigned char *pst, *c, *v;
        int i, j, ie, x, *lp, *lpx;

        pst = st;
        do
        {
            lp = (int *)(y0*gdd.p+gdd.f);
            for(j=1;j<256;y0++,lp=(int *)(((INT_PTR)lp)+gdd.p),j+=j)
                if ((unsigned)y0 < (unsigned)gdd.y)
                    for(c=pst,x=x0;c[0] && (c[0] != '\n');c++,x+=6)
                    {
                        v = ((int)(*c))*6 + ((unsigned char *)font6x8); lpx = &lp[x];
                        for(i=std::max(-x,0),ie=std::min(gdd.x-x,6);i<ie;i++) { if (v[i]&j) lpx[i] = fcol; else if (bcol < 0) lpx[i] = bcol; }
                        if ((*c) == 9) { if (bcol < 0) { for(i=std::max(-x,6),ie=std::min(gdd.x-x,18);i<ie;i++) lpx[i] = bcol; } x += 2*6; }
                    }
            if (!c[0]) return;
            pst = &c[1];
        } while (1);
    }
}

int oct_startdraw(window_state& window, point3d *ipos, point3d *irig, point3d *idow, point3d *ifor, double ghx, double ghy, double ghz)
{
    tiletype dd;
    INT_PTR i, j, k;

    if constexpr (oct_usegpu) {
        ((PFNGLMATRIXMODE)glfp[glMatrixMode])(GL_PROJECTION);
        ((PFNGLLOADIDENTITY)glfp[glLoadIdentity])();
        ((PFNGLFRUSTUM)glfp[glFrustum])(-gznear, gznear, -(float)window.yres / (float)window.xres * gznear, (float)window.yres / (float)window.xres * gznear, gznear, gzfar);
        ((PFNGLMATRIXMODE)glfp[glMatrixMode])(GL_MODELVIEW);
        ((PFNGLLOADIDENTITY)glfp[glLoadIdentity])();

        k = oct_fogcol;
        ((PFNGLCLEARCOLOR)glfp[glClearColor])(((k >> 16) & 255) / 255.0, ((k >> 8) & 255) / 255.0, (k & 255) / 255.0, 0.0);
        ((PFNGLCLEAR)glfp[glClear])(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
        ((PFNGLENABLE)glfp[glEnable])(GL_DEPTH_TEST);
        ((PFNGLENABLE)glfp[glEnable])(GL_BLEND);
        ((PFNGLBLENDFUNC)glfp[glBlendFunc])(GL_ONE, GL_NONE);
        gshadermode = 1;
    }
    else {
        if (!window.startdirectdraw(&dd.f,&dd.p,&dd.x,&dd.y)) return(-1);

        i = dd.p*dd.y;
        if (i > zbufmal) {
            zbufmal = i;
            zbuf = (int *)realloc(zbuf,zbufmal+256);
        }
    }

    i = ((window.xres + 127) & ~127);

    if ((i >> 3) * window.yres > gcbitmax) {
        gcbitmax = (i >> 3) * window.yres;
        gcbitpl = i;
        gcbitplo5 = (gcbitpl >> 5);
        gcbitplo3 = (gcbitpl >> 3);
        gcbpl = (gcbitpl >> 3);
        gcbitmal = (INT_PTR)realloc((void*)gcbitmal, gcbitmax + 256);
        gcbit = (unsigned int*)((gcbitmal + 15) & ~15);
    }

    if constexpr (!oct_usegpu) {
        //zbuffer aligns its memory to the same pixel boundaries as the screen!
        //WARNING: Pentium 4's L2 cache has severe slowdowns when 65536-64 <= (zbufoff&65535) < 64
        zbufoff = (((((INT_PTR)zbuf)-dd.f-128)+255)&~255)+128;
        for(i=0,j=dd.f+zbufoff;i<dd.y;i++,j+=dd.p) {
            memset16_safe((void *)j,0x7f7f7f7f,dd.x<<2); //Clear z-buffer
        }
    }

    oct_setcam(window.xres,window.yres,&dd,zbufoff,ipos,irig,idow,ifor,ghx,ghy,ghz);

    if ((oct_fogdist >= 1e32) && (gskyid)) {
        drawsky(window.xres,window.yres);
    }
    else if constexpr (!oct_usegpu) {
        drawrectfill(&dd,0,0,dd.x,dd.y,oct_fogcol);
    }

    if constexpr (oct_usegpu) {
        ((PFNGLBLENDFUNC)glfp[glBlendFunc])(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    }

    return 0;
}

int oct_startdraw(window_state& window, point3f *ipos, point3f *irig, point3f *idow, point3f *ifor, float ghx, float ghy, float ghz)
{
    point3d dpos, drig, ddow, dfor;
    dpos.x = ipos->x; dpos.y = ipos->y; dpos.z = ipos->z;
    drig.x = irig->x; drig.y = irig->y; drig.z = irig->z;
    ddow.x = idow->x; ddow.y = idow->y; ddow.z = idow->z;
    dfor.x = ifor->x; dfor.y = ifor->y; dfor.z = ifor->z;
    return(oct_startdraw(window, &dpos,&drig,&ddow,&dfor,ghx,ghy,ghz));
}

void oct_stopdraw(window_state& window) {
    if constexpr(oct_usegpu) {
        SwapBuffers(glhDC);
    }
    else {
        window.stopdirectdraw();
        window.nextpage();
    }
}
