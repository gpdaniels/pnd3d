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
    free((void *)((tiles_t *)ptr)->f);
    free((void *)ptr);
}

INT_PTR oct_loadtex (INT_PTR ptr, int xsiz, int ysiz, int flags)
{
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

static void drawsky(int xres, int yres)
{
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
        gxd.f = 0; gxd.p = 0; gxd.x = 0; gxd.y = 0;
    }
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

    gdd = (*dd);
    gzbufoff = zbufoff;
    gddz = gdd.f+zbufoff;

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

    memset16((void *)gcbit,-1,gcbpl*loct->yres); //Clear cover buffer

    gxmin = 0; gymin = 0; gxmax = loct->xres; gymax = loct->yres;
    memset4(glmost,         0,loct->yres*4);
    memset4(grmost,loct->xres,loct->yres*4);

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

    if (gdrawclose) oct_rendclosecubes(loct,&glorig);
    if (gdrawclose < 2)
    {
        htrun(oct_ftob2_dorect,loct,1,rectn,true);
    }

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

    ls = loct->lsid-1; ptr = &((octv_t *)loct->nod.buf)[loct->head]; j = 8-1;
    while (1)
    {
        i = (1<<j); if (!(ptr->chi&i)) goto tosibly;

        if (ls <= 0)
        {
            ind = popcount8(ptr->chi&(i-1)) + ptr->ind; psurf = &((surf_t *)loct->sur.buf)[ind];
            oct_surf_normalize(psurf);
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
    if (((unsigned)x >= (unsigned)gdd.x) || ((unsigned)y >= (unsigned)gdd.y)) {
        return;
    }
    *(int *)(gdd.p*y + (x<<2) + gdd.f) = col;
}

void oct_drawline (int xres, int yres, float x0, float y0, float x1, float y1, int col)
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

int oct_startdraw(window_state& window, point3d *ipos, point3d *irig, point3d *idow, point3d *ifor, double ghx, double ghy, double ghz)
{
    tiletype dd;
    INT_PTR i, j, k;

    if (!window.startdirectdraw(&dd.f,&dd.p,&dd.x,&dd.y)) return(-1);

    i = dd.p*dd.y;
    if (i > zbufmal) {
        zbufmal = i;
        zbuf = (int *)realloc(zbuf,zbufmal+256);
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

    //zbuffer aligns its memory to the same pixel boundaries as the screen!
    //WARNING: Pentium 4's L2 cache has severe slowdowns when 65536-64 <= (zbufoff&65535) < 64
    zbufoff = (((((INT_PTR)zbuf)-dd.f-128)+255)&~255)+128;
    for(i=0,j=dd.f+zbufoff;i<dd.y;i++,j+=dd.p) {
        memset16_safe((void *)j,0x7f7f7f7f,dd.x<<2); //Clear z-buffer
    }

    oct_setcam(window.xres,window.yres,&dd,zbufoff,ipos,irig,idow,ifor,ghx,ghy,ghz);

    if ((oct_fogdist >= 1e32) && (gskyid)) {
        drawsky(window.xres,window.yres);
    }

    drawrectfill(&dd,0,0,dd.x,dd.y,oct_fogcol);

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
    window.stopdirectdraw();
    window.nextpage();
}
