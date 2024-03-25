// This file has been modified from Ken Silverman's original release

#include "octree.hpp"

// TODO: REMOVE!
#include "octree_renderer.hpp"
#include "octree_modify.hpp"
#include "octree_physics.hpp"
// TODO: REMOVE!

#include "utilities/bits.hpp"
#include "utilities/gl.hpp"
#include "utilities/popcount.hpp"
#include "utilities/window.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <emmintrin.h>

#ifdef __GNUC__
#define stricmp strcasecmp
#endif

#define cvttss2si(f) _mm_cvtt_ss2si(_mm_set_ss(f))

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int oct_initonce (HWND ghwnd, float zfar)
{
    int x, y, xsiz, ysiz;
    char *cbuf, *cptr;

    gzfar = zfar;
    gznear = zfar * (1.0/8388608.0);
    shaderfunc = cpu_shader_texmap;

    return(0);
}

void oct_uninitonce(HWND ghwnd)
{
}

void oct_new (int xres, int yres, oct_t *loct, int los, INT_PTR tilid, int startnodes, int startsurfs, int hax4mark2spr)
{
    int i;
    octv_t *ptr;

    loct->head = 0;
    loct->lsid = los; loct->sid = (1<<los); loct->nsid = -loct->sid;
    loct->flags = 0;
    loct->edgeiswrap = 0;
    loct->edgeissol = 0;

    loct->recvoctfunc = 0;

    loct->xres = xres;
    loct->yres = yres;

    // octvn,64rndsph: model:
    // 1:      9       8*4^1  32
    // 2:     65       8*4^2  128
    // 3:    361       8*4^3  512
    // 4:    172       8*4^4  2048
    // 5:   3434       8*4^5  8192
    // 6:  20131       8*4^6  32768
    // 7: 103732       8*4^7  131072
    // 8: 474610       8*4^8  524288
    // 9:1841536       8*4^9  2097152
    //10:7417321       8*4^10 8388608

    if (startnodes)
        loct->nod.mal = startnodes;
    else if (los <= 10)
        loct->nod.mal = (8<<(los<<1));
    else
        loct->nod.mal = 16777216;

    loct->nod.siz = sizeof(octv_t);
    loct->nod.buf = malloc(loct->nod.mal*loct->nod.siz);
    ptr = (octv_t *)loct->nod.buf;
    ptr[0].chi = 0;
    ptr[0].sol = 0;
    ptr[0].mrk = 0;
    ptr[0].mrk2= 0;
    ptr[0].ind =-1;
    loct->nod.num = 1; loct->nod.ind = 1;
    if (!hax4mark2spr)
    {
        i = (((loct->nod.mal+63)>>5)<<2)+16;
        loct->nod.bit = (unsigned int *)malloc(i);
        memset4(loct->nod.bit,0,i); loct->nod.bit[0] = 1;
    }

    if (startsurfs) loct->sur.mal = startsurfs; else if (los <= 10) loct->sur.mal = (8<<(los<<1)); else loct->sur.mal = 16777216;
    loct->sur.siz = sizeof(surf_t);
    loct->sur.buf = (octv_t *)malloc(loct->sur.mal*loct->sur.siz);
    loct->sur.num = 1; loct->sur.ind = 1;
    if (!hax4mark2spr)
    {
        i = (((loct->sur.mal+63)>>5)<<2)+16;
        loct->sur.bit = (unsigned int *)malloc(i);
        memset4(loct->sur.bit,0,i);
        loct->sur.bit[0] = 1; //Can't use index 0 since it is used as background color in GLSL shaders.
    }
}

#if defined(ENABLE_LOADING_ROUTINES)
int oct_load (oct_t *loct, char *filnam, point3f *rpos, point3f *rrig, point3f *rdow, point3f *rfor)
{
    static_cast<void>(rpos);
    static_cast<void>(rrig);
    static_cast<void>(rdow);
    static_cast<void>(rfor);

    int i, ret;

    i = strlen(filnam);
    if (i >= 4)
    {
             if (!stricmp(&filnam[i-4],".kvo")) ret = loadkvo        (loct,filnam,rpos,rrig,rdow,rfor);
        else if (!stricmp(&filnam[i-4],".kvs")) ret = loadkvo        (loct,filnam,rpos,rrig,rdow,rfor);
        else if (!stricmp(&filnam[i-4],".kvn")) ret = loadkvo        (loct,filnam,rpos,rrig,rdow,rfor);
        else if (!stricmp(&filnam[i-4],".kv6")) ret = loadkv6_kvx_vox(loct,filnam,rpos,rrig,rdow,rfor,0);
        else if (!stricmp(&filnam[i-4],".kvx")) ret = loadkv6_kvx_vox(loct,filnam,rpos,rrig,rdow,rfor,1);
        else if (!stricmp(&filnam[i-4],".vox")) ret = loadkv6_kvx_vox(loct,filnam,rpos,rrig,rdow,rfor,2);
        else if (!stricmp(&filnam[i-4],".vxl")) ret = loadvxl        (loct,filnam,rpos,rrig,rdow,rfor);
        else if (!stricmp(&filnam[i-4],".png")) ret = loadpng        (loct,filnam,rpos,rrig,rdow,rfor);
        else ret = -1;
    } else ret = -1;

    if (ret < 0)
    {
        char tbuf[1024];
        fprintf(stderr,"Unable to load file: %s\n",filnam);
        oct_new(loct,6,0,0,0,0); //generate bogus model
    }

    if constexpr (oct_usegpu) {
        if constexpr (oct_usegpubo) {
            if (!loct->gsurf) loct->gsurf = (surf_t *)bo_begin(loct->bufid,0);
            memcpy(loct->gsurf,loct->sur.buf,loct->sur.mal*loct->sur.siz);
        }
        else {
            ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_2D,loct->octid);
            ((PFNGLTEXSUBIMAGE2D)glfp[glTexSubImage2D])(GL_TEXTURE_2D,0,0,0,loct->gxsid,(loct->sur.mal*2)>>loct->glxsid,GL_RGBA,GL_UNSIGNED_BYTE,(void *)loct->sur.buf);
        }
    }

    return(0);
}

int oct_load (oct_t *loct, char *filnam, point3d *rpos, point3d *rrig, point3d *rdow, point3d *rfor)
{
    int i;
    point3f fp, fr, fd, ff;
    i = oct_load(loct,filnam,&fp,&fr,&fd,&ff);
    rpos->x = fp.x; rpos->y = fp.y; rpos->z = fp.z;
    rrig->x = fr.x; rrig->y = fr.y; rrig->z = fr.z;
    rdow->x = fd.x; rdow->y = fd.y; rdow->z = fd.z;
    rfor->x = ff.x; rfor->y = ff.y; rfor->z = ff.z;
    return(i);
}
#endif

#if defined(ENABLE_SAVING_ROUTINES)
void oct_save (oct_t *loct, char *filnam, point3f *rpos, point3f *rrig, point3f *rdow, point3f *rfor)
{
    int i;

    i = strlen(filnam);
    if ((i >= 4) && (!stricmp(&filnam[i-4],".kvo"))) { savekvo(loct,filnam,rpos,rrig,rdow,rfor); return; }
    if ((i >= 4) && (!stricmp(&filnam[i-4],".kvs"))) { savekvo(loct,filnam,rpos,rrig,rdow,rfor); return; }
    if ((i >= 4) && (!stricmp(&filnam[i-4],".kvn"))) { savekvo(loct,filnam,rpos,rrig,rdow,rfor); return; }
    if ((i >= 4) && (!stricmp(&filnam[i-4],".kv6"))) { savekv6(loct,filnam,rpos,rrig,rdow,rfor); return; }
}

void oct_save (oct_t *loct, char *filnam, point3d *rpos, point3d *rrig, point3d *rdow, point3d *rfor)
{
    point3f fp, fr, fd, ff;
    fp.x = rpos->x; fp.y = rpos->y; fp.z = rpos->z;
    fr.x = rrig->x; fr.y = rrig->y; fr.z = rrig->z;
    fd.x = rdow->x; fd.y = rdow->y; fd.z = rdow->z;
    ff.x = rfor->x; ff.y = rfor->y; ff.z = rfor->z;
    oct_save(loct,filnam,&fp,&fr,&fd,&ff);
}
#endif

#if defined(ENABLE_UNUSED_ROUTINES)
void oct_dup (oct_t *ooct, oct_t *noct)
{
    int i;

    memcpy(noct,ooct,sizeof(oct_t));

    i = noct->nod.mal*noct->nod.siz;
    noct->nod.buf = (octv_t *)malloc(i);
    memcpy(noct->nod.buf,ooct->nod.buf,i);

    i = (((noct->nod.mal+63)>>5)<<2)+16;
    noct->nod.bit = (unsigned int *)malloc(i);
    memcpy(noct->nod.bit,ooct->nod.bit,i);


    i = noct->sur.mal*noct->sur.siz;
    noct->sur.buf = (octv_t *)malloc(i);
    memcpy(noct->sur.buf,ooct->sur.buf,i);

    i = (((noct->sur.mal+63)>>5)<<2)+16;
    noct->sur.bit = (unsigned int *)malloc(i);
    memcpy(noct->sur.bit,ooct->sur.bit,i);

    if constexpr (oct_usegpu) {
        ((PFNGLGENTEXTURES)glfp[glGenTextures])(1,&noct->octid);
        kglalloctex(noct->octid,0,noct->gxsid,noct->gysid,1,KGL_RGBA32+KGL_NEAREST); //only NEAREST makes sense here!
        if constexpr (oct_usegpubo) {
            noct->bufid = bo_init(noct->gxsid*noct->gysid*4);
            noct->gsurf = (surf_t *)bo_begin(noct->bufid,0);
            memcpy(noct->gsurf,noct->sur.buf,noct->sur.mal*noct->sur.siz);
        }
        else {
            ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_2D,noct->octid);
            ((PFNGLTEXSUBIMAGE2D)glfp[glTexSubImage2D])(GL_TEXTURE_2D,0,0,0,noct->gxsid,(noct->sur.mal*2)>>noct->glxsid,GL_RGBA,GL_UNSIGNED_BYTE,(void *)noct->sur.buf);
        }
    }
}
#endif

void oct_free (oct_t *loct)
{
    if (loct->sur.bit) free(loct->sur.bit);
    if (loct->sur.buf) free(loct->sur.buf);
    if (loct->nod.bit) free(loct->nod.bit);
    if (loct->nod.buf) free(loct->nod.buf);
    memset(loct,0,sizeof(oct_t));
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// THE REST OF THE CODE IN THIS FILE IS DISABLED

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//EQUIVEC code begins -----------------------------------------------------
#if defined(ENABLE_EQUIVEC_ROUTINES)
#define GOLDRAT 0.3819660112501052 //Golden Ratio: 1 - 1/((sqrt(5)+1)/2)
typedef struct
{
    float fibx[45], fiby[45];
    float azval[20], zmulk, zaddk;
    int fib[47], aztop, npoints;
    point3f *p;  //For memory version :/
    int pcur;
} equivectyp;
static equivectyp equivec;

static void equiind2vec (int i, float *x, float *y, float *z)
{
    float a, r;
    (*z) = (float)i*equivec.zmulk + equivec.zaddk; r = std::sqrt(1.f - (*z)*(*z));
    a = ((float)i)*(GOLDRAT*M_PI*2); (*x) = std::cos(a)*r; (*y) = std::sin(a)*r;
}

    //(Pass n=0 to free buffer)
static int equimemset (int n)
{
    int z;

    if (equivec.pcur == n) return(1); //Don't recalculate if same #
    if (equivec.p) { free(equivec.p); equivec.p = 0; }
    if (!n) return(1);

        //Init vector positions (equivec.p) for memory version
    equivec.p = (point3f *)malloc(((n+1)&~1)*sizeof(point3f));
    if (!equivec.p) {
        return(0);
    }

    equivec.pcur = n;
    equivec.zmulk = 2 / (float)n; equivec.zaddk = equivec.zmulk*.5 - 1.0;
    for(z=n-1;z>=0;z--)
        equiind2vec(z,&equivec.p[z].x,&equivec.p[z].y,&equivec.p[z].z);
    if (n&1) //Hack for when n=255 and want a <0,0,0> vector
        { equivec.p[n].x = equivec.p[n].y = equivec.p[n].z = 0; }
    return(1);
}

    //Very fast; good quality, requires equivec.p[] :/
static int equivec2indmem (float x, float y, float z)
{
    int b, i, j, k, bestc;
    float xy, zz, md, d;

    xy = atan2(y,x); //atan2 is 150 clock cycles!
    j = ((*(int *)&z)&0x7fffffff);
    bestc = equivec.aztop;
    do
    {
        if (j < *(int *)&equivec.azval[bestc]) break;
        bestc--;
    } while (bestc);

    zz = z + 1.f;
    i = (int)(equivec.fibx[bestc]*xy + equivec.fiby[bestc]*zz - .5);
    bestc++;
    j = (int)(equivec.fibx[bestc]*xy + equivec.fiby[bestc]*zz - .5);

    k = equivec.fib[bestc+2]*i + equivec.fib[bestc+1]*j;
    if ((unsigned)k < (unsigned)equivec.npoints)
    {
        md = equivec.p[k].x*x + equivec.p[k].y*y + equivec.p[k].z*z;
        j = k;
    } else md = -2.f;
    b = bestc+3;
    do
    {
        i = equivec.fib[b] + k;
        if ((unsigned)i < (unsigned)equivec.npoints)
        {
            d = equivec.p[i].x*x + equivec.p[i].y*y + equivec.p[i].z*z;
            if (*(int *)&d > *(int *)&md) { md = d; j = i; }
        }
        b--;
    } while (b != bestc);
    return(j);
}

static void equivecinit (int n)
{
    float t0, t1;
    int z;

        //Init constants for ind2vec
    equivec.npoints = n;
    equivec.zmulk = 2 / (float)n; equivec.zaddk = equivec.zmulk*.5 - 1.0;

    equimemset(n);

        //Init Fibonacci table
    equivec.fib[0] = 0; equivec.fib[1] = 1;
    for(z=2;z<47;z++) equivec.fib[z] = equivec.fib[z-2]+equivec.fib[z-1];

        //Init fibx/y LUT
    t0 = (float)(.5 / M_PI); t1 = (float)n * -.5f;
    for(z=0;z<45;z++)
    {
        t0 = -t0; equivec.fibx[z] = (float)equivec.fib[z+2]*t0;
        t1 = -t1; equivec.fiby[z] = ((float)equivec.fib[z+2]*GOLDRAT - (float)equivec.fib[z])*t1;
    }

    t0 = 1 / ((float)n * M_PI);
    for(equivec.aztop=0;equivec.aztop<20;equivec.aztop++)
    {
        t1 = 1 - (float)equivec.fib[(equivec.aztop<<1)+6]*t0; if (t1 < 0) break;
        equivec.azval[equivec.aztop+1] = sqrt(t1);
    }
}

static void equivecuninit() {
    equimemset(0);
}
#endif
//EQUIVEC code ends -------------------------------------------------------

//--------------------------------------------------------------------------------------------------
#if defined(ENABLE_SAVING_ROUTINES)
static int palhashead[4096];
typedef struct {
    int rgb, n;
} paltab_t;
static paltab_t paltab[4096];
static int paltabn;
static void palreset() {
    memset(palhashead,-1,sizeof(palhashead));
    paltabn = 0;
}
static int palget (int rgb)
{
    int i, hash;
    hash = ((((rgb>>22)^(rgb>>11))-rgb)&(sizeof(palhashead)/sizeof(palhashead[0])-1));
    for(i=palhashead[hash];i>=0;i=paltab[i].n) {
        if (paltab[i].rgb == rgb) {
            return(i);
        }
    }
    if (paltabn >= sizeof(paltab)/sizeof(paltab[0])) {
        return(-1);
    }
    paltab[paltabn].rgb = rgb;
    paltab[paltabn].n = palhashead[hash];
    palhashead[hash] = paltabn;
    paltabn++;
    return(paltabn-1);
}
#endif
//--------------------------------------------------------------------------------------------------

#if defined(ENABLE_SAVING_ROUTINES)
//KV6 format:
//   int sig; //0x6c78764b (Kvxl)
//   int xsiz, ysiz, zsiz;
//   float xpiv, ypiv, zpiv;
//   int numvoxs;
//   kv6data { char b, g, r, a; ushort z; char vis, dir; } [numvoxs];
//           { Nnnnnnnn--VvvvvvZzzzzzzzzzzzzzzz--------RrrrrrrrGgggggggBbbbbbbb }
//   int xlen[xsiz];
//   ushort ylen[xsiz][ysiz];
void savekv6 (oct_t *loct, char *filnam, point3f *rpos, point3f *rrig, point3f *rdow, point3f *rfor)
{
    static_cast<void>(rpos);
    static_cast<void>(rrig);
    static_cast<void>(rdow);
    static_cast<void>(rfor);

    FILE *fil;
    surf_t surf, *psurf;
    float f;
    int i, x, y, z, numvoxs, *xlen, x0, y0, z0, x1, y1, z1, vis, oy0, oz0, oy1, oz1;
    short *ylen;

    x0 = 0; y0 = 0; z0 = 0; x1 = loct->sid; y1 = loct->sid; z1 = loct->sid;
    oct_getsolbbox(loct,&x0,&y0,&z0,&x1,&y1,&z1);
    oy0 = y0; oz0 = z0; oy1 = y1; oz1 = z1;
    z0 = oy0; y0 = loct->sid-oz1;
    z1 = oy1; y1 = loct->sid-oz0;

    equivecinit(255);

    xlen = (int *)malloc((x1-x0)<<2); if (!xlen) return;
    ylen = (short *)malloc(((x1-x0)*(y1-y0))<<1); if (!ylen) { free(xlen); return; }

    fil = fopen(filnam,"wb"); if (!fil) { free(ylen); free(xlen); return; }
    i = 0x6c78764b; fwrite(&i,4,1,fil); //Kvxl
    i = x1-x0; fwrite(&i,4,1,fil);
    i = y1-y0; fwrite(&i,4,1,fil);
    i = z1-z0; fwrite(&i,4,1,fil);
    //f = (x1-x0)*.5; fwrite(&f,4,1,fil); //bbox of solid centered
    //f = (y1-y0)*.5; fwrite(&f,4,1,fil);
    //f = (z1-z0)*.5; fwrite(&f,4,1,fil);
    f = (float)loct->sid*.5 - x0; fwrite(&f,4,1,fil); //pivot preserved
    f = (float)loct->sid*.5 - y0; fwrite(&f,4,1,fil);
    f = (float)loct->sid*.5 - z0; fwrite(&f,4,1,fil);

    numvoxs = 0;
    memset(xlen,0,(x1-x0)<<2);
    memset(ylen,0,((x1-x0)*(y1-y0))<<1);

    fwrite(&numvoxs,4,1,fil); //dummy write

    for(x=x0;x<x1;x++)
        for(y=y0;y<y1;y++)
            for(z=oct_findsurfdowny(loct,x,-1,loct->sid-1-y,&psurf);z<loct->sid;z=oct_findsurfdowny(loct,x,z,loct->sid-1-y,&psurf))
            {
                surf = (*psurf);

                vis = 0;
                if (!oct_getsol(loct,x-1,z,loct->sid-1-(y))) vis |= (1<<0);
                if (!oct_getsol(loct,x+1,z,loct->sid-1-(y))) vis |= (1<<1);
                if (!oct_getsol(loct,x,z,loct->sid-1-(y-1))) vis |= (1<<2);
                if (!oct_getsol(loct,x,z,loct->sid-1-(y+1))) vis |= (1<<3);
                if (!oct_getsol(loct,x,z-1,loct->sid-1-(y))) vis |= (1<<4);
                if (!oct_getsol(loct,x,z+1,loct->sid-1-(y))) vis |= (1<<5);

                fputc(surf.b    ,fil);
                fputc(surf.g    ,fil);
                fputc(surf.r    ,fil);
                fputc(128       ,fil);
                fputc((z-z0)&255,fil);
                fputc((z-z0)>>8 ,fil);
                fputc(vis       ,fil);

                fputc(equivec2indmem(((float)surf.norm[0])*(+1.f/64.f),
                                            ((float)surf.norm[2])*(-1.f/64.f),
                                            ((float)surf.norm[1])*(+1.f/64.f)),fil);

                numvoxs++; xlen[x-x0]++; ylen[(x-x0)*(y1-y0)+y-y0]++;
            }

    fwrite(xlen,(x1-x0)<<2,1,fil);
    fwrite(ylen,((x1-x0)*(y1-y0))<<1,1,fil);
    fseek(fil,28,SEEK_SET); fwrite(&numvoxs,4,1,fil);
    fclose(fil);

    free(ylen); free(xlen);
}
#endif

#if 0
02/24/2012: KVO format, described in pseudo-C code. (assume: int=4 bytes; all fields Little Endian)

    fopen(..)

    int id;       //must be: 'KVOx' or 0x784f564b (Little Endian)
    char version; //must be 0x01
    char lsid;    //log(2) of octree side. sid=(1<<lsid) (Ex.: use 10 for 1024^3)

        //bit 0: pos&ori usage: 0:pivot point (for kvx/kv6 sprite), 1:camera point (for vxl world) //use edgeissol to determine!
        //bit 1: surface voxel color format:
        //    0:log2(paltabn)-bit palette (LUT std::max entries is 4096)
        //    1:24-bit RGB
        //bit 3-2: surface normal format (supports loadkvo only?)
        //  00:not stored
        //  01:equivec (see equivecbits for bit size; see KV6 format for derivation)
        //  10:24-bit signed char x, y, z;
        //bit 4: 1:store surface voxel alpha (8-bit) (supports loadkvo only?)
        //bit 5: 1:store surface voxel tex   (8-bit) (supports loadkvo only?)
        //bit 6: 1:use tail compression; store 32-bit ptr's (not strict bfs order) (not implemented)
        //bit 8-7: 0:store chi&~sol for ls>0, then sol for all. 1:store chi's for all, 2:store chi's for all, then sol (comp)
    short flags;
    char edgeiswrap; //bits 5-0: 0=air/sol, 1=wrap
    char edgeissol;  //bits 5-0: if (!edgeiswrap) { 0=air; 1=sol; } else reserved;
    char equivecbits; //# bits in equivec normal (if used)
    char reserved;

        //Coordinate system is: x=right, y=down, z=forward
    float px, py, pz; //pos of default view; voxel coords (i.e. 0..(1<<lsid)-1)
    float rx, ry, rz; //  right unit vector of default view
    float dx, dy, dz; //   down unit vector of default view
    float fx, fy, fz; //forward unit vector of default view
    float suggdia; //suggested diameter for rendering; usually: (float)(1<<lsid) but may differ for animations

    int octpren; //number of nodes in chi table
    int sur.num; //number of octree leaves (surface voxels); not needed but useful as guess for alloc
    int nod.num; //number of octree nodes after reconstruct; not needed but useful as guess for alloc

        //bits 8-7:description for save mode 0 here: (obsolete!)
    char chi[octpren];           //chi bytes,  BFS order, big to small (ls>0 only). 0=no child, 1=has child
    char sol[popcount8(~chi[*])]; //sol bits in BFS order, big to small, 1 bit per chi=0. sol/air (surfs calced later)
    //flush any extra bits to next byte boundary, i.e.: bitptr = (bitptr+7)&~7;

    if (!(flags&2)) //Palette (if applicable)
    {
        short paltabn; unsigned char bgrcol[paltabn][3]; //blue:[][0], green:[][1], red:[][2], range:{0-255}
    }

        //surface voxel data: colors and normals (optional)
    for(i=0;i<sur.num;i++) //visit all surface voxels (solid voxels next to air) in octree order..
    {
        if (!(flags&2)) { log2(paltabn)_bits palind; } //color:Paletted
                      else { unsigned char b, g, r;     } //color:24-bit RGB
        switch((flags>>2)&7) //normal:
        {
            case 0: /*ignore normal*/ break;
            case 1: equivecbits equivec; break;
            case 2: signed char norm[3]; break;
        }
        if (flags&16) unsigned char alpha; //alpha value
        if (flags&32) unsigned char tex //texture index
    }

    fclose(..);

#endif

#if defined(ENABLE_SAVING_ROUTINES)
    //See SVO_COMP.KC for derivation
static const char hsollut[116] = //log2(n),0,..,n-1
{
    0,0,
    1,0,1, 1,0,2, 1,0,3, 1,0,4, 1,0,5, 1,0,7, 1,0,8, 1,0,10, 1,0,11, 1,0,12, 1,0,13,
    1,0,14, 1,0,16, 1,0,17, 1,0,19, 1,0,21, 1,0,23, 1,0,32, 1,0,34, 1,0,35, 1,0,42, 1,0,43,
    2,0,1,8,9, 2,0,1,32,33, 2,0,2,4,6, 2,0,2,16,18, 2,0,4,16,20, 2,0,8,32,40,
    3,0,1,8,9,32,33,40,41, 3,0,2,4,6,16,18,20,22,
};
static const char hsolind[256] =
{
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 5, 0, 0,
    0, 0, 0, 0, 0, 0, 2, 0, 0, 5, 0, 0, 8, 5, 2, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0,11, 0,11, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 2, 0, 0,11,14,11, 0, 0, 2, 0,
    0,38, 0,38,0,38,41,38,50,107,47,88,44,83,41,38,
    0, 0, 0, 0, 0, 0, 2, 0,17,78,14,11, 8, 5, 2, 0,
    0, 0, 0, 0, 0, 0,20,20, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,53,53,65,62,98,93, 0,56,53,53,59,56,73,53,
    0, 0, 0, 0, 0,23,20,20, 0, 5, 0, 0, 0, 5, 0, 0,
    0, 0, 0, 0,26,23,68,20, 0, 5, 0, 0, 8, 5, 2, 0,
    0, 0, 0,29, 0, 0,20,20, 0,11, 0,11, 0, 0, 0, 0,
    0, 0,32,29, 0, 0,68,20, 0,11,14,11, 0, 0, 2, 0,
    0,35, 0,29, 0,23,20,20, 0,78, 0,11, 0, 5, 0, 0,
    0, 0, 0, 0, 0, 0, 2, 0, 0, 5, 0, 0, 0, 0, 0, 0,
};

    //------------------------------------------------------------
    //07/18/2010: NOTE source models were 5-mip KVX; converted from HRP by Devan
    //
    //                   kv6 kvx1mip    kvo kv6zip kvxzip kvozip
    //               ------- ------- ------ ------ ------ ------
    //0692_knife       99288   26176  14224  26075  11292   9171
    //1158_wiaderkob 1077364  376543 169895 224431 117772  66248
    //2282            748666  279360  60311 209869 116243  35647
    //4407_robotmouse 306708  117302  44310  93035  48478  27443
    //4562_doggy      305460  116592  42911  90414  49231  25366
    //922             307012  120842  33021  93255  48269  17102
    //caco            113804   48843  11630  23404  11195   5123
    //duke_stand      669762  220385 109975 245690 120125  60076
    //enforcer_walk1  987566  389378 157127 338132 196425  93589
    //octa_move1      814876  305706 134344 311672 178442  87857
    //pigtank1975    3854356 1376587 641405 885433 525100 240750
    //pig_shot1       758656  274802 123827 257446 143023  65518
    //reconcar1960    782328  326706 119738 183278 111579  49152
    //trooper_walk1  1079648  412565 174066 405456 230402 110240
    //
    //trike                   333399  56224         70839  21422
    //rocketgun2      254862          32234  55808         16286
    //untitled(vxl)      12425384  16765996    3288588   2297266

#define SAVEKVO_PRINTF 1
typedef struct
{
    FILE *fil;
    oct_t *loct;
    int val, cnt, ls, pass;
        //pass:
        //savmode 0: always requires slow surface search at load!
        // 0 chi&~sol:ls-1..1, bytes (recurse on interesting nodes (some sol + some air)
        // 1 sol/air :ls-1..1, comp (determine air/sol of nonleaf pure nodes)
        // 1 sol     :      0, bytes (uses same code as above because all chi's are 0 at ls==0)
        //savmode 1: ~7% bigger than savmode 0; loads ~5x faster! Hollow interior :/
        // 0 chi     :ls-1..1, bytes
        // 1 chi     :      0, bytes
        //savmode 2:
        // 0 chi     :ls-1..1, bytes
        // 1 chi     :      0, bytes
        // 2 sol     :ls-1..0, comp (nicely compressed using hsolind[]/hsollut[])
    int savmode;
} gsav_t;
    //Save obsolete KVO/KVS formats only:
static void save_recur (gsav_t *sav, octv_t *inode, int ls)
{
    octv_t *onode;
    int i, iup;

    if (!sav->savmode)
    {
        i = (~inode->sol)&inode->chi; //visit child w/ >0air & >0sol; skips pure air/sol
        if (ls == sav->ls)
        {
                //Nonleaf:write all child visit bits to file
            if (!sav->pass) { fputc(i,sav->fil); return; }

                //Leaf:write only if non surface
            for(i^=((1<<8)-1);i;i^=iup) { iup = (-i)&i; if (inode->sol&iup) { sav->val += (1<<sav->cnt); } sav->cnt++; }
            if (sav->cnt >= 8) { fputc(sav->val&255,sav->fil); sav->val >>= 8; sav->cnt -= 8; }
            return;
        }
    }
    else
    {
        if (ls == sav->ls) { fputc(inode->chi,sav->fil); return; }
        i = inode->chi;
    }
    onode = &((octv_t *)sav->loct->nod.buf)[inode->ind];
    for(;i;i^=iup) { iup = (-i)&i; save_recur(sav,&onode[popcount8(inode->chi&(iup-1))],ls-1); }
}

static int visitnodes_call (oct_t *loct, int lstop, int (*myfunc)(oct_t *loct, int ptr, int par, void *mydat), void *mydat)
{
    int ls, ptr, par = 0, stki[OCT_MAXLS], stkn[OCT_MAXLS];

    ls = loct->lsid; stki[ls] = loct->head; stkn[ls] = 1;
    while (1)
    {
        ptr = stki[ls]; stki[ls]++; stkn[ls]--; //2sibly
        if (ls == lstop) { if (myfunc(loct,ptr,par,mydat) < 0) return(-1); }
        if (ls >  lstop) { ls--; stki[ls] = ((octv_t *)loct->nod.buf)[ptr].ind; stkn[ls] = popcount8(((octv_t *)loct->nod.buf)[ptr].chi); par = ptr; } //2child
        while (stkn[ls] <= 0) { ls++; if (ls >= loct->lsid) return(0); } //2parent
    }
}

static int savegetcnt_cb (oct_t *loct, int ptr, int par, void *_) {
    static_cast<void>(loct);
    static_cast<void>(ptr);
    static_cast<void>(par);

    gsav_t *sav = (gsav_t *)_;
    sav->cnt++;
    return(0);
}

static int savechi_cb (oct_t *loct, int ptr, int par, void *_)
{
    static_cast<void>(loct);
    static_cast<void>(par);

    octv_t *inode;
    gsav_t *sav = (gsav_t *)_;

    inode = &((octv_t *)sav->loct->nod.buf)[ptr];

#if (SAVEKVO_PRINTF != 0)
    printf("ls=%d,chi:0x%02x\n",sav->ls,inode->chi);
#endif
    fputc(inode->chi,sav->fil);
    return(0);
}

static int savesol_cb (oct_t *loct, int ptr, int par, void *_)
{
    static_cast<void>(loct);
    static_cast<void>(par);

    octv_t *inode;
    int i, j, k, ln, n;
    gsav_t *sav = (gsav_t *)_;

    inode = &((octv_t *)sav->loct->nod.buf)[ptr];

    if (sav->ls)
    {
#if (SAVEKVO_PRINTF != 0)
        printf("ls=%d,sol:0x%02x\n",sav->ls,inode->sol);
#endif
        fputc(inode->sol,sav->fil); return(0);
    }

    i = hsolind[inode->chi]; ln = hsollut[i]; n = (1<<ln);
    for(k=n-1;k>=0;k--)
    {
        j = hsollut[i+k+1];
        if ((inode->sol&(~inode->chi)) ==  j     ) {                break; }
        if ((inode->sol|  inode->chi ) == (j^255)) { k ^= (n<<1)-1; break; }
    }
#if (SAVEKVO_PRINTF != 0)
    if (k < 0) { printf("Octree corrupt error: chi=0x%02x && sol =0x%02x\n",inode->chi,inode->sol); }
    printf("ls=%d,sol:0x%02x (%d/%d)\n",sav->ls,inode->sol,k,1<<(ln+1));
#endif
    sav->val += (k<<sav->cnt); sav->cnt += ln+1;
    if (sav->cnt >= 8) { fputc(sav->val&255,sav->fil); sav->val >>= 8; sav->cnt -= 8; }

    return(0);
}

//Callback functions for visitnodes_call() for grabbing palette, then saving colors in later pass
static int savegetpal_cb (oct_t *loct, int ptr, int par, void*) {
    static_cast<void>(par);

    surf_t *psurf = &((surf_t *)loct->sur.buf)[ptr];
    return(palget((*(int *)&psurf->b)&0xffffff));
}

typedef struct { int colmode, palval, palcnt, palbits; FILE *fil; } palwrite_t;
static int savecol_cb (oct_t *loct, int ptr, int par, void *_)
{
    static_cast<void>(par);

    int i, col;
    surf_t *psurf = &((surf_t *)loct->sur.buf)[ptr];
    palwrite_t *pw = (palwrite_t *)_;

    col = (*(int *)&psurf->b)&0xffffff;
        //dump on overflow or changing mode
    if (!pw->colmode)
    {
        i = palget(col);
#if (SAVEKVO_PRINTF != 0)
        printf("col:%d\n",i);
#endif
        pw->palval += (i<<pw->palcnt);
        pw->palcnt += pw->palbits;
        while (pw->palcnt >= 8) { fputc(pw->palval&255,pw->fil); pw->palval >>= 8; pw->palcnt -= 8; }
    }
    else
    {
#if (SAVEKVO_PRINTF != 0)
        printf("col:0x%06x\n",col);
#endif
        fwrite(&col,3,1,pw->fil);
    }

    return(0);
}

void savekvo (oct_t *loct, char *filnam, point3f *rpos, point3f *rrig, point3f *rdow, point3f *rfor)
{
    gsav_t sav;
    palwrite_t pw;
    FILE *fil;
    float f;
    int i, octpren;

    fil = fopen(filnam,"wb"); if (!fil) return;
    i = 0x784f564b; fwrite(&i,4,1,fil); //KVOx

    fputc(0x01,fil); //version 1
    fputc(loct->lsid,fil);
    i = 0; fwrite(&i,2,1,fil); //flags

    fputc(loct->edgeiswrap,fil); //edgeiswrap
    fputc(loct->edgeissol,fil); //edgeissol
    i = 0; fwrite(&i,2,1,fil); //reserved

    fwrite(rpos,sizeof(point3f),1,fil);
    fwrite(rrig,sizeof(point3f),1,fil);
    fwrite(rdow,sizeof(point3f),1,fil);
    fwrite(rfor,sizeof(point3f),1,fil);

    f = (float)(1<<loct->lsid);
    fwrite(&f,4,1,fil); //suggested diameter

    fwrite(&i,4,1,fil); //dummy for octpren
    fwrite(&loct->sur.num,4,1,fil);
    fwrite(&loct->nod.num,4,1,fil);

    i = toupper(filnam[std::max((int)strlen(filnam)-1, 0)]);
    sav.savmode = (i == 'S') + (i == 'N')*2;

    sav.fil = fil; sav.loct = loct; sav.cnt = 0; sav.val = 0;
    if (sav.savmode < 2)
    {
            //Write upper levels of octree, 1 byte at a time; octpren = # non-leaf nodes
        sav.pass = 0; i = ftell(fil);
        for(sav.ls=loct->lsid-1;sav.ls> 0;sav.ls--) save_recur(&sav,&((octv_t *)loct->nod.buf)[loct->head],loct->lsid-1); //saving breadth-first is annoying :/
        octpren = ftell(fil)-i;

        sav.pass = 1;
        if (!sav.savmode) { sav.ls = loct->lsid-1; } else { sav.ls = 0; }
        for(;sav.ls>=0;sav.ls--) save_recur(&sav,&((octv_t *)loct->nod.buf)[loct->head],loct->lsid-1); //saving breadth-first is annoying :/
    }
    else
    {
        for(sav.ls=loct->lsid-1;sav.ls> 0;sav.ls--) { visitnodes_call(loct,sav.ls+1,savegetcnt_cb,&sav); } octpren = sav.cnt; sav.cnt = 0;
        for(sav.ls=loct->lsid-1;sav.ls>=0;sav.ls--) { visitnodes_call(loct,sav.ls+1,savechi_cb   ,&sav); }
        for(sav.ls=loct->lsid-1;sav.ls>=0;sav.ls--) { visitnodes_call(loct,sav.ls+1,savesol_cb   ,&sav); }
    }
    while (sav.cnt > 0) { fputc(sav.val&255,fil); sav.val >>= 8; sav.cnt -= 8; } //flush bits to next byte boundary

         //if (<= 4096 unique colors) pw.colmode = 0; else pw.colmode = 1;
    palreset(); pw.colmode = 0;
    if (!visitnodes_call(loct,0,savegetpal_cb,0))
    {
        pw.colmode = 0;
        fwrite(&paltabn,2,1,fil);
        for(i=0;i<paltabn;i++)
        {
#if (SAVEKVO_PRINTF != 0)
            printf("pal(%d):0x%06x\n",i,paltab[i].rgb);
#endif
            fputc((paltab[i].rgb    )&255,fil);
            fputc((paltab[i].rgb>> 8)&255,fil);
            fputc((paltab[i].rgb>>16)&255,fil);
        }
        pw.palbits = bsr(std::max(paltabn-1,1))+1;
    } else pw.colmode = 1;
    pw.palval = 0; pw.palcnt = 0; pw.fil = fil;
    visitnodes_call(loct,0,savecol_cb,&pw);
    if (pw.palcnt) fputc(pw.palval&255,fil); //flush bits

    i = ftell(fil);
    fseek(fil, 6,SEEK_SET); loct->flags = (loct->flags&~2)|(pw.colmode<<1)|(sav.savmode<<7); fwrite(&loct->flags,2,1,fil); //write flags
    fseek(fil,64,SEEK_SET); fwrite(&octpren,4,1,fil);
    fseek(fil, i,SEEK_SET);

    fclose(fil);
}

#endif

//--------------------------------------------------------------------------------------------------

#if defined(ENABLE_LOADING_ROUTINES)
//tree node vs. kvo surf
typedef struct
{
    int  (*isins  )(brush_t *brush, int x0, int y0, int z0, int ls);
    void (*getsurf)(brush_t *brush, int x0, int y0, int z0, surf_t *surf);
    int mx0, my0, mz0, mx1, my1, mz1;
    int flags;

    int paltabn, palbits, palmask, bitcnt, equivecbits, equivecmask;
    float equiveczmulk, equiveczaddk;
    unsigned char *filptr;
    short filags;
} brush_kvo_t;
static void brush_loadkvo_getsurf (brush_t *brush, int x0, int y0, int z0, surf_t *surf)
{
    static_cast<void>(x0);
    static_cast<void>(y0);
    static_cast<void>(z0);

    brush_kvo_t *kvo = (brush_kvo_t *)brush;

        //Load color
    if (!(kvo->filags&2)) //Palette
    {
        *(int *)&surf->b = paltab[((*(int *)kvo->filptr)>>kvo->bitcnt)&kvo->palmask].rgb;
        kvo->bitcnt += kvo->palbits; kvo->filptr += (kvo->bitcnt>>3); kvo->bitcnt &= 7;
    }
    else
    {
        *(int *)&surf->b = (((*(int *)kvo->filptr)>>kvo->bitcnt)&0xffffff); kvo->filptr += 3;
    }

        //Load normal
    switch((kvo->filags>>2)&3)
    {
        case 0: surf->norm[0] = 0; surf->norm[1] = 0; surf->norm[2] = 0; break; //unused
        case 1:
        {
            #define GOLDRAT 0.3819660112501052 //Golden Ratio: 1 - 1/((sqrt(5)+1)/2)
            float fx, fy, fz, r, a;
            int i;

            i = ((*(int *)kvo->filptr)>>kvo->bitcnt)&(kvo->equivecmask);

            fz = (float)i*kvo->equiveczmulk + kvo->equiveczaddk; r = sqrt(1.f - fz*fz);
            a = ((float)i)*(GOLDRAT*M_PI*2); fx = cos(a)*r; fy = sin(a)*r;
            surf->norm[0] = (signed char)cvttss2si(fx*127.0);
            surf->norm[1] = (signed char)cvttss2si(fy*127.0);
            surf->norm[2] = (signed char)cvttss2si(fz*127.0);

            kvo->bitcnt += kvo->equivecbits; kvo->filptr += (kvo->bitcnt>>3); kvo->bitcnt &= 7;
            break;
        }
        case 2: *(int *)&surf->norm[0] = (((*(int *)kvo->filptr)>>kvo->bitcnt)&0xffffff); kvo->filptr += 3; break;
        case 3: surf->norm[0] = 0; surf->norm[1] = 0; surf->norm[2] = 0; break; //reserved
    }

        //Load alpha
    if (kvo->filags&16) { *(unsigned char *)&surf->a = (unsigned char)(((*(int *)kvo->filptr)>>kvo->bitcnt)&0xff); kvo->filptr++; }
                        else { surf->a = 255; }

        //Load tex
    if (kvo->filags&32) { *(unsigned char *)&surf->tex = (unsigned char)(((*(int *)kvo->filptr)>>kvo->bitcnt)&0xff); kvo->filptr++; }
                        else { surf->tex = 0; }
}

static int loadkvo (oct_t *loct, char *filnam, point3f *rpos, point3f *rrig, point3f *rdow, point3f *rfor)
{
    brush_kvo_t brush;
    int i, j, k, iup, nprenodes, nprenodes2, octcn, octvn, bitcnt, *iptr;
    unsigned char *filbuf;

    if (!kzopen(filnam)) return(-1);
    kzseek(0,SEEK_END); i = kztell(); kzseek(0,SEEK_SET);
    filbuf = (unsigned char *)malloc(i); if (!filbuf) { kzclose(); return(-1); }
    kzread(filbuf,i);
    kzclose();

    if (*(int *)&filbuf[0] != 0x784f564b) { free(filbuf); return(-1); } //KVOx

    if (filbuf[4] != 1) { free(filbuf); return(-1); } //version 1
    loct->lsid = (int)filbuf[5];
    brush.filags = *(short *)&filbuf[ 6];

    brush.filptr = &filbuf[12];
    memcpy(rpos,brush.filptr,sizeof(point3f)); brush.filptr += sizeof(point3f);
    memcpy(rrig,brush.filptr,sizeof(point3f)); brush.filptr += sizeof(point3f);
    memcpy(rdow,brush.filptr,sizeof(point3f)); brush.filptr += sizeof(point3f);
    memcpy(rfor,brush.filptr,sizeof(point3f)); brush.filptr += sizeof(point3f);
    brush.filptr += sizeof(float); //suggdia

    nprenodes = *(int *)brush.filptr; brush.filptr += 4;
    octcn     = *(int *)brush.filptr; brush.filptr += 4; //not needed; used only for guess at alloc
    octvn     = *(int *)brush.filptr; brush.filptr += 4; //not needed; used only for guess at alloc

    oct_new(loct,loct->lsid,loct->tilid,(octvn>>2)+octvn,(octcn>>2)+octcn,0);

      //NOTE:must copy edgeis* AFTER oct_new()!
    loct->edgeiswrap = filbuf[ 8]; //edgeiswrap
    loct->edgeissol  = filbuf[ 9]; //edgeissol
    if (loct->edgeissol) loct->flags |= 1;

    nprenodes2 = 1;
    for(i=0;i<nprenodes;i++)
    {
        ((octv_t *)loct->nod.buf)[i].chi = brush.filptr[i]; ((octv_t *)loct->nod.buf)[i].sol = 0; ((octv_t *)loct->nod.buf)[i].mrk = 0; ((octv_t *)loct->nod.buf)[i].mrk2 = 0;
        ((octv_t *)loct->nod.buf)[i].ind = nprenodes2; nprenodes2 += popcount8(((octv_t *)loct->nod.buf)[i].chi); //Generate ind's from chi's
    }

    if (!(brush.filags&384)) //if (non save surface format)
    {
        memset4(&((octv_t *)loct->nod.buf)[nprenodes],0,(nprenodes2-nprenodes)*loct->nod.siz);
        brush.filptr += nprenodes;

        bitcnt = 0; iptr = (int *)brush.filptr;
        for(j=0;j<nprenodes2;j++)
        {
            i = ((octv_t *)loct->nod.buf)[j].chi^((1<<8)-1);
            for(;i;i^=iup,bitcnt++) { iup = (-i)&i; if (iptr[bitcnt>>5]&(1<<bitcnt)) ((octv_t *)loct->nod.buf)[j].sol += iup; }
        }
        //NOTE: chi's of ls==0 calculated inside oct_updatesurfs()

        brush.filptr += ((bitcnt+7)>>3);
    }
    else
    {
        j = 1; //Must skip 0 for GPU shaders which use index 0 as discard
        for(i=nprenodes;i<nprenodes2;i++)
        {
            ((octv_t *)loct->nod.buf)[i].chi = brush.filptr[i]; ((octv_t *)loct->nod.buf)[i].sol = brush.filptr[i]; ((octv_t *)loct->nod.buf)[i].mrk = 0; ((octv_t *)loct->nod.buf)[i].mrk2 = 0;
            ((octv_t *)loct->nod.buf)[i].ind = j; j += popcount8(((octv_t *)loct->nod.buf)[i].chi); //Generate ind's from chi's
        }
        brush.filptr += nprenodes2;

        if ((brush.filags&384) == 256)
        {
            int m, n, chi, ln;

            bitcnt = 0;
            for(m=0;m<nprenodes2;m++)
            {
                if (m < nprenodes) { ((octv_t *)loct->nod.buf)[m].sol = brush.filptr[m]; bitcnt += 8; continue; }

                chi = ((octv_t *)loct->nod.buf)[m].chi;
                i = hsolind[chi]; ln = hsollut[i]; n = (1<<ln);

                k = (((*(int *)&brush.filptr[bitcnt>>3])>>(bitcnt&7)) & ((n<<1)-1));
                if (k < n) { k = hsollut[ k            +i+1]    ; }
                        else { k = hsollut[(k^((n<<1)-1))+i+1]^255; }
                ((octv_t *)loct->nod.buf)[m].sol = (k|chi);

                bitcnt += ln+1;
            }

            brush.filptr += ((bitcnt+7)>>3);
        }

        nprenodes2++; //Must skip 0 for GPU shaders which use index 0 as discard
    }

    loct->nod.num = nprenodes2; loct->nod.ind = nprenodes2;

         //loct->nod.bit: 0<=?<nprenodes2 are 1's; nprenodes2<=?<loct->nod.mal are 0's
    memset4(loct->nod.bit,-1,(nprenodes2>>5)<<2); loct->nod.bit[nprenodes2>>5] = pow2m1[nprenodes2&31];

    if (!(brush.filags&2)) //Palette
    {
        brush.paltabn = (int)(*(short *)brush.filptr); brush.filptr += 2;
        brush.palbits = bsr(std::max(brush.paltabn-1,1))+1; brush.palmask = (1<<brush.palbits)-1;
        for(i=0;i<brush.paltabn;i++) { paltab[i].rgb = *(int *)brush.filptr; brush.filptr += 3; }
    }
    else { brush.paltabn = 0; }
    brush.bitcnt = 0;

    brush.equivecbits = 0;
    if (brush.equivecbits)
    {
        brush.equivecmask = (1<<brush.equivecbits)-1;
        brush.equiveczmulk = 2.f / (float)(1<<brush.equivecbits); brush.equiveczaddk = brush.equiveczmulk*.5f - 1.f;
    }

    brush.mx0 = 0; brush.my0 = 0; brush.mz0 = 0; brush.mx1 = loct->sid; brush.my1 = loct->sid; brush.mz1 = loct->sid;
    brush.isins = 0;/*not used*/
    brush.getsurf = brush_loadkvo_getsurf;
    brush.flags = 0;
    if (!(brush.filags&384)) //if (non save surface format)
    {
        oct_updatesurfs(loct,brush.mx0,brush.my0,brush.mz0,brush.mx1,brush.my1,brush.mz1,(brush_t *)&brush,1);
    }
    else
    {
            //Must skip 0 for GPU shaders which use index 0 as discard
        for(i=1;i<j;i++) { brush_loadkvo_getsurf((brush_t *)&brush,0,0,0,&((surf_t *)loct->sur.buf)[i]); }

        loct->sur.num = j; loct->sur.ind = j;
        memset4(loct->sur.bit,-1,(j>>5)<<2); loct->sur.bit[j>>5] = pow2m1[j&31];
    }
    free(filbuf);

    return(0);
}
#endif

//--------------------------------------------------------------------------------------------------

#if defined(ENABLE_LOADING_ROUTINES)
static int lightvox (int i)
{
    int sh = (((unsigned)i)>>24);
    return((std::min((((i>>16)&255)*sh)>>7,255)<<16)+
             (std::min((((i>>8 )&255)*sh)>>7,255)<< 8)+
             (std::min((((i    )&255)*sh)>>7,255)    ));
}

    //tree node vs. kv6 format
typedef struct { unsigned char b, g, r, a; unsigned short z; char vis, dir; } kv6vox_t;
typedef struct
{
    int  (*isins  )(brush_t *, int x0, int y0, int z0, int ls);
    void (*getsurf)(brush_t *, int x0, int y0, int z0, surf_t *);
    int mx0, my0, mz0, mx1, my1, mz1;
    int flags;

    unsigned int *bitor_[OCT_MAXLS], *bitand_[OCT_MAXLS], lsid, xsiz, ysiz, zsiz, xof, yof, zof;
    void **vcolptr;
    int filtyp;
    char *ppal;
} brush_kv6_t;
    //Returns: 0: node doesn't intersect brush
    //         1: node partially  inside brush
    //         2: node fully      inside brush
static int brush_kv6_isins (brush_t *brush, int x0, int y0, int z0, int ls)
{
    int i, j;
    brush_kv6_t *kv6 = (brush_kv6_t *)brush;

    if (!ls) { i = (x0<<(kv6->lsid*2)) + (z0<<kv6->lsid) + y0; if (!(kv6->bitor_[0][i>>5]&(1<<i))) return(0); return(2); } //<-optimization only
    i = _lrotl(x0,kv6->lsid*2-ls*3) + _lrotl(z0,kv6->lsid-ls*2) + (y0>>ls); j = (1<<i); i >>= 5;
    if (!(kv6->bitor_ [ls][i]&j)) return(0);
    if (  kv6->bitand_[ls][i]&j ) return(2);
    return(1);
}
static void brush_kv6_getsurf (brush_t *brush, int x0, int y0, int z0, surf_t *surf)
{
    kv6vox_t *v;
    float fx, fy, fz;
    int i;
    brush_kv6_t *kv6 = (brush_kv6_t *)brush;

    //fixcnt_getsurf++;
    x0 -= kv6->xof;
    y0 -= kv6->zof;
    z0 -= kv6->yof;
    switch(kv6->filtyp)
    {
        case 0: //.KV6
            i = x0*kv6->ysiz + (kv6->ysiz-1-z0);
            for(v=(kv6vox_t *)kv6->vcolptr[i];v<(kv6vox_t *)kv6->vcolptr[i+1];v++)
            {
                if (v->z < y0) continue;
                if (v->z > y0) break;

                //i = lightvox(*(int *)&v->b);
                surf->b = v->b;
                surf->g = v->g;
                surf->r = v->r;
                surf->a = v->a;
                surf->tex = 0; //((x0^y0^z0)&63);

                equiind2vec(v->dir,&fx,&fy,&fz);
                surf->norm[0] = (signed char)(fx*+64.0);
                surf->norm[1] = (signed char)(fz*+64.0);
                surf->norm[2] = (signed char)(fy*-64.0);
                return;
            }
            break;
        case 1: //.KVX
            {
            unsigned char *v0, *v1;
            v0 = (unsigned char *)kv6->vcolptr[x0*kv6->ysiz + kv6->ysiz-1-z0];
            v1 = (unsigned char *)kv6->vcolptr[x0*kv6->ysiz + kv6->ysiz-1-z0 + 1];
            for(;v0<v1;v0+=((int)v0[1])+3)
            {
                //char slabztop, slabzleng, slabbackfacecullinfo, col[slabzleng];
                if ((y0 < ((int)v0[0])) || (y0 >= ((int)v0[0])+((int)v0[1]))) continue;

                i = v0[y0-((int)v0[0])+3];
                surf->b = kv6->ppal[i*3+2]*4;
                surf->g = kv6->ppal[i*3+1]*4;
                surf->r = kv6->ppal[i*3+0]*4;
                surf->a = 255;
                surf->tex = 0; //((x0^y0^z0)&63);
                surf->norm[0] = 0;
                surf->norm[1] = 0;
                surf->norm[2] = 0;
                return;
            }
            }
            break;
        case 2: //.VOX
            i = ((unsigned char *)kv6->vcolptr[x0*kv6->ysiz + kv6->ysiz-1-z0])[y0];
            surf->b = kv6->ppal[i*3+2]*4;
            surf->g = kv6->ppal[i*3+1]*4;
            surf->r = kv6->ppal[i*3+0]*4;
            surf->a = 255;
            surf->tex = 0; //((x0^y0^z0)&63);
            surf->norm[0] = 0;
            surf->norm[1] = 0;
            surf->norm[2] = 0;
            return;
    }

    *(int *)&surf->b = 0x80808080;
    surf->norm[0] = 0; surf->norm[1] = 0; surf->norm[2] = 0;
    surf->tex = -1;
}

static int loadkv6_kvx_vox (oct_t *loct, char *filnam, point3f *rpos, point3f *rrig, point3f *rdow, point3f *rfor, int filtyp)
{
    brush_kv6_t kv6;
    kv6vox_t *v, *vbuf;
    surf_t surf, *psurf;
    unsigned int *isanymal, *rptr, *wptr;
    unsigned short *ylen; //xsiz*ysiz*sizeof(short)
    float fx, fy, fz, xpiv, ypiv, zpiv;
    int i, j, s, x, y, z, z0, z1, col, xsiz, ysiz, zsiz, xyzsiz, numvoxs, datell;
    int numbytes;
    short *xyoffset;
    char *cptr, *voxdata, pal[768];

    if (!kzopen(filnam)) return(-1);
    switch(filtyp)
    {
        case 0: //.KV6
            kzread(&i,4); if (i != 0x6c78764b) { kzclose(); return(-1); } //Kvxl
            kzread(&xsiz,4); kzread(&ysiz,4); kzread(&zsiz,4);
            kzread(&xpiv,4); kzread(&ypiv,4); kzread(&zpiv,4);
            kzread(&numvoxs,4);
            break;
        case 1: //.KVX
            kzread(&numbytes,4);
            kzread(&xsiz,4); kzread(&ysiz,4); kzread(&zsiz,4);
            kzread(&i,4); xpiv = ((float)i)*(1.0/256.0);
            kzread(&i,4); ypiv = ((float)i)*(1.0/256.0);
            kzread(&i,4); zpiv = ((float)i)*(1.0/256.0);
            break;
        case 2: //.VOX
            kzread(&xsiz,4); kzread(&ysiz,4); kzread(&zsiz,4);
            xpiv = xsiz*.5; ypiv = ysiz*.5; zpiv = zsiz*.5;
            break;
    }

    i = std::max(std::max(xsiz,ysiz),zsiz); if (i <= 0) { kzclose(); return(-1); }
    oct_new(loct,bsr(i-1)+1,0,0,0,0);

    switch(filtyp)
    {
        case 0: //.KV6
            ylen = (unsigned short *)malloc(xsiz*ysiz*sizeof(short)); if (!ylen) { kzclose(); return(-1); }
            kzseek(32+numvoxs*sizeof(kv6vox_t)+(xsiz<<2),SEEK_SET);
            kzread(ylen,xsiz*ysiz*sizeof(short));
            kzseek(32,SEEK_SET);

            equivecinit(255);

            vbuf = (kv6vox_t *)malloc(numvoxs*sizeof(kv6vox_t)); if (!vbuf) { free(ylen); kzclose(); return(-1); }
            kzread(vbuf,numvoxs*sizeof(kv6vox_t));
            break;
        case 1: //.KVX
            xyoffset = (short *)malloc(xsiz*(ysiz+1)*2); if (!xyoffset) { kzclose(); return(-1); }
            voxdata = (char *)malloc(numbytes-24-(xsiz+1)*4-xsiz*(ysiz+1)*2); if (!voxdata) { free(xyoffset); kzclose(); return(-1); }

            kzseek((xsiz+1)*4,SEEK_CUR);
            kzread(xyoffset,xsiz*(ysiz+1)*2);
            kzread(voxdata,numbytes-24-(xsiz+1)*4-xsiz*(ysiz+1)*2);

            kzseek(-768,SEEK_END);
            kzread(pal,768); //0-63
            break;
        case 2: //.VOX
            voxdata = (char *)malloc(xsiz*ysiz*zsiz); //255=transparent
            if (!voxdata) { kzclose(); return(-1); }
            kzread(voxdata,xsiz*ysiz*zsiz);
            kzread(pal,768); //0-63
            break;
    }
    kzclose();

    kv6.xof = std::max(std::min((loct->sid>>1)-(int)xpiv,loct->sid-xsiz),0);
    kv6.yof = std::max(std::min((loct->sid>>1)-(int)ypiv,loct->sid-ysiz),0);
    kv6.zof = std::max(std::min((loct->sid>>1)-(int)zpiv,loct->sid-zsiz),0);

    i = (1<<(loct->lsid*3)); for(j=loct->lsid-1;j>=0;j--) i += (2<<(j*3));
    isanymal = (unsigned int *)malloc((i>>3)+256);
    if (!isanymal) { free(ylen); return(-1); }

    kv6.vcolptr = (void **)malloc((xsiz*ysiz+1)*sizeof(void *));
    if (!kv6.vcolptr) { free(isanymal); free(ylen); return(-1); }

        //Prepare bit buffer pointers
    kv6.bitor_[0] = (unsigned int *)((((int)isanymal)+15)&~15);
    j = (1<<(loct->lsid*3-3)); memset16(kv6.bitor_[0],0,j);

        //Generate bit buffer & vcolptr[][]
    switch(filtyp)
    {
        case 0: //.KV6
            v = vbuf;
            for(x=0;x<xsiz;x++)
                for(y=0;y<ysiz;y++)
                {
                    kv6.vcolptr[x*ysiz+y] = (void *)v;
                    for(i=ylen[x*ysiz+y];i>0;i--,v++)
                    {
                        if (v->vis&16) z = v->z; //air above
                        if (v->vis&32) //air below
                        {
                            j = ((x+kv6.xof)<<(loct->lsid*2)) + ((ysiz-1-y+kv6.yof)<<loct->lsid) + z+kv6.zof;
                            setzrange1(&kv6.bitor_[0][j>>5],j&31,(j&31)+(v->z+1)-z);
                        }
                    }
                }
            kv6.vcolptr[xsiz*ysiz] = (void *)v;
            free(ylen);
            break;
        case 1: //.KVX
            cptr = (char *)voxdata;
            for(x=0;x<xsiz;x++)
                for(y=0;y<ysiz;y++)
                {
                    kv6.vcolptr[x*ysiz+y] = (void *)cptr; z0 = 0;
                    for(i=xyoffset[x*(ysiz+1)+y+1]-xyoffset[x*(ysiz+1)+y];i>0;i-=j)
                    {
                        //char slabztop, slabzleng, slabbackfacecullinfo, col[slabzleng];
                        z1 = ((int)cptr[0])+((int)cptr[1]);
                        if (cptr[2]&16) z0 = ((int)cptr[0]);
                        j = ((x+kv6.xof)<<(loct->lsid*2)) + ((ysiz-1-y+kv6.yof)<<loct->lsid) + z0+kv6.zof;
                        setzrange1(&kv6.bitor_[0][j>>5],j&31,(j&31)+z1-z0); //((int)cptr[0]));
                        z0 = ((int)cptr[0])+((int)cptr[1]);
                        j = ((int)cptr[1])+3; cptr += j;
                    }
                }
            kv6.vcolptr[xsiz*ysiz] = (void *)cptr;
            free(xyoffset);
            break;
        case 2: //.VOX
            cptr = (char *)voxdata;
            for(x=0;x<xsiz;x++)
                for(y=0;y<ysiz;y++)
                {
                    kv6.vcolptr[x*ysiz+y] = (void *)cptr;
                    for(z=0;z<zsiz;z++)
                    {
                        if (cptr[z] == 255) continue;
                        j = ((x+kv6.xof)<<(loct->lsid*2)) + ((ysiz-1-y+kv6.yof)<<loct->lsid) + z+kv6.zof;
                        *(int *)&kv6.bitor_[0][j>>5] |= (1<<j);
                    }
                    cptr += zsiz;
                }
            kv6.vcolptr[xsiz*ysiz] = (void *)cptr;
            break;
    }

        //Generate mips (bitor_ & bitand_)
    kv6.bitand_[0] = kv6.bitor_[0]; //NOTE:bitand_[0] bitmap is same as bitor_[0]; saves 32MB! :)
    j = (1<<(loct->lsid*3-3));
    for(i=1;i<=loct->lsid;i++)
    {
        genmip_voxbits_t gm;

        xyzsiz = (1<<(loct->lsid-i));
        kv6.bitor_[i]  = &kv6.bitand_[i-1][(j+3)>>2]; j >>= 3;
        kv6.bitand_[i] = &kv6.bitor_ [i  ][(j+3)>>2];

        gm.ys = xyzsiz; gm.zs = xyzsiz;
        gm.rbit = kv6.bitor_ [i-1]; gm.wbit = kv6.bitor_ [i]; gm.isor = 1; htrun(genmip_voxbits,&gm,0,xyzsiz,(gm.zs>=16)); //NOTE:MT not safe in genmip_voxbits when gm.zs<16!
        gm.rbit = kv6.bitand_[i-1]; gm.wbit = kv6.bitand_[i]; gm.isor = 0; htrun(genmip_voxbits,&gm,0,xyzsiz,(gm.zs>=16));
    }

    kv6.isins   = brush_kv6_isins;
    kv6.getsurf = brush_kv6_getsurf;
    kv6.flags   = 0;
    kv6.lsid    = loct->lsid; kv6.xsiz = xsiz; kv6.ysiz = ysiz; kv6.zsiz = zsiz;
    kv6.filtyp = filtyp; kv6.ppal = pal;

    oct_mod(loct,(brush_t *)&kv6,1+2);

    free(kv6.vcolptr);
    free(isanymal);
    if (filtyp) free(voxdata);
    return(0);
}
#endif

//--------------------------------------------------------------------------------------------------
#if defined(ENABLE_LOADING_ROUTINES)
    //tree node vs. png format
typedef struct
{
    int  (*isins  )(brush_t *, int x0, int y0, int z0, int ls);
    void (*getsurf)(brush_t *, int x0, int y0, int z0, surf_t *);
    int mx0, my0, mz0, mx1, my1, mz1;
    int flags;

    tiletype pic;
    unsigned char *umax[OCT_MAXLS], *umin[OCT_MAXLS];
    int sid, yofs;
} brush_png_t;
    //Returns: 0: node doesn't intersect brush
    //         1: node partially  inside brush
    //         2: node fully      inside brush
static int brush_png_isins (brush_t *brush, int x0, int y0, int z0, int ls)
{
    brush_png_t *png = (brush_png_t *)brush;
    int i, xs, zs;

    if (ls >= 8) return(1);
    y0 -= png->yofs;
    if (y0&(-256)) { if (y0 < 0) return(0); else return(2); }
    z0 = png->sid-1-z0;
    x0 >>= ls; xs = ((png->pic.x+(1<<ls)-1)>>ls); if (x0 >= xs) return(0);
    z0 >>= ls; zs = ((png->pic.y+(1<<ls)-1)>>ls); if (z0 >= zs) return(0);
    i = z0*xs+x0;
    if (y0         >= png->umax[ls][i]) return(2);
    if (y0+(1<<ls) <= png->umin[ls][i]) return(0);
    return(1);
}
static void brush_png_getsurf (brush_t *brush, int x0, int y0, int z0, surf_t *surf)
{
    brush_png_t *png = (brush_png_t *)brush;
    int i;

    surf->norm[0] = 0; surf->norm[1] = 0; surf->norm[2] = 0;
    surf->a = 255; surf->tex = ((y0-png->yofs)%108);

    z0 = png->sid-1-z0;
    if ((x0 < png->pic.x) && (z0 < png->pic.y))
    {
        *(int *)&surf->b = *(int *)(png->pic.p*z0 + (x0<<2) + png->pic.f);
        return;
    }
    *(int *)&surf->b = 0x404040;
}

    //05/28/2012:     ri_2048: ri_4096:
    //---------------------------------
    //kpzload            208.1    791.3
    //oct_init             1.3      1.4
    //extract alpha        4.8     19.5
    //gen std::min/std::max mips    25.7     89.9
    //oct_mod()         2785.8   8272.3
    //

static int loadpng (oct_t *loct, char *filnam, point3f *rpos, point3f *rrig, point3f *rdow, point3f *rfor)
{
    brush_png_t png;
    int i, j, ls, x, y, xx, yy, xs, ys, nxs, nys, chi = 0, sol = 0, ind = 0;
    unsigned char cmin, cmax, *rptr, *wptr, xorhei;

    if (filnam[0] == '~') { xorhei =   0; filnam++; } //255 is lo (tomland.png)
                          else { xorhei = 255;           } //255 is hi (default/most others)
    kpzload(filnam,&png.pic.f,&png.pic.p,&png.pic.x,&png.pic.y);
    if (!png.pic.f) return(-1);

    oct_new(loct,bsr(std::max(png.pic.x,png.pic.y)-1)+1,0,0,0,0);

    loct->flags |= 1;
    loct->edgeissol = 61;

    rpos->x = loct->sid*.5; rpos->y = 64; rpos->z = loct->sid*.5;
    rrig->x = 1.0; rrig->y = 0.0; rrig->z = 0.0;
    rdow->x = 0.0; rdow->y = 1.0; rdow->z = 0.0;
    rfor->x = 0.0; rfor->y = 0.0; rfor->z = 1.0;

        //generate std::min/std::max mip tables..
    for(i=0,ls=0;ls<loct->lsid;ls++) { i += ((png.pic.x+(1<<ls)-1)>>ls) * ((png.pic.y+(1<<ls)-1)>>ls) * ((ls!=0)+1); }
    png.umin[0] = (unsigned char *)malloc(i); if (!png.umin[0]) { free((void *)png.pic.f); return(-1); }
    png.umax[0] = png.umin[0]; i = png.pic.x*png.pic.y;
    for(ls=1;ls<loct->lsid;ls++)
    {
        png.umin[ls] = png.umax[ls-1]+i; i = ((png.pic.x+(1<<ls)-1)>>ls) * ((png.pic.y+(1<<ls)-1)>>ls);
        png.umax[ls] = png.umin[ls  ]+i;
    }
    wptr = &png.umin[0][0]; rptr = (unsigned char *)(png.pic.f + 3);
    for(y=0;y<png.pic.y;y++,wptr+=png.pic.x,rptr+=png.pic.p)
        for(x=0;x<png.pic.x;x++)
            wptr[x] = rptr[x<<2]^xorhei;
    nxs = png.pic.x; nys = png.pic.y;
    for(ls=1;ls<loct->lsid;ls++)
    {
        xs = nxs; nxs = ((png.pic.x+(1<<ls)-1)>>ls);
        ys = nys; nys = ((png.pic.y+(1<<ls)-1)>>ls);
        for(y=0;y*2+1<ys;y++)
        {
            i = y*xs*2; j = y*nxs; xx = ((xs-1)>>1);
            for(x=0;x<xx;x++,i+=2)
            {
                rptr = &png.umin[ls-1][i]; png.umin[ls][j+x] = std::min(std::min(rptr[0],rptr[1]),std::min(rptr[xs],rptr[xs+1]));
                rptr = &png.umax[ls-1][i]; png.umax[ls][j+x] = std::max(std::max(rptr[0],rptr[1]),std::max(rptr[xs],rptr[xs+1]));
            }
            if (x < nxs)
            {
                rptr = &png.umin[ls-1][i]; png.umin[ls][j+x] = std::min(rptr[0],rptr[xs]);
                rptr = &png.umax[ls-1][i]; png.umax[ls][j+x] = std::max(rptr[0],rptr[xs]);
            }
        }
        if (y < nys)
        {
            for(x=0;x<nxs;x++)
            {
                cmin = 255; cmax = 0;
                for(xx=0;xx<2;xx++)
                {
                    if (x*2+xx >= xs) continue;
                    cmin = std::min(cmin,png.umin[ls-1][(y*2)*xs+(x*2+xx)]);
                    cmax = std::max(cmax,png.umax[ls-1][(y*2)*xs+(x*2+xx)]);
                }
                png.umin[ls][y*nxs+x] = cmin;
                png.umax[ls][y*nxs+x] = cmax;
            }
        }
    }

    png.isins   = brush_png_isins;
    png.getsurf = brush_png_getsurf;
    png.flags   = 0;
    png.sid     = loct->sid;
    png.yofs    = (std::max(loct->sid-256,0)>>1)&-256;
    oct_mod(loct,(brush_t *)&png,1+2);

    rpos->y += png.yofs;

    free(png.umin[0]);
    free((void *)png.pic.f);

    return(0);
}
#endif

//--------------------------------------------------------------------------------------------------

#if defined(ENABLE_LOADING_ROUTINES)
    //tree node vs. vxl format
typedef struct
{
    int  (*isins  )(brush_t *, int x0, int y0, int z0, int ls);
    void (*getsurf)(brush_t *, int x0, int y0, int z0, surf_t *);
    int mx0, my0, mz0, mx1, my1, mz1;
    int flags;

    unsigned int *bitor_[9], *bitand_[9];
    unsigned char **vcolptr;
    int vsid, lvsid;
} brush_vxl_t;
    //Returns: 0: node doesn't intersect brush
    //         1: node partially  inside brush
    //         2: node fully      inside brush
static int brush_vxl_isins (brush_t *brush, int x0, int y0, int z0, int ls)
{
    int i, j;
    brush_vxl_t *vxl = (brush_vxl_t *)brush;

        //0  128 256 384 512 640 768 896 1024
        //+---+---+---+---+---+---+---+---+
        //        | x x xx|xxxxxxxxxxxxxxxx
        //
        //0  128 256 384 512
        //+---+---+---+---+
        //        | x x xx|

    y0 -= 256; //arbitrary height adjustment
    if (ls >= 8) return(1);
    if (y0&(-256)) { if (y0 < 0) return(0); else return(2); }
    if (!ls) { i = (x0<<(vxl->lvsid+8)) + (z0<<8) + y0; if (!(vxl->bitor_[0][i>>5]&(1<<i))) return(0); return(2); } //<-optimization only
    i = _lrotl(x0,vxl->lvsid+8-ls*3) + _lrotl(z0,8-ls*2) + (y0>>ls); j = (1<<i); i >>= 5;
    if (!(vxl->bitor_ [ls][i]&j)) return(0);
    if (  vxl->bitand_[ls][i]&j ) return(2);
    return(1);
}
static void brush_vxl_getsurf (brush_t *brush, int x0, int y0, int z0, surf_t *surf)
{
    int z;
    unsigned char *v;
    brush_vxl_t *vxl = (brush_vxl_t *)brush;

    //fixcnt_getsurf++;
    surf->norm[0] = 0; surf->norm[1] = 0; surf->norm[2] = 0;
    surf->tex = ((rand()*108)>>15);

    *(int *)&surf->b = 0x40404040;
    y0 -= 256; if (y0&(-256)) return; //height adjustment
    v = vxl->vcolptr[(x0<<vxl->lvsid)+z0];
    while (1)
    {
        if (y0 <= v[2]) { *(int *)&surf->b = lightvox(*(int *)&v[(y0-v[1]+1)<<2]); break; }
        if (!v[0]) break; z = v[2]-v[1]-v[0]+2; v += v[0]*4;
        if (y0 < z+v[3]) break;
        if (y0 <  v[3]) { *(int *)&surf->b = lightvox(*(int *)&v[(y0-v[3]  )<<2]); break; }
    }
}

static void loadvxl_genbitbuf (int y, void *vptr)
{
    int x, z;
    unsigned char *v;
    brush_vxl_t *vxl = (brush_vxl_t *)vptr;

    memset4(&vxl->bitor_[0][y<<(vxl->lvsid+8-5)],-1,1<<(vxl->lvsid+8-3));
    for(x=0;x<vxl->vsid;x++)
    {
        v = vxl->vcolptr[(y<<vxl->lvsid)+x];
        for(z=0;1;v+=v[0]*4,z=v[3])
        {
            setzrange0(&vxl->bitor_[0][((y<<vxl->lvsid)+x)<<(8-5)],z,v[1]);
            if (!v[0]) break;
        }
        v += ((((int)v[2])-((int)v[1])+2)<<2);
    }
}

static int loadvxl (oct_t *loct, char *filnam, point3f *rpos, point3f *rrig, point3f *rdow, point3f *rfor)
{
    point3d dp;
    double d;
    brush_vxl_t vxl;
    genmip_voxbits_t gm;
    unsigned int *isanymal, *rptr, *wptr;
    unsigned char *v, *vbuf;
    int i, j, s, x, y, z;

    if (!kzopen(filnam)) return(-1);
    kzread(&i,4);
    if (i == 0x09072000)
    {
        kzread(&i,4);
        kzread(&j,4); if ((i != j) || (i > 16384)) { kzclose(); return(-1); }
        for(vxl.lvsid=0;(1<<vxl.lvsid)<i;vxl.lvsid++);
        vxl.vsid = (1<<vxl.lvsid); if (vxl.vsid != j) { kzclose(); return(-1); }
        kzread(&dp,24); rpos->x = dp.x; rpos->y = dp.z+256.0; rpos->z = vxl.vsid-dp.y; //camera position
        kzread(&dp,24); rrig->x = dp.x; rrig->y = dp.z; rrig->z = -dp.y; //unit right vector
        kzread(&dp,24); rdow->x = dp.x; rdow->y = dp.z; rdow->z = -dp.y; //unit down vector
        kzread(&dp,24); rfor->x = dp.x; rfor->y = dp.z; rfor->z = -dp.y; //unit forward vector
    }
    else //AoS format w/no header info :/
    {
        kzseek(0,SEEK_SET);
        vxl.lvsid = 9; vxl.vsid = (1<<vxl.lvsid); //AoS format has no header :/
        rpos->x = 256.0; rpos->y = 1024.0; rpos->z = +0.0; //camera position
        rrig->x = 1.0; rrig->y = 0.0; rrig->z = 0.0; //unit right vector
        rdow->x = 0.0; rdow->y = 0.0; rdow->z = 1.0; //unit down vector
        rfor->x = 0.0; rfor->y =-1.0; rfor->z = 0.0; //unit forward vector
    }

        //Allocate huge buffer and load rest of file into it...
    x = kztell(); kzseek(0,SEEK_END); i = kztell()-x; kzseek(x,SEEK_SET);
    vbuf = (unsigned char *)malloc(i); if (!vbuf) { kzclose(); return(-1); }
    kzread(vbuf,i);
    kzclose();

        //~41.1431MB for vsid=1024
    for(j=0,i=vxl.vsid,s=256;s>0;i>>=1,s>>=1)
    {
        if (!j) j += i*i*s; //don't need 2 buffers at full res
            else j += i*i*s*2;
    }
    j = (j>>2)+256;

    isanymal = (unsigned int *)malloc(j);
    if (!isanymal) return(-1);

        //4MB
    vxl.vcolptr = (unsigned char **)malloc(vxl.vsid*vxl.vsid*sizeof(unsigned char *));
    if (!vxl.vcolptr) { free(isanymal); return(-1); }

        //Generate vcolptr[][]
    v = vbuf;
    for(x=vxl.vsid-1;x>=0;x--) //NOTE:can't change for loop order
        for(y=0;y<vxl.vsid;y++) //NOTE:can't change for loop order
        {
            vxl.vcolptr[(y<<vxl.lvsid)+x] = v;
            for(z=0;v[0];v+=v[0]*4,z=v[3]);
            v += ((((int)v[2])-((int)v[1])+2)<<2);
        }

        //Prepare bit buffer pointers
    vxl.bitor_[0] = (unsigned int *)((((int)isanymal)+15)&~15);
    htrun(loadvxl_genbitbuf,&vxl,0,vxl.vsid,true); //Generate bit buffer using vcolptr's

        //Generate mips (bitor_ & bitand_)
    j = ((vxl.vsid*vxl.vsid*256)>>3);
    vxl.bitand_[0] = vxl.bitor_[0]; //NOTE:bitand_[0] bitmap is same as bitor_[0]; saves 32MB! :)
    for(i=1;i<=8;i++)
    {
        vxl.bitor_[i]  = &vxl.bitand_[i-1][(j+3)>>2]; j >>= 3;
        vxl.bitand_[i] = &vxl.bitor_ [i  ][(j+3)>>2];

        gm.ys = (vxl.vsid>>i); gm.zs = (256>>i);
        gm.rbit = vxl.bitor_ [i-1]; gm.wbit = vxl.bitor_ [i]; gm.isor = 1; htrun(genmip_voxbits,&gm,0,gm.ys,(gm.zs>=16)); //NOTE:MT not safe in genmip_voxbits when gm.zs<16!
        gm.rbit = vxl.bitand_[i-1]; gm.wbit = vxl.bitand_[i]; gm.isor = 0; htrun(genmip_voxbits,&gm,0,gm.ys,(gm.zs>=16));
    }

    vxl.isins   = brush_vxl_isins;
    vxl.getsurf = brush_vxl_getsurf;
    vxl.flags = 1;

    oct_new(loct,vxl.lvsid,0,0,0,0);
    loct->flags |= 1;
    loct->edgeissol = 61;
    oct_mod(loct,(brush_t *)&vxl,1+2);

    free(vxl.vcolptr);
    free(isanymal);

    free(vbuf);
    return(0);
}
#endif
