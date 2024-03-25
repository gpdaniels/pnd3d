// This file has been modified from Ken Silverman's original release

#include "octree_modify.hpp"

// TODO: REMOVE!
#include "octree_renderer.hpp"
// TODO: REMOVE!

#include "octree.hpp"
#include "octree_physics.hpp"
#include "brushes.hpp"

#include "types/parallelepiped.hpp"

#include "utilities/bits.hpp"
#include "utilities/gl.hpp"
#include "utilities/macros.hpp"
#include "utilities/morton.hpp"
#include "utilities/popcount.hpp"

#include <algorithm>
#include <malloc.h>
#include <emmintrin.h>

#define cvttss2si(f) _mm_cvtt_ss2si(_mm_set_ss(f))

//(mode&1)==0:air, !=0:sol
//(mode&2)!=0:do oct_updatesurfs() here (helpful because bbox calculated inside oct_mod_recur)
//(mode&4)!=0:do hover check, ==0:not
//(mode&8)!=0:do normal refresh, ==0:not

void oct_mod (oct_t *loct, brush_t *brush, int mode)
{
    //temp needed in case bitalloc() realloc's nod.buf
    octv_t nnode;

    if constexpr (oct_usegpu && oct_usegpubo) {
        if (!loct->gsurf) {
            loct->gsurf = (surf_t *)bo_begin(loct->bufid,0);
        }
    }

    brush->mx0 = 0x7fffffff; brush->mx1 = 0x80000000;
    brush->my0 = 0x7fffffff; brush->my1 = 0x80000000;
    brush->mz0 = 0x7fffffff; brush->mz1 = 0x80000000;
    oct_mod_recur(loct,loct->head,0,0,0,loct->lsid-1,&nnode,brush,(mode&1)*((1<<8)-1));
    ((octv_t *)loct->nod.buf)[loct->head] = nnode; //can't write node directly to loct->nod.buf[loct->head] because of possible realloc inside

    if (mode&2) { oct_updatesurfs(loct,brush->mx0,brush->my0,brush->mz0,brush->mx1,brush->my1,brush->mz1,brush,mode); }
    if (mode&4) { oct_hover_check(loct,brush->mx0,brush->my0,brush->mz0,brush->mx1,brush->my1,brush->mz1,loct->recvoctfunc); }
    if (mode&8) { int i = (mode&1)^1; oct_refreshnorms(loct,2,brush->mx0-i,brush->my0-i,brush->mz0-i,brush->mx1+i,brush->my1+i,brush->mz1+i); }

    oct_checkreducesizes(loct);

    if constexpr (oct_usegpu && !oct_usegpubo) {
        ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_2D,loct->octid);
        ((PFNGLTEXSUBIMAGE2D)glfp[glTexSubImage2D])(GL_TEXTURE_2D,0,0,0,loct->gxsid,(loct->sur.mal*2)>>loct->glxsid,GL_RGBA,GL_UNSIGNED_BYTE,(void *)loct->sur.buf);
    }
}

void oct_paint (oct_t *loct, brush_t *brush)
{
    typedef struct { octv_t *ptr; int x, y, z, j; } stk_t;
    stk_t stk[OCT_MAXLS];
    octv_t *ptr;
    int i, j, ls, s, x, y, z, nx, ny, nz, ind;

    if constexpr (oct_usegpu) {
        if constexpr (oct_usegpubo) {
            if (!loct->gsurf) loct->gsurf = (surf_t *)bo_begin(loct->bufid,0);
        }
        else {
            ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_2D,loct->octid);
        }
    }

    ls = loct->lsid-1; s = (1<<ls); ptr = &((octv_t *)loct->nod.buf)[loct->head];
    x = 0; y = 0; z = 0; j = 8-1;
    while (1)
    {
        i = (1<<j); if (!(ptr->chi&i)) goto tosibly;

        nx = (( j    &1)<<ls)+x;
        ny = (((j>>1)&1)<<ls)+y;
        nz = (((j>>2)&1)<<ls)+z;
        if (!brush->isins(brush,nx,ny,nz,ls)) goto tosibly;

        if (ls > 0)
        {
            stk[ls].ptr = ptr; stk[ls].x = x; stk[ls].y = y; stk[ls].z = z; stk[ls].j = j; ls--; s >>= 1; //2child
            ptr = &((octv_t *)loct->nod.buf)[popcount8(ptr->chi&(i-1)) + ptr->ind];
            x = nx; y = ny; z = nz; j = 8-1; continue;
        }

        ind = popcount8(ptr->chi&(i-1)) + ptr->ind;

        brush->getsurf(brush,nx,ny,nz,&((surf_t *)loct->sur.buf)[ind]);
        if constexpr (oct_usegpu) {
            if constexpr (oct_usegpubo) {
                memcpy(&loct->gsurf[ind],&((surf_t *)loct->sur.buf)[ind],loct->sur.siz);
            }
            else {
                ((PFNGLTEXSUBIMAGE2D)glfp[glTexSubImage2D])(GL_TEXTURE_2D,0,(ind&((loct->gxsid>>1)-1))<<1,ind>>(loct->glysid-1),2,1,GL_RGBA,GL_UNSIGNED_BYTE,(void *)&((surf_t *)loct->sur.buf)[ind]);
            }
        }

tosibly:;
        j--; if (j >= 0) continue;
        do { ls++; s <<= 1; if (ls >= loct->lsid) goto break2; j = stk[ls].j-1; } while (j < 0); //2parent
        ptr = stk[ls].ptr; x = stk[ls].x; y = stk[ls].y; z = stk[ls].z;
    }

break2:;
}

void oct_writesurf (oct_t *loct, int ind, surf_t *psurf)
{
    memcpy(&((surf_t *)loct->sur.buf)[ind],psurf,loct->sur.siz);
    if constexpr (oct_usegpu) {
        if constexpr (oct_usegpubo) {
            if (!loct->gsurf) {
                loct->gsurf = (surf_t *)bo_begin(loct->bufid,0);
            }
            memcpy(&loct->gsurf[ind],psurf,loct->sur.siz);
        }
        else {
            ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_2D,loct->octid);
            ((PFNGLTEXSUBIMAGE2D)glfp[glTexSubImage2D])(GL_TEXTURE_2D,0,(ind&((loct->gxsid>>1)-1))<<1,ind>>(loct->glysid-1),2,1,GL_RGBA,GL_UNSIGNED_BYTE,(void *)&((surf_t *)loct->sur.buf)[ind]);
        }
    }
}

void oct_copysurfs (oct_t *loct)
{
    if constexpr (oct_usegpu && oct_usegpubo)
    {
        if (!loct->gsurf) loct->gsurf = (surf_t *)bo_begin(loct->bufid,0);
        memcpy(loct->gsurf,loct->sur.buf,loct->sur.mal*loct->sur.siz);
    }
}

//--------------------------------------------------------------------------------------------------
//tree node vs. other octree
//Returns: 0: node doesn't intersect brush
//         1: node partially  inside brush
//         2: node fully      inside brush
typedef struct
{
    BRUSH_HEADER;
    oct_t *oct;
    point3f p[OCT_MAXLS], r[OCT_MAXLS], d[OCT_MAXLS], f[OCT_MAXLS];
    pgram3d_t bi[OCT_MAXLS*2];
} brush_oct_t;
static int brush_oct_isins (brush_t *brush, int x0, int y0, int z0, int ls)
{
    typedef struct { octv_t *ptr; int x, y, z, j; } stk_t;
    stk_t stk[OCT_MAXLS];
    oct_t *loct;
    octv_t *ptr;
    float dx, dy, dz, ex, ey, ez;
    int i, j, k, s, ls2, x, y, z, x1, y1, z1, nx, ny, nz;
    brush_oct_t *boct = (brush_oct_t *)brush;
    loct = boct->oct;

    dx = x0*boct->r[1].x + y0*boct->d[1].x + z0*boct->f[1].x + boct->p[ls].x;
    dy = x0*boct->r[1].y + y0*boct->d[1].y + z0*boct->f[1].y + boct->p[ls].y;
    dz = x0*boct->r[1].z + y0*boct->d[1].z + z0*boct->f[1].z + boct->p[ls].z;

        //Bounded box init:
    ex = fabs(boct->r[ls].x) + fabs(boct->d[ls].x) + fabs(boct->f[ls].x);
    ey = fabs(boct->r[ls].y) + fabs(boct->d[ls].y) + fabs(boct->f[ls].y);
    ez = fabs(boct->r[ls].z) + fabs(boct->d[ls].z) + fabs(boct->f[ls].z);
    x0 = std::max(cvttss2si(dx-ex),0); x1 = std::min(cvttss2si(dx+ex),loct->sid);
    y0 = std::max(cvttss2si(dy-ey),0); y1 = std::min(cvttss2si(dy+ey),loct->sid);
    z0 = std::max(cvttss2si(dz-ez),0); z1 = std::min(cvttss2si(dz+ez),loct->sid);

    k = 0;
    i = pgram3d_isint(&boct->bi[OCT_MAXLS-loct->lsid+ls],loct->sid-dx*2,loct->sid-dy*2,loct->sid-dz*2);
    if (i == 0) return(0);
    if (i == 1) k = 1;

    ls2 = loct->lsid-1; s = (1<<ls2); ptr = &((octv_t *)loct->nod.buf)[loct->head];
    x = 0; y = 0; z = 0; j = 8-1;
    while (1)
    {
        nx = (( j    &1)<<ls2)+x; if ((nx > x1) || (nx+s <= x0)) goto tosibly;
        ny = (((j>>1)&1)<<ls2)+y; if ((ny > y1) || (ny+s <= y0)) goto tosibly;
        nz = (((j>>2)&1)<<ls2)+z; if ((nz > z1) || (nz+s <= z0)) goto tosibly;

        i = (1<<j);
        if (ptr->chi&(~ptr->sol)&i) //mixed air&sol
        {
            stk[ls2].ptr = ptr; stk[ls2].x = x; stk[ls2].y = y; stk[ls2].z = z; stk[ls2].j = j; ls2--; s >>= 1; //2child
            ptr = &((octv_t *)loct->nod.buf)[popcount8(ptr->chi&(i-1)) + ptr->ind]; x = nx; y = ny; z = nz; j = 8-1;
            continue;
        }

            //Do slow exact check for leaf (pure air/sol) nodes
        if (!pgram3d_isint(&boct->bi[OCT_MAXLS-ls2+ls],(nx-dx)*2+s,(ny-dy)*2+s,(nz-dz)*2+s)) goto tosibly;

        if (ptr->sol&i) k |= 2; else k |= 1;
        if (k == 3) return(1);

tosibly:;
        j--; if (j >= 0) continue;
        do { ls2++; s <<= 1; if (ls2 >= loct->lsid) return(k&2); j = stk[ls2].j-1; } while (j < 0); //2parent
        ptr = stk[ls2].ptr; x = stk[ls2].x; y = stk[ls2].y; z = stk[ls2].z;
    }
}
static void brush_oct_getsurf (brush_t *brush, int x0, int y0, int z0, surf_t *surf)
{
    static const int dirx[27] = {0, -1,+1, 0, 0, 0, 0, -1,+1,-1,+1, 0, 0, 0, 0,-1,-1,+1,+1, -1,+1,-1,+1,-1,+1,-1,+1};
    static const int diry[27] = {0,  0, 0,-1,+1, 0, 0, -1,-1,+1,+1,-1,+1,-1,+1, 0, 0, 0, 0, -1,-1,+1,+1,-1,-1,+1,+1};
    static const int dirz[27] = {0,  0, 0, 0, 0,-1,+1,  0, 0, 0, 0,-1,-1,+1,+1,-1,+1,-1,+1, -1,-1,-1,-1,+1,+1,+1,+1};
    int i, j, x, y, z;
    brush_oct_t *boct = (brush_oct_t *)brush;

    x = cvttss2si(x0*boct->r[1].x + y0*boct->d[1].x + z0*boct->f[1].x + boct->p[0].x);
    y = cvttss2si(x0*boct->r[1].y + y0*boct->d[1].y + z0*boct->f[1].y + boct->p[0].y);
    z = cvttss2si(x0*boct->r[1].z + y0*boct->d[1].z + z0*boct->f[1].z + boct->p[0].z);

    for(i=0;i<27;i++)
    {
        j = oct_getsurf(boct->oct,dirx[i]+x,diry[i]+y,dirz[i]+z); if (j == -1) continue;
        (*surf) = ((surf_t *)boct->oct->sur.buf)[j]; return;
    }
      //:/
    surf->b = 255; surf->g = 0; surf->r = 255; surf->a = 255;
    surf->norm[0] = 0; surf->norm[1] = 0; surf->norm[2] = 0; surf->tex = 0;
}
static void brush_oct_init (brush_oct_t *boct, oct_t *loct0, point3d *p0, point3d *r0, point3d *d0, point3d *f0, oct_t *loct1, point3d *p1, point3d *r1, point3d *d1, point3d *f1)
{
    point3d cr, cd, cf, dr, dd, df;
    point3f ap, ar, ad, af;
    double d, mat[9];
    int i;

    boct->isins   = brush_oct_isins;
    boct->getsurf = brush_oct_getsurf;
    boct->oct = loct1;

    //FIXFIXFIXFIX:can remove this inverse now that pgram3d handles it!
    //transform loct0 to make loct1 identity (p1=<0,0,0>, r1=<1,0,0>, d1=<0,1,0>, f1=<0,0,1>)
    //[r,d,f,p] = loct1^-1 * loct0
    invert3x3(r1,d1,f1,mat);
    ar.x = r0->x*mat[0] + r0->y*mat[1] + r0->z*mat[2];
    ar.y = r0->x*mat[3] + r0->y*mat[4] + r0->z*mat[5];
    ar.z = r0->x*mat[6] + r0->y*mat[7] + r0->z*mat[8];
    ad.x = d0->x*mat[0] + d0->y*mat[1] + d0->z*mat[2];
    ad.y = d0->x*mat[3] + d0->y*mat[4] + d0->z*mat[5];
    ad.z = d0->x*mat[6] + d0->y*mat[7] + d0->z*mat[8];
    af.x = f0->x*mat[0] + f0->y*mat[1] + f0->z*mat[2];
    af.y = f0->x*mat[3] + f0->y*mat[4] + f0->z*mat[5];
    af.z = f0->x*mat[6] + f0->y*mat[7] + f0->z*mat[8];
    dr.x = p0->x-p1->x; dr.y = p0->y-p1->y; dr.z = p0->z-p1->z;
    ap.x = dr.x*mat[0] + dr.y*mat[1] + dr.z*mat[2];
    ap.y = dr.x*mat[3] + dr.y*mat[4] + dr.z*mat[5];
    ap.z = dr.x*mat[6] + dr.y*mat[7] + dr.z*mat[8];

    for(i=loct0->lsid-1;i>=0;i--)
    {
        d = pow(2.0,(double)i)*0.5;
        boct->r[i].x = ar.x*d; boct->d[i].x = ad.x*d; boct->f[i].x = af.x*d; boct->p[i].x = boct->r[i].x + boct->d[i].x + boct->f[i].x + ap.x;
        boct->r[i].y = ar.y*d; boct->d[i].y = ad.y*d; boct->f[i].y = af.y*d; boct->p[i].y = boct->r[i].y + boct->d[i].y + boct->f[i].y + ap.y;
        boct->r[i].z = ar.z*d; boct->d[i].z = ad.z*d; boct->f[i].z = af.z*d; boct->p[i].z = boct->r[i].z + boct->d[i].z + boct->f[i].z + ap.z;
    }

    cr.x = 1.0; cr.y = 0.0; cr.z = 0.0;
    cd.x = 0.0; cd.y = 1.0; cd.z = 0.0;
    cf.x = 0.0; cf.y = 0.0; cf.z = 1.0;
    for(i=-loct1->lsid;i<loct0->lsid;i++)
    {
        d = (double)(1<<(i+17));
        dr.x = ar.x*d; dr.y = ar.y*d; dr.z = ar.z*d;
        dd.x = ad.x*d; dd.y = ad.y*d; dd.z = ad.z*d;
        df.x = af.x*d; df.y = af.y*d; df.z = af.z*d;
        pgram3d_init(&boct->bi[i+OCT_MAXLS],&cr,&cd,&cf,&dr,&dd,&df);
    }
}

//--------------------------------------------------------------------------------------------------
    //NOTE:volume of solid must be < 65536 for brush_bmp() to work! (example: 40^3 is always safe)
static void brush_bmp_calcsum (unsigned short *boxsum, char *bmp, int valeq, int xs, int ys, int zs)
{
    unsigned short *nboxsum;
    int i, j, x, y, z, xmul, ymul, zmul;

    zmul = 1; ymul = xs+1; xmul = (ys+1)*ymul;
    for(i=0;i<=xs;i++) boxsum[i*xmul] = 0;
    for(i=0;i<=ys;i++) boxsum[i*ymul] = 0;
    for(i=0;i<=zs;i++) boxsum[i*zmul] = 0;

    for(i=0,z=0;z<zs;z++)
        for(y=0;y<ys;y++)
        {
            nboxsum = &boxsum[y*ymul + z*zmul];
            for(j=0,x=0;x<xs;x++,i++,nboxsum+=xmul)
            {
                j += (bmp[i] == valeq);
                nboxsum[xmul+ymul+zmul] = nboxsum[xmul+ymul] + nboxsum[xmul+zmul] - nboxsum[xmul] + j;
            }
        }
}

//returns: 0:air, 1:surf (c valid), 2:interior
int oct_getvox (oct_t *loct, int x, int y, int z, surf_t **surf)
{
    octv_t *ptr;
    int i, s;

    if ((x|y|z)&loct->nsid) return(0);
    for(i=loct->head,s=(loct->sid>>1);1;s>>=1)
    {
        ptr = &((octv_t *)loct->nod.buf)[i]; //2child
        i = 1;
        if (x&s) i <<= 1;
        if (y&s) i <<= 2;
        if (z&s) i <<= 4;
        if (!(ptr->chi&i)) { if (ptr->sol&i) return(2); return(0); }
        i = popcount8(ptr->chi&(i-1)) + ptr->ind; //2child
        if (s <= 1) { if (surf) (*surf) = &((surf_t *)loct->sur.buf)[i]; return(1); }
    }
}

//faster than oct_getvox() if only sol vs. air is needed (doesn't find surface pointer)
//returns: 0:air, 1:sol (surface or interior)
int oct_getsol (oct_t *loct, int x, int y, int z)
{
    octv_t *ptr;
    int i, s;

    if ((x|y|z)&loct->nsid) return(0);
    for(i=loct->head,s=(loct->sid>>1);1;s>>=1)
    {
        ptr = &((octv_t *)loct->nod.buf)[i]; //2child
        i = 1;
        if (x&s) i <<= 1;
        if (y&s) i <<= 2;
        if (z&s) i <<= 4;
        if (ptr->sol&i) return(1);
        if ((!(ptr->chi&i)) || (s <= 1)) return(0);
        i = popcount8(ptr->chi&(i-1)) + ptr->ind; //2child
    }
}

//returns: -1:invalid, {0..sur.mal-1}:surface index
int oct_getsurf (oct_t *loct, int x, int y, int z)
{
    octv_t *ptr;
    int i, s;

    if ((x|y|z)&loct->nsid) return(-1);
    for(i=loct->head,s=(loct->sid>>1);1;s>>=1)
    {
        ptr = &((octv_t *)loct->nod.buf)[i]; //2child
        i = 1;
        if (x&s) i <<= 1;
        if (y&s) i <<= 2;
        if (z&s) i <<= 4;
        if (!(ptr->chi&i)) return(-1);
        i = popcount8(ptr->chi&(i-1)) + ptr->ind; //2child
        if (s <= 1) return(i);
    }
}

#if 0
static int oct_isboxanyair (oct_t *loct, int x0, int y0, int z0, int s, octv_t *ptr, int bx, int by, int bz, int bs)
{
    octv_t *octv;
    int x, y, z, v;

    v = ptr->sol^((1<<8)-1);
    x = x0+s-bx; if (x <= 0) v &=~0x55; else if (x >= bs) v &=~0xaa;
    y = y0+s-by; if (y <= 0) v &=~0x33; else if (y >= bs) v &=~0xcc;
    z = z0+s-bz; if (z <= 0) v &=~0x0f; else if (z >= bs) v &=~0xf0;
    if (v&~ptr->chi) return(1);
    if (s == 1) return(0);
    octv = (octv_t *)loct->nod.buf;
    if ((v&(1<<0)) && (oct_isboxanyair(loct,x0  ,y0  ,z0  ,s>>1,&octv[popcount8(ptr->chi&((1<<0)-1))+ptr->ind],bx,by,bz,bs))) return(1);
    if ((v&(1<<1)) && (oct_isboxanyair(loct,x0+s,y0  ,z0  ,s>>1,&octv[popcount8(ptr->chi&((1<<1)-1))+ptr->ind],bx,by,bz,bs))) return(1);
    if ((v&(1<<2)) && (oct_isboxanyair(loct,x0  ,y0+s,z0  ,s>>1,&octv[popcount8(ptr->chi&((1<<2)-1))+ptr->ind],bx,by,bz,bs))) return(1);
    if ((v&(1<<3)) && (oct_isboxanyair(loct,x0+s,y0+s,z0  ,s>>1,&octv[popcount8(ptr->chi&((1<<3)-1))+ptr->ind],bx,by,bz,bs))) return(1);
    if ((v&(1<<4)) && (oct_isboxanyair(loct,x0  ,y0  ,z0+s,s>>1,&octv[popcount8(ptr->chi&((1<<4)-1))+ptr->ind],bx,by,bz,bs))) return(1);
    if ((v&(1<<5)) && (oct_isboxanyair(loct,x0+s,y0  ,z0+s,s>>1,&octv[popcount8(ptr->chi&((1<<5)-1))+ptr->ind],bx,by,bz,bs))) return(1);
    if ((v&(1<<6)) && (oct_isboxanyair(loct,x0  ,y0+s,z0+s,s>>1,&octv[popcount8(ptr->chi&((1<<6)-1))+ptr->ind],bx,by,bz,bs))) return(1);
    if ((v&(1<<7)) && (oct_isboxanyair(loct,x0+s,y0+s,z0+s,s>>1,&octv[popcount8(ptr->chi&((1<<7)-1))+ptr->ind],bx,by,bz,bs))) return(1);
    return(0);
}
#endif

//Function assumes (x,y,z,s) already known to not be pure air
//returns 1 if any voxels in cube (x,y,z,s) or its 1 voxel neighbors are air
static int oct_issurf (oct_t *loct, int x, int y, int z, int ls, oct_getvox_hint_t *och)
{
    int xx, yy, s, e;

    s = pow2[ls]; e = loct->sid-s;

    if ((!(loct->edgeissol& 1)) && (x == 0)) return(1);
    if ((!(loct->edgeissol& 2)) && (y == 0)) return(1);
    if ((!(loct->edgeissol& 4)) && (z == 0)) return(1);
    if ((!(loct->edgeissol& 8)) && (x == e)) return(1);
    if ((!(loct->edgeissol&16)) && (y == e)) return(1);
    if ((!(loct->edgeissol&32)) && (z == e)) return(1);

        //test cube; if interior's not pure solid, it must contain surfs
    if ((!oct_getsol_hint(loct,x,y,z,och)) || (och->minls < ls)) return(1); //NOTE:och->mins must be compared AFTER oct_getsol_hint() call!

        //Fast&elegant algo! :)
    if (x != 0) { xx = 0; yy = 0; do { if (!oct_getsol_hint(loct,x- 1,y+xx,z+yy,och)) return(1); mort2add(och->mins,xx,yy); } while (xx < s); }
    if (x != e) { xx = 0; yy = 0; do { if (!oct_getsol_hint(loct,x+ s,y+xx,z+yy,och)) return(1); mort2add(och->mins,xx,yy); } while (xx < s); }
    if (y != 0) { xx = 0; yy = 0; do { if (!oct_getsol_hint(loct,x+xx,y- 1,z+yy,och)) return(1); mort2add(och->mins,xx,yy); } while (xx < s); }
    if (y != e) { xx = 0; yy = 0; do { if (!oct_getsol_hint(loct,x+xx,y+ s,z+yy,och)) return(1); mort2add(och->mins,xx,yy); } while (xx < s); }
    if (z != 0) { xx = 0; yy = 0; do { if (!oct_getsol_hint(loct,x+xx,y+yy,z- 1,och)) return(1); mort2add(och->mins,xx,yy); } while (xx < s); }
    if (z != e) { xx = 0; yy = 0; do { if (!oct_getsol_hint(loct,x+xx,y+yy,z+ s,och)) return(1); mort2add(och->mins,xx,yy); } while (xx < s); }

    return(0);
}

int oct_findsurfdowny (oct_t *loct, int ox, int oy, int oz, surf_t **surf)
{
    typedef struct { octv_t *ptr; int x, y, z, j; } stk_t;
    stk_t stk[OCT_MAXLS];
    octv_t *ptr;
    int i, j, ls, x, y, z, nx, ny, nz;

    ls = loct->lsid-1; ptr = &((octv_t *)loct->nod.buf)[loct->head]; x = 0; y = 0; z = 0; j = 0;
    while (1)
    {
        i = (1<<j); if (!(ptr->chi&i)) goto tosibly;

        nx = (( j    &1)<<ls)+x; if ((nx+(1<<ls)-1 < ox) || (nx > ox)) goto tosibly;
        ny = (((j>>1)&1)<<ls)+y; if (ny+(1<<ls)-1 <= oy) goto tosibly;
        nz = (((j>>2)&1)<<ls)+z; if ((nz+(1<<ls)-1 < oz) || (nz > oz)) goto tosibly;

        if (ls <= 0) { (*surf) = &((surf_t *)loct->sur.buf)[popcount8(ptr->chi&(i-1)) + ptr->ind]; return(ny); }

        stk[ls].ptr = ptr; stk[ls].x = x; stk[ls].y = y; stk[ls].z = z; stk[ls].j = j; ls--; //2child
        ptr = &((octv_t *)loct->nod.buf)[popcount8(ptr->chi&(i-1)) + ptr->ind]; x = nx; y = ny; z = nz; j = 0;
        continue;

tosibly:;
        j++; if (j < 8) continue;
        do { ls++; if (ls >= loct->lsid) return(loct->sid); j = stk[ls].j+1; } while (j >= 8); //2parent
        ptr = stk[ls].ptr; x = stk[ls].x; y = stk[ls].y; z = stk[ls].z;
    }
}

void oct_setvox (oct_t *loct, int x, int y, int z, surf_t *surf, int mode)
{
    brush_vox_t vox;
    brush_vox_init(&vox,x,y,z,surf);
    oct_mod(loct,(brush_t *)&vox,mode);
}

//Generates bit cube with bx0,by0,bz0 as top-left-front
//lsdax: least significant dimension axis (-=mirrored): 1=x, 2=y, 3=z, -1=-x, -2=-y, -3=-z
//NOTE! It is assumed that bitvis allocation x size is rounded up to next multiple of 32
void oct_sol2bit (oct_t *loct, unsigned int *bitvis, int bx0, int by0, int bz0, int dx, int dy, int dz, int lsdax)
{
    typedef struct { octv_t *ptr; int j; } stk_t;
    stk_t stk[OCT_MAXLS];
    octv_t *ptr;
    int i, j, ls, s, x, y, z, nx, ny, nz, x0, y0, z0, x1, y1, z1, xx, yy, zz, pit, isneg;

        //obtain bit array to calculate visibility info quickly
    ls = loct->lsid-1; s = (1<<ls); ptr = &((octv_t *)loct->nod.buf)[loct->head];
    x = -bx0; y = -by0; z = -bz0; j = 8-1;
    if (lsdax < 0) { isneg = 1; lsdax = -lsdax; } else isneg = 0;
    switch(lsdax)
    {
        case 1: pit = ((dx+31)>>5); memset(bitvis,0,pit*dy*dz*sizeof(bitvis[0])); break;
        case 2: pit = ((dy+31)>>5); memset(bitvis,0,dx*pit*dz*sizeof(bitvis[0])); break;
        case 3: pit = ((dz+31)>>5); memset(bitvis,0,dx*dy*pit*sizeof(bitvis[0])); break;
    }
    while (1)
    {
        nx = x; if (j&1) nx += s; if ((nx >= dx) || (nx+s <= 0)) goto tosibly;
        ny = y; if (j&2) ny += s; if ((ny >= dy) || (ny+s <= 0)) goto tosibly;
        nz = z; if (j&4) nz += s; if ((nz >= dz) || (nz+s <= 0)) goto tosibly;

        i = (1<<j);
        if (ptr->sol&i)
        {
            switch (lsdax)
            {
                case 1:
                    if (isneg) nx = dx-s-nx;
                    if (dx <= 32)
                    {
                        if (!ls) bitvis[nz*dy+ny] += (1<<nx);
                        else
                        {
                            x0 = std::max(nx,0); x1 = std::min(nx+s,dx); x0 = (-1<<x0)&~(-2<<(x1-1)); //NOTE:..&~(-1<<x1) gives wrong result when x1==32
                            y0 = std::max(ny,0); y1 = std::min(ny+s,dy);
                            z0 = std::max(nz,0); z1 = std::min(nz+s,dz);
                            for(zz=z0;zz<z1;zz++) for(yy=y0;yy<y1;yy++) bitvis[zz*dy+yy] += x0;
                        }
                    }
                    else
                    {
                        if (!ls) bitvis[(nz*dy+ny)*pit + (nx>>5)] += (1<<nx);
                        else
                        {
                            x0 = std::max(nx,0); x1 = std::min(nx+s,dx);
                            y0 = std::max(ny,0); y1 = std::min(ny+s,dy);
                            z0 = std::max(nz,0); z1 = std::min(nz+s,dz);
                            for(zz=z0;zz<z1;zz++) for(yy=y0;yy<y1;yy++) setzrange1(&bitvis[(zz*dy+yy)*pit],x0,x1);
                        }
                    }
                    break;
                case 2:
                    if (isneg) ny = dy-s-ny;
                    if (dy <= 32)
                    {
                        if (!ls) bitvis[nz*dx+nx] += (1<<ny);
                        else
                        {
                            x0 = std::max(nx,0); x1 = std::min(nx+s,dx);
                            y0 = std::max(ny,0); y1 = std::min(ny+s,dy); y0 = (-1<<y0)&~(-2<<(y1-1)); //NOTE:..&~(-1<<y1) gives wrong result when y1==32
                            z0 = std::max(nz,0); z1 = std::min(nz+s,dz);
                            for(zz=z0;zz<z1;zz++) for(xx=x0;xx<x1;xx++) bitvis[zz*dx+xx] += y0;
                        }
                    }
                    else
                    {
                        if (!ls) bitvis[(nz*dx+nx)*pit + (ny>>5)] += (1<<ny);
                        else
                        {
                            x0 = std::max(nx,0); x1 = std::min(nx+s,dx);
                            y0 = std::max(ny,0); y1 = std::min(ny+s,dy);
                            z0 = std::max(nz,0); z1 = std::min(nz+s,dz);
                            for(zz=z0;zz<z1;zz++) for(xx=x0;xx<x1;xx++) setzrange1(&bitvis[(zz*dx+xx)*pit],y0,y1);
                        }
                    }
                    break;
                case 3:
                    if (isneg) nz = dz-s-nz;
                    if (dz <= 32)
                    {
                        if (!ls) bitvis[ny*dx+nx] += (1<<nz);
                        else
                        {
                            x0 = std::max(nx,0); x1 = std::min(nx+s,dx);
                            y0 = std::max(ny,0); y1 = std::min(ny+s,dy);
                            z0 = std::max(nz,0); z1 = std::min(nz+s,dz); z0 = (-1<<z0)&~(-2<<(z1-1)); //NOTE:..&~(-1<<z1) gives wrong result when z1==32
                            for(yy=y0;yy<y1;yy++) for(xx=x0;xx<x1;xx++) bitvis[yy*dx+xx] += z0;
                        }
                    }
                    else
                    {
                        if (!ls) bitvis[(ny*dx+nx)*pit + (nz>>5)] += (1<<nz);
                        else
                        {
                            x0 = std::max(nx,0); x1 = std::min(nx+s,dx);
                            y0 = std::max(ny,0); y1 = std::min(ny+s,dy);
                            z0 = std::max(nz,0); z1 = std::min(nz+s,dz);
                            for(yy=y0;yy<y1;yy++) for(xx=x0;xx<x1;xx++) setzrange1(&bitvis[(yy*dx+xx)*pit],z0,z1);
                        }
                    }
                    break;
            }
            goto tosibly;
        }

        if (ptr->chi&i) //Recurse only if mixture of air&sol
        {
            stk[ls].ptr = ptr; stk[ls].j = j; ls--; s >>= 1; //2child
            ptr = &((octv_t *)loct->nod.buf)[popcount8(ptr->chi&(i-1)) + ptr->ind]; x = nx; y = ny; z = nz; j = 8-1;
            continue;
        }

tosibly:;
        j--; if (j >= 0) continue;
        do { ls++; s <<= 1; if (ls >= loct->lsid) return; j = stk[ls].j-1; } while (j < 0); //2parent
        ptr = stk[ls].ptr; i = s*-2;
        x = ((x+bx0)&i)-bx0;
        y = ((y+by0)&i)-by0;
        z = ((z+bz0)&i)-bz0;
    }
}

#if defined(ENABLE_UNUSED_ROUTINES)
static void oct_surf2bit (oct_t *loct, unsigned int *bitsurf, int bx0, int by0, int bz0, int dx, int dy, int dz)
{
    typedef struct { octv_t *ptr; int x, y, z, j; } stk_t;
    stk_t stk[OCT_MAXLS];
    octv_t *ptr;
    int i, j, ls, s, x, y, z, nx, ny, nz, pit;

        //obtain bit array to calculate visibility info quickly
    ls = loct->lsid-1; s = (1<<ls); ptr = &((octv_t *)loct->nod.buf)[loct->head];
    x = -bx0; y = -by0; z = -bz0; j = 8-1;
    pit = ((dx+31)>>5); memset(bitsurf,0,pit*dy*dz*sizeof(bitsurf[0]));
    while (1)
    {
        nx = (( j    &1)<<ls)+x; if ((nx >= dx) || (nx+s <= 0)) goto tosibly;
        ny = (((j>>1)&1)<<ls)+y; if ((ny >= dy) || (ny+s <= 0)) goto tosibly;
        nz = (((j>>2)&1)<<ls)+z; if ((nz >= dz) || (nz+s <= 0)) goto tosibly;

        i = (1<<j);

        if (ptr->chi&i) //Recurse only if surface
        {
            if (ls)
            {
                stk[ls].ptr = ptr; stk[ls].x = x; stk[ls].y = y; stk[ls].z = z; stk[ls].j = j; ls--; s >>= 1; //2child
                ptr = &((octv_t *)loct->nod.buf)[popcount8(ptr->chi&(i-1)) + ptr->ind]; x = nx; y = ny; z = nz; j = 8-1;
                continue;
            }
            bitsurf[(nz*dy + ny)*pit + (nx>>5)] |= (1<<nx);
        }

tosibly:;
        j--; if (j >= 0) continue;
        do { ls++; s <<= 1; if (ls >= loct->lsid) return; j = stk[ls].j-1; } while (j < 0); //2parent
        ptr = stk[ls].ptr; x = stk[ls].x; y = stk[ls].y; z = stk[ls].z;
    }
}
#endif

static int oct_touch_brush_rec (oct_t *loct, octv_t *inode, int x0, int y0, int z0, int ls, brush_t *brush)
{
    int i, j, x, y, z, s, ind;

    ls--; s = (1<<ls); ind = inode->ind;
    for(i=1,z=0;z<=s;z+=s)
        for(y=0;y<=s;y+=s)
            for(x=0;x<=s;x+=s,i<<=1)
            {
                if (!((inode->chi|inode->sol)&i)) continue;
                j = brush->isins(brush,x+x0,y+y0,z+z0,ls);
                if (inode->chi&i)       //FIX:this works because at least 1 node approaches small size.
                //if (!(inode->sol&i))  //FIX:this fails due to approximation of bbox vs. bbox. (see quadtre3.c:see oct_touch_oct_recur() for potential solution)
                {
                    if (j == 2) return(1);
                    if (j == 1)
                    {
                        if (s <= 1) return(1); //NOTE:for safety only - a properly written brush->isins function should never return 1 for ls==0
                        if (oct_touch_brush_rec(loct,&((octv_t *)loct->nod.buf)[ind],x+x0,y+y0,z+z0,ls,brush)) return(1);
                    }
                    ind++;
                }
                else if (j) return(1);
            }
    return(0);
}
int oct_touch_brush (oct_t *loct, brush_t *brush) {
    return(oct_touch_brush_rec(loct,&((octv_t *)loct->nod.buf)[loct->head],0,0,0,loct->lsid,brush));
}

#define COMPACT_LBLKSIZ 5 //# ints per block (4 r's: (3:1.43ms, 4:1.63ms, 5:1.94ms, 6:2.51ms, 7:3.68ms, 8:6.15ms, 9:10.69ms, 10:20.04ms))
static void oct_compact (oct_t *loct, int issur) //issur==0:nod, issur==1:sur
{
    bitmal_t *bm;
    int i, j, k, ls, ptr, optr, num0, *popcountlut, stkptr[OCT_MAXLS], stknum[OCT_MAXLS];

    if (!((octv_t *)loct->nod.buf)[loct->head].chi) return; //don't compact empty oct_t (needed to prevent crash in loop later)
    if (!issur) bm = &loct->nod; else bm = &loct->sur;

    popcountlut = (int *)malloc(((bm->mal>>(COMPACT_LBLKSIZ+5))+1)*sizeof(int));
    for(i=0,k=0;i<(bm->mal>>5);i+=(1<<COMPACT_LBLKSIZ))
    {
        popcountlut[i>>COMPACT_LBLKSIZ] = k;
        //for(j=(1<<COMPACT_LBLKSIZ)-1;j>=0;j--) k += popcount32(bm->bit[i+j]);
        for(j=(1<<COMPACT_LBLKSIZ)-4;j>=0;j-=4) k += popcount128(&bm->bit[i+j]);
    }

        //Subtract appropriate amount (# holes preceding address) from all ptr's
    ls = loct->lsid-1; stkptr[ls] = loct->head; stknum[ls] = 1;
    while (1)
    {
        ptr = stkptr[ls]; stkptr[ls]++; stknum[ls]--; //2sibly
        if (ls >= 0)
        {
            optr = ((octv_t *)loct->nod.buf)[ptr].ind;
            if ((!issur) == (ls != 0))
            {
                    //popcount8; uses LUT value every (1<<COMPACT_LBLKSIZ) ints
                i = optr>>5; j = i&(-1<<COMPACT_LBLKSIZ); num0 = popcount32(bm->bit[i]|(-1<<optr)) + popcountlut[j>>COMPACT_LBLKSIZ];
                while (i&3) { i -= 1; num0 += popcount32 ( bm->bit[i]); } //Even faster; requires SSE2
                while (i>j) { i -= 4; num0 += popcount128(&bm->bit[i]); }
                num0 = ((optr&~31)+32) - num0;

                ((octv_t *)loct->nod.buf)[ptr].ind = optr-num0;
            }
            if (ls > 0) { ls--; stkptr[ls] = optr; stknum[ls] = popcount8(((octv_t *)loct->nod.buf)[ptr].chi); } //2child
        }
        while (stknum[ls] <= 0) { ls++; if (ls >= loct->lsid) goto break2; } //2parent
    }
break2:;

        //shift all nodes towards index 0
    if (!issur)
    {
        for(k=0,i=dntil1(bm->bit,0,bm->mal);i<bm->mal;k+=j-i,i=dntil1(bm->bit,j+1,bm->mal)) //faster
            { j = dntil0(bm->bit,i+1,bm->mal); memmove(&((octv_t *)bm->buf)[k],&((octv_t *)bm->buf)[i],(j-i)*bm->siz); }
    }
    else
    {
        for(k=0,i=dntil1(bm->bit,0,bm->mal);i<bm->mal;k+=j-i,i=dntil1(bm->bit,j+1,bm->mal)) //faster
            { j = dntil0(bm->bit,i+1,bm->mal); memmove(&((surf_t *)bm->buf)[k],&((surf_t *)bm->buf)[i],(j-i)*bm->siz); }
    }

    setzrange1(bm->bit,      0,bm->num);
    setzrange0(bm->bit,bm->num,bm->mal);
    bm->ind = bm->num;

    free(popcountlut);

    if constexpr (oct_usegpu) {
        if (issur) {
            if constexpr (oct_usegpubo) {
                if (!loct->gsurf) {
                    loct->gsurf = (surf_t *)bo_begin(loct->bufid,0);
                }
                memcpy(loct->gsurf,loct->sur.buf,loct->sur.num*loct->sur.siz);
            }
            else {
                ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_2D,loct->octid);
                ((PFNGLTEXSUBIMAGE2D)glfp[glTexSubImage2D])(GL_TEXTURE_2D,0,0,0,loct->gxsid,(loct->sur.num*2+loct->gxsid-1)>>loct->glxsid,GL_RGBA,GL_UNSIGNED_BYTE,(void *)loct->sur.buf);
            }
        }
    }
}

void oct_checkreducesizes (oct_t *loct)
{
    if ((loct->nod.num < (loct->nod.mal>>2)) && (loct->nod.mal > 256))
    {
        oct_compact(loct,0);
        loct->nod.mal >>= 1; //halve space
        loct->nod.buf = (octv_t       *)realloc(loct->nod.buf,   (int64_t)loct->nod.mal*loct->nod.siz);
        loct->nod.bit = (unsigned int *)realloc(loct->nod.bit,((((int64_t)loct->nod.mal+63)>>5)<<2)+16);
    }
    if ((loct->sur.num < (loct->sur.mal>>2)) && (loct->sur.mal > 256))
    {
        oct_compact(loct,1);
        loct->sur.mal >>= 1; //halve space
        loct->sur.buf = (octv_t       *)realloc(loct->sur.buf,   (int64_t)loct->sur.mal*loct->sur.siz);
        loct->sur.bit = (unsigned int *)realloc(loct->sur.bit,((((int64_t)loct->sur.mal+63)>>5)<<2)+16);
        if constexpr (oct_usegpu)
        {
            if constexpr (oct_usegpubo) {
                if (loct->gsurf)
                {
                    loct->gsurf = 0;
                    ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_2D,loct->octid);
                    bo_end(loct->bufid,0,0,0,0,GL_RGBA,GL_UNSIGNED_BYTE,0);
                }
                ((PFNGLDELETEBUFFERS)glfp[glDeleteBuffers])(1,&loct->bufid);
            }

            if (loct->glxsid > loct->glysid) loct->glxsid--; else loct->glysid--;
            loct->gxsid = (1<<loct->glxsid);
            loct->gysid = (1<<loct->glysid);
            loct->sur.mal = loct->gxsid*loct->gysid/(sizeof(surf_t)>>2);
            kglalloctex(loct->octid,0,loct->gxsid,loct->gysid,1,KGL_RGBA32+KGL_NEAREST); //only NEAREST makes sense here!

            if constexpr (oct_usegpubo) {
                loct->bufid = bo_init(loct->gxsid*loct->gysid*4);
                loct->gsurf = (surf_t *)bo_begin(loct->bufid,0);
                memcpy(loct->gsurf,loct->sur.buf,loct->sur.num*loct->sur.siz);
            }
            else {
                ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_2D,loct->octid);
                ((PFNGLTEXSUBIMAGE2D)glfp[glTexSubImage2D])(GL_TEXTURE_2D,0,0,0,loct->gxsid,(loct->sur.num*2+loct->gxsid-1)>>loct->glxsid,GL_RGBA,GL_UNSIGNED_BYTE,(void *)loct->sur.buf);
            }
        }
    }
}

//Handle case at edges of defined octree space
static int isins_func (oct_t *loct, brush_t *brush, int x0, int y0, int z0, int ls)
{
if ((x0|y0|z0)&loct->nsid) return(0);
return(brush->isins(brush,x0,y0,z0,ls));
}

static void getsurf_func (oct_t *loct, brush_t *brush, int x0, int y0, int z0, surf_t *surf)
{
brush->getsurf(brush,x0,y0,z0,surf);
if (!(brush->flags&1)) return;
if (x0 ==           0) { surf->norm[0] =+64; surf->norm[1] =  0; surf->norm[2] =  0; }
if (y0 ==           0) { surf->norm[0] =  0; surf->norm[1] =+64; surf->norm[2] =  0; }
if (z0 ==           0) { surf->norm[0] =  0; surf->norm[1] =  0; surf->norm[2] =+64; }
if (x0 == loct->sid-1) { surf->norm[0] =-64; surf->norm[1] =  0; surf->norm[2] =  0; }
if (y0 == loct->sid-1) { surf->norm[0] =  0; surf->norm[1] =-64; surf->norm[2] =  0; }
if (z0 == loct->sid-1) { surf->norm[0] =  0; surf->norm[1] =  0; surf->norm[2] =-64; }
}

static void oct_dealloctree (oct_t *loct, int ls, int ind, int chi)
{
int nbits;
nbits = popcount8(chi);
if (ls <= 0) { xorzrangesmall(loct->sur.bit,ind,nbits); loct->sur.num -= nbits; return; }
xorzrangesmall(loct->nod.bit,ind,nbits); loct->nod.num -= nbits;
for(nbits--;nbits>=0;nbits--,ind++) oct_dealloctree(loct,ls-1,((octv_t *)loct->nod.buf)[ind].ind,((octv_t *)loct->nod.buf)[ind].chi);
}

//--------------------------------------------------------------------------------------------------
//02/22/2012 floodfill algo:
//for (each unmarked solid in bbox)
//{
//   while (floodfill)
//   {
//      if ((mrk2) || (hit bottom)) { oct_copymrk2mrk2(); exit_early; }
//      if (mrk) continue;
//      mrk = 1;
//      write_fif(x,y,z,ls);
//   }
//   if (didn't exit_early)
//   {
//      zot or gen sprite;
//      oct_removemrks();
//   }
//}
//oct_clearmrk2s();

typedef struct { unsigned short x, y, z; char ls, skipdir; } ffif_t;
#define FFIFMAX 65536 //FIX:make dynamic!
static ffif_t ffif[FFIFMAX];

//returns: 0:not hover (connected to ground/out of mem/invalid nuc), 1:is hover
static int oct_floodfill (oct_t *loct, int nx, int ny, int nz)
{
    oct_getvox_hint_t och;
    octv_t *ptr;
    static const int wwadd[3][3] = {2,0,1, 0,2,1, 0,1,2};
    const int *pwwadd;
    int i, c, s, ls, x, y, z, skipdir, ww[3], ffifn;
    ffif_t *pfif;

    if ((nx|ny|nz)&loct->nsid) return(0); //out of bounds
    oct_getvox_hint_init(loct,&och);
    ffifn = 0; ww[0] = 0; ww[1] = 0; s = 0; c = -1; skipdir = -1; goto in2it_cpp;
    do
    {
        ffifn--; pfif = &ffif[ffifn]; x = pfif->x; y = pfif->y; z = pfif->z; s = pow2[pfif->ls]; skipdir = pfif->skipdir;

        for(c=6-1;c>=0;c--) //visit all neighbors of octree node (x,y,s)
        {
            if (skipdir == c) continue;
            ww[0] = 0; ww[1] = 0; ww[2] = ((c&1)-1)|s; pwwadd = &wwadd[c>>1][0];
            nx = ww[pwwadd[0]]+x;
            ny = ww[pwwadd[1]]+y;
            nz = ww[pwwadd[2]]+z; if ((nx|ny|nz)&loct->nsid) continue;
            do
            {
in2it_cpp:;     i = oct_getsol_hint(loct,nx,ny,nz,&och);
                if (i) //(nx&(-och.mins),ny&(-och.mins),nz&(-och.mins),och.mins) is neighbor
                {
                    ptr = &((octv_t *)loct->nod.buf)[och.stkind[och.minls]];
                    if (ptr->mrk2&i) return(0); //if hit previous grounded floodfill, exit
                    if (!(ptr->mrk&i))
                    {
                        if (ny+och.mins >= loct->sid) return(0); //hit bottom; exit early

                            //add neighbor to local fif
                        if (ffifn >= FFIFMAX) return(0); //out of mem :/
                        pfif = &ffif[ffifn]; pfif->x = nx&(-och.mins); pfif->y = ny&(-och.mins); pfif->z = nz&(-och.mins);
                        pfif->ls = och.minls; pfif->skipdir = (c^1)|((s-och.mins)>>31); ffifn++;

                        ptr->mrk |= i; //mark node
                        for(ls=och.minls+1;ls<loct->lsid;ls++)
                        {
                            ptr = &((octv_t *)loct->nod.buf)[och.stkind[ls]];
                            i = 1;
                            if (nx&(1<<ls)) i <<= 1;
                            if (ny&(1<<ls)) i <<= 2;
                            if (nz&(1<<ls)) i <<= 4;
                            if (ptr->mrk&i) break; ptr->mrk |= i;
                        }
                    }
                }

                mort2add(och.mins,ww[0],ww[1]); if (ww[0] >= s) break;
                nx = ww[pwwadd[0]]+x;
                ny = ww[pwwadd[1]]+y;
                nz = ww[pwwadd[2]]+z;
            } while (1);
        }
    } while (ffifn > 0);
    return(1);
}

static void oct_copymrk2mrk2 (oct_t *loct, int inode, int ls) //copy mrk to mrk2
{
    octv_t *ptr;
    int iup, v;

    ptr = &((octv_t *)loct->nod.buf)[inode];
    if (ls)
    {
        for(v=ptr->chi&ptr->mrk;v;v^=iup)
        {
            iup = (-v)&v;
            oct_copymrk2mrk2(loct,popcount8(ptr->chi&(iup-1))+ptr->ind,ls-1);
        }
    }
    ptr->mrk2 |= ptr->mrk; ptr->mrk = 0;
}

static void oct_removemrks (oct_t *loct, int inode, int ls, octv_t *roct)
{
    surf_t surf[8];
    octv_t *ooct, noct[8];
    int iup, n, o, v;

    ooct = &((octv_t *)loct->nod.buf)[inode];
    roct->sol  = ooct->sol &~ooct->mrk;
    roct->chi  = ooct->chi &~ooct->mrk;
    roct->mrk  = 0;
    roct->mrk2 = ooct->mrk2&~ooct->mrk;
    o = ooct->ind; n = 0;
    if (!ls)
    {
        for(v=ooct->chi;v;v^=iup,o++) //visit only nodes that may differ from brush color
        {
            iup = (-v)&v;
            if (!(ooct->mrk&iup)) { memcpy(&surf[n],&((surf_t *)loct->sur.buf)[o],loct->sur.siz); n++; continue; } //no intersect:copy
        }
        o = popcount8(ooct->chi); roct->ind = ooct->ind;
        xorzrangesmall(loct->sur.bit,ooct->ind+n,o-n); //shorten node
        memcpy(&((surf_t *)loct->sur.buf)[ooct->ind],surf,n*loct->sur.siz);
        loct->sur.num += n-o;
    }
    else
    {
        for(v=ooct->chi;v;v^=iup,o++) //visit only nodes that may differ from brush color
        {
            iup = (-v)&v;
            if (!(ooct->mrk&iup)) { noct[n] = ((octv_t *)loct->nod.buf)[o]; n++; continue; } //no intersect:copy
            if (ooct->sol&iup) { oct_dealloctree(loct,ls-1,((octv_t *)loct->nod.buf)[o].ind,((octv_t *)loct->nod.buf)[o].chi); continue; } //clear all
            oct_removemrks(loct,o,ls-1,&noct[n]); //intersects partially:recurse
            if (noct[n].chi || ((unsigned)(noct[n].sol-1) < (unsigned)((1<<8)-2))) { roct->chi += iup; roct->mrk2 += (ooct->mrk2&iup); n++; }
        }
        o = popcount8(ooct->chi); roct->ind = ooct->ind;
        xorzrangesmall(loct->nod.bit,ooct->ind+n,o-n); //shorten node
        memcpy(&((octv_t *)loct->nod.buf)[ooct->ind],noct,n*loct->nod.siz);
        loct->nod.num += n-o;
    }
}

//copies marked (mrk) sections of loct to newoct; newoct must be fresh from an oct_new()
static void oct_mark2spr (oct_t *loct, oct_t *newoct, int inode, int ls, octv_t *roct, int mskor)
{
    surf_t surf[8];
    octv_t *ooct, noct[8];
    int iup, n, o, v, omrk;

    ooct = &((octv_t *)loct->nod.buf)[inode]; omrk = ooct->mrk|mskor;
    roct->sol = ooct->sol&omrk; roct->chi = 0; roct->mrk = 0; roct->mrk2 = 0;
    n = 0;
    if (!ls)
    {
        for(v=(ooct->chi&omrk);v;v^=iup) //visit only nodes that may differ from brush color
        {
            iup = (-v)&v; o = popcount8(ooct->chi&(iup-1)) + ooct->ind;
            memcpy(&surf[n],&((surf_t *)loct->sur.buf)[o],loct->sur.siz);
            roct->chi += iup; n++;
        }

        if (!n) {
            roct->ind = -1;
            return;
        }

        if (newoct->sur.num+n > newoct->sur.mal) {
            if constexpr (oct_usegpu) {
                //grow space by 100%
                if (newoct->glxsid < newoct->glysid) {
                    newoct->glxsid++; newoct->gxsid <<= 1;
                }
                else {
                    newoct->glysid++; newoct->gysid <<= 1;
                }
                newoct->sur.mal <<= 1;
            }
            else {
                newoct->sur.mal = std::max(((1<<bsr(newoct->sur.mal))>>2) + newoct->sur.mal,newoct->sur.num+n); //grow space by ~25%
            }

            newoct->sur.buf = (octv_t       *)realloc(newoct->sur.buf,   (int64_t)newoct->sur.mal*newoct->sur.siz);
        }
        roct->ind = newoct->sur.num; newoct->sur.num += n; //simple allocator (don't need bitalloc since no deallocation)
        memcpy(&((surf_t *)newoct->sur.buf)[roct->ind],surf,n*newoct->sur.siz);
    }
    else
    {
        for(v=(ooct->chi&omrk);v;v^=iup) //visit only nodes that may differ from brush color
        {
            iup = (-v)&v; o = popcount8(ooct->chi&(iup-1)) + ooct->ind;
            oct_mark2spr(loct,newoct,o,ls-1,&noct[n],(-(ooct->sol&iup))>>8); //intersects partially:recurse
            if (noct[n].chi || ((unsigned)(noct[n].sol-1) < (unsigned)((1<<8)-2))) { roct->chi += iup; n++; }
        }

        if (!n) {
            roct->ind = -1;
            return;
        }

        if (newoct->nod.num+n > newoct->nod.mal)
        {
            newoct->nod.mal = std::max(((1<<bsr(newoct->nod.mal))>>2) + newoct->nod.mal,newoct->nod.num+n); //grow space by ~25%
            newoct->nod.buf = (octv_t       *)realloc(newoct->nod.buf,   (int64_t)newoct->nod.mal*newoct->nod.siz);
        }

        roct->ind = newoct->nod.num; newoct->nod.num += n; //simple allocator (don't need bitalloc since no deallocation)
        memcpy(&((octv_t *)newoct->nod.buf)[roct->ind],noct,n*newoct->nod.siz);
    }
}

static void oct_clearmrk2s (oct_t *loct, int inode, int ls)
{
    octv_t *ptr;
    int iup, v;

    ptr = &((octv_t *)loct->nod.buf)[inode];
    if (ls)
    {
        for(v=ptr->chi&ptr->mrk2;v;v^=iup)
        {
            iup = (-v)&v;
            oct_clearmrk2s(loct,popcount8(ptr->chi&(iup-1))+ptr->ind,ls-1);
        }
    }
    ptr->mrk2 = 0;
}

void oct_hover_check (oct_t *loct, int x0, int y0, int z0, int x1, int y1, int z1, void (*recvoctfunc)(oct_t *ooct, oct_t *noct))
{
    typedef struct { int ind, j; } stk_t;
    stk_t stk[OCT_MAXLS];
    oct_t newoct;
    octv_t nnode, *ptr;
    int i, j, ls, s, x, y, z, nx, ny, nz, didx, didy, didz, ind;

    didx = 0; didy = 0; didz = 0;
ohc_restart:;
    ls = loct->lsid-1; s = (1<<ls); ind = loct->head; ptr = &((octv_t *)loct->nod.buf)[ind]; x = 0; y = 0; z = 0; j = 0; //WARNING:j must ascend for morton comparison to work!
    while (1)
    {
        i = (1<<j);
        nx = ((-((j   )&1))&s)+x; if ((nx > x1) || (nx+s < x0)) goto tosibly;
        ny = ((-((j>>1)&1))&s)+y; if ((ny > y1) || (ny+s < y0)) goto tosibly;
        nz = ((-((j>>2)&1))&s)+z; if ((nz > z1) || (nz+s < z0)) goto tosibly;

        if (ptr->sol&i)
        {
            if ((ptr->mrk|ptr->mrk2)&i) goto tosibly;
            if (mortcmp_fast(nx,ny,nz,didx,didy,didz) < 0) goto tosibly;

            if (!oct_floodfill(loct,nx,ny,nz)) //returns:0=hit bottom/out of mem, 1=is hover
            {
                oct_copymrk2mrk2(loct,loct->head,loct->lsid-1);
            }
            else
            {
                if (recvoctfunc)
                {
                    oct_new(loct->xres, loct->yres, &newoct,loct->lsid,loct->tilid,256,256,1);

                    oct_mark2spr(loct,&newoct,loct->head,loct->lsid-1,&nnode,0);
                    ((octv_t *)newoct.nod.buf)[newoct.head] = nnode; //can't write node directly to loct->nod.buf[loct->head] because of possible realloc inside

                    //if ((newoct.nod.num <= 1) && (newoct.sur.num <= 1)) MessageBox(ghwnd,"Blank sprite :/",prognam,MB_OK);

                    newoct.nod.bit = (unsigned int *)malloc(((((int64_t)newoct.nod.mal+63)>>5)<<2)+16);
                    newoct.nod.ind = newoct.nod.num;
                    setzrange1(newoct.nod.bit,             0,newoct.nod.num);
                    setzrange0(newoct.nod.bit,newoct.nod.num,newoct.nod.mal);

                    newoct.sur.bit = (unsigned int *)malloc(((((int64_t)newoct.sur.mal+63)>>5)<<2)+16);
                    newoct.sur.ind = newoct.sur.num;
                    setzrange1(newoct.sur.bit,             0,newoct.sur.num);
                    setzrange0(newoct.sur.bit,newoct.sur.num,newoct.sur.mal);

                    if constexpr (oct_usegpu)
                    {
                        kglalloctex(newoct.octid,0,newoct.gxsid,newoct.gysid,1,KGL_RGBA32+KGL_NEAREST); //only NEAREST makes sense here!

                        if constexpr (oct_usegpubo) {
                            newoct.bufid = bo_init(newoct.gxsid*newoct.gysid*4);
                            newoct.gsurf = (surf_t *)bo_begin(newoct.bufid,0);
                            memcpy(newoct.gsurf,newoct.sur.buf,newoct.sur.num*newoct.sur.siz);
                        }
                        else {
                            ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_2D,newoct.octid);
                            ((PFNGLTEXSUBIMAGE2D)glfp[glTexSubImage2D])(GL_TEXTURE_2D,0,0,0,newoct.gxsid,(newoct.sur.num*2+newoct.gxsid-1)>>newoct.glxsid,GL_RGBA,GL_UNSIGNED_BYTE,(void *)newoct.sur.buf);
                        }
                    }
                    recvoctfunc(loct,&newoct);
                }
                oct_removemrks(loct,loct->head,loct->lsid-1,&nnode);
                ((octv_t *)loct->nod.buf)[loct->head] = nnode; //can't write node directly to loct->nod.buf[loct->head] because of possible realloc inside

                didx = nx; didy = ny; didz = nz;
                goto ohc_restart;
            }
            goto tosibly;
        }
        if (ptr->chi&i)
        {
            stk[ls].ind = ind; stk[ls].j = j; ls--; s >>= 1; //2child
            ind = popcount8(ptr->chi&(i-1)) + ptr->ind; ptr = &((octv_t *)loct->nod.buf)[ind]; x = nx; y = ny; z = nz; j = 0;
            continue;
        }

tosibly:;
        j++; if (j < 8) continue;
        do { s <<= 1; ls++; if (ls >= loct->lsid) goto endit; j = stk[ls].j+1; } while (j >= 8); //2parent
        ind = stk[ls].ind; ptr = &((octv_t *)loct->nod.buf)[ind]; i = -(s<<1); x &= i; y &= i; z &= i;
    }
endit:;
    oct_clearmrk2s(loct,loct->head,loct->lsid-1);
    oct_checkreducesizes(loct);
}

//--------------------------------------------------------------------------------------------------
//      n: # structs to alloc
//returns: bit index or crash if realloc fails
static int bitalloc (oct_t *loct, bitmal_t *bm, int n)
{
    int i, oi, ie, i0, i1, cnt;

    #if defined(__x86_64__) || defined(_WIN64)
        int64_t j, k;
        uint64_t *bitbuf = (uint64_t *)bm->bit;
    #elif defined(__i386__) || defined(_WIN32)
        int j, k;
        unsigned int *bitbuf = bm->bit;
    #else
        #error
    #endif

    i = bm->ind; oi = i; ie = bm->mal-n;
    for(cnt=1;1;cnt--)
    {
        switch(n)
        {
            #if defined(__x86_64__) || defined(_WIN64)
                //NOTE:changing to 64-bit had absolutely no effect on initboard() speed :/
                case 1: ie >>= 6; for(i>>=6;i<ie;i++) { j = bitbuf[i];                                                                              if (j != -1ll) goto found; } break;
                case 2: ie >>= 6; for(i>>=6;i<ie;i++) { j = bitbuf[i]; j |= (j>>1);                           k = bitbuf[i+1]; j |= (      k <<63); if (j != -1ll) goto found; } break;
                case 3: ie >>= 6; for(i>>=6;i<ie;i++) { j = bitbuf[i]; j |= (j>>1)|(j>>2);                    k = bitbuf[i+1]; j |= (((-k)|k)<<62); if (j != -1ll) goto found; } break;
                case 4: ie >>= 6; for(i>>=6;i<ie;i++) { j = bitbuf[i]; j |= (j>>1); j |= (j>>2);              k = bitbuf[i+1]; j |= (((-k)|k)<<61); if (j != -1ll) goto found; } break;
                case 5: ie >>= 6; for(i>>=6;i<ie;i++) { j = bitbuf[i]; j |= (j>>1); j |= (j>>2)|(j>>3);       k = bitbuf[i+1]; j |= (((-k)|k)<<60); if (j != -1ll) goto found; } break;
                case 6: ie >>= 6; for(i>>=6;i<ie;i++) { j = bitbuf[i]; j |= (j>>1); j |= (j>>2)|(j>>4);       k = bitbuf[i+1]; j |= (((-k)|k)<<59); if (j != -1ll) goto found; } break;
                case 7: ie >>= 6; for(i>>=6;i<ie;i++) { j = bitbuf[i]; j |= (j>>1); j |= (j>>2); j |= (j>>3); k = bitbuf[i+1]; j |= (((-k)|k)<<58); if (j != -1ll) goto found; } break;
                case 8: ie >>= 6; for(i>>=6;i<ie;i++) { j = bitbuf[i]; j |= (j>>1); j |= (j>>2); j |= (j>>4); k = bitbuf[i+1]; j |= (((-k)|k)<<57); if (j != -1ll) goto found; } break;
            #elif defined(__i386__) || defined(_WIN32)
                case 1: ie >>= 5; for(i>>=5;i<ie;i++) { j = bitbuf[i];                                                                              if (j != -1) goto found; } break;
                case 2: ie >>= 5; for(i>>=5;i<ie;i++) { j = bitbuf[i]; j |= (j>>1);                           k = bitbuf[i+1]; j |= (      k <<31); if (j != -1) goto found; } break;
                case 3: ie >>= 5; for(i>>=5;i<ie;i++) { j = bitbuf[i]; j |= (j>>1)|(j>>2);                    k = bitbuf[i+1]; j |= (((-k)|k)<<30); if (j != -1) goto found; } break;
                case 4: ie >>= 5; for(i>>=5;i<ie;i++) { j = bitbuf[i]; j |= (j>>1); j |= (j>>2);              k = bitbuf[i+1]; j |= (((-k)|k)<<29); if (j != -1) goto found; } break;
                case 5: ie >>= 5; for(i>>=5;i<ie;i++) { j = bitbuf[i]; j |= (j>>1); j |= (j>>2)|(j>>3);       k = bitbuf[i+1]; j |= (((-k)|k)<<28); if (j != -1) goto found; } break;
                case 6: ie >>= 5; for(i>>=5;i<ie;i++) { j = bitbuf[i]; j |= (j>>1); j |= (j>>2)|(j>>4);       k = bitbuf[i+1]; j |= (((-k)|k)<<27); if (j != -1) goto found; } break;
                case 7: ie >>= 5; for(i>>=5;i<ie;i++) { j = bitbuf[i]; j |= (j>>1); j |= (j>>2); j |= (j>>3); k = bitbuf[i+1]; j |= (((-k)|k)<<26); if (j != -1) goto found; } break;
                case 8: ie >>= 5; for(i>>=5;i<ie;i++) { j = bitbuf[i]; j |= (j>>1); j |= (j>>2); j |= (j>>4); k = bitbuf[i+1]; j |= (((-k)|k)<<25); if (j != -1) goto found; } break;
            #else
                #error
            #endif
            default:
                for(;i<ie;i=i1+1) //NOTE:this seems to be faster than the above cases
                {
                    i0 = dntil0((unsigned int *)bitbuf,i   ,ie  ); if (i0 >= ie) break;
                    i1 = dntil1((unsigned int *)bitbuf,i0+1,i0+n); if (i1-i0 < n) continue;
                    setzrange1(bitbuf,i0,i0+n); bm->ind = i0+n; return(i0);
                }
                break;
        }
        cnt--;
        if (cnt < 0) break;
        i = 0;
        ie = std::min(oi, (int)(bm->mal - n));
    }

    i = bm->mal;

    //Only surf needs to do GPU stuff :P
    if ((oct_usegpu) && (&loct->sur == bm)) {
        //grow space by 100%
        if constexpr (oct_usegpubo) {
            if (loct->gsurf)
            {
                loct->gsurf = 0;
                ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_2D,loct->octid);
                bo_end(loct->bufid,0,0,0,0,GL_RGBA,GL_UNSIGNED_BYTE,0);
            }
            ((PFNGLDELETEBUFFERS)glfp[glDeleteBuffers])(1,&loct->bufid);
        }

        if (loct->glxsid < loct->glysid) loct->glxsid++; else loct->glysid++;
        loct->gxsid = (1<<loct->glxsid);
        loct->gysid = (1<<loct->glysid);
        bm->mal = loct->gxsid*loct->gysid/(sizeof(surf_t)>>2);
        //only NEAREST makes sense here!
        kglalloctex(loct->octid,0,loct->gxsid,loct->gysid,1,KGL_RGBA32+KGL_NEAREST);

        if constexpr (oct_usegpubo) {
            loct->bufid = bo_init(loct->gxsid*loct->gysid*4);
            loct->gsurf = (surf_t *)bo_begin(loct->bufid,0);
            memcpy(loct->gsurf,bm->buf,i*bm->siz);
        }
        else {
            ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_2D,loct->octid);
            ((PFNGLTEXSUBIMAGE2D)glfp[glTexSubImage2D])(GL_TEXTURE_2D,0,0,0,loct->gxsid,(i*2)>>loct->glxsid,GL_RGBA,GL_UNSIGNED_BYTE,(void *)bm->buf);
        }
    }
    else {
        //grow space by ~25%
        bm->mal = std::max(((1<<bsr(bm->mal))>>2) + bm->mal,bm->mal+n);
    }

    bm->buf = (octv_t       *)realloc(bm->buf,   (int64_t)bm->mal*bm->siz);        if (!bm->buf) { fprintf(stderr, "realloc(bm->buf,%lld) failed\n",   (long long)bm->mal*bm->siz);        }
    bm->bit = (unsigned int *)realloc(bm->bit,((((int64_t)bm->mal+63)>>5)<<2)+16); if (!bm->bit) { fprintf(stderr, "realloc(bm->bit,%lld) failed\n",((((long long)bm->mal+63)>>5)<<2)+16); }

    setzrange1(bm->bit,i  ,i+n);
    setzrange0(bm->bit,i+n,bm->mal);
    bm->ind = i+n;
    return(i);

found:
    #if defined(__x86_64__) || defined(_WIN64)
        i = (i<<6)+bsf64(~j);
    #elif defined(__i386__) || defined(_WIN32)
        i = (i<<5)+bsf(~j);
    #else
        #error
    #endif

    xorzrangesmall(bitbuf,i,n); bm->ind = i; return(i);
}

//Input: o=n(old), n={0..8}, ind=octree index, noct[0..n-1]
static int oct_refreshnode (oct_t *loct, int ind, int o, int n, octv_t *noct)
{
    bitmal_t *bm;
    bm = &loct->nod;

    if (n < o)
    {
        xorzrangesmall(bm->bit,ind+n,o-n); //shorten node
        if (!n) ind = -1;
        bm->num += n-o;
    }
    else if (n > o)
    {
        if (o) xorzrangesmall(bm->bit,ind,o); //delete node (if(o) prevents possible read fault if !o)
        ind = bitalloc(loct,bm,n); //alloc node
        bm->num += n-o;
    }
    if (n)
    {
        memcpy(&((octv_t *)bm->buf)[ind],noct,n*bm->siz);
    }
    return(ind);
}

//Input: o=n(old), n={0..8}, ind=octree index, surf[0..n-1]
static int oct_refreshsurf (oct_t *loct, int ind, int o, int n, surf_t *surf)
{
    bitmal_t *bm;
    bm = &loct->sur;

    if (n < o)
    {
        xorzrangesmall(bm->bit,ind+n,o-n); //shorten node
        if (!n) ind = -1;
        bm->num += n-o;
    }
    else if (n > o)
    {
        if (o) xorzrangesmall(bm->bit,ind,o); //delete node (if(o) prevents possible read fault if !o)
        ind = bitalloc(loct,bm,n); //alloc node
        bm->num += n-o;
    }
    if (n)
    {
        //memcpy(&((octv_t *)bm->buf)[ind],surf,n*bm->siz);
        memcpy(&((surf_t *)bm->buf)[ind],surf,n*bm->siz); //FIXFIXFIXFIX::? ????
        if constexpr (oct_usegpu && oct_usegpubo) {
            memcpy(&loct->gsurf[ind],surf,n*bm->siz);
        }
    }
    return(ind);
}

static void oct_modnew_copy_recur (oct_t *loct, int inode, int ls, octv_t *roct, oct_t *newoct)
{
    octv_t noct[8];
    int i, n;
    (*roct) = ((octv_t *)loct->nod.buf)[inode]; n = popcount8(roct->chi);
    if (!ls) { roct->ind = oct_refreshsurf(newoct,0,0,n,&((surf_t *)loct->sur.buf)[roct->ind]); return; }
    for(i=0;i<n;i++) oct_modnew_copy_recur(loct,roct->ind+i,ls-1,&noct[i],newoct);
    roct->ind = oct_refreshnode(newoct,0,0,n,noct);
}
    //(mode&1)==0:air, !=0:sol
    //(mode&2)==0:nand, 1:and
static void oct_modnew_recur (oct_t *loct, int inode, int x0, int y0, int z0, int ls, octv_t *roct, brush_t *brush, oct_t *newoct, int mode)
{
    surf_t surf[8];
    octv_t ooct, noct[8];
    int i, iup, n, o, s, x, y, z, xsol, issol;

    issol = (mode&1)*255;
    if (inode >= 0) { ooct = ((octv_t *)loct->nod.buf)[inode]; xsol = ooct.sol^issol; }
                  else { xsol = (1<<8)-1; ooct.chi = 0; ooct.sol = issol^xsol; ooct.ind = -1; ooct.mrk = 0; ooct.mrk2 = 0; }

    if (!ls)
    {
        roct->sol = ooct.sol; roct->mrk = 0; roct->mrk2 = 0;
        for(;xsol;xsol^=iup)
        {
            iup = (-xsol)&xsol;
            i = isins_func(loct,brush,(((iup&0xaa)+0xff)>>8)+x0,(((iup&0xcc)+0xff)>>8)+y0,(((iup&0xf0)+0xff)>>8)+z0,0);
            if (!(mode&2)) i = 2-i;
            if (i) { roct->sol ^= iup; }
        }

        roct->chi = ooct.chi&roct->sol;
        for(n=0,xsol=roct->chi;xsol;xsol^=iup,n++)
        {
            iup = (-xsol)&xsol;
            memcpy(&surf[n],&((surf_t *)loct->sur.buf)[popcount8((iup-1)&ooct.chi)+ooct.ind],loct->sur.siz); //copy existing surf
        }
        roct->ind = oct_refreshsurf(newoct,0,0,n,surf);

        brush->mx0 = std::min(brush->mx0,x0); brush->mx1 = std::max(brush->mx1,x0+2);
        brush->my0 = std::min(brush->my0,y0); brush->my1 = std::max(brush->my1,y0+2);
        brush->mz0 = std::min(brush->mz0,z0); brush->mz1 = std::max(brush->mz1,z0+2);
        return;
    }

    xsol |= ooct.chi; roct->chi = 0; roct->sol = ~xsol&issol; roct->mrk = 0; roct->mrk2 = 0; //copy solid if already brush color
    o = ooct.ind; n = 0; s = pow2[ls];
    for(;xsol;xsol^=iup) //visit only nodes that may differ from brush color
    {
        iup = (-xsol)&xsol;
        x = x0; if (iup&0xaa) x += s;
        y = y0; if (iup&0xcc) y += s;
        z = z0; if (iup&0xf0) z += s;
        i = isins_func(loct,brush,x,y,z,ls);
        if (!(mode&2)) i = 2-i;
        switch(i)
        {
            case 0: //octree node doesn't intersect brush:copy old tree
                roct->sol += (ooct.sol&iup);
                if (ooct.chi&iup) { roct->chi += iup; oct_modnew_copy_recur(loct,o,ls-1,&noct[n],newoct); o++; n++; }
                break;
            case 1:
                if (ooct.chi&iup) { i = o; o++; } //octree node intersects brush partially:must recurse
                                 else { i =-1;      } //leaf; split pure node
                oct_modnew_recur(loct,i,x,y,z,ls-1,&noct[n],brush,newoct,mode);
                if (noct[n].sol == (1<<8)-1) roct->sol += iup;
                if (noct[n].chi || ((unsigned)(noct[n].sol-1) < (unsigned)((1<<8)-2))) { roct->chi += iup; n++; }
                break;
            case 2: //brush fully covers octree node:all brush
                roct->sol += (issol&iup);
                if (ooct.chi&iup) o++;
                brush->mx0 = std::min(brush->mx0,x); brush->mx1 = std::max(brush->mx1,x+s);
                brush->my0 = std::min(brush->my0,y); brush->my1 = std::max(brush->my1,y+s);
                brush->mz0 = std::min(brush->mz0,z); brush->mz1 = std::max(brush->mz1,z+s);
                break;
            default:
                // The default can't be reached
                std::abort();
        }
    }

    roct->ind = oct_refreshnode(newoct,0,0,n,noct);
}

//mode==0: newoct =  oct & brush          //most useful - extracts piece without requiring full copy
//mode==1: newoct = (oct & brush)|~brush  //not useful - very unusual bool op
//mode==2: newoct =  oct &~brush          //useful, but copy&oct_mod() also works
//mode==3: newoct =  oct | brush          //useful, but copy&oct_mod() also works
void oct_modnew (oct_t *loct, brush_t *brush, int mode)
{
    oct_t newoct;
    octv_t nnode;

    oct_new(loct->xres, loct->yres, &newoct,loct->lsid,loct->tilid,256,256,0);
    if constexpr (oct_usegpu && oct_usegpubo) {
        if (!newoct.gsurf) newoct.gsurf = (surf_t *)bo_begin(newoct.bufid,0);
    }

    brush->mx0 = 0x7fffffff; brush->mx1 = 0x80000000;
    brush->my0 = 0x7fffffff; brush->my1 = 0x80000000;
    brush->mz0 = 0x7fffffff; brush->mz1 = 0x80000000;
    oct_modnew_recur(loct,loct->head,0,0,0,loct->lsid-1,&nnode,brush,&newoct,mode);
    ((octv_t *)newoct.nod.buf)[newoct.head] = nnode; //can't write node directly because of possible realloc inside

    oct_updatesurfs(&newoct,0,0,0,newoct.sid,newoct.sid,newoct.sid,brush,mode);
    oct_checkreducesizes(&newoct);

    loct->recvoctfunc(loct,&newoct);
}

void oct_mod_recur(oct_t *loct, int inode, int x0, int y0, int z0, int ls, octv_t *roct, brush_t *brush, int issol)
{
    surf_t surf[8];
    octv_t ooct, noct[8];
    int i, iup, n, o, s, x, y, z, xsol;

    if (inode >= 0) {
        ooct = ((octv_t *)loct->nod.buf)[inode];
        xsol = ooct.sol^issol;
    }
    else {
        xsol = (1<<8)-1;
        ooct.chi = 0;
        ooct.sol = issol^xsol;
        ooct.ind = -1;
        ooct.mrk = 0;
        ooct.mrk2 = 0;
    }

    if (!ls)
    {
        roct->sol = ooct.sol;
        roct->mrk = 0;
        roct->mrk2 = 0;

        if constexpr (false) {
            for(;xsol;xsol^=iup) {
                iup = (-xsol)&xsol;
                if (isins_func(loct,brush,(((iup&0xaa)+0xff)>>8)+x0,(((iup&0xcc)+0xff)>>8)+y0,(((iup&0xf0)+0xff)>>8)+z0,0)) { roct->sol ^= iup; } //roct->mrk ^= iup; }
            }
        }
        else {
            if ((xsol&(1<<0)) && (isins_func(loct,brush,x0  ,y0  ,z0  ,0))) { roct->sol ^= (1<<0); } //roct->mrk ^= (1<<0); }
            if ((xsol&(1<<1)) && (isins_func(loct,brush,x0+1,y0  ,z0  ,0))) { roct->sol ^= (1<<1); } //roct->mrk ^= (1<<1); }
            if ((xsol&(1<<2)) && (isins_func(loct,brush,x0  ,y0+1,z0  ,0))) { roct->sol ^= (1<<2); } //roct->mrk ^= (1<<2); }
            if ((xsol&(1<<3)) && (isins_func(loct,brush,x0+1,y0+1,z0  ,0))) { roct->sol ^= (1<<3); } //roct->mrk ^= (1<<3); }
            if ((xsol&(1<<4)) && (isins_func(loct,brush,x0  ,y0  ,z0+1,0))) { roct->sol ^= (1<<4); } //roct->mrk ^= (1<<4); }
            if ((xsol&(1<<5)) && (isins_func(loct,brush,x0+1,y0  ,z0+1,0))) { roct->sol ^= (1<<5); } //roct->mrk ^= (1<<5); }
            if ((xsol&(1<<6)) && (isins_func(loct,brush,x0  ,y0+1,z0+1,0))) { roct->sol ^= (1<<6); } //roct->mrk ^= (1<<6); }
            if ((xsol&(1<<7)) && (isins_func(loct,brush,x0+1,y0+1,z0+1,0))) { roct->sol ^= (1<<7); } //roct->mrk ^= (1<<7); }
        }

        if (ooct.chi&~roct->sol)
        {
            roct->chi = ooct.chi&roct->sol;
            for(n=0,xsol=roct->chi;xsol;xsol^=iup,n++)
            {
                iup = (-xsol)&xsol;
                memcpy(&surf[n],&((surf_t *)loct->sur.buf)[popcount8((iup-1)&ooct.chi)+ooct.ind],loct->sur.siz); //copy existing surf
            }
            roct->ind = oct_refreshsurf(loct,ooct.ind,popcount8(ooct.chi),n,surf);
        } else { roct->chi = ooct.chi; roct->ind = ooct.ind; }

        brush->mx0 = std::min(brush->mx0,x0); brush->mx1 = std::max(brush->mx1,x0+2);
        brush->my0 = std::min(brush->my0,y0); brush->my1 = std::max(brush->my1,y0+2);
        brush->mz0 = std::min(brush->mz0,z0); brush->mz1 = std::max(brush->mz1,z0+2);
        return;
    }

    //copy solid if already brush color
    xsol |= ooct.chi;
    roct->chi = 0;
    roct->sol = ~xsol&issol;
    roct->mrk = 0;
    roct->mrk2 = 0;

    o = ooct.ind;
    n = 0;
    s = pow2[ls];

    for(;xsol;xsol^=iup) //visit only nodes that may differ from brush color
    {
        iup = (-xsol)&xsol;
        x = x0; if (iup&0xaa) x += s;
        y = y0; if (iup&0xcc) y += s;
        z = z0; if (iup&0xf0) z += s;
        switch(isins_func(loct,brush,x,y,z,ls))
        {
            case 0: //octree node doesn't intersect brush:copy old tree
                roct->sol += (ooct.sol&iup);
                if (ooct.chi&iup) { roct->chi += iup; noct[n] = ((octv_t *)loct->nod.buf)[o]; o++; n++; }
                break;
            case 1:
                if (ooct.chi&iup) { i = o; o++; } //octree node intersects brush partially:must recurse
                                 else { i =-1;      } //leaf; split pure node
                oct_mod_recur(loct,i,x,y,z,ls-1,&noct[n],brush,issol);
                //if (noct[n].mrk) roct->mrk |= iup;
                if (noct[n].sol == (1<<8)-1) roct->sol += iup;
                if (noct[n].chi || ((unsigned)(noct[n].sol-1) < (unsigned)((1<<8)-2))) { roct->chi += iup; n++; }
                break;
            case 2: //brush fully covers octree node:all brush
                roct->sol += (issol&iup);
                if (ooct.chi&iup) { oct_dealloctree(loct,ls-1,((octv_t *)loct->nod.buf)[o].ind,((octv_t *)loct->nod.buf)[o].chi); o++; }
                brush->mx0 = std::min(brush->mx0,x); brush->mx1 = std::max(brush->mx1,x+s);
                brush->my0 = std::min(brush->my0,y); brush->my1 = std::max(brush->my1,y+s);
                brush->mz0 = std::min(brush->mz0,z); brush->mz1 = std::max(brush->mz1,z+s);
                //roct->mrk |= iup;
                break;
            default:
                // The default can't be reached
                std::abort();
        }
    }
    roct->ind = oct_refreshnode(loct,ooct.ind,o-ooct.ind,n,noct);
}

typedef struct { oct_t *loct; int mode, mx0, my0, mz0, mx1, my1, mz1, bitmask; oct_getvox_hint_t och; } oct_updatesurfs_t;
static void oct_updatesurfs_recur (oct_updatesurfs_t *ous, int inode, int x0, int y0, int z0, int ls, octv_t *roct, brush_t *brush)
{
    surf_t surf[8];
    octv_t ooct, noct[8];
    int i, iup, n, o, s, x, y, z, xsol;

    if (inode >= 0) { ooct = ((octv_t *)ous->loct->nod.buf)[inode]; }
                  else { ooct.chi = 0; ooct.sol = (1<<8)-1; ooct.ind = -1; ooct.mrk = 0; ooct.mrk2 = 0; }
    roct->sol = ooct.sol;
    roct->mrk = 0; roct->mrk2 = 0;

    if (!ls)
    {
        if (0)
        {
              //Slow&brute:
            roct->chi = 0;
            for(xsol=ooct.sol;xsol;xsol^=iup)
            {
                iup = (-xsol)&xsol;
                x = (((iup&0xaa)+0xff)>>8) + x0;
                y = (((iup&0xcc)+0xff)>>8) + y0;
                z = (((iup&0xf0)+0xff)>>8) + z0;
                if (oct_issurf(ous->loct,x,y,z,0,&ous->och)) roct->chi += iup;
            }
        }
        else
        {
                //above:      top:       bot:     below:
                //           18 19      22 23
                //+-----+    +---+      +---+    +-----+
                //|28 29|   9|0 1|8   13|4 5|12  |24 25|
                //|30 31|  11|2 3|10  15|6 7|14  |26 27|
                //+-----+    +---+      +---+    +-----+
                //           16 17      20 21
            i = ooct.sol; s = ous->loct->sid-2;
            if (x0 == 0) i += 0x0000aa00; else if (x0 == s) i += 0x00005500;
            if (y0 == 0) i += 0x00cc0000; else if (y0 == s) i += 0x00330000;
            if (z0 == 0) i += 0xf0000000; else if (z0 == s) i += 0x0f000000;
            i &= ous->bitmask;

            if (((i|~((1<<0)+(1<<1)+(1<<2)+(1<<4)                )) == -1) ||
                 ((i|~((1<<2)+(1<<0)+(1<<3)+(1<<6)                )) == -1) ||
                 ((i|~((1<<4)+(1<<0)+(1<<5)+(1<<6)                )) == -1) ||
                 ((i|~((1<<6)+(1<<2)+(1<<4)+(1<<7)                )) == -1)) i += ((oct_getsol_hint_2x2x2(ous->loct,x0-2,y0  ,z0  ,&ous->och)&0xaa)<< 8);
            if (((i|~((1<<1)+(1<<0)+(1<<3)+(1<<5)                )) == -1) ||
                 ((i|~((1<<3)+(1<<1)+(1<<2)+(1<<7)                )) == -1) ||
                 ((i|~((1<<5)+(1<<1)+(1<<4)+(1<<7)                )) == -1) ||
                 ((i|~((1<<7)+(1<<3)+(1<<5)+(1<<6)                )) == -1)) i += ((oct_getsol_hint_2x2x2(ous->loct,x0+2,y0  ,z0  ,&ous->och)&0x55)<< 8);
            if (((i|~((1<<0)+(1<<1)+(1<<2)+(1<<4)+(1<< 9)        )) == -1) ||
                 ((i|~((1<<1)+(1<<0)+(1<<3)+(1<<5)+(1<< 8)        )) == -1) ||
                 ((i|~((1<<4)+(1<<0)+(1<<5)+(1<<6)+(1<<13)        )) == -1) ||
                 ((i|~((1<<5)+(1<<1)+(1<<4)+(1<<7)+(1<<12)        )) == -1)) i += ((oct_getsol_hint_2x2x2(ous->loct,x0  ,y0-2,z0  ,&ous->och)&0xcc)<<16);
            if (((i|~((1<<2)+(1<<0)+(1<<3)+(1<<6)+(1<<11)        )) == -1) ||
                 ((i|~((1<<3)+(1<<1)+(1<<2)+(1<<7)+(1<<10)        )) == -1) ||
                 ((i|~((1<<6)+(1<<2)+(1<<4)+(1<<7)+(1<<15)        )) == -1) ||
                 ((i|~((1<<7)+(1<<3)+(1<<5)+(1<<6)+(1<<14)        )) == -1)) i += ((oct_getsol_hint_2x2x2(ous->loct,x0  ,y0+2,z0  ,&ous->och)&0x33)<<16);
            if (((i|~((1<<0)+(1<<1)+(1<<2)+(1<<4)+(1<< 9)+(1<<18))) == -1) ||
                 ((i|~((1<<1)+(1<<0)+(1<<3)+(1<<5)+(1<< 8)+(1<<19))) == -1) ||
                 ((i|~((1<<2)+(1<<0)+(1<<3)+(1<<6)+(1<<11)+(1<<16))) == -1) ||
                 ((i|~((1<<3)+(1<<1)+(1<<2)+(1<<7)+(1<<10)+(1<<17))) == -1)) i += ((oct_getsol_hint_2x2x2(ous->loct,x0  ,y0  ,z0-2,&ous->och)&0xf0)<<24);
            if (((i|~((1<<4)+(1<<0)+(1<<5)+(1<<6)+(1<<13)+(1<<22))) == -1) ||
                 ((i|~((1<<5)+(1<<1)+(1<<4)+(1<<7)+(1<<12)+(1<<23))) == -1) ||
                 ((i|~((1<<6)+(1<<2)+(1<<4)+(1<<7)+(1<<15)+(1<<20))) == -1) ||
                 ((i|~((1<<7)+(1<<3)+(1<<5)+(1<<6)+(1<<14)+(1<<21))) == -1)) i += ((oct_getsol_hint_2x2x2(ous->loct,x0  ,y0  ,z0+2,&ous->och)&0x0f)<<24);

            if constexpr (false) {
                roct->chi = i; i = ~i;
                if (!(i&((1<<1)+(1<< 9) + (1<<2)+(1<<18) + (1<<4)+(1<<28)))) roct->chi &=~(1<<0);
                if (!(i&((1<<0)+(1<< 8) + (1<<3)+(1<<19) + (1<<5)+(1<<29)))) roct->chi &=~(1<<1);
                if (!(i&((1<<3)+(1<<11) + (1<<0)+(1<<16) + (1<<6)+(1<<30)))) roct->chi &=~(1<<2);
                if (!(i&((1<<2)+(1<<10) + (1<<1)+(1<<17) + (1<<7)+(1<<31)))) roct->chi &=~(1<<3);
                if (!(i&((1<<5)+(1<<13) + (1<<6)+(1<<22) + (1<<0)+(1<<24)))) roct->chi &=~(1<<4);
                if (!(i&((1<<4)+(1<<12) + (1<<7)+(1<<23) + (1<<1)+(1<<25)))) roct->chi &=~(1<<5);
                if (!(i&((1<<7)+(1<<15) + (1<<4)+(1<<20) + (1<<2)+(1<<26)))) roct->chi &=~(1<<6);
                if (!(i&((1<<6)+(1<<14) + (1<<5)+(1<<21) + (1<<3)+(1<<27)))) roct->chi &=~(1<<7);
            }
            else {
                o = (i>> 8)&i; s  = ((o&0x55)<<1) + ((o&0xaa)>>1);
                o = (i>>16)&i; s &= ((o&0x33)<<2) + ((o&0xcc)>>2);
                o = (i>>24)&i; s &= ((o&0x0f)<<4) + ((o&0xf0)>>4);
                roct->chi = (~s)&i;
            }
        }

        for(n=0,xsol=roct->chi;xsol;xsol^=iup,n++)
        {
            iup = (-xsol)&xsol;
            if (ooct.chi&iup) { memcpy(&surf[n],&((surf_t *)ous->loct->sur.buf)[popcount8((iup-1)&ooct.chi)+ooct.ind],ous->loct->sur.siz); continue; } //copy existing surf
            x = (((iup&0xaa)+0xff)>>8) + x0;
            y = (((iup&0xcc)+0xff)>>8) + y0;
            z = (((iup&0xf0)+0xff)>>8) + z0;
            getsurf_func(ous->loct,brush,x,y,z,&surf[n]);
        }
        roct->ind = oct_refreshsurf(ous->loct,ooct.ind,popcount8(ooct.chi),n,surf);
        //NOTE:don't need to update ind&chi of parent here because it is surface (not geometry) info
        return;
    }

    roct->chi = ooct.chi; xsol = ooct.sol|ooct.chi; //skip pure air
    o = ooct.ind; n = 0; s = pow2[ls];
    for(;xsol;xsol^=iup) //visit nodes that are not pure air
    {
        iup = (-xsol)&xsol;
        x = x0; if (iup&0xaa) x += s;
        y = y0; if (iup&0xcc) y += s;
        z = z0; if (iup&0xf0) z += s;

        //if (!(ooct.mrk&iup)) //FIX:can't use here
        if ((x+s < ous->mx0) || (x > ous->mx1) || (y+s < ous->my0) || (y > ous->my1) || (z+s < ous->mz0) || (z > ous->mz1)) //outside update region: copy tree
        {
            if (ooct.chi&iup) { noct[n] = ((octv_t *)ous->loct->nod.buf)[o]; o++; n++; }
            continue;
        }

        if (!oct_issurf(ous->loct,x,y,z,ls,&ous->och)) //no surfs inside: remove tree
        {
            if (ooct.chi&iup)
            {
                oct_dealloctree(ous->loct,ls-1,((octv_t *)ous->loct->nod.buf)[o].ind,((octv_t *)ous->loct->nod.buf)[o].chi); o++;
                roct->chi ^= iup;
            }
            roct->sol |= iup;
            continue;
        }

        if (ooct.chi&iup) { i = o; o++; } //octree node intersects brush partially:must recurse
                         else { i =-1;      } //crack solid.. insert 1 node at all lower levels

        oct_updatesurfs_recur(ous,i,x,y,z,ls-1,&noct[n],brush); //recurse (crack solid if o<0)
        roct->chi |= iup;
        n++;
    }
    roct->ind = oct_refreshnode(ous->loct,ooct.ind,o-ooct.ind,n,noct);
    if (inode >= 0) ((octv_t *)ous->loct->nod.buf)[inode] = *roct; //NOTE! must keep tree&hint cache valid because of oct_getsol..() during mod
}

void oct_updatesurfs (oct_t *loct, int mx0, int my0, int mz0, int mx1, int my1, int mz1, brush_t *brush, int mode)
{
    oct_updatesurfs_t ous;
    octv_t nnode;

    if constexpr (oct_usegpu && oct_usegpubo) {
        if (!loct->gsurf) {
            loct->gsurf = (surf_t *)bo_begin(loct->bufid,0);
        }
    }

    ous.loct = loct; ous.mode = mode;
    ous.mx0 = std::max(mx0,0); ous.mx1 = std::min(mx1,loct->sid);
    ous.my0 = std::max(my0,0); ous.my1 = std::min(my1,loct->sid);
    ous.mz0 = std::max(mz0,0); ous.mz1 = std::min(mz1,loct->sid);

    ous.bitmask = 0x000000ff;
    if (loct->edgeissol& 1) ous.bitmask += 0x0000aa00;
    if (loct->edgeissol& 2) ous.bitmask += 0x00cc0000;
    if (loct->edgeissol& 4) ous.bitmask += 0xf0000000;
    if (loct->edgeissol& 8) ous.bitmask += 0x00005500;
    if (loct->edgeissol&16) ous.bitmask += 0x00330000;
    if (loct->edgeissol&32) ous.bitmask += 0x0f000000;

    oct_getvox_hint_init(loct,&ous.och);
    oct_updatesurfs_recur(&ous,loct->head,0,0,0,loct->lsid-1,&nnode,brush);
    ((octv_t *)loct->nod.buf)[loct->head] = nnode;

    oct_checkreducesizes(loct);
}

void oct_getvox_hint_init (oct_t *loct, oct_getvox_hint_t *och) {
    och->stkind[loct->lsid-1] = loct->head;
    och->minls = loct->lsid-1;
    och->mins = (loct->sid>>1);
    och->ox = 0;
    och->oy = 0;
    och->oz = 0;
}

// WARNING:assumes x,y,z inside grid
int oct_getsol_hint (oct_t *loct, int x, int y, int z, oct_getvox_hint_t *och)
{
    octv_t *ptr;
    int i, ls, s;

    ls = bsr((och->ox^x)|(och->oy^y)|(och->oz^z)|och->mins); s = pow2[ls];
    for(;1;s>>=1,ls--)
    {
        ptr = &((octv_t *)loct->nod.buf)[och->stkind[ls]];
        i = 1;
        if (x&s) i <<= 1;
        if (y&s) i <<= 2;
        if (z&s) i <<= 4;
        if (ptr->sol&i) break;
        if ((!(ptr->chi&i)) || (!ls)) { i = 0; break; }
        och->stkind[ls-1] = popcount8(ptr->chi&(i-1)) + ptr->ind; //2child
    }
    och->minls = ls; och->mins = s; och->ox = x; och->oy = y; och->oz = z; return(i);
}

// returns sol mask of 2x2x2 cell (optimization for smallest node)
// WARNING:assumes x,y,z even #'s
__forceinline
int oct_getsol_hint_2x2x2(oct_t *loct, int x, int y, int z, oct_getvox_hint_t *och) //WARNING:assumes x,y,z even #'s
{
    octv_t *ptr;
    int i, ls, s;

    ls = (och->ox^x)|(och->oy^y)|(och->oz^z); if (ls&loct->nsid) return(0);
    ls = bsr(ls|och->mins); s = pow2[ls];
    for(;1;s>>=1,ls--)
    {
        ptr = &((octv_t *)loct->nod.buf)[och->stkind[ls]];
        if (!ls) { i = ptr->sol; break; }
        i = 1;
        if (x&s) i <<= 1;
        if (y&s) i <<= 2;
        if (z&s) i <<= 4;
        if (ptr->sol&i) { i = (1<<8)-1; break; }
        if (!(ptr->chi&i)) { i = 0; break; }
        och->stkind[ls-1] = popcount8(ptr->chi&(i-1)) + ptr->ind; //2child
    }
    och->minls = ls; och->mins = s; och->ox = x; och->oy = y; och->oz = z; return(i);
}

//returns: -1:invalid, {0..sur.mal-1}:surface index
int oct_getsurf_hint (oct_t *loct, int x, int y, int z, oct_getvox_hint_t *och)
{
    octv_t *ptr;
    int i, ls, s;

    if ((x|y|z)&loct->nsid) return(-1);
    ls = bsr((och->ox^x)|(och->oy^y)|(och->oz^z)|och->mins); s = pow2[ls];
    for(;1;s>>=1,ls--)
    {
        ptr = &((octv_t *)loct->nod.buf)[och->stkind[ls]];
        i = 1;
        if (x&s) i <<= 1;
        if (y&s) i <<= 2;
        if (z&s) i <<= 4;
        if (!(ptr->chi&i)) { i = -1; break; }
        i = popcount8(ptr->chi&(i-1)) + ptr->ind;
        if (!ls) break;
        och->stkind[ls-1] = i; //2child
    }
    och->minls = ls; och->mins = s; och->ox = x; och->oy = y; och->oz = z; return(i);
}

typedef struct { oct_t *loct; int remap[8]; char invrebits[256]; } osr_t;
static void oct_swizzle_recur (osr_t *osr, int ind, int ls)
{
    surf_t surf[8];
    octv_t node[8], *ptr;
    int i, j, n, ochi;

    ptr = &((octv_t *)osr->loct->nod.buf)[ind];
    ochi = ptr->chi;
    ptr->chi = osr->invrebits[ptr->chi];
    ptr->sol = osr->invrebits[ptr->sol];
    n = 0;
    for(i=0;i<8;i++)
    {
        j = osr->remap[i];
        if (!(ochi&(1<<j))) continue;
        if (ls) memcpy(&node[n],&((surf_t *)osr->loct->nod.buf)[popcount8(((1<<j)-1)&ochi)+ptr->ind],osr->loct->nod.siz);
            else memcpy(&surf[n],&((surf_t *)osr->loct->sur.buf)[popcount8(((1<<j)-1)&ochi)+ptr->ind],osr->loct->sur.siz);
        n++;
    }

    if (ls) memcpy(&((surf_t *)osr->loct->nod.buf)[ptr->ind],node,osr->loct->nod.siz*n);
        else memcpy(&((surf_t *)osr->loct->sur.buf)[ptr->ind],surf,osr->loct->sur.siz*n);

    if (ls) for(i=0;i<n;i++) oct_swizzle_recur(osr,ptr->ind+i,ls-1);
}
//+x=+1  +y=+2  +z=+3
//-x=-1  -y=-2  -z=-3
void oct_swizzle (oct_t *loct, int ax0, int ax1, int ax2)
{
    osr_t osr;
    int i, j, k, m, invremap[8];

    if ((ax0 == 0) || (ax1 == 0) || (ax2 == 0)) return;
    if ((labs(ax0) > 3) || (labs(ax1) > 3) || (labs(ax2) > 3)) return;
    if ((labs(ax0) == labs(ax1)) || (labs(ax0) == labs(ax2)) || (labs(ax1) == labs(ax2))) return;

    if constexpr (oct_usegpu && oct_usegpubo) {
        if (!loct->gsurf) loct->gsurf = (surf_t *)bo_begin(loct->bufid,0);
    }

    osr.loct = loct;
    switch(labs(ax0))
    {
        case 1: if (labs(ax1) == 2) j = 0x01234567; else j = 0x01452367; break;
        case 2: if (labs(ax1) == 1) j = 0x02134657; else j = 0x02461357; break;
        case 3: if (labs(ax1) == 1) j = 0x04152637; else j = 0x04261537; break;
    }
    k = ((ax0<0)<<(labs(ax0)-1)) + ((ax1<0)<<(labs(ax1)-1)) + ((ax2<0)<<(labs(ax2)-1));
    for(i=8-1;i>=0;i--) osr.remap[i] = ((j>>(28-(i<<2)))&15)^k;

    for(i=8-1;i>=0;i--) invremap[osr.remap[i]] = i;
    memset(osr.invrebits,0,sizeof(osr.invrebits));
    for(j=8-1;j>=0;j--)
    {
        k = (1<<j); m = (1<<invremap[j]);
        for(i=k;i<256;i=((i+1)|k)) osr.invrebits[i] += m;
    }

    oct_swizzle_recur(&osr,loct->head,loct->lsid-1);

    if constexpr (oct_usegpu) {
        if constexpr (oct_usegpubo) {
            memcpy(loct->gsurf,loct->sur.buf,loct->sur.mal*loct->sur.siz);
        }
        else {
            ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_2D,loct->octid);
            ((PFNGLTEXSUBIMAGE2D)glfp[glTexSubImage2D])(GL_TEXTURE_2D,0,0,0,loct->gxsid,(loct->sur.mal*2)>>loct->glxsid,GL_RGBA,GL_UNSIGNED_BYTE,(void *)loct->sur.buf);
        }
    }
}

//                                      nsol:
//chi sol   ls=0   ls>0            | ls=0  ls>0
// 0   0    air    pure air        |   1     1
// 1   0    N/A    hybrid          |  (0)    0
// 0   1    sol    pure sol        |   0     0
// 1   1    surf   pure sol w/surf |   1   &OfChildren
static void oct_swap_sol_air_recur (oct_t *loct, int ind, int ls)
{
    octv_t *ptr;
    int o, v, iup;

    ptr = &((octv_t *)loct->nod.buf)[ind];
    if (!ls)
    {
        ptr->sol = (~ptr->sol)|ptr->chi;
        return;
    }
    ptr->sol = ~(ptr->chi|ptr->sol); o = ptr->ind;
    for(v=ptr->chi;v;v^=iup)
    {
        iup = (-v)&v;
        oct_swap_sol_air_recur(loct,o,ls-1);
        if (((octv_t *)loct->nod.buf)[o].sol == (1<<8)-1) ptr->sol |= iup;
        o++;
    }
}
void oct_swap_sol_air (oct_t *loct)
{
    oct_swap_sol_air_recur(loct,loct->head,loct->lsid-1);
    loct->edgeissol ^= 63;
}

//Ensures box is defined in octree (i.e. box is inside 0..sid-1)
//Algo is: grow .. translate .. shrink
int oct_rebox (oct_t *loct, int bx0, int by0, int bz0, int bx1, int by1, int bz1, int *dx, int *dy, int *dz)
{
    static const char wshlut[3][8] = {-1,1,7,9,15,17,23,25, -1,6,2,10,14,22,18,26, -1,4,12,20,4,12,20,28};
    static const char wshr[3][8] = {0,-1,8,-1,16,-1,24,-1, 0,8,-1,-1,16,24,-1,-1, 0,8,16,24,-1,-1,-1,-1};
    static const int windex[3][4] = {0,2,4,6, 0,1,4,5, 0,1,2,3};
    octv_t nod0, nod1, onod[32];
    int nchi0, nsol0, nind0, nchi1, nsol1, nind1;
    int w, w0[3], w1[3];
    int i, j, k, m, n, ind, got, anygot;

    if (dx) (*dx) = 0;
    if (dy) (*dy) = 0;
    if (dz) (*dz) = 0;
    anygot = 0;

        //force bbox to be outside defined solid
    for(i=0;i<3;i++) { w0[i] = 0; w1[i] = loct->sid; }
    oct_getsolbbox(loct,&w0[0],&w0[1],&w0[2],&w1[0],&w1[1],&w1[2]);
    w0[0] = std::min(w0[0],bx0); w1[0] = std::max(w1[0],bx1);
    w0[1] = std::min(w0[1],by0); w1[1] = std::max(w1[1],by1);
    w0[2] = std::min(w0[2],bz0); w1[2] = std::max(w1[2],bz1);

        //grow (double)
    while ((loct->lsid < OCT_MAXLS) && (((w1[0]-1)|(w1[1]-1)|(w1[2]-1)|w0[0]|w0[1]|w0[2])&loct->nsid))
    {
        ind = bitalloc(loct,&loct->nod,1); loct->nod.num++;
        ((octv_t *)loct->nod.buf)[ind] = ((octv_t *)loct->nod.buf)[loct->head];

        i = 1;
        if (w0[0]+w1[0] < loct->sid) i <<= 1;
        if (w0[1]+w1[1] < loct->sid) i <<= 2;
        if (w0[2]+w1[2] < loct->sid) i <<= 4;
        ((octv_t *)loct->nod.buf)[loct->head].chi = i;

        if (i&0xaa) { w0[0] += loct->sid; w1[0] += loct->sid; if (dx) (*dx) += loct->sid; }
        if (i&0xcc) { w0[1] += loct->sid; w1[1] += loct->sid; if (dy) (*dy) += loct->sid; }
        if (i&0xf0) { w0[2] += loct->sid; w1[2] += loct->sid; if (dz) (*dz) += loct->sid; }
        anygot = 1;

        ((octv_t *)loct->nod.buf)[loct->head].sol = 0;
        ((octv_t *)loct->nod.buf)[loct->head].ind = ind;
        loct->lsid++; loct->sid <<= 1; loct->nsid = -loct->sid;
    }

        //translate (needed to allow model with tiny cluster in center to shrink)

        //algo guarantees bbox region no smaller than sid/4+1 on the largest sid
    do
    {
        got = 0;

            //translate by sid/2 if possible
        i = (loct->sid>>1); ind = loct->head;
        for(w=0;w<3;w++)
        {
            if (w0[w] < i) continue;
            w0[w] -= i; w1[w] -= i; got = 1;
            if ((w == 0) && dx) (*dx) -= i;
            if ((w == 1) && dy) (*dy) -= i;
            if ((w == 2) && dz) (*dz) -= i;

            ((octv_t *)loct->nod.buf)[ind].chi >>= (1<<w); ((octv_t *)loct->nod.buf)[ind].sol >>= (1<<w);
        }

            //translate by sid/4 if possible
            //
            //   //Before:                                                       //After:
            //buf[ 0].chi = 0xff; buf[ 0].sol = 0x00; buf[ 0].ind = 10;       buf[  0].chi = 0x55; buf[  0].sol = 0x55; buf[  0].ind = 100;
            //..                                                              ..
            //buf[10].chi = 0xaa; buf[10].sol = 0xaa; buf[10].ind = 20;       buf[100].chi = 0xff; buf[100].sol = 0xff; buf[100].ind = 110;
            //buf[11].chi = 0x55; buf[11].sol = 0x55; buf[11].ind = 30;       buf[101].chi = 0xff; buf[101].sol = 0xff; buf[101].ind = 120;
            //buf[12].chi = 0xaa; buf[12].sol = 0xaa; buf[12].ind = 40;       buf[102].chi = 0xff; buf[102].sol = 0xff; buf[102].ind = 130;
            //buf[13].chi = 0x55; buf[13].sol = 0x55; buf[13].ind = 50;       buf[103].chi = 0xff; buf[103].sol = 0xff; buf[103].ind = 140;
            //buf[14].chi = 0xaa; buf[14].sol = 0xaa; buf[14].ind = 60;       ..
            //buf[15].chi = 0x55; buf[15].sol = 0x55; buf[15].ind = 70;       buf[110].chi = 0x7f; buf[110].sol = 0xff; buf[110].ind =   A;
            //buf[16].chi = 0xaa; buf[16].sol = 0xaa; buf[16].ind = 80;       buf[111].chi = 0xbf; buf[111].sol = 0xff; buf[111].ind =   E;
            //buf[17].chi = 0x55; buf[17].sol = 0x55; buf[17].ind = 90;       buf[112].chi = 0x5f; buf[112].sol = 0xff; buf[112].ind =   B;
            //..                                                              buf[113].chi = 0xaf; buf[113].sol = 0xff; buf[113].ind =   F;
            //buf[20].chi = 0x7f; buf[20].sol = 0xff; buf[20].ind =  A;       buf[114].chi = 0x77; buf[114].sol = 0xff; buf[114].ind =   C;
            //buf[21].chi = 0x5f; buf[21].sol = 0xff; buf[21].ind =  B;       buf[115].chi = 0xbb; buf[115].sol = 0xff; buf[115].ind =   G;
            //buf[22].chi = 0x77; buf[22].sol = 0xff; buf[22].ind =  C;       buf[116].chi = 0x55; buf[116].sol = 0xff; buf[116].ind =   D;
            //buf[23].chi = 0x55; buf[23].sol = 0xff; buf[23].ind =  D;       buf[117].chi = 0xaa; buf[117].sol = 0xff; buf[117].ind =   H;
            //..                                                              ..
            //buf[30].chi = 0xbf; buf[30].sol = 0xff; buf[30].ind =  E;
            //buf[31].chi = 0xaf; buf[31].sol = 0xff; buf[31].ind =  F;
            //buf[32].chi = 0xbb; buf[32].sol = 0xff; buf[32].ind =  G;
            //buf[33].chi = 0xaa; buf[33].sol = 0xff; buf[33].ind =  H;
            //..
            //
            //X: A--+--B--+--+  +--+--+--+--+  E--+--F--+--+  +--+--+--+--+       j=0, 1, 2, 3, 4, 5, 6, 7
            //   |  | 0| 1|  |  |  | 4| 5|  |  |  |16|17|  |  |  |20|21|  |   i=0      0     2     4     6
            //   +--+--+--+--+  +--+--+--+--+  +--+--+--+--+  +--+--+--+--+   i=1   1     3     5     7
            //   |  | 2| 3|  |  |  | 6| 7|  |  |  |18|19|  |  |  |22|23|  |   i=2      8    10    12    14
            //   C--+--D--+--+  +--+--+--+--+  G--+--H--+--+  +--+--+--+--+   i=3   9    11    13    15
            //   |  | 8| 9|  |  |  |12|13|  |  |  |24|25|  |  |  |28|29|  |   i=4     16    18    20    22
            //   +--+--+--+--+  +--+--+--+--+  +--+--+--+--+  +--+--+--+--+   i=5  17    19    21    23
            //   |  |10|11|  |  |  |14|15|  |  |  |26|27|  |  |  |30|31|  |   i=6     24    26    28    30
            //   +--+--+--+--+  +--+--+--+--+  +--+--+--+--+  +--+--+--+--+   i=7  25    27    29    31
            //
            //Y: A--+--B--+--+  +--+--+--+--+  E--+--F--+--+  +--+--+--+--+       j=0, 1, 2, 3, 4, 5, 6, 7
            //   |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |   i=0         0  1        4  5
            //   +--+--+--+--+  +--+--+--+--+  +--+--+--+--+  +--+--+--+--+   i=1         8  9       12 13
            //   | 0| 1| 8| 9|  | 4| 5|12|13|  |16|17|24|25|  |20|21|28|29|   i=2   2  3        6  7
            //   C--+--D--+--+  +--+--+--+--+  G--+--H--+--+  +--+--+--+--+   i=3  10 11       14 15
            //   | 2| 3|10|11|  | 6| 7|14|15|  |18|19|26|27|  |22|23|30|31|   i=4        16 17       20 21
            //   +--+--+--+--+  +--+--+--+--+  +--+--+--+--+  +--+--+--+--+   i=5        24 25       28 29
            //   |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |   i=6  18 19       22 23
            //   +--+--+--+--+  +--+--+--+--+  +--+--+--+--+  +--+--+--+--+   i=7  26 27       30 31
            //
            //Z: A--+--B--+--+  +--+--+--+--+  E--+--F--+--+  +--+--+--+--+       j=0, 1, 2, 3, 4, 5, 6, 7
            //   |  |  |  |  |  | 0| 1| 8| 9|  | 4| 5|12|13|  |  |  |  |  |   i=0               0  1  2  3
            //   +--+--+--+--+  +--+--+--+--+  +--+--+--+--+  +--+--+--+--+   i=1               8  9 10 11
            //   |  |  |  |  |  | 2| 3|10|11|  | 6| 7|14|15|  |  |  |  |  |   i=2              16 17 18 19
            //   C--+--D--+--+  +--+--+--+--+  G--+--H--+--+  +--+--+--+--+   i=3              24 25 25 26
            //   |  |  |  |  |  |16|17|24|25|  |20|21|28|29|  |  |  |  |  |   i=4   4  5  6  7
            //   +--+--+--+--+  +--+--+--+--+  +--+--+--+--+  +--+--+--+--+   i=5  12 13 14 15
            //   |  |  |  |  |  |18|19|26|27|  |22|23|30|31|  |  |  |  |  |   i=6  20 21 22 23
            //   +--+--+--+--+  +--+--+--+--+  +--+--+--+--+  +--+--+--+--+   i=7  28 29 30 31
        for(w=0;w<3;w++)
        {
            if ((w0[w] < (loct->sid>>2)) || (w1[w] >= ((loct->sid*3)>>2))) continue;
            w0[w] -= (loct->sid>>2); w1[w] -= (loct->sid>>2); got = 1;
            if ((w == 0) && dx) (*dx) -= (loct->sid>>2);
            if ((w == 1) && dy) (*dy) -= (loct->sid>>2);
            if ((w == 2) && dz) (*dz) -= (loct->sid>>2);

            nchi1 = 0; nsol1 = 0;
            nod0 = ((octv_t *)loct->nod.buf)[loct->head];
            n = popcount8(nod0.chi); xorzrangesmall(loct->nod.bit,nod0.ind,n); loct->nod.num -= n; //remove 2nd level nodes (8 std::max)
            for(i=0;i<8;i++)
            {
                if (!(nod0.chi&(1<<i))) continue;
                nod1 = ((octv_t *)loct->nod.buf)[popcount8(nod0.chi&((1<<i)-1)) + nod0.ind];
                n = popcount8(nod1.chi); xorzrangesmall(loct->nod.bit,nod1.ind,n); loct->nod.num -= n; //remove 3rd level nodes (64 std::max)
                if (!i) //can't really encode shift left&right in same lut without platform incompatibilities :/
                {
                    nchi1 += ((int)nod1.chi>>(1<<w));
                    nsol1 += ((int)nod1.sol>>(1<<w));
                }
                else
                {
                    nchi1 += ((int)nod1.chi<<wshlut[w][i]);
                    nsol1 += ((int)nod1.sol<<wshlut[w][i]);
                }
                for(j=0;j<8;j++)
                {
                    if (!(nod1.chi&(1<<j))) continue;
                          if (w == 0) onod[((i&6)<<2)              + (j^1)] = ((octv_t *)loct->nod.buf)[popcount8(nod1.chi&((1<<j)-1)) + nod1.ind];
                    else if (w == 1) onod[((i&4)<<2) + ((i&1)<<3) + (j^2)] = ((octv_t *)loct->nod.buf)[popcount8(nod1.chi&((1<<j)-1)) + nod1.ind];
                    else             onod[((i&3)<<3)              + (j^4)] = ((octv_t *)loct->nod.buf)[popcount8(nod1.chi&((1<<j)-1)) + nod1.ind];
                }
            }
            nchi0 = 0; nsol0 = 0;
            for(i=0;i<4;i++)
            {
                nchi0 += (((nchi1& (255<<(i<<3)))!= 0)<<windex[w][i]);
                nsol0 += (((nsol1|~(255<<(i<<3)))==-1)<<windex[w][i]);
            }
            ((octv_t *)loct->nod.buf)[loct->head].chi = nchi0;
            ((octv_t *)loct->nod.buf)[loct->head].sol = nsol0;

            n = popcount8(nchi0);
            nind0 = bitalloc(loct,&loct->nod,n); loct->nod.num += n;
            ((octv_t *)loct->nod.buf)[loct->head].ind = nind0;
            for(j=0,i=0;i<8;i++)
            {
                if (!(nchi0&(1<<i))) continue;
                ((octv_t *)loct->nod.buf)[nind0+j].chi = (nchi1>>wshr[w][i]);
                ((octv_t *)loct->nod.buf)[nind0+j].sol = (nsol1>>wshr[w][i]);
                n = popcount8(((octv_t *)loct->nod.buf)[nind0+j].chi);
                nind1 = bitalloc(loct,&loct->nod,n); loct->nod.num += n;
                ((octv_t *)loct->nod.buf)[nind0+j].ind = nind1;
                for(m=0,k=0;k<8;k++)
                {
                    if (!(((octv_t *)loct->nod.buf)[nind0+j].chi&(1<<k))) continue;
                          if (w == 0) { ((octv_t *)loct->nod.buf)[nind1+m] = onod[((i&6)<<2)              + k]; }
                    else if (w == 1) { ((octv_t *)loct->nod.buf)[nind1+m] = onod[((i&4)<<2) + ((i&1)<<3) + k]; }
                    else             { ((octv_t *)loct->nod.buf)[nind1+m] = onod[((i&3)<<3)              + k]; }
                    m++;
                }
                j++;
            }
        }

            //shrink (halve) if possible
        while ((loct->lsid > 2) && (!(((w1[0]-1)|(w1[1]-1)|(w1[2]-1))&(loct->nsid>>1))))
        {
            if (popcount8(((octv_t *)loct->nod.buf)[loct->head].chi) != 1) break;
            got = 1;
            i = ((octv_t *)loct->nod.buf)[loct->head].ind;
            xorzrangesmall(loct->nod.bit,i,1); loct->nod.num--; //remove node
            ((octv_t *)loct->nod.buf)[loct->head] = ((octv_t *)loct->nod.buf)[i];
            loct->lsid--; loct->sid >>= 1; loct->nsid = -loct->sid;
        }
        if (got) anygot = 1;
    } while (got);
    return(anygot);
}
