// This file has been modified from Ken Silverman's original release

#include "brushes.hpp"

#include "utilities/maths.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <emmintrin.h>

#define cvttss2si(f) _mm_cvtt_ss2si(_mm_set_ss(f))

__forceinline
int block_radius(int x, int y, int z) {
    int i = x*x + y*y + z*z;

    #if 1
        i = (int)(-64.0*65536.0 / std::sqrt((double)i));
    #else
        constexpr static const float kmul = -64.0*65536.0;
        _asm
        {
            cvtsi2ss xmm0, i
            rsqrtss xmm0, xmm0
            mulss xmm0, kmul
            cvtss2si eax, xmm0
            mov i, eax
        }
    #endif

    return i;
}


//tree node vs. sphere
//Returns: 0: node doesn't intersect brush
//         1: node partially  inside brush
//         2: node fully      inside brush
static int brush_sph_isins (brush_t *brush, int x0, int y0, int z0, int ls)
{
    brush_sph_t *sph;
    int x, y, z, sm1, dx, dy, dz, nx, ny, nz;

    sph = (brush_sph_t *)brush;
    sm1 = (1<<ls)-1;
    dx = sph->x-x0; nx = sm1-dx;
    dy = sph->y-y0; ny = sm1-dy;
    dz = sph->z-z0; nz = sm1-dz;
    x = std::max(    dx   , nx); y = std::max(    dy   , ny); z = std::max(    dz   , nz); if (x*x + y*y + z*z <  sph->r2) return(2); //farthest point inside?
    x = std::max(std::min(dx,0),-nx); y = std::max(std::min(dy,0),-ny); z = std::max(std::min(dz,0),-nz); if (x*x + y*y + z*z >= sph->r2) return(0); //nearest point outside?
    return (1);
}

static void brush_sph_getsurf (brush_t *brush, int x0, int y0, int z0, surf_t *surf)
{
    brush_sph_t *sph;
    int i, dx, dy, dz;

    //fixcnt_getsurf++;
    sph = (brush_sph_t *)brush;
    i = (x0+13)*(y0+17)*(z0+19)*0x25d3eb;
    surf->b = std::min(std::max(( sph->col     &255) + /*((i>> 8)&15)*/ - 8,0),255);
    surf->g = std::min(std::max(((sph->col>> 8)&255) + /*((i>>12)&15)*/ - 8,0),255);
    surf->r = std::min(std::max(((sph->col>>16)&255) + /*((i>>16)&15)*/ - 8,0),255);
    surf->a = 255;
    surf->tex = ((x0^y0^z0)&63);
    dx = x0-sph->x;
    dy = y0-sph->y;
    dz = z0-sph->z;

    i = block_radius(dx, dy, dz);

    surf->norm[0] = (signed char)((dx*i)>>16);
    surf->norm[1] = (signed char)((dy*i)>>16);
    surf->norm[2] = (signed char)((dz*i)>>16);
}

void brush_sph_init (brush_sph_t *sph, int x, int y, int z, int r, int issol)
{
    static_cast<void>(issol);
    sph->isins = brush_sph_isins;
    sph->getsurf = brush_sph_getsurf;
    sph->flags = 1;
    sph->x = x; sph->y = y; sph->z = z; sph->r2 = r*r;
    sph->col = (((x>>2)&255)<<16) + (((y>>2)&255)<< 8) + (((z>>2)&255)<< 0) + 0x404040;
}

//--------------------------------------------------------------------------------------------------
//tree node vs. box

//Returns: 0: node doesn't intersect brush
//         1: node partially  inside brush
//         2: node fully      inside brush
static int brush_box_isins (brush_t *brush, int x0, int y0, int z0, int ls)
{
    brush_box_t *box;
    int s;

    box = (brush_box_t *)brush;

    s = (1<<ls);

    if ((x0 >= box->x0) && (x0+s <= box->x1) &&
         (y0 >= box->y0) && (y0+s <= box->y1) &&
         (z0 >= box->z0) && (z0+s <= box->z1)) return(2);
    if ((x0+s <= box->x0) || (x0 >= box->x1) ||
         (y0+s <= box->y0) || (y0 >= box->y1) ||
         (z0+s <= box->z0) || (z0 >= box->z1)) return(0);
    return(1);
}

static void brush_box_getsurf (brush_t *brush, int x0, int y0, int z0, surf_t *surf)
{
    brush_box_t *box;
    int i, dx, dy, dz;

    box = (brush_box_t *)brush;
    i = (x0+13)*(y0+17)*(z0+19)*0x25d3eb;
    surf->b = std::min(std::max(( box->col     &255)/* + ((i>> 8)&15)*/ - 8,0),255);
    surf->g = std::min(std::max(((box->col>> 8)&255)/* + ((i>>12)&15)*/ - 8,0),255);
    surf->r = std::min(std::max(((box->col>>16)&255)/* + ((i>>16)&15)*/ - 8,0),255);
    surf->a = 255;
    surf->tex = ((x0^y0^z0)&63);

    dx = 0; dy = 0; dz = 0;
    if (x0 == box->x0) dx--; else if (x0 == box->x1-1) dx++;
    if (y0 == box->y0) dy--; else if (y0 == box->y1-1) dy++;
    if (z0 == box->z0) dz--; else if (z0 == box->z1-1) dz++;

    i = block_radius(dx, dy, dz);

    surf->norm[0] = (signed char)((dx*i)>>16);
    surf->norm[1] = (signed char)((dy*i)>>16);
    surf->norm[2] = (signed char)((dz*i)>>16);
}

void brush_box_init (brush_box_t *box, float x0, float y0, float z0, float x1, float y1, float z1, int issol)
{
    static_cast<void>(issol);
    box->isins = brush_box_isins;
    box->getsurf = brush_box_getsurf;
    box->flags = 1;
    box->x0 = cvttss2si(x0);
    box->y0 = cvttss2si(y0);
    box->z0 = cvttss2si(z0);
    box->x1 = cvttss2si(x1);
    box->y1 = cvttss2si(y1);
    box->z1 = cvttss2si(z1);
    box->col = ((((int)cvttss2si(x0+x1)>>3)&255)<<16) + ((((int)cvttss2si(y0+y1)>>3)&255)<< 8) + ((((int)cvttss2si(z0+z1)>>3)&255)<< 0) + 0x404040;
}

void brush_box_init (brush_box_t *box, int x0, int y0, int z0, int x1, int y1, int z1, int issol)
{
    static_cast<void>(issol);
    box->isins = brush_box_isins;
    box->getsurf = brush_box_getsurf;
    box->flags = 1;
    box->x0 = x0;
    box->y0 = y0;
    box->z0 = z0;
    box->x1 = x1;
    box->y1 = y1;
    box->z1 = z1;
    box->col = ((((int)cvttss2si(x0+x1)>>3)&255)<<16) + ((((int)cvttss2si(y0+y1)>>3)&255)<< 8) + ((((int)cvttss2si(z0+z1)>>3)&255)<< 0) + 0x404040;
}

//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
    //tree node vs. sphere, solid surf - for oct_paint()
static void brush_sph_sol_getsurf (brush_t *brush, int x0, int y0, int z0, surf_t *surf)
{
    static_cast<void>(x0);
    static_cast<void>(y0);
    static_cast<void>(z0);
    brush_sph_t *sph;
    sph = (brush_sph_t *)brush;
    std::memcpy(surf,sph->surf,sizeof(surf_t));
}
void brush_sph_sol_init (brush_sph_t *sph, int x, int y, int z, int r, surf_t *surf)
{
    sph->isins = brush_sph_isins;
    sph->getsurf = brush_sph_sol_getsurf;
    sph->flags = 1;
    sph->x = x; sph->y = y; sph->z = z; sph->r2 = r*r;
    sph->surf = surf;
}
//--------------------------------------------------------------------------------------------------
//FIXFIXFIX

//tree node vs. cone

//See CONE_INTERSECT.KC for derivation
static void indrawcone3d_init (brush_cone_t *ic)
{
    float x0, y0, z0, r0, x1, y1, z1, r1;
    float t, dx, dy, dz, dr, d2, f, k, r;

    x0 = ic->x0; y0 = ic->y0; z0 = ic->z0; r0 = ic->r0;
    x1 = ic->x1; y1 = ic->y1; z1 = ic->z1; r1 = ic->r1;
    if (fabs(r0) > fabs(r1)) //avoid bad case of 1 fully inside 0
    {
        t = x0; x0 = x1; x1 = t;
        t = y0; y0 = y1; y1 = t;
        t = z0; z0 = z1; z1 = t;
        t = r0; r0 = r1; r1 = t;
    }
    ic->x0 = x0; ic->y0 = y0; ic->z0 = z0; ic->r02 = std::max(r0,0.f); ic->r02 *= ic->r02; r0 = fabs(r0);
    ic->x1 = x1; ic->y1 = y1; ic->z1 = z1; ic->r12 = std::max(r1,0.f); ic->r12 *= ic->r12; r1 = fabs(r1);
    dx = x1-x0; dy = y1-y0; dz = z1-z0; dr = r1-r0; d2 = dx*dx + dy*dy + dz*dz - dr*dr; f = 1/sqrt(d2);
    if (fabs(dr) < r0*1e-6) { ic->hak = r0*r0; ic->dr = 0.f;  ic->cr = r0;  k =      f; r =   0.f; r0 = 0.f; r1 = d2; }
                             else { ic->hak = 0.f;   ic->dr = dr*f; ic->cr = 0.f; k = 1.f/dr; r = -r0*k; k *= f*d2; }
    ic->cx = dx*r + x0; ic->dx = dx*f; ic->k0 = r0*k;
    ic->cy = dy*r + y0; ic->dy = dy*f; ic->k1 = r1*k;
    ic->cz = dz*r + z0; ic->dz = dz*f;
}

static int indrawcone3d (brush_cone_t *ic, float x, float y, float z)
{
    float nx, ny, nz, d;
    nx = x-ic->cx; ny = y-ic->cy; nz = z-ic->cz; d = nx*ic->dx + ny*ic->dy + nz*ic->dz;
    if (d <= ic->k0) return((x-ic->x0)*(x-ic->x0) + (y-ic->y0)*(y-ic->y0) + (z-ic->z0)*(z-ic->z0) < ic->r02);
    if (d >= ic->k1) return((x-ic->x1)*(x-ic->x1) + (y-ic->y1)*(y-ic->y1) + (z-ic->z1)*(z-ic->z1) < ic->r12);
    return(nx*nx + ny*ny + nz*nz < d*d + ic->hak);
}

//find point on line segment that minimizes cone angle
//Derivation:
//   x = (x1-x0)*t + x0
//   y = (y1-y0)*t + y0
//   z = (z1-z0)*t + z0
//   ((x-cx)*dx + (y-cy)*dy + (z-cz)*dz) = sqrt((x-cx)^2 + (y-cy)^2 + (z-cz)^2)*sqrt(dx^2 + dy^2 + dz^2)*cos(ang)
//To solve, plug in x/y/z, solve derivate w.r.t. t (quotient rule), then solve numerator=0 (quadratic equation)
static float nearestptline2cone (float dx, float dy, float dz, float dr, float cx, float cy, float cz, float cr, float x0, float y0, float z0, float x1, float y1, float z1)
{
    static_cast<void>(dr);
    static_cast<void>(cr);
    float k0, k1, k2, k3, k4, k5, k6, k7, Za, Zb, Zc, insqr, t, s;
    k0 = (x1-x0)*dx + (y1-y0)*dy + (z1-z0)*dz;
    k1 = (x0-cx)*dx + (y0-cy)*dy + (z0-cz)*dz;
    k2 = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0);
    k3 =((x1-x0)*(x0-cx) + (y1-y0)*(y0-cy) + (z1-z0)*(z0-cz))*2;
    k4 = (x0-cx)*(x0-cx) + (y0-cy)*(y0-cy) + (z0-cz)*(z0-cz);
    k5 = k0*k0;
    k6 = k0*k1*2;
    k7 = k1*k1;
    Za = k3*k5 - k2*k6;
    Zb =(k4*k5 - k2*k7)*2;
    Zc = k4*k6 - k3*k7;
    insqr = Zb*Zb - 4.f*Za*Zc; if (insqr < 0.f) return(-1);
    for(s=-1.f;s<=1.f;s+=2.f)
    {
        t = (-Zb + s*sqrt(insqr))/(Za+Za);
        if ((t > 0.f) && (t < 1.f)) return(t);
    }
    return(-1.f);
}

static int box_sph_isins (float bx, float by, float bz, float bs, float cx, float cy, float cz, float cr)
{
    float dx, dy, dz, nx, ny, nz, x, y, z;
    dx = cx-bx; nx = bs-dx;
    dy = cy-by; ny = bs-dy;
    dz = cz-bz; nz = bs-dz;
    x = std::max(    dx     , nx); y = std::max(    dy     , ny); z = std::max(    dz,      nz); if (x*x + y*y + z*z <  cr*cr) return(2); //far  pt in?
    x = std::max(std::min(dx,0.f),-nx); y = std::max(std::min(dy,0.f),-ny); z = std::max(std::min(dz,0.f),-nz); if (x*x + y*y + z*z >= cr*cr) return(0); //near pt out?
    return(1);
}

static int indrawcone3d_intbox (brush_cone_t *ic, float bx, float by, float bz, float bs)
{
    static const signed char edge[12][6] =
    {
        0,0,0,1,0,0, 0,1,0,1,1,0, 0,0,1,1,0,1, 0,1,1,1,1,1,
        0,0,0,0,1,0, 1,0,0,1,1,0, 0,0,1,0,1,1, 1,0,1,1,1,1,
        0,0,0,0,0,1, 1,0,0,1,0,1, 0,1,0,0,1,1, 1,1,0,1,1,1,
    };
    float t, x, y, z, r, nx, ny, nz, nx0, ny0, nz0, nx1, ny1, nz1;
    int i;

    for(i=8-1;i>=0;i--)
    {
        nx = bx+((i   )&1)*bs;
        ny = by+((i>>1)&1)*bs;
        nz = bz+((i>>2)&1)*bs;
        if (!indrawcone3d(ic,nx,ny,nz)) break;
    }
    if (i < 0) return(2);
    for(i=0;i<12;i++)
    {
        nx0 = bx + (float)edge[i][0]*bs; ny0 = by + (float)edge[i][1]*bs; nz0 = bz + (float)edge[i][2]*bs;
        nx1 = bx + (float)edge[i][3]*bs; ny1 = by + (float)edge[i][4]*bs; nz1 = bz + (float)edge[i][5]*bs;
        t = nearestptline2cone(ic->dx,ic->dy,ic->dz,ic->dr, ic->cx,ic->cy,ic->cz,ic->cr, nx0,ny0,nz0, nx1,ny1,nz1);
        if (t < 0.f) continue;

        nx = (nx1-nx0)*t + nx0;
        ny = (ny1-ny0)*t + ny0;
        nz = (nz1-nz0)*t + nz0;
        t = std::min(std::max((nx-ic->cx)*ic->dx + (ny-ic->cy)*ic->dy + (nz-ic->cz)*ic->dz,ic->k0),ic->k1);
        x = t*ic->dx + ic->cx;
        y = t*ic->dy + ic->cy;
        z = t*ic->dz + ic->cz;
        r = t*ic->dr + ic->cr;
        if (box_sph_isins(bx,by,bz,bs,x,y,z,r)) return(1);
    }
    for(i=0;i<8;i++)
    {
        nx = bx+((i   )&1)*bs;
        ny = by+((i>>1)&1)*bs;
        nz = bz+((i>>2)&1)*bs;
        t = std::min(std::max((nx-ic->cx)*ic->dx + (ny-ic->cy)*ic->dy + (nz-ic->cz)*ic->dz,ic->k0),ic->k1);
        x = t*ic->dx + ic->cx;
        y = t*ic->dy + ic->cy;
        z = t*ic->dz + ic->cz;
        r = t*ic->dr + ic->cr;
        if (box_sph_isins(bx,by,bz,bs,x,y,z,r)) return(1);
    }
    return(0);
}

//Returns: 0: node doesn't intersect brush
//         1: node partially  inside brush
//         2: node fully      inside brush
static int brush_cone_isins (brush_t *brush, int x0, int y0, int z0, int ls)
{
    brush_cone_t *c;
    float d, dx, dy, dz, dx2, dy2, dz2;
    int i, x, y, z, x1, y1, z1;

    c = (brush_cone_t *)brush;
    //x1 = x0 + (1<<ls);
    //y1 = y0 + (1<<ls);
    //z1 = z0 + (1<<ls);

    //cube: (x0,y0,z0) - (x1,y1,z1)
    //   vs.
    //cone: (c->x0,c->y0,c->z0,c->r0) - (c->x1,c->y1,c->z1,c->r1)

    return(indrawcone3d_intbox(c,(float)x0,(float)y0,(float)z0,(float)(1<<ls)));

    if constexpr (false) {
        //if all 8 points of cube are inside cone, return(2);
        for(i=8-1;i>=0;i--)
        {
            if (i&1) x = x0; else x = x1;
            if (i&2) y = y0; else y = y1;
            if (i&4) z = z0; else z = z1;

            dx = ((float)x)-c->x1;
            dy = ((float)y)-c->y1;
            dz = ((float)z)-c->z1;
            if (dx*dx + dy*dy + dz*dz >= c->r1*c->r1) break;

            dx = ((float)x)-c->x0;
            dy = ((float)y)-c->y0;
            dz = ((float)z)-c->z0;
            if (dx*dx + dy*dy + dz*dz >= c->r0*c->r0) break;

            dx2 = c->x1-c->x0;
            dy2 = c->y1-c->y0;
            dz2 = c->z1-c->z0;
            d = (dx*dx2 + dy*dy2 + dz*dz2) / (dx2*dx2 + dy2*dy2 + dz2*dz2);
            //????//FIXFIXFIX
        }
        if (i < 0) {
            return(2);
        }


        //if (x*x + y*y + z*z >= c->r1) {
        //    return(0); //nearest point outside?
        //}
    }

    return 0;
}

static void brush_cone_getsurf (brush_t *brush, int x0, int y0, int z0, surf_t *surf)
{
    brush_cone_t *cone;
    int i, dx, dy, dz;

    //fixcnt_getsurf++;
    cone = (brush_cone_t *)brush;
    *(int *)&surf->b = cone->col;
    surf->tex = ((x0^y0^z0)&63);
    dx = x0-cvttss2si((cone->x0+cone->x1)*.5f);
    dy = y0-cvttss2si((cone->y0+cone->y1)*.5f);
    dz = z0-cvttss2si((cone->z0+cone->z1)*.5f);

    i = block_radius(dx, dy, dz);

    surf->norm[0] = (signed char)((dx*i)>>16);
    surf->norm[1] = (signed char)((dy*i)>>16);
    surf->norm[2] = (signed char)((dz*i)>>16);
}

void brush_cone_init (brush_cone_t *cone, float x0, float y0, float z0, float r0, float x1, float y1, float z1, float r1)
{
    cone->isins = brush_cone_isins;
    cone->getsurf = brush_cone_getsurf;
    cone->flags = 1;
    cone->x0 = x0; cone->y0 = y0; cone->z0 = z0; cone->r0 = r0;
    cone->x1 = x1; cone->y1 = y1; cone->z1 = z1; cone->r1 = r1;
    indrawcone3d_init(cone);
    cone->col = ((((int)x0>>2)&255)<<16) + ((((int)y0>>2)&255)<< 8) + ((((int)z0>>2)&255)<< 0) + 0x404040;
}

//tree node vs. single voxel

//Returns: 0: node doesn't intersect brush
//         1: node partially  inside brush
//         2: node fully      inside brush
static int brush_vox_isins (brush_t *brush, int x0, int y0, int z0, int ls)
{
    brush_vox_t *vox;
    int s;

    vox = (brush_vox_t *)brush;
    s = (1<<ls);
    if ((vox->x < x0) || (vox->x >= x0+s)) {
        return(0);
    }
    if ((vox->y < y0) || (vox->y >= y0+s)) {
        return(0);
    }
    if ((vox->z < z0) || (vox->z >= z0+s)) {
        return(0);
    }
    return((!ls)+1);
}

static void brush_vox_getsurf (brush_t *brush, int x0, int y0, int z0, surf_t *surf)
{
    static_cast<void>(x0);
    static_cast<void>(y0);
    static_cast<void>(z0);
    brush_vox_t *vox = (brush_vox_t *)brush;
    //fixcnt_getsurf++;
    (*surf) = vox->surf;
}

void brush_vox_init (brush_vox_t *vox, int x, int y, int z, surf_t *surf)
{
    vox->isins = brush_vox_isins;
    vox->getsurf = brush_vox_getsurf;
    vox->flags = 1;
    vox->x = x;
    vox->y = y;
    vox->z = z;
    vox->surf = (*surf);
}

//tree node vs. bitmap
//Returns: 0: node doesn't intersect brush
//         1: node partially  inside brush
//         2: node fully      inside brush
static int brush_bmp_isins (brush_t *brush, int x0, int y0, int z0, int ls)
{
    int i, ox0, oy0, oz0, ox1, oy1, oz1, nx0, ny0, nz0, nx1, ny1, nz1;
    brush_bmp_t *bmp = (brush_bmp_t *)brush;

    i = bmp->iux*x0 + bmp->ivx*y0 + bmp->iwx*z0; ox0 = ((bmp->iox0[ls]+i)>>16); ox1 = ((bmp->iox1[ls]+i)>>16); nx0 = std::max(ox0,0); nx1 = std::min(ox1,bmp->xs); if (nx0 >= nx1) return(0);
    i = bmp->iuy*x0 + bmp->ivy*y0 + bmp->iwy*z0; oy0 = ((bmp->ioy0[ls]+i)>>16); oy1 = ((bmp->ioy1[ls]+i)>>16); ny0 = std::max(oy0,0); ny1 = std::min(oy1,bmp->ys); if (ny0 >= ny1) return(0);
    i = bmp->iuz*x0 + bmp->ivz*y0 + bmp->iwz*z0; oz0 = ((bmp->ioz0[ls]+i)>>16); oz1 = ((bmp->ioz1[ls]+i)>>16); nz0 = std::max(oz0,0); nz1 = std::min(oz1,bmp->zs); if (nz0 >= nz1) return(0);

    i =  (bmp->xs+1); ny0 *= i; ny1 *= i;
    i *= (bmp->ys+1); nz0 *= i; nz1 *= i;
    i  = (bmp->boxsum[nz1+ny0+nx0] - bmp->boxsum[nz1+ny1+nx0] - bmp->boxsum[nz1+ny0+nx1] + bmp->boxsum[nz1+ny1+nx1]);
    i -= (bmp->boxsum[nz0+ny0+nx0] - bmp->boxsum[nz0+ny1+nx0] - bmp->boxsum[nz0+ny0+nx1] + bmp->boxsum[nz0+ny1+nx1]);
    if (i ==                             0) return(0);
    if (i == (ox1-ox0)*(oy1-oy0)*(oz1-oz0)) return(2);
    return(1);
}

static void brush_bmp_getsurf (brush_t *brush, int x0, int y0, int z0, surf_t *surf)
{
    static_cast<void>(brush);
    static_cast<void>(x0);
    static_cast<void>(y0);
    static_cast<void>(z0);

    //brush_bmp_t *bmp = (brush_bmp_t *)brush;
    //int x = ((bmp->iux*x0 + bmp->ivx*y0 + bmp->iwx*z0 + bmp->iox0[0])>>16);
    //int y = ((bmp->iuy*x0 + bmp->ivy*y0 + bmp->iwy*z0 + bmp->ioy0[0])>>16);
    //int z = ((bmp->iuz*x0 + bmp->ivz*y0 + bmp->iwz*z0 + bmp->ioz0[0])>>16);

    surf->b = 128;
    surf->g = 128;
    surf->r = 192;
    surf->a = 0;
    surf->norm[0] = 0;
    surf->norm[1] = 0;
    surf->norm[2] = 0;
    surf->tex = 0;
}

//pp,pr,pd,pf is bmp's pos&ori in oct's coords
void brush_bmp_init (brush_bmp_t *bmp, unsigned short *boxsum, int xs, int ys, int zs, point3f *pp, point3f *pr, point3f *pd, point3f *pf)
{
    float mat[9];
    int i, ls, iux, iuy, iuz, ivx, ivy, ivz, iwx, iwy, iwz, iox, ioy, ioz;

    bmp->isins   = brush_bmp_isins;
    bmp->getsurf = brush_bmp_getsurf;
    bmp->boxsum = boxsum; bmp->xs = xs; bmp->ys = ys; bmp->zs = zs;

    //convert bmp to loct coords using:
    //x2 = x*pr->x + y*pd->x + z*pf->x + pp->x;
    //y2 = x*pr->y + y*pd->y + z*pf->y + pp->y;
    //z2 = x*pr->z + y*pd->z + z*pf->z + pp->z;
    //
    //convert loct to bmp coords using:
    //x = (x2-pp->x)*mat[0] + (y2-pp->y)*mat[1] + (z2-pp->z)*mat[2]
    //y = (x2-pp->x)*mat[3] + (y2-pp->y)*mat[4] + (z2-pp->z)*mat[5]
    //z = (x2-pp->x)*mat[6] + (y2-pp->y)*mat[7] + (z2-pp->z)*mat[8]

    invert3x3(pr,pd,pf,mat);
    bmp->iux = cvttss2si(mat[0]*65536.0); bmp->ivx = cvttss2si(mat[1]*65536.0); bmp->iwx = cvttss2si(mat[2]*65536.0);
    bmp->iuy = cvttss2si(mat[3]*65536.0); bmp->ivy = cvttss2si(mat[4]*65536.0); bmp->iwy = cvttss2si(mat[5]*65536.0);
    bmp->iuz = cvttss2si(mat[6]*65536.0); bmp->ivz = cvttss2si(mat[7]*65536.0); bmp->iwz = cvttss2si(mat[8]*65536.0);
    iox = cvttss2si((pp->x*mat[0] + pp->y*mat[1] + pp->z*mat[2])*-65536);
    ioy = cvttss2si((pp->x*mat[3] + pp->y*mat[4] + pp->z*mat[5])*-65536);
    ioz = cvttss2si((pp->x*mat[6] + pp->y*mat[7] + pp->z*mat[8])*-65536);
    for(ls=0;ls<OCT_MAXLS;ls++)
    {
        i = (1<<ls)-1;
        iux = bmp->iux*i; ivx = bmp->ivx*i; iwx = bmp->iwx*i;
        iuy = bmp->iuy*i; ivy = bmp->ivy*i; iwy = bmp->iwy*i;
        iuz = bmp->iuz*i; ivz = bmp->ivz*i; iwz = bmp->iwz*i;
        bmp->iox0[ls] = std::min(iux,0) + std::min(ivx,0) + std::min(iwx,0) + iox;
        bmp->iox1[ls] = std::max(iux,0) + std::max(ivx,0) + std::max(iwx,0) + iox + 65536;
        bmp->ioy0[ls] = std::min(iuy,0) + std::min(ivy,0) + std::min(iwy,0) + ioy;
        bmp->ioy1[ls] = std::max(iuy,0) + std::max(ivy,0) + std::max(iwy,0) + ioy + 65536;
        bmp->ioz0[ls] = std::min(iuz,0) + std::min(ivz,0) + std::min(iwz,0) + ioz;
        bmp->ioz1[ls] = std::max(iuz,0) + std::max(ivz,0) + std::max(iwz,0) + ioz + 65536;
    }
}
