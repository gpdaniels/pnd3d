// This file has been modified from Ken Silverman's original release

#pragma once
#ifndef TYPES_PARALLELEPIPED_HPP
#define TYPES_PARALLELEPIPED_HPP

#include "utilities/macros.hpp"
#include "utilities/maths.hpp"
#include "types/point.hpp"

#include <algorithm>

//--------------------------------------------------------------------------------------------------
//Test intersection between 2 parallelepipeds.
//      bi: internal structure to be filled
//r0,d0,f0: axes of base parallelepiped (vectors specify side lengths)
//r1,d1,f1: axes of test parallelepiped (vectors specify side lengths)
struct pgram3d_t {
    double m0[9], m1[9], n[9], m[9], k[9];
    point3d r0, d0, f0, r1, d1, f1;
};

constexpr static const bool MAT0ISIDENTITY = false;

static void pgram3d_init (pgram3d_t *bi, point3d *r0, point3d *d0, point3d *f0, point3d *r1, point3d *d1, point3d *f1)
{
    int i, j0, j1;

    bi->r0 = (*r0); bi->d0 = (*d0); bi->f0 = (*f0);
    bi->r1 = (*r1); bi->d1 = (*d1); bi->f1 = (*f1);

    invert3x3(r1,d1,f1,bi->m1);
    if constexpr (MAT0ISIDENTITY) {
        bi->m[0] = fabs(bi->m1[0]) + fabs(bi->m1[1]) + fabs(bi->m1[2]) + 1.0;
        bi->m[1] = fabs(bi->m1[3]) + fabs(bi->m1[4]) + fabs(bi->m1[5]) + 1.0;
        bi->m[2] = fabs(bi->m1[6]) + fabs(bi->m1[7]) + fabs(bi->m1[8]) + 1.0;

        bi->n[0] = r1->x; bi->n[1] = r1->y; bi->n[2] = r1->z;
        bi->n[3] = d1->x; bi->n[4] = d1->y; bi->n[5] = d1->z;
        bi->n[6] = f1->x; bi->n[7] = f1->y; bi->n[8] = f1->z;
    }
    else {
        bi->m[0] = fabs(r0->x*bi->m1[0] + r0->y*bi->m1[1] + r0->z*bi->m1[2]) +
                      fabs(d0->x*bi->m1[0] + d0->y*bi->m1[1] + d0->z*bi->m1[2]) +
                      fabs(f0->x*bi->m1[0] + f0->y*bi->m1[1] + f0->z*bi->m1[2]) + 1.0;
        bi->m[1] = fabs(r0->x*bi->m1[3] + r0->y*bi->m1[4] + r0->z*bi->m1[5]) +
                      fabs(d0->x*bi->m1[3] + d0->y*bi->m1[4] + d0->z*bi->m1[5]) +
                      fabs(f0->x*bi->m1[3] + f0->y*bi->m1[4] + f0->z*bi->m1[5]) + 1.0;
        bi->m[2] = fabs(r0->x*bi->m1[6] + r0->y*bi->m1[7] + r0->z*bi->m1[8]) +
                      fabs(d0->x*bi->m1[6] + d0->y*bi->m1[7] + d0->z*bi->m1[8]) +
                      fabs(f0->x*bi->m1[6] + f0->y*bi->m1[7] + f0->z*bi->m1[8]) + 1.0;

        invert3x3(r0,d0,f0,bi->m0);
        bi->n[0] = r1->x*bi->m0[0] + r1->y*bi->m0[1] + r1->z*bi->m0[2];
        bi->n[1] = r1->x*bi->m0[3] + r1->y*bi->m0[4] + r1->z*bi->m0[5];
        bi->n[2] = r1->x*bi->m0[6] + r1->y*bi->m0[7] + r1->z*bi->m0[8];
        bi->n[3] = d1->x*bi->m0[0] + d1->y*bi->m0[1] + d1->z*bi->m0[2];
        bi->n[4] = d1->x*bi->m0[3] + d1->y*bi->m0[4] + d1->z*bi->m0[5];
        bi->n[5] = d1->x*bi->m0[6] + d1->y*bi->m0[7] + d1->z*bi->m0[8];
        bi->n[6] = f1->x*bi->m0[0] + f1->y*bi->m0[1] + f1->z*bi->m0[2];
        bi->n[7] = f1->x*bi->m0[3] + f1->y*bi->m0[4] + f1->z*bi->m0[5];
        bi->n[8] = f1->x*bi->m0[6] + f1->y*bi->m0[7] + f1->z*bi->m0[8];
    }

    for(i=0;i<3;i++) { bi->m[i+3] = fabs(bi->n[i]) + fabs(bi->n[i+3]) + fabs(bi->n[i+6]) + 1.0; bi->m[i+6] = 2.0 - bi->m[i+3]; }
    j0 = 2; j1 = 0;
    do
    {
        //NOTE:VC6 gave internal compiler error using: for(i0=..,i1=..,i2..)
        bi->k[j0*3+0] = fabs(bi->n[0+j0]*bi->n[3+j1] - bi->n[0+j1]*bi->n[3+j0])
                          + fabs(bi->n[0+j0]*bi->n[6+j1] - bi->n[0+j1]*bi->n[6+j0])
                          + fabs(bi->n[0+j0])       + fabs(bi->n[0+j1]);
        bi->k[j0*3+1] = fabs(bi->n[3+j0]*bi->n[6+j1] - bi->n[3+j1]*bi->n[6+j0])
                          + fabs(bi->n[3+j0]*bi->n[0+j1] - bi->n[3+j1]*bi->n[0+j0])
                          + fabs(bi->n[3+j0])       + fabs(bi->n[3+j1]);
        bi->k[j0*3+2] = fabs(bi->n[6+j0]*bi->n[0+j1] - bi->n[6+j1]*bi->n[0+j0])
                          + fabs(bi->n[6+j0]*bi->n[3+j1] - bi->n[6+j1]*bi->n[3+j0])
                          + fabs(bi->n[6+j0])       + fabs(bi->n[6+j1]);
        j0 = j1; j1++;
    } while (j1 < 3);
}

static int pgram3d_isint (pgram3d_t *bi, double dx, double dy, double dz)
{
    double x, y, z, ax, ay, az;

    if (fabs(dx*bi->m1[0] + dy*bi->m1[1] + dz*bi->m1[2]) >= bi->m[0]) return(0); //blue outside red face?
    if (fabs(dx*bi->m1[3] + dy*bi->m1[4] + dz*bi->m1[5]) >= bi->m[1]) return(0);
    if (fabs(dx*bi->m1[6] + dy*bi->m1[7] + dz*bi->m1[8]) >= bi->m[2]) return(0);

    if constexpr (MAT0ISIDENTITY) {
        x = dx; ax = fabs(dx); if (ax >= bi->m[3]) return(0); //red outside a blue face?
        y = dy; ay = fabs(dy); if (ay >= bi->m[4]) return(0);
        z = dz; az = fabs(dz); if (az >= bi->m[5]) return(0);
    }
    else {
        x = dx*bi->m0[0] + dy*bi->m0[1] + dz*bi->m0[2]; ax = fabs(x); if (ax >= bi->m[3]) return(0); //red outside blue face?
        y = dx*bi->m0[3] + dy*bi->m0[4] + dz*bi->m0[5]; ay = fabs(y); if (ay >= bi->m[4]) return(0);
        z = dx*bi->m0[6] + dy*bi->m0[7] + dz*bi->m0[8]; az = fabs(z); if (az >= bi->m[5]) return(0);
    }

    if ((ax <= bi->m[6]) && (ay <= bi->m[7]) && (az <= bi->m[8])) return(2); //blue inside red?
    if (fabs(x*bi->n[1] - y*bi->n[0]) > bi->k[0]) return(0); //blue outside red, 2D xy plane
    if (fabs(x*bi->n[4] - y*bi->n[3]) > bi->k[1]) return(0);
    if (fabs(x*bi->n[7] - y*bi->n[6]) > bi->k[2]) return(0);
    if (fabs(y*bi->n[2] - z*bi->n[1]) > bi->k[3]) return(0); //blue outside red, 2D yz plane
    if (fabs(y*bi->n[5] - z*bi->n[4]) > bi->k[4]) return(0);
    if (fabs(y*bi->n[8] - z*bi->n[7]) > bi->k[5]) return(0);
    if (fabs(z*bi->n[0] - x*bi->n[2]) > bi->k[6]) return(0); //blue outside red, 2D zx plane
    if (fabs(z*bi->n[3] - x*bi->n[5]) > bi->k[7]) return(0);
    if (fabs(z*bi->n[6] - x*bi->n[8]) > bi->k[8]) return(0);
    return(1);
}

static int pgram3d_gethit (pgram3d_t *bi, point3d *c0, point3d *c1, point3d *dpos, point3d *hit, point3d *norm)
{
    point3d e0, ev, f0, fv;
    double fk, f, s, t, x, y, z, dx, dy, dz, bx, by, bz, cx, cy, cz, bt, d2, d2min, rb;
    int i, j, besti, k0, k1, k2, k3;

    dx = c1->x-c0->x;
    dy = c1->y-c0->y;
    dz = c1->z-c0->z; d2min = 1e32;
    for(i=0;i<15;i++)
    {
              if (i < 3) { j =  i   *3; x = bi->m1[j]; y = bi->m1[j+1]; z = bi->m1[j+2]; fk = bi->m[i]; }
        else if (i < 6) { j = (i-3)*3; x = bi->m0[j]; y = bi->m0[j+1]; z = bi->m0[j+2]; fk = bi->m[i]; }
        else
        {
            j = i-6;
            k0 = (j/3)*3; k1 = (k0+3)%9; k2 = (((j/3)+1)%3) + (j%3)*3; k3 = (j/3)+(j%3)*3;
            x = bi->m0[k0  ]*bi->n[k2] - bi->m0[k1  ]*bi->n[k3];
            y = bi->m0[k0+1]*bi->n[k2] - bi->m0[k1+1]*bi->n[k3];
            z = bi->m0[k0+2]*bi->n[k2] - bi->m0[k1+2]*bi->n[k3];
            fk = bi->k[j];
        }
        f = x*x + y*y + z*z; t = dx*x + dy*y + dz*z; if (fabs(t) < 1e-16) continue;
        if (t >= 0) s = fk; else s = -fk;
        t = (s-t) / f;
        d2 = t*t*f; if (d2 < d2min) { d2min = d2; bx = x; by = y; bz = z; bt = t; besti = i; }
    }
    if (d2min >= 1e32) return(0);

    bx *= bt; dpos->x = bx;
    by *= bt; dpos->y = by;
    bz *= bt; dpos->z = bz;

    if (besti < 3)
    {
        cx = c0->x; cy = c0->y; cz = c0->z;
        s = bx*bi->r0.x + by*bi->r0.y + bz*bi->r0.z; s = (s>=0.0)*2-1.0; cx += bi->r0.x*s; cy += bi->r0.y*s; cz += bi->r0.z*s;
        s = bx*bi->d0.x + by*bi->d0.y + bz*bi->d0.z; s = (s>=0.0)*2-1.0; cx += bi->d0.x*s; cy += bi->d0.y*s; cz += bi->d0.z*s;
        s = bx*bi->f0.x + by*bi->f0.y + bz*bi->f0.z; s = (s>=0.0)*2-1.0; cx += bi->f0.x*s; cy += bi->f0.y*s; cz += bi->f0.z*s;
        bx *= -1; by *= -1; bz *= -1;

            //Make sure (cx,cy,cz) is inside both paralellepipeds..
        x = (cx-c1->x)*bi->m1[0] + (cy-c1->y)*bi->m1[1] + (cz-c1->z)*bi->m1[2]; x = std::min(std::max(x,-1.0),1.0);
        y = (cx-c1->x)*bi->m1[3] + (cy-c1->y)*bi->m1[4] + (cz-c1->z)*bi->m1[5]; y = std::min(std::max(y,-1.0),1.0);
        z = (cx-c1->x)*bi->m1[6] + (cy-c1->y)*bi->m1[7] + (cz-c1->z)*bi->m1[8]; z = std::min(std::max(z,-1.0),1.0);
        cx = bi->r1.x*x + bi->d1.x*y + bi->f1.x*z + c1->x;
        cy = bi->r1.y*x + bi->d1.y*y + bi->f1.y*z + c1->y;
        cz = bi->r1.z*x + bi->d1.z*y + bi->f1.z*z + c1->z;
    }
    else if (besti < 6)
    {
        cx = c1->x; cy = c1->y; cz = c1->z;
        s = bx*bi->r1.x + by*bi->r1.y + bz*bi->r1.z; s = (s>=0.0)*2-1.0; cx -= bi->r1.x*s; cy -= bi->r1.y*s; cz -= bi->r1.z*s;
        s = bx*bi->d1.x + by*bi->d1.y + bz*bi->d1.z; s = (s>=0.0)*2-1.0; cx -= bi->d1.x*s; cy -= bi->d1.y*s; cz -= bi->d1.z*s;
        s = bx*bi->f1.x + by*bi->f1.y + bz*bi->f1.z; s = (s>=0.0)*2-1.0; cx -= bi->f1.x*s; cy -= bi->f1.y*s; cz -= bi->f1.z*s;

            //Make sure (cx,cy,cz) is inside both paralellepipeds..
        x = (cx-c0->x)*bi->m0[0] + (cy-c0->y)*bi->m0[1] + (cz-c0->z)*bi->m0[2]; x = std::min(std::max(x,-1.0),1.0);
        y = (cx-c0->x)*bi->m0[3] + (cy-c0->y)*bi->m0[4] + (cz-c0->z)*bi->m0[5]; y = std::min(std::max(y,-1.0),1.0);
        z = (cx-c0->x)*bi->m0[6] + (cy-c0->y)*bi->m0[7] + (cz-c0->z)*bi->m0[8]; z = std::min(std::max(z,-1.0),1.0);
        cx = bi->r0.x*x + bi->d0.x*y + bi->f0.x*z + c0->x;
        cy = bi->r0.y*x + bi->d0.y*y + bi->f0.y*z + c0->y;
        cz = bi->r0.z*x + bi->d0.z*y + bi->f0.z*z + c0->z;
    }
    else
    {
        i = (besti%3); e0 = (*c1);
        if (i == 0) ev = bi->r1; else { s = bx*bi->r1.x + by*bi->r1.y + bz*bi->r1.z; s = (s>=0.0)*2-1.0; e0.x -= bi->r1.x*s; e0.y -= bi->r1.y*s; e0.z -= bi->r1.z*s; }
        if (i == 1) ev = bi->d1; else { s = bx*bi->d1.x + by*bi->d1.y + bz*bi->d1.z; s = (s>=0.0)*2-1.0; e0.x -= bi->d1.x*s; e0.y -= bi->d1.y*s; e0.z -= bi->d1.z*s; }
        if (i == 2) ev = bi->f1; else { s = bx*bi->f1.x + by*bi->f1.y + bz*bi->f1.z; s = (s>=0.0)*2-1.0; e0.x -= bi->f1.x*s; e0.y -= bi->f1.y*s; e0.z -= bi->f1.z*s; }
        i = (((besti-6)/3)+2)%3; f0 = (*c0);
        if (i == 0) fv = bi->r0; else { s = bx*bi->r0.x + by*bi->r0.y + bz*bi->r0.z; s = (s>=0.0)*2-1.0; f0.x += bi->r0.x*s; f0.y += bi->r0.y*s; f0.z += bi->r0.z*s; }
        if (i == 1) fv = bi->d0; else { s = bx*bi->d0.x + by*bi->d0.y + bz*bi->d0.z; s = (s>=0.0)*2-1.0; f0.x += bi->d0.x*s; f0.y += bi->d0.y*s; f0.z += bi->d0.z*s; }
        if (i == 2) fv = bi->f0; else { s = bx*bi->f0.x + by*bi->f0.y + bz*bi->f0.z; s = (s>=0.0)*2-1.0; f0.x += bi->f0.x*s; f0.y += bi->f0.y*s; f0.z += bi->f0.z*s; }
        dx = by*ev.z - bz*ev.y; //f0.x + fv.x*i + bx*b = e0.x + ev.x*e  |  (fv.x)*i + (-ev.x)*e + (bx)*b = (e0.x-f0.x)
        dy = bz*ev.x - bx*ev.z; //f0.y + fv.y*i + by*b = e0.y + ev.y*e  |  (fv.y)*i + (-ev.y)*e + (by)*b = (e0.y-f0.y)
        dz = bx*ev.y - by*ev.x; //f0.z + fv.z*i + bz*b = e0.z + ev.z*e  |  (fv.z)*i + (-ev.z)*e + (bz)*b = (e0.z-f0.z)
        f = ((e0.x-f0.x)*dx + (e0.y-f0.y)*dy + (e0.z-f0.z)*dz) / (fv.x*dx + fv.y*dy + fv.z*dz);
        cx = fv.x*f + f0.x; bx *= -1.0;
        cy = fv.y*f + f0.y; by *= -1.0;
        cz = fv.z*f + f0.z; bz *= -1.0;
    }

    hit->x = bx*.5 + cx;
    hit->y = by*.5 + cy;
    hit->z = bz*.5 + cz;
    rb = bx*bx + by*by + bz*bz;
    if (rb == 0.0) { norm->x = 0.0; norm->y = 1.0; norm->z = 0.0; return(1); } //FIXFIXFIXFIX:wrong hack!
    rb = 1.0/sqrt(rb);
    norm->x = bx*rb;
    norm->y = by*rb;
    norm->z = bz*rb;
    return(1);
}

#endif // TYPES_PARALLELEPIPED_HPP
