// This file has been modified from Ken Silverman's original release

#pragma once
#ifndef MATHS_HPP
#define MATHS_HPP

#include "macros.hpp"
#include "types/point.hpp"

#include <cmath>

#if 1
    __forceinline
    float klog (float f) {
        // Super fast natural log, error = +/-0.02983002966476566
        // ORIGINAL: return static_cast<float>(((float)((*(int*)&f) - 0x3f7a7dcf)) * 8.262958294867817E-8);
        return static_cast<float>(static_cast<double>(static_cast<float>((*reinterpret_cast<int*>(&f)) - 0x3F7A7DCF)) * 8.262958294867817E-8);
    }
#else
    __forceinline
    float klog (float f) {
        // Fast natural log, error = +/-0.003423966993
        // ORIGINAL: float g = ((float)(((*(int*)&f) & 8388607) - 4074142)) * 5.828231702537851e-8;
        // ORIGINAL: return static_cast<float>(((float)(*(int*)&f)) * 8.262958294867817E-8 - g * g - 87.96988524938206);
        float g = static_cast<float>(static_cast<double>(static_cast<float>(((*reinterpret_cast<int*>(&f)) & 8388607) - 4074142)) * 5.828231702537851E-8);
        return static_cast<float>(static_cast<double>(static_cast<float>(*reinterpret_cast<int*>(&f))) * 8.262958294867817E-8 - static_cast<double>(g * g) - 87.96988524938206);
    }
#endif

__forceinline
float klog2(float f) {
    // ORIGINAL: return static_cast<float>(klog(f) * 1.442695040888963);
    return static_cast<float>(klog(f) * 1.442695040888963f);
}

__forceinline
float klog10(float f) {
    // ORIGINAL: return static_cast<float>(klog(f) * 0.4342944819032518);
    return static_cast<float>(klog(f) * 0.4342944819032518f);
}

__forceinline
float klog2up7 (float f) {
    // Same as the super fast natural log, but with the up7 multiplication done first. Returns: klog2(f) * 128.0
    // ORIGINAL: return static_cast<float>(((float)((*(int *)&f) - 0x3F7A7DCF)) * (8.262958294867817E-8 * 1.442695040888963 * 128.0));
    return static_cast<float>(static_cast<double>(static_cast<float>((*reinterpret_cast<int*>(&f)) - 0x3F7A7DCF)) * (8.262958294867817E-8 * 1.442695040888963 * 128.0));
}

template <typename type>
inline void invert3x3(point<type, 3> *r, point<type, 3> *d, point<type, 3> *f, type *mat)
{
    type rdet;

    mat[0] = d->y * f->z - d->z * f->y;
    mat[1] = d->z * f->x - d->x * f->z;
    mat[2] = d->x * f->y - d->y * f->x;
    rdet = 1.0 / (mat[0] * r->x + mat[1] * r->y + mat[2] * r->z);
    mat[0] *= rdet;
    mat[1] *= rdet;
    mat[2] *= rdet;
    mat[3] = (f->y * r->z - f->z * r->y) * rdet;
    mat[4] = (f->z * r->x - f->x * r->z) * rdet;
    mat[5] = (f->x * r->y - f->y * r->x) * rdet;
    mat[6] = (r->y * d->z - r->z * d->y) * rdet;
    mat[7] = (r->z * d->x - r->x * d->z) * rdet;
    mat[8] = (r->x * d->y - r->y * d->x) * rdet;
}

//optimization for (r,d,f) symmetrix matrix
template <typename type>
inline void invert3x3sym(point<type, 3> *r, point<type, 3> *d, point<type, 3> *f, type *mat)
{
    type rdet;

    mat[0] = d->y * f->z - d->z * d->z;
    mat[1] = r->z * d->z - r->y * f->z;
    mat[2] = r->y * d->z - r->z * d->y;
    rdet = 1.0 / (mat[0] * r->x + mat[1] * r->y + mat[2] * r->z);
    mat[0] *= rdet;
    mat[1] *= rdet;
    mat[2] *= rdet;
    mat[5] = (r->y * r->z - r->x * d->z) * rdet;
    mat[4] = (r->x * f->z - r->z * r->z) * rdet;
    mat[8] = (r->x * d->y - r->y * r->y) * rdet;
    mat[3] = mat[1];
    mat[6] = mat[2];
    mat[7] = mat[5];
}

// Rotate vectors a & b around their common plane, by ang
template <typename type>
void rotvex(type ang, point3<type> *a, point3<type> *b) {
    type c, s, f;
    c = cos(ang); s = sin(ang);
    f = a->x; a->x = f*c + b->x*s; b->x = b->x*c - f*s;
    f = a->y; a->y = f*c + b->y*s; b->y = b->y*c - f*s;
    f = a->z; a->z = f*c + b->z*s; b->z = b->z*c - f*s;
}

inline void vecadd   (point3d *c, point3d *a, point3d *b)
{
    c->x = a->x+b->x;
    c->y = a->y+b->y;
    c->z = a->z+b->z;
}

inline void vecsub   (point3d *c, point3d *a, point3d *b)
{
    c->x = a->x-b->x;
    c->y = a->y-b->y;
    c->z = a->z-b->z;
}

inline void vecscale (point3d *c, point3d *a, double sc)
{
    c->x = a->x*sc;
    c->y = a->y*sc;
    c->z = a->z*sc;
}

inline double vecdot (point3d *a, point3d *b)
{
    return(a->x*b->x + a->y*b->y + a->z*b->z);
}

inline void veccross (point3d *c, point3d *a, point3d *b) //C = A x B
{
    point3d nc;
    nc.x = a->y*b->z - a->z*b->y;
    nc.y = a->z*b->x - a->x*b->z;
    nc.z = a->x*b->y - a->y*b->x;
    (*c) = nc;
}

inline void matvecmul (point3d *c, double *a, point3d *b)
{
    point3d nc;
    nc.x = a[0]*b->x + a[1]*b->y + a[2]*b->z;
    nc.y = a[3]*b->x + a[4]*b->y + a[5]*b->z;
    nc.z = a[6]*b->x + a[7]*b->y + a[8]*b->z;
    (*c) = nc;
}

inline void vecmatmul (point3d *c, point3d *a, double *b)
{
    point3d nc;
    nc.x = a->x*b[0] + a->y*b[3] + a->z*b[6];
    nc.y = a->x*b[1] + a->y*b[4] + a->z*b[7];
    nc.z = a->x*b[2] + a->y*b[5] + a->z*b[8];
    (*c) = nc;
}

inline void simxform (double *b, double *r, double *a) //B = R*A*R^T
{
    double t[9];

        //[t0 t1 t2]   [r0 r1 r2][a0 a1 a2]
        //[t3 t4 t5] = [r3 r4 r5][a3 a4 a5]
        //[t6 t7 t8]   [r6 r7 r8][a6 a7 a8]
    t[0] = r[0]*a[0] + r[1]*a[3] + r[2]*a[6];
    t[1] = r[0]*a[1] + r[1]*a[4] + r[2]*a[7];
    t[2] = r[0]*a[2] + r[1]*a[5] + r[2]*a[8];
    t[3] = r[3]*a[0] + r[4]*a[3] + r[5]*a[6];
    t[4] = r[3]*a[1] + r[4]*a[4] + r[5]*a[7];
    t[5] = r[3]*a[2] + r[4]*a[5] + r[5]*a[8];
    t[6] = r[6]*a[0] + r[7]*a[3] + r[8]*a[6];
    t[7] = r[6]*a[1] + r[7]*a[4] + r[8]*a[7];
    t[8] = r[6]*a[2] + r[7]*a[5] + r[8]*a[8];

        //[b0 b1 b2]   [t0 t1 t2][r0 r3 r6]
        //[b3 b4 b5] = [t3 t4 t5][r1 r4 r7]
        //[b6 b7 b8]   [t6 t7 t8][r2 r5 r8]
    b[0] = t[0]*r[0] + t[1]*r[1] + t[2]*r[2];
    b[1] = t[0]*r[3] + t[1]*r[4] + t[2]*r[5];
    b[2] = t[0]*r[6] + t[1]*r[7] + t[2]*r[8];
    b[3] = t[3]*r[0] + t[4]*r[1] + t[5]*r[2];
    b[4] = t[3]*r[3] + t[4]*r[4] + t[5]*r[5];
    b[5] = t[3]*r[6] + t[4]*r[7] + t[5]*r[8];
    b[6] = t[6]*r[0] + t[7]*r[1] + t[8]*r[2];
    b[7] = t[6]*r[3] + t[7]*r[4] + t[8]*r[5];
    b[8] = t[6]*r[6] + t[7]*r[7] + t[8]*r[8];
}

#endif // MATHS_HPP
