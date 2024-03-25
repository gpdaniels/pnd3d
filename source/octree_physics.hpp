// This file has been modified from Ken Silverman's original release

#pragma once
#ifndef OCTREE_PHYSICS_HPP
#define OCTREE_PHYSICS_HPP

#include "octree.hpp"
#include "brushes.hpp"
#include "types/point.hpp"
#include "types/surface.hpp"

// ?
int oct_hitscan(oct_t *oct, point3f *p, point3f *padd, point3i *hit, int *rhitdir, float  *fracwent);

// ?
int oct_hitscan(oct_t *oct, point3d *p, point3d *padd, point3i *hit, int *rhitdir, double *fracwent);

// ?
surf_t *oct_sphtrace(oct_t *oct, point3f *p, point3f *padd, float rad, point3i *hit, point3f *pgoal, point3f *hitnorm);

// ?
surf_t *oct_sphtrace(oct_t *oct, point3d *p, point3d *padd, double rad, point3i *hit, point3d *pgoal, point3d *hitnorm);

// ?
double oct_balloonrad(oct_t *oct,  point3f *p, double cr, point3i *hit, surf_t **hitsurf);

// ?
double oct_balloonrad(oct_t *oct, point3d *p, double cr, point3i *hit, surf_t **hitsurf);

// ?
int oct_slidemove(oct_t *oct, point3f *p, point3f *padd, double fat, point3f *pgoal);

// ?
int oct_slidemove(oct_t *oct, point3d *p, point3d *padd, double fat, point3d *pgoal);

#define USENEWCAC 0
#if(USENEWCAC == 0)
    // cache structure for repeated use of estnorm
    #pragma pack(push,1)
    struct oct_bitcac_t {
        point3i pt;
        unsigned int buf[32*32];
    };
    #pragma pack(pop)
#else
    // cache structure for repeated use of estnorm
    #define CACHASHN 32
    #define CACN 32
    #pragma pack(push,1)
    struct oct_bitcac_t {
        point3i pt[CACN];
        int cur;
        int ptn[CACN];
        int n;
        int hashead[CACHASHN];
        unsigned int buf[CACN][32*32];
    };
    #pragma pack(pop)
#endif

// ?
void oct_bitcac_reset(oct_bitcac_t *);

// ?
void oct_estnorm(oct_t *oct, int x, int y, int z, int r,  point3f *norm, oct_bitcac_t *bitcac);

// ?
void oct_estnorm(oct_t *oct, int x, int y, int z, int r, point3d *norm, oct_bitcac_t *bitcac);

// ?
void oct_refreshnorms(oct_t *oct, int dist, int x0, int y0, int z0, int x1, int y1, int z1);

// ?
void oct_getmoi(oct_t *oct, float *mas, float *cx, float *cy, float *cz, float *ixx, float *iyy, float *izz, float *ixy, float *ixz, float *iyz);

// ?
void oct_getmoi(oct_t *oct, double *mas, double *cx, double *cy, double *cz, double *ixx, double *iyy, double *izz, double *ixy, double *ixz, double *iyz);

// returns: 0:pure air, 1:some air&sol, 2:pure sol
int oct_getboxstate(oct_t *oct, int bx0, int by0, int bz0, int bx1, int by1, int bz1);

// find bounding box of solid inside user bbox (params are in&out)
int oct_getsolbbox(oct_t *oct, int *x0, int *y0, int *z0, int *x1, int *y1, int *z1);

// ?
void oct_sol2bit(oct_t *oct, unsigned int *bitvis, int bx0, int by0, int bz0, int dx, int dy, int dz, int lsdax);

// ?
int oct_touch_brush(oct_t *oct, brush_t *brush);

// octree node hit(ls is log(2) of node side;0=1x1x1, 1=2x2x2, 2=4x4x4, ..)
#pragma pack(push,1)
struct oct_hit_t {
    int x, y, z, ls;
};
#pragma pack(pop)

// hit = 1st intersecting node per octree(array of 2!)
int oct_touch_oct(oct_t *oct0, point3f *p0, point3f *r0, point3f *d0, point3f *f0, oct_t *oct1, point3f *p1, point3f *r1, point3f *d1, point3f *f1, oct_hit_t *hit);

// hit = 1st intersecting node per octree(array of 2!)
int oct_touch_oct(oct_t *oct0, point3d *p0, point3d *r0, point3d *d0, point3d *f0, oct_t *oct1, point3d *p1, point3d *r1, point3d *d1, point3d *f1, oct_hit_t *hit);

#endif // OCTREE_PHYSICS_HPP
