// This file has been modified from Ken Silverman's original release

#pragma once
#ifndef BRUSHES_HPP
#define BRUSHES_HPP

#include "types/point.hpp"
#include "types/surface.hpp"
#include "octree.hpp"

//---------------------------------------------- Brush ----------------------------------------------
// isins() must return: 0:no int, 1:partial int, 2:total int
// getsurf() must write surf structure
// flags bit0:1=override normals at voxel boundaries

#define BRUSH_HEADER                                                        \
    int  (*isins  )(brush_t *brush, int x0, int y0, int z0, int log2sid);   \
    void (*getsurf)(brush_t *brush, int x0, int y0, int z0, surf_t *surf);  \
    int mx0, my0, mz0, mx1, my1, mz1, flags

// generic brush used by mod functions
#pragma pack(push,1)
struct brush_t {
    BRUSH_HEADER;
    /*add user data here*/
};
#pragma pack(pop)

#pragma pack(push,1)
struct brush_vox_t {
    BRUSH_HEADER;
    int x, y, z;
    surf_t surf;
};
#pragma pack(pop)

void brush_vox_init (brush_vox_t *vox, int x, int y, int z, surf_t *surf);

#pragma pack(push,1)
struct brush_sph_t {
    BRUSH_HEADER;
    int x, y, z, r2, col;
    surf_t *surf;
};
#pragma pack(pop)

void brush_sph_init (brush_sph_t *sph, int x, int y, int z, int r, int issol);

#pragma pack(push,1)
struct brush_box_t {
    BRUSH_HEADER;
    int x0, y0, z0, x1, y1, z1, col;
};
#pragma pack(pop)

void brush_box_init(brush_box_t* box, float x0, float y0, float z0, float x1, float y1, float z1, int issol);

#pragma pack(push,1)
struct brush_cone_t {
    BRUSH_HEADER;
    float x0, y0, z0, r0, x1, y1, z1, r1, cx, cy, cz, cr, dx, dy, dz, dr, r02, r12, k0, k1, hak;
    int col;
};
#pragma pack(pop)

void brush_cone_init (brush_cone_t *cone, float x0, float y0, float z0, float r0, float x1, float y1, float z1, float r1);

#pragma pack(push,1)
struct brush_bmp_t {
    BRUSH_HEADER;
    unsigned short *boxsum;
    int xs, ys, zs, iux, iuy, iuz, ivx, ivy, ivz, iwx, iwy, iwz;
    int iox0[OCT_MAXLS], ioy0[OCT_MAXLS], ioz0[OCT_MAXLS], iox1[OCT_MAXLS], ioy1[OCT_MAXLS], ioz1[OCT_MAXLS];
};
#pragma pack(pop)

void brush_bmp_init (brush_bmp_t *bmp, unsigned short *boxsum, int xs, int ys, int zs, point3f *pp, point3f *pr, point3f *pd, point3f *pf);

#endif
