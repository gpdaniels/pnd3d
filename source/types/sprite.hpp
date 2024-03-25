// This file has been modified from Ken Silverman's original release

#pragma once
#ifndef TYPES_SPRITE_HPP
#define TYPES_SPRITE_HPP

#include "octree.hpp"
#include "types/point.hpp"
#include "utilities/maths.hpp"

#pragma pack(push,1)
struct spr_t {
    oct_t *oct;
    point3d p, v, ax;   //xlat&inst.ang vel
    double ori[9];      //unit orthonormal rotation matrix, [br,bd,bf]*ori = [r,d,f] body->world space
    point3d br, bd, bf; //pgram, body (moi) space
    point3d r, d, f;    //pgram, world space
    point3d cen;
    double mas, rmas, moi[9], rmoi[9];

    float mixval;
    int imulcol;
    int cnt;
    double tim;
};
#pragma pack(pop)

// transform pos
template <typename type>
void oct_vox2worldpos (spr_t *sp, point<type, 3> *pvox, point<type, 3> *pworld) {
    point<type, 3> opvox = *pvox;
    pworld->x = opvox.x*sp->r.x + opvox.y*sp->d.x + opvox.z*sp->f.x + sp->p.x;
    pworld->y = opvox.x*sp->r.y + opvox.y*sp->d.y + opvox.z*sp->f.y + sp->p.y;
    pworld->z = opvox.x*sp->r.z + opvox.y*sp->d.z + opvox.z*sp->f.z + sp->p.z;
}

//transform vec (ignore origin)
template <typename type>
void oct_vox2worlddir (spr_t *sp, point<type, 3> *pvox, point<type, 3> *pworld) {
    point<type, 3> opvox = *pvox;
    pworld->x = opvox.x*sp->r.x + opvox.y*sp->d.x + opvox.z*sp->f.x;
    pworld->y = opvox.x*sp->r.y + opvox.y*sp->d.y + opvox.z*sp->f.y;
    pworld->z = opvox.x*sp->r.z + opvox.y*sp->d.z + opvox.z*sp->f.z;
}

// transform pos
template <typename type>
void oct_world2voxpos (spr_t *sp, point<type, 3> *pworld, point<type, 3> *pvox) {
    type fx, fy, fz, mat[9];

    invert3x3(&sp->r,&sp->d,&sp->f,mat);

    fx = pworld->x - sp->p.x;
    fy = pworld->y - sp->p.y;
    fz = pworld->z - sp->p.z;

    pvox->x = fx*mat[0] + fy*mat[1] + fz*mat[2];
    pvox->y = fx*mat[3] + fy*mat[4] + fz*mat[5];
    pvox->z = fx*mat[6] + fy*mat[7] + fz*mat[8];
}

// transform vec (ignore origin)
template <typename type>
void oct_world2voxdir (spr_t *sp, point<type, 3> *pworld, point<type, 3> *pvox)
{
    type fx, fy, fz, mat[9];

    invert3x3(&sp->r,&sp->d,&sp->f,mat);

    fx = pworld->x;
    fy = pworld->y;
    fz = pworld->z;

    pvox->x = fx*mat[0] + fy*mat[1] + fz*mat[2];
    pvox->y = fx*mat[3] + fy*mat[4] + fz*mat[5];
    pvox->z = fx*mat[6] + fy*mat[7] + fz*mat[8];
}

#endif // TYPES_SPRITE_HPP
