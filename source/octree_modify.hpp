// This file has been modified from Ken Silverman's original release

#pragma once
#ifndef OCTREE_MODIFY_HPP
#define OCTREE_MODIFY_HPP

#include "octree.hpp"
#include "brushes.hpp"

// mode info:
// (mode&1)==0:brush is air, !=0:brush is solid
// (mode&2)!=0:do oct_updatesurfs() here (nice because bounding box calculated inside)
// (mode&4)!=0:hover check

// mods existing oct; faster than oct_setvox() for large volumes
void oct_mod(oct_t *oct, brush_t *brush, int mode);

// like oct_mod but modify surfaces only
void oct_paint(oct_t *oct, brush_t *brush);

// mode==0: newoct =  oct & brush          //most useful - extracts piece without requiring full copy
// mode==1: newoct = (oct & brush)|~brush  //not useful - very unusual bool op
// mode==2: newoct =  oct &~brush          //useful, but copy&oct_mod() also works
// mode==3: newoct =  oct | brush          //useful, but copy&oct_mod() also works
// doesn't touch existing oct; uses oct->recvoctfunc to generate new oct_t
void oct_modnew(oct_t *oct, brush_t *brush, int mode);

void oct_mod_recur(oct_t *loct, int inode, int x0, int y0, int z0, int ls, octv_t *roct, brush_t *brush, int issol);

void oct_updatesurfs(oct_t* oct, int x0, int y0, int z0, int x1, int y1, int z1, brush_t* brush, int mode);

// returns: 0:air, 1:surface (surf written), 2:interior
int oct_getvox(oct_t *oct, int x, int y, int z, surf_t **surf);

// faster than oct_getvox() if only sol vs. air needed; returns: 0:air, 1:solid
int oct_getsol(oct_t *oct, int x, int y, int z);

// returns: -1:invalid, {0..sur.mal-1}:surface index
int oct_getsurf(oct_t *oct, int x, int y, int z);

// equivalent to using brush_vox_t with oct_mod()
void oct_setvox(oct_t *oct, int x, int y, int z, surf_t  *surf, int mode);

// proper way to write a surface - this function updates GPU memory
void oct_writesurf(oct_t *oct, int ind, surf_t *psurf);

// copy all surfs from CPU to GPU memory
void oct_copysurfs(oct_t *oct);

// searches on increasing y, returns y and surf, or oct->sid if hit nothing
int oct_findsurfdowny(oct_t *oct, int x, int y, int z, surf_t **surf);

// ?
void oct_hover_check(oct_t *oct, int x0, int y0, int z0, int x1, int y1, int z1, void (*recvoctfunc)(oct_t *ooct, oct_t *noct));

// change octree size in-place
int oct_rebox(oct_t *oct, int x0, int y0, int z0, int x1, int y1, int z1, int *dx, int *dy, int *dz);

// re-orient axes; x=1, y=2, z=3, -x=-1, -y=-2, -z=-3
void oct_swizzle(oct_t *oct, int ax0, int ax1, int ax2);

// ?
void oct_checkreducesizes(oct_t *loct);

#pragma pack(push,1)
struct oct_getvox_hint_t {
    int stkind[OCT_MAXLS];
    int minls;
    int mins;
    int ox;
    int oy;
    int oz;
};
#pragma pack(pop)

void oct_getvox_hint_init(oct_t *oct, oct_getvox_hint_t *och);

//WARNING:assumes x,y,z inside grid
int oct_getsol_hint(oct_t *oct, int x, int y, int z, oct_getvox_hint_t *och);

// WARNING:assumes x,y,z even #'s
int oct_getsol_hint_2x2x2(oct_t *loct, int x, int y, int z, oct_getvox_hint_t *och);

//returns: -1:invalid, {0..sur.mal-1}:surface index
int oct_getsurf_hint(oct_t *oct, int x, int y, int z, oct_getvox_hint_t *och);

#endif // OCTREE_MODIFY_HPP
