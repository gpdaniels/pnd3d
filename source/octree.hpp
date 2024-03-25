// This file has been modified from Ken Silverman's original release

#pragma once
#ifndef OCTREE_HPP
#define OCTREE_HPP

#include "types/point.hpp"
#include "types/surface.hpp"
#include "utilities/macros.hpp"
#include "utilities/window.hpp"

#define LOCT_MAXLS 4
#define OCT_MAXLS (1<<LOCT_MAXLS) //should be >= loctsid

#pragma pack(push,1)
struct octv_t {
    unsigned char chi;
    unsigned char sol;
    unsigned char mrk;
    unsigned char mrk2;
    int ind;
};
#pragma pack(pop)

//buf: cast to: octv_t* or surf_t*
//bit: 1 bit per sizeof(buf[0]); 0=free, 1=occupied
#pragma pack(push,1)
struct bitmal_t {
    void *buf;
    unsigned int mal;
    unsigned int *bit;
    unsigned int ind;
    unsigned int num;
    unsigned int siz;
};
#pragma pack(pop)

//internal structure; most fields not useful to caller
#pragma pack(push,1)
struct oct_t {
    int head;
    int lsid;
    int sid;
    int nsid;
    bitmal_t nod;
    bitmal_t sur;

    //assume neighboring voxels at boundary for oct_mod()
    char edgeiswrap; //bits 5-0: 0=air/sol, 1=wrap
    char edgeissol;  //bits 5-0: if (!edgeiswrap) { 0=air; 1=sol; } else reserved;

    //pnd3d private (don't modify these vars)
    unsigned int octid;  //gl texture indices
    unsigned int tilid;  //gl texture indices
    unsigned int bufid;  //gl buffer indices
    int glxsid;
    int glysid;
    int gxsid;
    int gysid;
    surf_t *gsurf;

    void (*recvoctfunc)(oct_t *ooct, oct_t *noct);
    int flags; //optional vars

    int xres;
    int yres;
};
#pragma pack(pop)

//use ~256.0 for 1 world unit = sprite; use ~65536.0 for 1 world unit = 1 voxel
int  oct_initonce(HWND ghwnd, float zfar);

// ?
void oct_uninitonce(HWND ghwnd);

// ? ex: loctsid = 8 for 256^3
void oct_new(int xres, int yres, oct_t *oct, int loctsid, INT_PTR tilid, int startnodes, int startsurfs, int hax4mark2spr);

#if defined(ENABLE_LOADING_ROUTINES)
// ? supports:KVO,KV6,KVX,VOX,VXL,PNG
int  oct_load(oct_t *oct, char *filnam, point3f *p, point3f *r, point3f *d, point3f *f);

// ?
int  oct_load(oct_t *oct, char *filnam, point3d *p, point3d *r, point3d *d, point3d *f);
#endif

#if defined(ENABLE_SAVING_ROUTINES)
// ? supports:KVO,KV6
void oct_save(oct_t *oct, char *filnam, point3f *p, point3f *r, point3f *d, point3f *f);

// ?
void oct_save(oct_t *oct, char *filnam, point3d *p, point3d *r, point3d *d, point3d *f);
#endif

#if defined(ENABLE_UNUSED_ROUTINES)
// ?
void oct_dup(oct_t *oct, oct_t *newoct);
#endif

// ?
void oct_free(oct_t *oct);

#if defined(ENABLE_SAVING_ROUTINES)
// SAVING ROUTINES:
void savekv6(oct_t *loct, char *filnam, point3f *rpos, point3f *rrig, point3f *rdow, point3f *rfor);
void savekvo(oct_t *loct, char *filnam, point3f *rpos, point3f *rrig, point3f *rdow, point3f *rfor);
#endif

#endif // OCTREE_HPP
