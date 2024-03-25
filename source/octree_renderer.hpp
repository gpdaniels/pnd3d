// This file has been modified from Ken Silverman's original release

#pragma once
#ifndef OCTREE_RENDERER_HPP
#define OCTREE_RENDERER_HPP

#include "brushes.hpp"
#include "octree.hpp"
#include "types/point.hpp"
#include "types/surface.hpp"
#include "utilities/macros.hpp"
#include "utilities/window.hpp"

#include <cstdint>

// use world coords
extern float oct_fogdist;
extern unsigned int oct_fogcol;
extern int oct_sideshade[6];

#if defined(_WIN32)
    #define NOMINMAX
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
#endif

extern HDC glhDC;
extern HGLRC glhRC;

extern unsigned int shadcur;
#define SHADNUM 5
//[0]=oct_*, [1]=drawcone, [2]=drawsky, [3]=drawtext, [4]=drawpol(need shader to brighten :/)
extern unsigned int shadprog[SHADNUM];
extern unsigned int shadvert[SHADNUM];
extern unsigned int shadfrag[SHADNUM];

extern void (*shaderfunc)(int,void*);
extern int swapinterval; //0=max, 1=60fps, 2=30fps, ..

extern double gznear;
extern double gzfar;

extern char vshad_drawoct[];
extern char fshad_drawoct[];
extern char vshad_drawsky[];
extern char fshad_drawsky[];
extern char vshad_drawtext[];
extern char fshad_drawtext[];

extern char vshadasm_drawoct[];
extern char fshadasm_drawoct[];
extern char vshadasm_drawsky[];
extern char fshadasm_drawsky[];
extern char vshadasm_drawtext[];
extern char fshadasm_drawtext[];

#define PIXBUFBYPP 4

//glBindTexture   id's (see also oct_t.octid, oct_t.tilid)
extern unsigned int gpixtexid;

//INT_PTR to support CPU as pointer & GPU as index for drawpol()
extern INT_PTR gskyid;

extern unsigned int gfont6x8id;
extern unsigned int gtextbufid;

//glBufferTexture id's (see also oct_t.gbufid)
extern unsigned int gpixbufid;

extern char *gpixbuf;
extern int gpixbufmal;
extern int gpixxdim;
extern int gpixydim;

void cpu_shader_solcol(int sy, void *_);
void cpu_shader_znotex(int sy, void *_);
void cpu_shader_texmap(int sy, void *_);
void cpu_shader_texmap_mc(int sy, void *_);

#define TXTBUFSIZ 1024
extern const uint64_t font6x8[];

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// All rendering must go between start & stop per frame
int oct_startdraw(window_state& window, point3f *ipos, point3f *irig, point3f *idow, point3f *ifor, float ghx, float ghy, float ghz);
int oct_startdraw(window_state& window, point3d *ipos, point3d *irig, point3d *idow, point3d *ifor, double ghx, double ghy, double ghz);
void oct_stopdraw(window_state& window);

// Simple 2D stuff
void oct_drawpix(int xres, int yres, int x, int y, int col);
void oct_drawline(int xres, int yres, float x0, float y0, float x1, float y1, int col);
void oct_drawtext6x8(int xres, int yres, int x0, int y0, int fcol, int bcol, const char *fmt, ...);

// Texture-mapped poly
// bitmap is copied & may be destroyed after call
// same as oct_loadtex() but with hacks for tile buffer
INT_PTR oct_loadtex(INT_PTR ptr, int xsiz, int ysiz, int flags);
void oct_freetex(INT_PTR ptr);

// Octree
void oct_drawoct(oct_t *oct, point3f *p, point3f *r, point3f *d, point3f *f, float mixval, int imulcol);
void oct_drawoct(oct_t *oct, point3d *p, point3d *r, point3d *d, point3d *f, float mixval, int imulcol);

#endif // OCTREE_RENDERER_HPP
