// This file has been modified from Ken Silverman's original release

#pragma once
#ifndef GRAPHICS_HPP
#define GRAPHICS_HPP

#include "utilities/macros.hpp"

struct tiletype {
    // memory?
    union {
        INT_PTR f;
        char* framebuffer;
    };

    // bytes per line (aka for RGBA it's 4 x width)
    union {
        int p;
        int stride;
    };

    union {
        int x;
        int width;
    };

    union {
        int y;
        int height;
    };
};

void print6x8(tiletype *dd, int ox, int y, int colour_foreground, int colour_background, const char* format, ...);
void drawpix(tiletype *dd, int x, int y, int colour);
void drawhlin(tiletype *dd, int x0, int x1, int y, int colour);
void drawline(tiletype *dd, float x0, float y0, float x1, float y1, int colour);
void drawcirc(tiletype *dd, int xc, int yc, int r, int colour);
void drawrectfill(tiletype *dd, int x0, int y0, int x1, int y1, int colour);
void drawcircfill(tiletype *dd, int xc, int yc, int r, int colour);

#endif // GRAPHICS_HPP
