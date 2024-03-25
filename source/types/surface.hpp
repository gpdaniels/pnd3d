// This file has been modified from Ken Silverman's original release

#pragma once
#ifndef TYPES_SURFACE_HPP
#define TYPES_SURFACE_HPP

#pragma pack(push,1)
struct surf_t {
    unsigned char b;
    unsigned char g;
    unsigned char r;
    unsigned char a;
    signed char norm[3];
    signed char tex;
};
#pragma pack(pop)

#endif // TYPES_SURFACE_HPP
