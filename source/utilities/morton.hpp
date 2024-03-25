// This file has been modified from Ken Silverman's original release

#pragma once
#ifndef MORTON_HPP
#define MORTON_HPP

#include "utilities/macros.hpp"

//02/10/2012: Adds round values to x&y as if they form a morton code; tracks carry bits using tricks to the extreme! :)
//* s is side length of square (area) to add (must be pow2)
//* assumes: s is pow2, (x%s) == 0, (y%s) == 0
//
//Algo: find next 0 bit, starting from s, and using morton order on mx&my
//   Performs ripple carry manually. x/y with 1st 0 is summed with s; other gets bits cleared up to that point.
//
//Ex#1:                                  | Ex#2:
//   000000000100 s                      |    000000000100 s
//   011010111100 x  ;clear bits 5-0     |    0110101o1100 x  ;add s
//   001111o11100 y  ;add s              |    001111011100 y  ;clear bits 3-0
//                                       |
//   011011000000 x+s                    |    011010110000 x+s
//   001111100000 y+s                    |    001111100000 y+s
//   100101000000 -(x+s)                 |    100101010000 -(x+s)
//   110000100000 -(y+s)                 |    110000100000 -(y+s)
//
//   111111000000 (-(x+s))|(x+s)         |    111111110000 (-(x+s))|(x+s)
//   111111000000 (-(y+s))^(y+s)         |    111111000000 (-(y+s))^(y+s)
//
#define mort2add(s,mx,my)\
{\
    int _x, _y;\
    (mx) += (s); _x = (-(mx))|(mx);\
    (my) += (s); _y = (-(my))^(my);\
    if (_y >= _x) (mx) += _y; else (my) += _x;\
}

__forceinline
static int mortcmp_fast (int x0, int y0, int z0, int x1, int y1, int z1)
{
    #if !defined(_WIN32)
        //7.83ns
        int i0, i1, i2, i3;
        i0 = x0^x1;
        i1 = y0^y1;
        i2 = z0^z1;
        i3 = i0|i1;
        return ( ((((z0-z1)>>31)&-4) - (((z1-z0)>>31)&-4)) & ((((i3^i2)&i3)-i2)>>31) ) +
               ( ((((y0-y1)>>31)&-2) - (((y1-y0)>>31)&-2)) & ((((i0^i1)&i0)-i1)>>31) ) +
               ( ((((x0-x1)>>31)   ) - (((x1-x0)>>31)   ))                           );
    #else
        //5.92ns
        alignas(16) static const int dqzero[4] = {0,0,0,0};
        _asm
        {
            movd xmm0, x0
            movd xmm1, x1
            movd xmm2, y0
            movd xmm3, y1
            punpckldq xmm0, xmm2
            punpckldq xmm1, xmm3
            movd xmm2, z0
            movd xmm3, z1
            movlhps xmm0, xmm2      ;xmm0:[ 0 z0 y0 x0]
            movlhps xmm1, xmm3      ;xmm1:[ 0 z1 y1 x1]

            movaps xmm2, xmm0
            pxor xmm2, xmm1         ;xmm2:[ 0 i2 i1 i0]

            pshufd xmm3, xmm2, 0x3f ;xmm3:[i0  0  0  0]
            pshufd xmm2, xmm2, 0x64 ;xmm3:[i1 i2 i1 i0]
            por xmm2, xmm3          ;xmm2:[i3 i2 i1 i0]

            pshufd xmm3, xmm2, 0xb1 ;xmm3:[i2 i3 i0 i1]
            movaps xmm4, xmm3       ;xmm4:[i2 i3 i0 i1]
            pxor xmm3, xmm2         ;xmm3:[i2^i3 i2^i3 i0^i1 i0^i1]
            pand xmm3, xmm4         ;xmm3:[(i2^i3)&i2 (i2^i3)&i3 (i0^i1)&i0 (i0^i1)&i1]
            pcmpgtd xmm2, xmm3      ;xmm2:[((i2^i3)&i2)<i3 ((i3^i2)&i3)<i2 ((i0^i1)&i0)<i1 ((i0^i1)&i1)<i0]

            psubd xmm1, xmm0        ;xmm1:[0 z1-z0 y1-y0 x1-x0]
            pand xmm1, xmm2         ;xmm1:[0 (z1-z0)&? (y1-y0)&? (x1-x0)&?]
            movmskps eax, xmm1
            pcmpgtd xmm1, dqzero    ;xmm1:[0 (z0-z1)&? (y0-y1)&? (x0-x1)&?]
            movmskps edx, xmm1

            sub eax, edx
        }
    #endif
}

#endif // MORTON_HPP
