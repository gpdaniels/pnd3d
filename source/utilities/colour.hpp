// This file has been modified from Ken Silverman's original release

#pragma once
#ifndef COLOUR_HPP
#define COLOUR_HPP

#include "macros.hpp"

// Multiply colours.
//returns ARGB; 0x404040 in i or j is no change
__forceinline
static int mulcol (int i, int j)
{
#if 0
    _asm
    {
        pxor xmm7, xmm7         ; xmm7 is empty
        movd xmm0, i            ; move i into xmm0
        punpcklbw xmm0, xmm7    ; interleave low-order bytes from xmm0 (i) and xmm7 (zeros) into xmm0. = 0xAA00RR00GG00BB
        movd xmm1, j            ; move i into xmm1
        punpcklbw xmm1, xmm7    ; interleave low-order bytes from xmm1 (i) and xmm7 (zeros) into xmm1. = 0xAA00RR00GG00BB
        pmullw xmm0, xmm1       ; multiply the packed signed word integers in xmm0 register and xmm1, and store the low 16 bits of the results in xmm0.
        psrlw xmm0, 6           ; shift words in xmm0 right by 6 while shifting in 0s.
        packuswb xmm0, xmm0     ; converts 4 signed word integers from xmm0 and 4 signed word integers from xmm0 into 8 unsigned byte integers in mm using unsigned saturation.
        movd eax, xmm0          ; return lowest 32 bits of xmm0
    }
#else
    const int extracted_i[4] = {
        (i >> 24) & 0xFF,
        (i >> 16) & 0xFF,
        (i >> 8) & 0xFF,
        (i >> 0) & 0xFF
    };
    const int extracted_j[4] = {
        (j >> 24) & 0xFF,
        (j >> 16) & 0xFF,
        (j >> 8) & 0xFF,
        (j >> 0) & 0xFF
    };
    const int expanded_i[4] = {
        (extracted_i[0] << 8) | 0,
        (extracted_i[1] << 8) | 0,
        (extracted_i[2] << 8) | 0,
        (extracted_i[3] << 8) | 0
    };
    const int expanded_j[4] = {
        (extracted_j[0] << 8) | 0,
        (extracted_j[1] << 8) | 0,
        (extracted_j[2] << 8) | 0,
        (extracted_j[3] << 8) | 0
    };
    const int multiplied[4] = {
        expanded_i[0] * expanded_j[0],
        expanded_i[1] * expanded_j[1],
        expanded_i[2] * expanded_j[2],
        expanded_i[3] * expanded_j[3]
    };
    const int divided[4] = {
        multiplied[0] >> 16,
        multiplied[1] >> 16,
        multiplied[2] >> 16,
        multiplied[3] >> 16
    };
    const int shifted[4] = {
        divided[0] >> 6,
        divided[1] >> 6,
        divided[2] >> 6,
        divided[3] >> 6
    };
    const int saturated[4] = {
        shifted[0] > 0xFF ? 0xFF : shifted[0],
        shifted[1] > 0xFF ? 0xFF : shifted[1],
        shifted[2] > 0xFF ? 0xFF : shifted[2],
        shifted[3] > 0xFF ? 0xFF : shifted[3]
    };
    const int packed = ((saturated[0] & 0xFF) << 24) | ((saturated[1] & 0xFF) << 16) | ((saturated[2] & 0xFF) << 8) | ((saturated[3] & 0xFF) << 0);
    return packed;
#endif
}

// Scale colour.
//i:ARGB, j:scale (256 is no change), returns ARGB
__forceinline
static int mulsc (int i, int j)
{
#if 0
    _asm
    {
        movd xmm0, i
        punpcklbw xmm0, xmm0
        movd xmm1, j
        pshuflw xmm1, xmm1, 0
        pmulhuw xmm0, xmm1
        packuswb xmm0, xmm0
        movd eax, xmm0
    }
#else
    const int extracted_i[4] = {
        (i >> 24) & 0xFF,
        (i >> 16) & 0xFF,
        (i >> 8) & 0xFF,
        (i >> 0) & 0xFF
    };
    const int expanded_i[4] = {
        (extracted_i[0] << 8) | extracted_i[0],
        (extracted_i[1] << 8) | extracted_i[1],
        (extracted_i[2] << 8) | extracted_i[2],
        (extracted_i[3] << 8) | extracted_i[3]
    };
    const int scaled[4] = {
        expanded_i[0] * j,
        expanded_i[1] * j,
        expanded_i[2] * j,
        expanded_i[3] * j
    };
    const int divided[4] = {
        scaled[0] >> 16,
        scaled[1] >> 16,
        scaled[2] >> 16,
        scaled[3] >> 16
    };
    const int saturated[4] = {
        divided[0] > 0xFF ? 0xFF : divided[0],
        divided[1] > 0xFF ? 0xFF : divided[1],
        divided[2] > 0xFF ? 0xFF : divided[2],
        divided[3] > 0xFF ? 0xFF : divided[3]
    };
    const int packed = ((saturated[0] & 0xFF) << 24) | ((saturated[1] & 0xFF) << 16) | ((saturated[2] & 0xFF) << 8) | ((saturated[3] & 0xFF) << 0);
    return packed;
#endif
}

// Add colours. (Currently unused)
__forceinline
static int addcol (int i, int j)
{
#if 0
    _asm
    {
        movd xmm0, i
        movd xmm1, j
        paddusb xmm0, xmm1
        movd eax, xmm0
    }
#else
    const int expanded_i[4] = {
        (i >> 24) & 0xFF,
        (i >> 16) & 0xFF,
        (i >> 8) & 0xFF,
        (i >> 0) & 0xFF
    };
    const int expanded_j[4] = {
        (j >> 24) & 0xFF,
        (j >> 16) & 0xFF,
        (j >> 8) & 0xFF,
        (j >> 0) & 0xFF
    };
    const int added[4] = {
        expanded_i[0] + expanded_j[0],
        expanded_i[1] + expanded_j[1],
        expanded_i[2] + expanded_j[2],
        expanded_i[3] + expanded_j[3]
    };
    const int saturated[4] = {
        added[0] > 0xFF ? 0xFF : added[0],
        added[1] > 0xFF ? 0xFF : added[1],
        added[2] > 0xFF ? 0xFF : added[2],
        added[3] > 0xFF ? 0xFF : added[3]
    };
    const int packed = ((saturated[0] & 0xFF) << 24) | ((saturated[1] & 0xFF) << 16) | ((saturated[2] & 0xFF) << 8) | ((saturated[3] & 0xFF) << 0);
    return packed;
#endif
}

// Subtract colours. (Currently unused)
__forceinline
static int subcol (int i, int j)
{
#if 0
    _asm
    {
        movd xmm0, i
        movd xmm1, j
        psubusb xmm0, xmm1
        movd eax, xmm0
    }
#else
    const int expanded_i[4] = {
        (i >> 24) & 0xFF,
        (i >> 16) & 0xFF,
        (i >> 8) & 0xFF,
        (i >> 0) & 0xFF
    };
    const int expanded_j[4] = {
        (j >> 24) & 0xFF,
        (j >> 16) & 0xFF,
        (j >> 8) & 0xFF,
        (j >> 0) & 0xFF
    };
    const int subtracted[4] = {
        expanded_i[0] - expanded_j[0],
        expanded_i[1] - expanded_j[1],
        expanded_i[2] - expanded_j[2],
        expanded_i[3] - expanded_j[3]
    };
    const int saturated[4] = {
        subtracted[0] < 0 ? 0 : subtracted[0],
        subtracted[1] < 0 ? 0 : subtracted[1],
        subtracted[2] < 0 ? 0 : subtracted[2],
        subtracted[3] < 0 ? 0 : subtracted[3]
    };
    const int packed = ((saturated[0] & 0xFF) << 24) | ((saturated[1] & 0xFF) << 16) | ((saturated[2] & 0xFF) << 8) | ((saturated[3] & 0xFF) << 0);
    return packed;
#endif
}

#endif // COLOUR_HPP
