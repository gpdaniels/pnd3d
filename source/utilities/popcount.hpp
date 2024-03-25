// This file has been modified from Ken Silverman's original release

#pragma once
#ifndef POPCOUNT_HPP
#define POPCOUNT_HPP

#include "macros.hpp"

__forceinline
static int popcount8(unsigned char i) {
#if 1
    //http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetTable
    #define B2(n)    n ,    n+1 ,    n+1 ,    n+2
    #define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
    #define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
    constexpr static const int data[256] = {
        B6(0), B6(1), B6(1), B6(2)
    };
    #undef B6
    #undef B4
    #undef B2
    return data[i];
#else
    // Don't think this works
    constexpr static const unsigned char m1 = 0x55;
    constexpr static const unsigned char m2 = 0x33;
    constexpr static const unsigned char m4 = 0x0F;
    constexpr static const unsigned char h01 = 0x01;

    i -= ((i >> 1) & m1);
    i = (i & m2) + ((i >> 2) & m2);
    i = (i + (i >> 4)) & m4;

    return (i * h01) >> 6;
#endif
}

__forceinline
static int popcount16(unsigned short i) {
    // No idea if this works.
    constexpr static const unsigned short m1 = 0x5555;
    constexpr static const unsigned short m2 = 0x3333;
    constexpr static const unsigned short m4 = 0x0F0F;
    constexpr static const unsigned short h01 = 0x0101;

    i -= ((i >> 1) & m1);
    i = (i & m2) + ((i >> 2) & m2);
    i = (i + (i >> 4)) & m4;

    return (i * h01) >> 12;
}

__forceinline
static int popcount32(unsigned int i) {
    constexpr static const unsigned int m1 = 0x55555555;
    constexpr static const unsigned int m2 = 0x33333333;
    constexpr static const unsigned int m4 = 0x0F0F0F0F;
    constexpr static const unsigned int h01 = 0x01010101;

    i -= ((i >> 1) & m1);
    i = (i & m2) + ((i >> 2) & m2);
    i = (i + (i >> 4)) & m4;

    return (i * h01) >> 24;
}

__forceinline
static int popcount64(unsigned long long int i)
{
    constexpr static const unsigned long long int m1 = 0x5555555555555555ll;
    constexpr static const unsigned long long int m2 = 0x3333333333333333ll;
    constexpr static const unsigned long long int m4 = 0x0F0F0F0F0F0F0F0Fll;
    constexpr static const unsigned long long int h01 = 0x0101010101010101ll;

    i -= (i >> 1) & m1;
    i = (i & m2) + ((i >> 2) & m2);
    i = (i + (i >> 4)) & m4;

    return (i * h01) >> 56;
}

__forceinline
static int popcount128(void *buf)
{
#if 0
    alignas(16) static const int dqzeros[4] = {0,0,0,0};
    alignas(16) static const int dq5555[4] = {0x55555555, 0x55555555, 0x55555555, 0x55555555};
    alignas(16) static const int dq3333[4] = {0x33333333, 0x33333333, 0x33333333, 0x33333333};
    alignas(16) static const int dq0f0f[4] = {0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F};

    _asm
    {
        ;input: xmm0
        mov eax, buf
        movdqu xmm0, [eax]
        movaps xmm1, xmm0
        psrld xmm0, 1
        pand xmm0, dq5555
        psubd xmm1, xmm0
        movaps xmm0, xmm1
        psrld xmm1, 2
        pand xmm0, dq3333
        pand xmm1, dq3333
        paddd xmm0, xmm1
        movaps xmm1, xmm0
        psrld xmm0, 4
        paddd xmm0, xmm1
        pand xmm0, dq0f0f
        psadbw xmm0, dqzeros

        ;output: xmm0
        movhlps xmm1, xmm0
        paddd xmm0, xmm1
        movd eax, xmm0
    }
#else
#if defined(__x86_64__) || defined(_WIN64)
    return popcount64(static_cast<unsigned long long int*>(buf)[0]) +
           popcount64(static_cast<unsigned long long int*>(buf)[1]);
#elif defined(__i386__) || defined(_WIN32)
    return popcount32(static_cast<unsigned int*>(buf)[0]) +
           popcount32(static_cast<unsigned int*>(buf)[1]) +
           popcount32(static_cast<unsigned int*>(buf)[2]) +
           popcount32(static_cast<unsigned int*>(buf)[3]);
#else
    #error
#endif
#endif
}

#endif // POPCOUNT_HPP
