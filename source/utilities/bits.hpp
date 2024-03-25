// This file has been modified from Ken Silverman's original release

#pragma once
#ifndef BITS_HPP
#define BITS_HPP

#include "macros.hpp"

#include <cstdint>
#include <cstring>

static void memset4(void* d, unsigned int v, int64_t n) {
    int64_t i;
    unsigned int* iptr = static_cast<unsigned int*>(d);
    n >>= 2;
    for (i = n - 1; i >= 0; i--) {
        iptr[i] = v;
    }
}

//Compatible with memset except:
//   1. All 32 bits of v are used and expected to be filled
//   2. Writes std::max((n+15)&~15,16) bytes
//   3. Assumes d is aligned on 16 byte boundary
#if 0
NAKED_ATTRIBUTE
static void memset16_safe(void* d, int v, int n) {
    (void)d;
    (void)v;
    (void)n;
    _asm {
        mov edx, [esp+4] ;d
        mov ecx, [esp+12] ;n
        movd xmm0, [esp+8] ;v
        pshufd xmm0, xmm0, 0
        add edx, ecx
        neg ecx
        test edx, 15
        jz short memset16aligned
memset16unaligned:
        movdqu [edx+ecx], xmm0
        add ecx, 16
        jl short memset16unaligned
        ret
memset16aligned:
        movaps [edx+ecx], xmm0
        add ecx, 16
        jl short memset16aligned
        ret
    }
}
#else
// Marginally slower (20.3FPS vs 20.8FPS)
static void memset16_safe(void* d, int v, int n) {
    unsigned char* buf = reinterpret_cast<unsigned char*>(d);
    int blob[4] = {v, v, v, v};
    buf += n;
    n = -n;
    while(n != 0) {
        std::memcpy(buf + n, reinterpret_cast<unsigned char*>(&blob[0]), 4 * 4);
        n += 16;
    }
}
#endif

//Compatible with memset except:
//   1. All 32 bits of v are used and expected to be filled
//   2. Writes std::max((n+15)&~15,16) bytes
//   3. Assumes d is aligned on 16 byte boundary
#if 0
NAKED_ATTRIBUTE
static void memset16(void* d, int v, int n) {
    (void)d;
    (void)v;
    (void)n;
    _asm {
        mov edx, [esp+4] ;d
        mov ecx, [esp+12] ;n
        movd xmm0, [esp+8] ;v
        pshufd xmm0, xmm0, 0
        add edx, ecx
        neg ecx
memset16aligned:
        movaps [edx+ecx], xmm0
        add ecx, 16
        jl short memset16aligned
        ret
    }
}
#else
// Marginally slower (20.3FPS vs 20.8FPS)
static void memset16(void* d, int v, int n) {
    unsigned char* buf = reinterpret_cast<unsigned char*>(d);
    int blob[4] = {v, v, v, v};
    buf += n;
    n = -n;
    while(n != 0) {
        std::memcpy(buf + n, reinterpret_cast<unsigned char*>(&blob[0]), 4 * 4);
        n += 16;
    }
}
#endif

static const int pow2[32] =
{
    0x00000001,0x00000002,0x00000004,0x00000008,0x00000010,0x00000020,0x00000040,0x00000080,
    0x00000100,0x00000200,0x00000400,0x00000800,0x00001000,0x00002000,0x00004000,0x00008000,
    0x00010000,0x00020000,0x00040000,0x00080000,0x00100000,0x00200000,0x00400000,0x00800000,
    0x01000000,0x02000000,0x04000000,0x08000000,0x10000000,0x20000000,0x40000000,static_cast<int>(0x80000000),
};

static const int pow2m1[32] =
{
    0x00000000,0x00000001,0x00000003,0x00000007,0x0000000f,0x0000001f,0x0000003f,0x0000007f,
    0x000000ff,0x000001ff,0x000003ff,0x000007ff,0x00000fff,0x00001fff,0x00003fff,0x00007fff,
    0x0000ffff,0x0001ffff,0x0003ffff,0x0007ffff,0x000fffff,0x001fffff,0x003fffff,0x007fffff,
    0x00ffffff,0x01ffffff,0x03ffffff,0x07ffffff,0x0fffffff,0x1fffffff,0x3fffffff,0x7fffffff,
};

static const int npow2[32] =
{
    static_cast<int>(0xffffffff),static_cast<int>(0xfffffffe),static_cast<int>(0xfffffffc),static_cast<int>(0xfffffff8),static_cast<int>(0xfffffff0),static_cast<int>(0xffffffe0),static_cast<int>(0xffffffc0),static_cast<int>(0xffffff80),
    static_cast<int>(0xffffff00),static_cast<int>(0xfffffe00),static_cast<int>(0xfffffc00),static_cast<int>(0xfffff800),static_cast<int>(0xfffff000),static_cast<int>(0xffffe000),static_cast<int>(0xffffc000),static_cast<int>(0xffff8000),
    static_cast<int>(0xffff0000),static_cast<int>(0xfffe0000),static_cast<int>(0xfffc0000),static_cast<int>(0xfff80000),static_cast<int>(0xfff00000),static_cast<int>(0xffe00000),static_cast<int>(0xffc00000),static_cast<int>(0xff800000),
    static_cast<int>(0xff000000),static_cast<int>(0xfe000000),static_cast<int>(0xfc000000),static_cast<int>(0xf8000000),static_cast<int>(0xf0000000),static_cast<int>(0xe0000000),static_cast<int>(0xc0000000),static_cast<int>(0x80000000),
};


#ifdef _MSC_VER

#if !defined(_WIN64)

static __forceinline unsigned int bsf(unsigned int a) {
    _asm bsf eax, a
}
static __forceinline unsigned int bsr(unsigned int a) {
    _asm bsr eax, a
}

#else

static __forceinline unsigned int bsf(unsigned int a) {
    unsigned int r;
    _BitScanForward(&r, a);
    return (r);
}
static __forceinline unsigned int bsr(unsigned int a) {
    unsigned int r;
    _BitScanReverse(&r, a);
    return (r);
}
static __forceinline unsigned int bsf64(unsigned long long int a) {
    unsigned int r;
    _BitScanForward64(&r, a);
    return (r);
}
static __forceinline unsigned int bsr64(unsigned long long int a) {
    unsigned int r;
    _BitScanReverse64(&r, a);
    return (r);
}

#endif

#else

//#define bsf(x) __builtin_clz(x) //Count Leading  Zeros; compiler generates function call!
//#define bsr(x) __builtin_ctz(x) //Count Trailing Zeros; compiler generates function call! GCC bug: same as "_clz!
//#define bsf64(x) __builtin_clzl(x) //Count Leading  Zeros 64-bit (untested)
//#define bsr64(x) __builtin_ctzl(x) //Count Trailing Zeros 64-bit (untested)
#define bsf(r)   ({ unsigned int           __r=(r); __asm__ __volatile__ ("bsf %0, %0;": "=a" (__r): "a" (__r): "cc"); __r; })
#define bsr(r)   ({ unsigned int           __r=(r); __asm__ __volatile__ ("bsr %0, %0;": "=a" (__r): "a" (__r): "cc"); __r; })
#define bsf64(r) ({ unsigned long long int __r=(r); __asm__ __volatile__ ("bsf %0, %0;": "=a" (__r): "a" (__r): "cc"); __r; })
#define bsr64(r) ({ unsigned long long int __r=(r); __asm__ __volatile__ ("bsr %0, %0;": "=a" (__r): "a" (__r): "cc"); __r; })

#endif

__forceinline static int dntil0(unsigned int* lptr, int z, int zsiz) {
#if 0
    //This line does the same thing (but slow & brute force)
    while ((z < zsiz) && (lptr[z >> 5] & (1 << (z & 31)))) {
        z++;
    }
    return z;
#else
    // WARNING: zsiz must be multiple of 32!
    int i = (lptr[z >> 5] | ((1 << (z & 31)) - 1));
    z &= ~31;
    while (i == 0xffffffff) {
        z += 32;
        if (z >= zsiz) {
            return zsiz;
        }
        i = lptr[z >> 5];
    }
    return bsf(~i) + z;
#endif
}

__forceinline static int uptil0(unsigned int* lptr, int z) {
#if 0
    //This line does the same thing (but slow & brute force)
    while ((z > 0) && (lptr[(z - 1) >> 5] & (1 << ((z - 1) & 31)))) {
        z--;
    }
    return z;
#else
    //Prevent possible crash
    if (!z) {
        return 0;
    }

    int i = (lptr[(z - 1) >> 5] | (-1 << (z & 31)));
    z &= ~31;
    while (i == 0xffffffff) {
        z -= 32;
        if (z < 0) {
            return 0;
        }
        i = lptr[z >> 5];
    }
    return bsr(~i) + z + 1;
#endif
}

__forceinline static int dntil1(unsigned int* lptr, int z, int zsiz) {
#if 0
    //This line does the same thing (but slow & brute force)
    while ((z < zsiz) && (!(lptr[z >> 5] & (1 << (z & 31))))) {
        z++;
    }
    return z;
#else
    //WARNING: zsiz must be multiple of 32!
    int i = (lptr[z >> 5] & (-1 << (z & 31)));
    z &= ~31;
    while (!i) {
        z += 32;
        if (z >= zsiz) {
            return zsiz;
        }
        i = lptr[z >> 5];
    }
    return bsf(i) + z;
#endif
}

__forceinline static int uptil1(unsigned int* lptr, int z) {
#if 0
    //This line does the same thing (but slow & brute force)
    while ((z > 0) && (!(lptr[(z - 1) >> 5] & (1 << ((z - 1) & 31))))) {
        z--;
    }
    return z;
#else
    //Prevent possible crash
    if (!z) {
        return 0;
    }

    // ADDED:
    // Seems to fix an invalid read.
    // Sometimes z is -1
    // Therefore ((-1 - 1) >> 5) => -2 >> 5 => 0xFFFFFFFF (as sign bit is extended).
    // lptr[0xFFFFFFFF] is definitely invalid.
    if (z < 0) {
        return z;
    }

    int i = (lptr[(z - 1) >> 5] & ((1 << (z & 31)) - 1));
    z &= ~31;
    while (!i) {
        z -= 32;
        if (z < 0) {
            return 0;
        }
        i = lptr[z >> 5];
    }
    return bsr(i) + z + 1;
#endif
}

// Swap bits in range {i0<=?<i0+n}; n must be <= 25
__forceinline static void xorzrangesmall(void* iptr, int i0, int n) {
    *(int*)((char*)iptr + (i0 >> 3)) ^= (pow2m1[n] << (i0 & 7));
}

#if defined(__x86_64__) || defined(_WIN64)

// Set all bits in vbit from (x,y,z0) to (x,y,z1-1) to 0's
__forceinline static void setzrange0(void* vptr, int64_t z0, int64_t z1) {
    uint64_t z, ze, *iptr = (uint64_t*)vptr;
    if (!((z0 ^ z1) & ~63)) {
        iptr[z0 >> 6] &= ((~(-1ll << z0)) | (-1ll << z1));
        return;
    }
    z = (z0 >> 6);
    ze = (z1 >> 6);
    iptr[z] &= ~(-1ll << z0);
    for (z++; z < ze; z++) {
        iptr[z] = 0ll;
    }
    iptr[z] &= (-1ll << z1);
}

//Set all bits in vbit from (x,y,z0) to (x,y,z1-1) to 1's
__forceinline static void setzrange1(void* vptr, int64_t z0, int64_t z1) {
    uint64_t z, ze, *iptr = (uint64_t*)vptr;
    if (!((z0 ^ z1) & ~63)) {
        iptr[z0 >> 6] |= ((~(-1ll << z1)) & (-1ll << z0));
        return;
    }
    z = (z0 >> 6);
    ze = (z1 >> 6);
    iptr[z] |= (-1ll << z0);
    for (z++; z < ze; z++) {
        iptr[z] = -1ll;
    }
    iptr[z] |= ~(-1ll << z1);
}

#elif defined(__i386__) || defined(_WIN32)

// Set all bits in vbit from (x,y,z0) to (x,y,z1-1) to 0's
__forceinline static void setzrange0(void* vptr, int z0, int z1) {
    int z, ze, *iptr = (int*)vptr;
    if (!((z0 ^ z1) & ~31)) {
        iptr[z0 >> 5] &= ((~(-1 << z0)) | (-1 << z1));
        return;
    }
    z = (z0 >> 5);
    ze = (z1 >> 5);
    iptr[z] &= ~(-1 << z0);
    for (z++; z < ze; z++) {
        iptr[z] = 0;
    }
    iptr[z] &= (-1 << z1);
}

//Set all bits in vbit from (x,y,z0) to (x,y,z1-1) to 1's
__forceinline static void setzrange1(void* vptr, int z0, int z1) {
    int z, ze, *iptr = (int*)vptr;
    if (!((z0 ^ z1) & ~31)) {
        iptr[z0 >> 5] |= ((~(-1 << z1)) & (-1 << z0));
        return;
    }
    z = (z0 >> 5);
    ze = (z1 >> 5);
    iptr[z] |= (-1 << z0);
    for (z++; z < ze; z++) {
        iptr[z] = -1;
    }
    iptr[z] |= ~(-1 << z1);
}
#endif

//Returns 1 if any bits in vbit from (x,y,z0) to (x,y,z1-1) are 0
__forceinline
int anyzrange0(unsigned int* lptr, int z0, int z1) {
    int z, ze;
    if (!((z0 ^ z1) & ~31)) {
        return ((lptr[z0 >> 5] | ((~(-1 << (z0 & 31))) | (-1 << (z1 & 31)))) != 0xffffffff);
    }
    z = (z0 >> 5);
    ze = (z1 >> 5);
    if ((lptr[z] | ~(-1 << (z0 & 31))) != 0xffffffff) {
        return 1;
    }
    for (z++; z < ze; z++) {
        if (lptr[z] != 0xffffffff) {
            return 1;
        }
    }
    if ((lptr[z] | (-1 << (z1 & 31))) != 0xffffffff) {
        return 1;
    }
    return 0;
}

#if !defined(_WIN32)

__forceinline
unsigned long _lrotl(unsigned long value, int shift) {
    return (value << shift) | (value >> (sizeof (unsigned long) * 8 - shift));
}

#endif

#endif // BITS_HPP
