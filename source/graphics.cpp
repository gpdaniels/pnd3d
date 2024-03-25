// This file has been modified from Ken Silverman's original release

#include "graphics.hpp"

#include <algorithm>
#include <cstdarg>
#include <cstdio>
#include <cmath>

void print6x8(tiletype* dd, int ox, int y, int colour_foreground, int colour_background, const char* format, ...) {
    //NOTE: font is stored vertically first!
    //256 DOS chars, from: DOSAPP.FON (tab blank)
    constexpr static const unsigned int font6x8[384]
        = {0x00000000, 0x3E000000, 0x3E455145, 0x6F6B3E00, 0x1C003E6B, 0x1C3E7C3E, 0x7E3C1800, 0x3000183C, 0x30367F36, 0x7E5C1800, 0x0000185C, 0x00001818, 0xE7E7FFFF, 0x0000FFFF, 0x00000000,
           0xDBDBC3FF, 0x3000FFC3, 0x0E364A48, 0x79290600, 0x60000629, 0x04023F70, 0x0A7E6000, 0x2A003F35, 0x2A1C361C, 0x3E7F0000, 0x0800081C, 0x007F3E1C, 0x7F361400, 0x00001436, 0x005F005F,
           0x7F090600, 0x22007F01, 0x2259554D, 0x60606000, 0x14000060, 0x14B6FFB6, 0x7F060400, 0x10000406, 0x10307F30, 0x3E080800, 0x0800081C, 0x08083E1C, 0x40407800, 0x08004040, 0x083E083E,
           0x3F3C3000, 0x0300303C, 0x030F3F0F, 0x00000000, 0x00000000, 0x00065F06, 0x00030700, 0x24000307, 0x247E247E, 0x6A2B2400, 0x63000012, 0x63640813, 0x56493600, 0x00005020, 0x00000307,
           0x413E0000, 0x00000000, 0x00003E41, 0x1C3E0800, 0x0800083E, 0x08083E08, 0x60E00000, 0x08000000, 0x08080808, 0x60600000, 0x20000000, 0x02040810, 0x49513E00, 0x00003E45, 0x00407F42,
           0x49516200, 0x22004649, 0x36494949, 0x12141800, 0x2F00107F, 0x31494949, 0x494A3C00, 0x01003049, 0x03050971, 0x49493600, 0x06003649, 0x1E294949, 0x6C6C0000, 0x00000000, 0x00006CEC,
           0x22140800, 0x24000041, 0x24242424, 0x22410000, 0x02000814, 0x06095901, 0x5D413E00, 0x7E001E55, 0x7E111111, 0x49497F00, 0x3E003649, 0x22414141, 0x41417F00, 0x7F003E41, 0x41494949,
           0x09097F00, 0x3E000109, 0x7A494941, 0x08087F00, 0x00007F08, 0x00417F41, 0x40403000, 0x7F003F40, 0x41221408, 0x40407F00, 0x7F004040, 0x7F020402, 0x04027F00, 0x3E007F08, 0x3E414141,
           0x09097F00, 0x3E000609, 0x5E215141, 0x09097F00, 0x26006619, 0x32494949, 0x7F010100, 0x3F000101, 0x3F404040, 0x40201F00, 0x3F001F20, 0x3F403C40, 0x08146300, 0x07006314, 0x07087008,
           0x45497100, 0x00000043, 0x0041417F, 0x08040200, 0x00002010, 0x007F4141, 0x01020400, 0x80800402, 0x80808080, 0x07030000, 0x20000000, 0x78545454, 0x44447F00, 0x38003844, 0x28444444,
           0x44443800, 0x38007F44, 0x08545454, 0x097E0800, 0x18000009, 0x7CA4A4A4, 0x04047F00, 0x00000078, 0x00407D00, 0x84804000, 0x7F00007D, 0x00442810, 0x7F000000, 0x7C000040, 0x78041804,
           0x04047C00, 0x38000078, 0x38444444, 0x4444FC00, 0x38003844, 0xFC444444, 0x44784400, 0x08000804, 0x20545454, 0x443E0400, 0x3C000024, 0x007C2040, 0x40201C00, 0x3C001C20, 0x3C603060,
           0x10106C00, 0x9C00006C, 0x003C60A0, 0x54546400, 0x0800004C, 0x0041413E, 0x77000000, 0x00000000, 0x083E4141, 0x02010200, 0x3C000001, 0x3C262326, 0xE1A11E00, 0x3D001221, 0x007D2040,
           0x54543800, 0x20000955, 0x78555555, 0x54552000, 0x20007855, 0x78545555, 0x55572000, 0x1C007857, 0x1422E2A2, 0x55553800, 0x38000855, 0x08555455, 0x55553800, 0x00000854, 0x00417C01,
           0x79020000, 0x00000042, 0x00407C01, 0x24297000, 0x78007029, 0x782F252F, 0x54547C00, 0x34004555, 0x58547C54, 0x7F097E00, 0x38004949, 0x00394545, 0x44453800, 0x39000039, 0x00384445,
           0x21413C00, 0x3D00007D, 0x007C2041, 0x60A19C00, 0x3D00003D, 0x003D4242, 0x40413C00, 0x1800003D, 0x00246624, 0x493E4800, 0x29006249, 0x292A7C2A, 0x16097F00, 0x40001078, 0x02097E88,
           0x55542000, 0x00007855, 0x00417D00, 0x45443800, 0x3C000039, 0x007D2140, 0x0A097A00, 0x7A000071, 0x00792211, 0x55550800, 0x4E005E55, 0x004E5151, 0x4D483000, 0x3C002040, 0x04040404,
           0x04040404, 0x17001C04, 0x506A4C08, 0x34081700, 0x0000782A, 0x00307D30, 0x00140800, 0x14001408, 0x08140008, 0x11441144, 0x55AA1144, 0x55AA55AA, 0xEEBBEEBB, 0x0000EEBB, 0x0000FF00,
           0xFF080808, 0x0A0A0000, 0x0000FF0A, 0xFF00FF08, 0xF8080000, 0x0000F808, 0xFE0A0A0A, 0xFB0A0000, 0x0000FF00, 0xFF00FF00, 0xFA0A0000, 0x0000FE02, 0x0F080B0A, 0x0F080000, 0x00000F08,
           0x0F0A0A0A, 0x08080000, 0x0000F808, 0x0F000000, 0x08080808, 0x08080F08, 0xF8080808, 0x00000808, 0x0808FF00, 0x08080808, 0x08080808, 0x0808FF08, 0xFF000000, 0xFF000A0A, 0x0808FF00,
           0x0B080F00, 0xFE000A0A, 0x0A0AFA02, 0x0B080B0A, 0xFA0A0A0A, 0x0A0AFA02, 0xFB00FF00, 0x0A0A0A0A, 0x0A0A0A0A, 0xFB00FB0A, 0x0A0A0A0A, 0x0A0A0B0A, 0x0F080F08, 0x0A0A0808, 0x0A0AFA0A,
           0xF808F808, 0x0F000808, 0x08080F08, 0x0F000000, 0x00000A0A, 0x0A0AFE00, 0xF808F800, 0xFF080808, 0x0808FF00, 0xFB0A0A0A, 0x08080A0A, 0x00000F08, 0xF8000000, 0xFFFF0808, 0xFFFFFFFF,
           0xF0F0F0F0, 0xFFFFF0F0, 0x000000FF, 0xFF000000, 0x0F0FFFFF, 0x0F0F0F0F, 0x24241800, 0xFE002418, 0x00344A4A, 0x01017F00, 0x02000003, 0x027E027E, 0x49556300, 0x18000063, 0x041C2424,
           0x2020FC00, 0x0800001C, 0x00047804, 0x77550800, 0x3E000855, 0x003E4949, 0x02724C00, 0x22004C72, 0x00305955, 0x18241800, 0x18001824, 0x18247E24, 0x2A2A1C00, 0x3C00002A, 0x003C0202,
           0x2A2A2A00, 0x0000002A, 0x00242E24, 0x4A4A5100, 0x44000044, 0x00514A4A, 0xFC000000, 0x20000402, 0x00003F40, 0x2A080800, 0x24000808, 0x00122412, 0x09090600, 0x00000006, 0x00001818,
           0x00080000, 0x30000000, 0x02023E40, 0x010E0100, 0x0900000E, 0x00000A0D, 0x3C3C3C00, 0x0000003C, 0x00000000};

    if (!format) {
        return;
    }

    #if defined(__GNUC__)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wformat-nonliteral"
    #endif

    char string[1024];
    va_list arglist;
    va_start(arglist, format);
    if (_vsnprintf(&string[0], sizeof(string) - 1, format, arglist)) {
        string[sizeof(string) - 1] = 0;
    }
    va_end(arglist);

    #if defined(__GNUC__)
    #pragma GCC diagnostic pop
    #endif

    int i, ie;
    int* lpx;

    int* lp = reinterpret_cast<int*>(y * dd->p + dd->f);

    for (int j = 1; j < 256; y++, lp = reinterpret_cast<int*>((reinterpret_cast<INT_PTR>(lp)) + dd->p), j += j) {
        if (static_cast<unsigned int>(y) < static_cast<unsigned int>(dd->y)) {
            char* character = string;
            int x = ox;
            for (; *character; character++, x += 6) {
                const char* v = (static_cast<int>(*character)) * 6 + (reinterpret_cast<const char*>(font6x8));
                lpx = &lp[x];
                for (i = std::max(-x, 0), ie = std::min(dd->x - x, 6); i < ie; i++) {
                    if (v[i] & j) {
                        lpx[i] = colour_foreground;
                    }
                    else if (colour_background >= 0) {
                        lpx[i] = colour_background;
                    }
                }
                if ((*character) == 9) {
                    if (colour_background >= 0) {
                        for (i = std::max(-x, 6), ie = std::min(dd->x - x, 18); i < ie; i++) {
                            lpx[i] = colour_background;
                        }
                    }
                    x += 2 * 6;
                }
            }
        }
    }
}

void drawpix(tiletype* dd, int x, int y, int colour) {
    if ((static_cast<unsigned int>(x) < static_cast<unsigned int>(dd->x)) && (static_cast<unsigned int>(y) < static_cast<unsigned int>(dd->y))) {
        int* lp = reinterpret_cast<int*>(dd->f + dd->p * y + (x << 2));
        *lp = colour;
    }
}

void drawhlin(tiletype* dd, int x0, int x1, int y, int colour) {
    int x, *iptr;

    if (static_cast<unsigned int>(y) >= static_cast<unsigned int>(dd->y)) {
        return;
    }
    x0 = std::max(x0, 0);
    x1 = std::min(x1, dd->x);
    if (x0 >= x1) {
        return;
    }
    iptr = reinterpret_cast<int*>(dd->p * y + dd->f);
    for (x = x0; x < x1; x++) {
        iptr[x] = colour;
    }
}

void drawrectfill(tiletype* dd, int x0, int y0, int x1, int y1, int colour) {
    x0 = std::max(x0, 0);
    x1 = std::min(x1, dd->x);
    if (x0 >= x1)
        return;
    y0 = std::max(y0, 0);
    y1 = std::min(y1, dd->y);
    if (y0 >= y1)
        return;

    int* iptr = reinterpret_cast<int*>(dd->f + dd->p * y0);
    for (int y = y0; y < y1; y++, iptr = reinterpret_cast<int*>((reinterpret_cast<INT_PTR>(iptr)) + dd->p)) {
        for (int x = x0; x < x1; x++) {
            iptr[x] = colour;
        }
    }
}

void drawline(tiletype* dd, float x0, float y0, float x1, float y1, int colour) {
    if (x0 < 0.f) {
        if (x1 < 0.f)
            return;
        y0 = (0.f - x0) * (y1 - y0) / (x1 - x0) + y0;
        x0 = 0.f;
    }
    else if (x0 > dd->x) {
        if (x1 > dd->x)
            return;
        y0 = static_cast<float>(dd->x - x0) * (y1 - y0) / (x1 - x0) + y0;
        x0 = static_cast<float>(dd->x);
    }

    if (y0 < 0.f) {
        if (y1 < 0.f)
            return;
        x0 = (0.f - y0) * (x1 - x0) / (y1 - y0) + x0;
        y0 = 0.f;
    }
    else if (y0 > dd->y) {
        if (y1 > dd->y)
            return;
        x0 = static_cast<float>(dd->y - y0) * (x1 - x0) / (y1 - y0) + x0;
        y0 = static_cast<float>(dd->y);
    }

    if (x1 < 0.f) {
        y1 = (0.f - x1) * (y1 - y0) / (x1 - x0) + y1;
        x1 = 0.f;
    }
    else if (x1 > dd->x) {
        y1 = static_cast<float>(dd->x - x1) * (y1 - y0) / (x1 - x0) + y1;
        x1 = static_cast<float>(dd->x);
    }

    if (y1 < 0.f) {
        x1 = (0.f - y1) * (x1 - x0) / (y1 - y0) + x1;
        y1 = 0.f;
    }
    else if (y1 > dd->y) {
        x1 = static_cast<float>(dd->y - y1) * (x1 - x0) / (y1 - y0) + x1;
        y1 = static_cast<float>(dd->y);
    }

    x1 -= x0;
    y1 -= y0;
    int i = static_cast<int>(std::max(std::fabs(x1), std::fabs(y1)));
    /*FUK:32-bit!*/
    if (!(i & 0x7fffffff))
        return;

    float f = 65536.0f / static_cast<float>(i);

    int ipx0 = static_cast<int>(x0 * 65536.0f);
    int ipx1 = static_cast<int>(x1 * f);
    int ipy0 = static_cast<int>(y0 * 65536.0f);
    int ipy1 = static_cast<int>(y1 * f);

    for (; i > 0; i--) {
        int x = (ipx0 >> 16);
        int y = (ipy0 >> 16);
        if ((static_cast<unsigned int>(x) < static_cast<unsigned int>(dd->x)) && (static_cast<unsigned int>(y) < static_cast<unsigned int>(dd->y))) {
            int* lp = reinterpret_cast<int*>(dd->f + dd->p * y + (x << 2));
            *lp = colour;
        }
        ipx0 += ipx1;
        ipy0 += ipy1;
    }
}

void drawcirc(tiletype* dd, int xc, int yc, int r, int colour) {
    for (int x = 0, y = r, d = 1 - r; 1; x++) {
        drawpix(dd, xc + x, yc + y, colour);
        drawpix(dd, xc + y, yc + x, colour);
        drawpix(dd, xc + y, yc - x, colour);
        drawpix(dd, xc + x, yc - y, colour);
        drawpix(dd, xc - x, yc + y, colour);
        drawpix(dd, xc - y, yc + x, colour);
        drawpix(dd, xc - y, yc - x, colour);
        drawpix(dd, xc - x, yc - y, colour);
        if (x >= y) {
            break;
        }
        if (d < 0) {
            d += (x << 1) + 3;
        }
        else {
            d += ((x - y) << 1) + 5;
            y--;
        }
    }
}

void drawcircfill(tiletype* dd, int xc, int yc, int r, int colour) {
    for (int y = yc - r; y <= yc + r; y++) {
        int i = static_cast<int>(std::sqrt(static_cast<double>(r * r - (yc - y) * (yc - y))));
        drawhlin(dd, xc - i, xc + i, y, colour);
    }
}
