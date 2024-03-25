// This file has been modified from Ken Silverman's original release

#include "window.hpp"

#if defined(_WIN32)

#include "timestamp.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <cmath>

//--------------------------------------------------------------------------------------------------

typedef struct
{
    HBITMAP hBitmap;
    BITMAPINFO* pbmi;
    INT_PTR ptr, bpl, xsiz, ysiz;
} bitmap_t;

static bitmap_t gscr;

static HMODULE huser32lib;

typedef struct
{
    USHORT usUsagePage, usUsage;
    DWORD dwFlags;
    HWND hwndTarget;
} KRAWINPUTDEVICE;

typedef struct
{
    DWORD dwType, dwSize;
    HANDLE hDevice;
    WPARAM wParam;
} KRAWINPUTHEADER;

static BOOL(__stdcall* KRegisterRawInputDevices)(KRAWINPUTDEVICE* pRawInputDevices, UINT uiNumDevices, UINT cbSize);

static UINT(__stdcall* KGetRawInputData)(HANDLE /*HRAWINPUT*/ hRawInput, UINT uiCommand, LPVOID pData, PUINT pcbSize, UINT cbSizeHeader);

static void wminput_uninit(void) {
    if (huser32lib) {
        FreeLibrary(huser32lib);
        huser32lib = 0;
    }
}

static int wminput_init(HWND hwnd) {
    KRAWINPUTDEVICE rid;
    char tbuf[MAX_PATH];

    GetSystemDirectory(tbuf, sizeof(tbuf) - 12);
    strcat(tbuf, "\\user32.dll");
    huser32lib = LoadLibrary(tbuf);
    if (!huser32lib) {
        MessageBox(0, "LoadLibrary failed", ":/", MB_OK);
        return (-1);
    }
    KRegisterRawInputDevices = (BOOL(_stdcall*)(KRAWINPUTDEVICE*, UINT, UINT))GetProcAddress(huser32lib, "RegisterRawInputDevices");
    if (!KRegisterRawInputDevices)
        return (-1);
    KGetRawInputData = (UINT(_stdcall*)(HANDLE /*HRAWINPUT*/, UINT, LPVOID, PUINT, UINT))GetProcAddress(huser32lib, "GetRawInputData");
    if (!KGetRawInputData)
        return (-1);

    rid.usUsagePage = 1;      //HID_USAGE_PAGE_GENERIC;
    rid.usUsage = 2;          //HID_USAGE_GENERIC_MOUSE;
    rid.dwFlags = 0x00000100; //RIDEV_INPUTSINK;
    rid.hwndTarget = hwnd;
    KRegisterRawInputDevices(&rid, 1, sizeof(KRAWINPUTDEVICE));
    return (0);
}

int window_state::startdirectdraw(INT_PTR* f, int* p, int* x, int* y) {
    HDC hdc;
    static int gcxres, gcyres = -1;

    if ((xres != gcxres) || (yres != gcyres) || (gscr.hBitmap == INVALID_HANDLE_VALUE)) {
        hdc = GetDC(ghwnd);
        if (gscr.hBitmap != INVALID_HANDLE_VALUE) {
            DeleteObject(gscr.hBitmap);
        }
        gcxres = xres;
        gcyres = yres;
        gscr.xsiz = xres;
        gscr.ysiz = yres;
        gscr.pbmi = (BITMAPINFO*)LocalAlloc(LMEM_FIXED, sizeof(BITMAPINFO));
        gscr.pbmi->bmiHeader.biSize = sizeof(BITMAPINFO);
        gscr.pbmi->bmiHeader.biWidth = xres;
        gscr.pbmi->bmiHeader.biHeight = -yres;
        gscr.pbmi->bmiHeader.biPlanes = 1;
        gscr.pbmi->bmiHeader.biBitCount = 32;
        gscr.pbmi->bmiHeader.biCompression = BI_RGB;
        gscr.pbmi->bmiHeader.biSizeImage = 0;
        gscr.pbmi->bmiHeader.biXPelsPerMeter = xres;
        gscr.pbmi->bmiHeader.biYPelsPerMeter = yres;
        gscr.pbmi->bmiHeader.biClrUsed = 0;
        gscr.pbmi->bmiHeader.biClrImportant = 0;
        gscr.hBitmap = CreateDIBSection(hdc, gscr.pbmi, DIB_RGB_COLORS, (void**)&gscr.ptr, 0, 0);
        gscr.bpl = xres * 4;
        ReleaseDC(0, hdc);
    }
    (*f) = gscr.ptr;
    (*p) = gscr.bpl;
    (*x) = gscr.xsiz;
    (*y) = gscr.ysiz;
    return (1);
}

void window_state::stopdirectdraw() {
}

void window_state::nextpage() {
    HDC hdc;
    hdc = GetDC(ghwnd);
    StretchDIBits(hdc, 0, 0, gscr.xsiz, gscr.ysiz, 0, 0, gscr.xsiz, gscr.ysiz, (void*)gscr.ptr, gscr.pbmi, DIB_RGB_COLORS, SRCCOPY);
    ReleaseDC(0, hdc);
}

static LRESULT CALLBACK WindowProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    window_state* me = reinterpret_cast<window_state*>(GetWindowLongPtr(hwnd, GWLP_USERDATA));
    if (me) return me->WindowProc(hwnd, msg, wParam, lParam);
    return DefWindowProc(hwnd, msg, wParam, lParam);
}

LRESULT window_state::WindowProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    switch (msg) {
    case WM_SIZE:
        if (hwnd == ghwnd) {
            xres = (lParam & 65535);
            yres = (lParam >> 16);
        }
        break;
    case WM_MOUSEMOVE:
        mousx = (lParam & 65535);
        mousy = (lParam >> 16);
        break;
    case 0x020a:
        dmousz = std::min(std::max(dmousz + ((signed short)HIWORD(wParam)), (signed)0x80010000), (signed)0x7fff0000);
        break;   //WM_MOUSEWHEEL
    case 0x00ff: //WM_INPUT
    {
        unsigned char tbuf[64];
        unsigned int u = sizeof(tbuf);
        KGetRawInputData((HANDLE /*HRAWINPUT*/)lParam, 0x10000003 /*RID_INPUT*/, tbuf, &u, sizeof(KRAWINPUTHEADER));
        if (*(int*)&tbuf /*raw->header.dwType*/ != 0 /*RIM_TYPEMOUSE*/)
            break;

        //printf("Mouse data is: size:%d device:%I64d wparam:%d\n",*(int *)&tbuf[4],*(__int64 *)&tbuf[8],*(int *)&tbuf[12]);
        if (!*(__int64*)&tbuf[8])
            break; //Ignore touch screen messages

        u = sizeof(INT_PTR) * 2 + 12;
        dmousx += *(int*)&tbuf[u + 8];  //raw->data.mouse.lLastX;
        dmousy += *(int*)&tbuf[u + 12]; //raw->data.mouse.lLastY;
        if (tbuf[u] & 1)
            bstatus |= 1;
        if (tbuf[u] & 2)
            bstatus &= ~1;
        if (tbuf[u] & 4)
            bstatus |= 2;
        if (tbuf[u] & 8)
            bstatus &= ~2;
        if (tbuf[u] & 16)
            bstatus |= 4;
        if (tbuf[u] & 32)
            bstatus &= ~4;
    } break;

    case WM_KEYDOWN:
    case WM_SYSKEYDOWN:
        keystatus[((lParam >> 16) & 127) + ((lParam >> 17) & 128)] |= 1;
        if ((wParam & 0xff) == 0xff)
            break; //Fixes SHIFT+[ext key] on XP
        switch (lParam & 0x17f0000) {
        case 0x02a0000:
            shkeystatus |= (1 << 16);
            break; //0x2a
        case 0x0360000:
            shkeystatus |= (1 << 17);
            break; //0x36
        case 0x01d0000:
            shkeystatus |= (1 << 18);
            break; //0x1d
        case 0x11d0000:
            shkeystatus |= (1 << 19);
            break; //0x9d
        case 0x0380000:
            shkeystatus |= (1 << 20);
            break; //0x38
        case 0x1380000:
            shkeystatus |= (1 << 21);
            break; //0xb8
        default: {
            int i = ((keybufw2 + 1) & (KEYBUFSIZ - 1));
            keybuf[keybufw2 & (KEYBUFSIZ - 1)] = (((lParam >> 8) & 0x7f00) + ((lParam >> 9) & 0x8000)) + shkeystatus;
            if (i == keybufr)
                keybufr = ((keybufr + 1) & (KEYBUFSIZ - 1)); //lose old keystrokes to prevent fifo overlap
            keybufw2 = i;
        }
        }
        return (0);
    case WM_KEYUP:
    case WM_SYSKEYUP:
        keystatus[((lParam >> 16) & 127) + ((lParam >> 17) & 128)] &= ~1;
        if ((wParam & 0xff) == 0xff)
            break; //Fixes SHIFT+[ext key] on XP
        switch (lParam & 0x17f0000) {
        case 0x02a0000:
            shkeystatus &= ~(3 << 16);
            break; //0x2a
        case 0x0360000:
            shkeystatus &= ~(3 << 16);
            break; //0x36
        case 0x01d0000:
            shkeystatus &= ~(1 << 18);
            break; //0x1d
        case 0x11d0000:
            shkeystatus &= ~(1 << 19);
            break; //0x9d
        case 0x0380000:
            shkeystatus &= ~(1 << 20);
            break; //0x38
        case 0x1380000:
            shkeystatus &= ~(1 << 21);
            break; //0xb8
        }
        return (0);
    case WM_SYSCHAR:
    case WM_CHAR:
        if (keybufw2 != keybufr) //stick ASCII code in last FIFO value
            keybuf[(keybufw2 - 1) & (KEYBUFSIZ - 1)] |= (wParam & 255);
        return (0);

    case WM_ACTIVATEAPP:
        ActiveApp = (BOOL)wParam;
        shkeystatus = 0;
        break;
    case WM_DESTROY:
        PostQuitMessage(0);
        break;
    }
    return (DefWindowProc(hwnd, msg, wParam, lParam));
}

static WNDCLASS wc;
int window_state::initapp(HINSTANCE hinst) {
    RECT rw;
    wc.style = CS_HREDRAW | CS_VREDRAW;
    wc.lpfnWndProc = ::WindowProc;
    wc.cbClsExtra = wc.cbWndExtra = 0;
    wc.hInstance = hinst;
    wc.hIcon = LoadIcon(0, IDI_APPLICATION);
    wc.hbrBackground = (HBRUSH__*)GetStockObject(BLACK_BRUSH);
    wc.lpszMenuName = wc.lpszClassName = prognam;
    if (!RegisterClass(&wc))
        return (0);

    SystemParametersInfo(SPI_GETWORKAREA, 0, &rw, 0);
    ghwnd = CreateWindowEx(0,
                           prognam,
                           prognam,
                           progwndflags,
                           ((rw.right - rw.left - (xres + GetSystemMetrics(SM_CXSIZEFRAME) * 2)) >> 1) + rw.left,
                           ((rw.bottom - rw.top - (yres + GetSystemMetrics(SM_CYCAPTION) + GetSystemMetrics(SM_CYSIZEFRAME) * 2)) >> 1) + rw.top,
                           xres + GetSystemMetrics(SM_CXSIZEFRAME) * 2,
                           yres + GetSystemMetrics(SM_CYCAPTION) + GetSystemMetrics(SM_CYSIZEFRAME) * 2,
                           0,
                           0,
                           hinst,
                           reinterpret_cast<LPVOID>(this));
    if (!ghwnd)
        return (0);

    SetWindowLongPtr(ghwnd, GWLP_USERDATA, (long)this);

    ShowWindow(ghwnd, showwindmode); //SW_NORMAL/SW_MAXIMIZE/etc..
    UpdateWindow(ghwnd);

    memset(keystatus, 0, sizeof(keystatus));

    if (wminput_init(ghwnd) < 0) { /*bad*/
        return (0);
    }

    return (1);
}
void window_state::uninitapp(void) {
    DeleteObject(gscr.hBitmap);
    gscr.hBitmap = (HBITMAP)INVALID_HANDLE_VALUE;
    wminput_uninit();
    DestroyWindow(ghwnd);
    UnregisterClass(wc.lpszClassName, wc.hInstance);
}

int window_state::breath(void) {
    MSG msg;
    while (PeekMessage(&msg, 0, 0, 0, PM_REMOVE)) {
        if (msg.message == WM_QUIT)
            return (-1);
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    keybufw = keybufw2; //to be safe with multithreads
    return (0);
}

void window_state::quitloop(void) {
    PostMessage(ghwnd, WM_CLOSE, 0, 0);
}

int window_state::keyread(void) {
    int i;
    if (keybufr == keybufw)
        return (0);
    i = keybuf[keybufr];
    keybufr = ((keybufr + 1) & (KEYBUFSIZ - 1));
    return (i);
}

#else

#include "gtl/io/window"

static gtl::window main_window;

int window_state::initapp(HINSTANCE hinst) {
    main_window.open(0, 0, xres, yres, prognam);
    main_window.set_visible(true);
    main_window.prepare();
    return 1;
}

void window_state::uninitapp() {
    main_window.close();
}

int window_state::breath() {
    // Have to present the previous frame here, as there's no other easy way of modifying the current code.
    main_window.present();

    gtl::window::event_type event;
    while (main_window.process(event)) {
        //...
    }
    keybufw = keybufw2; //to be safe with multithreads

    // Have to prepare the next frame here, as there's no other easy way of modifying the current code.
    main_window.prepare();

    return !main_window.is_open();
}

void window_state::quitloop() {
    main_window.close();
}

static int* framebuffer;

int window_state::startdirectdraw(INT_PTR *f, int *p, int *x, int *y) {

    static bool initialised = false;
    if (!initialised) {
        initialised = true;
        framebuffer = new int[xres * yres];
    }

    (*f) = (INT_PTR)framebuffer;
    (*p) = xres * 4;
    (*x) = xres;
    (*y) = yres;
    return 1;
}

void window_state::stopdirectdraw() {
    static GLuint texture = 0;
    if (texture == 0) {
        glGenTextures(1, &texture);
    }
    if (texture == -1) {
        exit(1);
    }

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, xres, yres, 0, GL_RGBA, GL_UNSIGNED_BYTE, framebuffer);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, texture);
    glDisable(GL_LIGHTING);

    // Set up ortographic projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, xres, 0, yres, -1, 1);

    // Render a quad
    glBegin(GL_QUADS);
    glTexCoord2f(0,1); glVertex2f(0,0);
    glTexCoord2f(0,0); glVertex2f(0,yres);
    glTexCoord2f(1,0); glVertex2f(xres, yres);
    glTexCoord2f(1,1); glVertex2f(xres,0);
    glEnd();

    // Reset Projection Matrix
    glDisable(GL_TEXTURE_2D);
    glEnable(GL_LIGHTING);
}

void window_state::nextpage() {
}

LRESULT window_state::WindowProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    return 0;
}

int window_state::keyread() {
    int i;
    if (keybufr == keybufw)
        return (0);
    i = keybuf[keybufr];
    keybufr = ((keybufr + 1) & (KEYBUFSIZ - 1));
    return (i);
}

#endif
