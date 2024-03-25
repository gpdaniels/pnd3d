// This file has been modified from Ken Silverman's original release

#pragma once
#ifndef MACROS_HPP
#define MACROS_HPP

#if defined(_WIN32)
    #if defined(_WIN64)
        #define INT_PTR long long int
        using INT_PTR = __int64;
    #elif defined(_WIN32)
        using INT_PTR = int;
    #endif
#else
    #if defined(__x86_64__)
        using INT_PTR = long long int;
    #elif defined(__i386__)
        using INT_PTR = int;
    #endif

    using HANDLE = void*;
    using HMODULE = void*;
    using HINSTANCE = void*;
    using HTASK = void*;
    using HKEY = void*;
    using HDESK = void*;
    using HMF = void*;
    using HEMF = void*;
    using HPEN = void*;
    using HRSRC = void*;
    using HSTR = void*;
    using HWINSTA = void*;
    using HKL = void*;
    using HGDIOBJ = void*;

    using HDWP = HANDLE;
    using LPHANDLE = HANDLE*;

    using HWND = void*;
    using HMENU = void*;
    using HACCEL = void*;
    using HBRUSH = void*;
    using HFONT = void*;
    using HDC = void*;
    using HICON = void*;
    using HRGN = void*;
    using HMONITOR = void*;

    using HGLRC = void*;

    using UINT = unsigned int;
    using WPARAM = unsigned int;
    using LPARAM = long long int;

    using LRESULT = long long int;

    #define SW_NORMAL 1
    #define WS_OVERLAPPED 1
    #define WS_CAPTION 1
    #define WS_SYSMENU 1
    #define WS_MINIMIZEBOX 1
    #define WS_MAXIMIZEBOX 1
    #define WS_THICKFRAME 1

    #define WINAPI
    #define LPSTR char*
    struct RECT { int right, left, top, bottom; };
    struct POINT { int x, y; };

    #define CALLBACK
    #define MAX_PATH 256

    inline void DestroyWindow(HWND){}
    inline void SwapBuffers(HGLRC){}

    inline void ClipCursor(RECT*) {}
    inline void GetClipCursor(RECT*) {}
    inline void ClientToScreen(HWND, POINT*) {}

    inline void ShowCursor(int) {}
#endif

#if defined(_WIN32)
    #define NAKED_ATTRIBUTE __declspec(naked)
#else
    #define NAKED_ATTRIBUTE __attribute__((naked))
#endif

#ifndef __inline__
    #define __inline__ __inline
#endif

#ifndef __forceinline
    #define __forceinline __inline__
#endif

#if !defined(_WIN32)
    #define _vsnprintf vsnprintf
#endif

#endif // MACROS_HPP
