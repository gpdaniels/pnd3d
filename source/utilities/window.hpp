// This file has been modified from Ken Silverman's original release

#pragma once
#ifndef WINDOW_HPP
#define WINDOW_HPP

#include "macros.hpp"

#if defined(_WIN32)
    #define NOMINMAX
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
#endif

struct window_state {
    //Program flow
    HWND ghwnd = nullptr;
    int xres = 0;
    int yres = 0;
    char* prognam = const_cast<char*>("window_state");
    int shkeystatus = 0;
    int showwindmode = SW_NORMAL;
    int progwndflags = (WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU | WS_MINIMIZEBOX | WS_MAXIMIZEBOX | WS_THICKFRAME);

    //Keyboard
    #define KEYBUFSIZ 256
    int keybuf[KEYBUFSIZ] = {};
    int keybufr = 0;
    int keybufw = 0;
    int keybufw2 = 0;
    char keystatus[256] = {};

    //Mouse
    int mousx = 0;
    int mousy = 0;
    int bstatus = 0;
    int dmousx = 0;
    int dmousy = 0;
    int dmousz = 0;
    int ActiveApp = 1;

    //Program flow
    int initapp(HINSTANCE hinst);
    void uninitapp();
    int breath();
    void quitloop();

    //Screen
    int startdirectdraw (INT_PTR *f, int *p, int *x, int *y);
    void stopdirectdraw ();
    void nextpage();
    LRESULT WindowProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam);

    //Keyboard
    int keyread();
};

#endif // WINDOW_HPP
