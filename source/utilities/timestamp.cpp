// This file has been modified from Ken Silverman's original release

#include "timestamp.hpp"

#include <chrono>

#if 0

#if defined(_WIN32)
    #define NOMINMAX
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
#endif

float time_now_seconds(void)
{
    static float gdrqper = [](){
        __int64 gqper;
        QueryPerformanceFrequency(reinterpret_cast<LARGE_INTEGER*>(&gqper));
        return 1.0f / static_cast<float>(gqper);
    }();

    static __int64 gq0 = [](){
        __int64 gq0_temp;
        QueryPerformanceCounter(reinterpret_cast<LARGE_INTEGER*>(&gq0_temp));
        return gq0_temp;
    }();

    __int64 q;
    QueryPerformanceCounter((LARGE_INTEGER*)&q);
    return static_cast<float>(q - gq0) * gdrqper;
}

#else

float time_now_seconds()
{
    static std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<float>(std::chrono::high_resolution_clock::now() - start).count();
}

#endif
