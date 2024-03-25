// This file has been modified from Ken Silverman's original release

#include "thread_pool.hpp"

#define THREAD_POOL_ENABLED 1

#if defined(_WIN32) && defined(THREAD_POOL_ENABLED) && THREAD_POOL_ENABLED

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <process.h>
#include <algorithm>

#include "utilities/macros.hpp"

#define MAXCPU 64

static HANDLE gthand[MAXCPU-1];
static HANDLE ghevent[2][MAXCPU-1]; //WARNING: WaitForMultipleObjects has a limit of 64 threads
static int gnthreadcnt = 0, glincnt[2];

#define getnextindexsynced() InterlockedExchangeAdd((long *)&glincnt[0],1)

static void (*ghtcallfunc)(int, void *);
static void *ghtdata;

static void htuninit()
{
    int i, j;
    for(i=gnthreadcnt-1;i>=0;i--)
    {
        TerminateThread(gthand[i],0); //Dangerous if using system resources in thread, but we are not
        for(j=1;j>=0;j--)
        { if (ghevent[j][i] != (HANDLE)-1) { CloseHandle(ghevent[j][i]); ghevent[j][i] = (HANDLE)-1; } }
    }
}

static unsigned int _stdcall ghtfunc (void *_)
{
    void (*lcallfunc)(int, void *);
    void *luserdata;
    int i, thnum = (int)_, endcnt;

    while (1)
    {
        WaitForSingleObject(ghevent[0][thnum],INFINITE);
        lcallfunc = ghtcallfunc; luserdata = ghtdata; endcnt = glincnt[1];
        while ((i = getnextindexsynced()) < endcnt) lcallfunc(i,luserdata);
        SetEvent(ghevent[1][thnum]);
    }
}

void htrun (void (*lcallfunc)(int, void *), void *luserdata, int v0, int v1, bool should_parallelise)
{
    static int oct_numcpu = [](){
        SYSTEM_INFO si;
        GetSystemInfo(&si);
        return std::min(si.dwNumberOfProcessors, static_cast<DWORD>(MAXCPU));
    }();

    int lnumcpu = oct_numcpu * should_parallelise;

    int i, threaduse; unsigned win98requiresme;

    //1 CPU requested; execute here and quit.
    if (lnumcpu <= 1) { for(i=v0;i<v1;i++) lcallfunc(i,luserdata); return; }

    //Initialize new threads if necessary
    threaduse = std::min(lnumcpu,MAXCPU)-1;
    while (gnthreadcnt < threaduse)
    {
        if (!gnthreadcnt) { SetThreadAffinityMask(GetCurrentThread(),1<<(lnumcpu-1)); atexit(htuninit); }
        for(i=0;i<2;i++) ghevent[i][gnthreadcnt] = CreateEvent(0,0,0,0);
        gthand[gnthreadcnt] = (HANDLE)_beginthreadex(0,0,ghtfunc,(void *)gnthreadcnt,0,&win98requiresme);
        SetThreadAffinityMask(gthand[gnthreadcnt],1<<gnthreadcnt);
        gnthreadcnt++;
    }

    //Start other threads
    ghtcallfunc = lcallfunc; ghtdata = luserdata; glincnt[0] = v0; glincnt[1] = v1;
    for(i=threaduse-1;i>=0;i--) SetEvent(ghevent[0][i]);

    //Do some processing in this thread too :)
    while ((i = getnextindexsynced()) < v1) lcallfunc(i,luserdata);

    //Wait for all other threads to finish (for safety reasons)
    WaitForMultipleObjects(threaduse,&ghevent[1][0],1,INFINITE);
}

#elif 0

#include "gtl/execution/thread_pool"

// Non-Windows use the thread_pool.
void htrun(void (*lcallfunc)(int, void *), void *luserdata, int v0, int v1, bool should_parallelise)
{

    if (!should_parallelise) {
        for(int i = v0; i < v1; ++i) {
            lcallfunc(i, luserdata);
        }
    }
    else {
        gtl::thread_pool pool;
        gtl::thread_pool::queue queue(pool);
        for(int i = v0; i < v1; ++i) {
            queue.push([=](){
                lcallfunc(i, luserdata);
            });
        }
        queue.drain();
        pool.join();
    }
}

#else

// Non-Windows fallback to the single threaded.
void htrun(void (*lcallfunc)(int, void *), void *luserdata, int v0, int v1, bool should_parallelise)
{
    static_cast<void>(should_parallelise);

    for(int i = v0; i < v1; ++i) {
        lcallfunc(i, luserdata);
    }
}

#endif
