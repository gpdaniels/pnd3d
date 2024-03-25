// This file has been modified from Ken Silverman's original release

#pragma once
#ifndef SELF_MODIFYING_CODE_HPP
#define SELF_MODIFYING_CODE_HPP

#if defined(_WIN32)

#if 0
#define SELFMODVAL(dalab,daval) \
{  void *_daddr; unsigned int _oprot; \
    _asm { mov _daddr, offset dalab } \
    VirtualProtect(_daddr,sizeof(daval),PAGE_EXECUTE_READWRITE,&_oprot); \
    switch(sizeof(daval)) \
    {  case 1: *(char *)_daddr = daval; break; \
        case 2: *(short *)_daddr = daval; break; \
        case 4: *(int *)_daddr = daval; break; \
        case 8: *(int64_t *)_daddr = daval; break; \
    } \
    VirtualProtect(_daddr,sizeof(daval),PAGE_EXECUTE,&_oprot); \
    /*FlushInstructionCache(GetCurrentProcess(),_daddr,sizeof(daval));*/ \
}
#else
#define SELFMODVAL(init,dalab,daval) \
{  void *_daddr; unsigned long _oprot; \
    _asm { mov _daddr, offset dalab } \
    if (init) VirtualProtect(_daddr,sizeof(daval),PAGE_EXECUTE_READWRITE,&_oprot); \
    switch(sizeof(daval)) \
    {  case 1: *(char *)_daddr = daval; break; \
        case 2: *(short *)_daddr = daval; break; \
        case 4: *(int *)_daddr = daval; break; \
        case 8: *(int64_t *)_daddr = daval; break; \
    } \
}
#endif

#else

#define SELFMODVAL(init,dalab,daval) /*TODO*/

#endif

#endif // SELF_MODIFYING_CODE_HPP
