// This file has been modified from Ken Silverman's original release

#pragma once
#ifndef GL_HPP
#define GL_HPP

#ifdef __linux__
//#   include <GL/gl.h>
#   define GLAPIENTRY
#endif

#if defined(_WIN32)
#   pragma comment(lib, "opengl32.lib")
#   define GLAPIENTRY __stdcall
#endif

#include <cstddef>
#include <cstdio>
#include "settings.hpp"

#define GL_MIRRORED_REPEAT                  0x8370
#define GL_TEXTURE_WRAP_R                   0x8072
#define GL_CLAMP_TO_EDGE                    0x812F
#define GL_TEXTURE_CUBE_MAP                 0x8513
#define GL_TEXTURE_CUBE_MAP_POSITIVE_X      0x8515
#define GL_TEXTURE_CUBE_MAP_NEGATIVE_X      0x8516
#define GL_TEXTURE_CUBE_MAP_POSITIVE_Y      0x8517
#define GL_TEXTURE_CUBE_MAP_NEGATIVE_Y      0x8518
#define GL_TEXTURE_CUBE_MAP_POSITIVE_Z      0x8519
#define GL_TEXTURE_CUBE_MAP_NEGATIVE_Z      0x851A
#define GL_TEXTURE_MAX_ANISOTROPY_EXT       0x84FE
#define GL_TEXTURE_3D                       0x806F
#define GL_FRAGMENT_SHADER                  0x8B30
#define GL_VERTEX_SHADER                    0x8B31
#define GL_COMPILE_STATUS                   0x8B81
#define GL_LINK_STATUS                      0x8B82
#define GL_VERTEX_PROGRAM_ARB               0x8620
#define GL_FRAGMENT_PROGRAM_ARB             0x8804
#define GL_PROGRAM_ERROR_STRING_ARB         0x8874
#define GL_PROGRAM_FORMAT_ASCII_ARB         0x8875
#define GL_TEXTURE0                         0x84c0
#define GL_FRAMEBUFFER_EXT                  0x8D40
#define GL_UNSIGNED_BYTE                    0x1401
#define GL_FLOAT                            0x1406
#define GL_TEXTURE_RECTANGLE_ARB            0x84F5
#define GL_RGBA32F_ARB                      0x8814
#define GL_LUMINANCE32F_ARB                 0x8818

#define GL_COLOR_ATTACHMENT0_EXT            0x8CE0
#define GL_PIXEL_UNPACK_BUFFER              0x88EC
#define GL_STREAM_DRAW                      0x88E0
#define GL_STREAM_READ                      0x88E1
#define GL_STREAM_COPY                      0x88E2
#define GL_STATIC_DRAW                      0x88E4
#define GL_STATIC_READ                      0x88E5
#define GL_STATIC_COPY                      0x88E6
#define GL_DYNAMIC_DRAW                     0x88E8
#define GL_DYNAMIC_READ                     0x88E9
#define GL_DYNAMIC_COPY                     0x88EA
#define GL_WRITE_ONLY                       0x88B9

#define GL_NO_ERROR                         0x0000
#define GL_NONE                             0x0000
#define GL_ONE                              0x0001
#define GL_TEXTURE_1D                       0x0DE0
#define GL_TEXTURE_2D                       0x0DE1
#define GL_TEXTURE_WRAP_S                   0x2802
#define GL_TEXTURE_WRAP_T                   0x2803
#define GL_TEXTURE_MAG_FILTER               0x2800
#define GL_TEXTURE_MIN_FILTER               0x2801
#define GL_CLAMP                            0x2900
#define GL_REPEAT                           0x2901
#define GL_NEAREST                          0x2600
#define GL_LINEAR                           0x2601
#define GL_NEAREST_MIPMAP_NEAREST           0x2700
#define GL_LINEAR_MIPMAP_NEAREST            0x2701
#define GL_NEAREST_MIPMAP_LINEAR            0x2702
#define GL_LINEAR_MIPMAP_LINEAR             0x2703
#define GL_TEXTURE_WRAP_S                   0x2802
#define GL_TEXTURE_WRAP_T                   0x2803
#define GL_BGRA_EXT                         0x80E1
#define GL_RGBA                             0x1908
#define GL_LUMINANCE                        0x1909
#define GL_LUMINANCE8                       0x8040
#define GL_LUMINANCE16                      0x8042
#define GL_UNSIGNED_SHORT                   0x1403
#define GL_MODELVIEW                        0x1700
#define GL_PROJECTION                       0x1701
#define GL_DEPTH_TEST                       0x0B71
#define GL_POINTS                           0x0000
#define GL_LINES                            0x0001
#define GL_TRIANGLE_FAN                     0x0006
#define GL_QUADS                            0x0007
#define GL_SRC_ALPHA                        0x0302
#define GL_ONE_MINUS_SRC_ALPHA              0x0303
#define GL_EXTENSIONS                       0x1F03
#define GL_BLEND                            0x0BE2
#define GL_COLOR_BUFFER_BIT                 0x00004000
#define GL_DEPTH_BUFFER_BIT                 0x00000100
#define GL_STENCIL_BUFFER_BIT               0x00000400

using GLvoid                                = void;
using GLboolean                             = unsigned char;
using GLchar                                = char;
using GLbyte                                = signed char;
using GLubyte                               = unsigned char;
using GLshort                               = short;
using GLushort                              = unsigned short;
using GLint                                 = int;
using GLuint                                = unsigned int;
// using GLfixed                              = ???;
using GLint64                               = long long int;
using GLuint64                              = unsigned long long int;
using GLsizei                               = int;
using GLenum                                = unsigned int;
using GLintptr                              = ptrdiff_t;
using GLsizeiptr                            = ptrdiff_t;
using GLsync                                = ptrdiff_t;
using GLbitfield                            = unsigned int;
//using GLhalf                                = ???;
using GLfloat                               = float;
using GLclampf                              = float;
using GLdouble                              = double;
using GLclampd                              = double;
using GLcharARB                             = char;
using GLhandleARB                           = unsigned int;

using PFNGLGETERROR                         = GLenum (GLAPIENTRY*)();
using PFNGLGETSTRING                        = const GLubyte* (GLAPIENTRY*)(GLenum);
using PFNGLTEXIMAGE1D                       = void (GLAPIENTRY*)(GLenum target, GLint level, GLint internalformat, GLsizei width, GLint border, GLenum format, GLenum type, const void* data);
using PFNGLTEXIMAGE2D                       = void (GLAPIENTRY*)(GLenum target, GLint level, GLint internalformat, GLsizei width, GLsizei height, GLint border, GLenum format, GLenum type, const void* data);
using PFNGLTEXSUBIMAGE1D                    = void (GLAPIENTRY*)(GLenum target, GLint level, GLint xoffset, GLsizei width, GLenum format, GLenum type, const GLvoid *pixels);
using PFNGLTEXSUBIMAGE2D                    = void (GLAPIENTRY*)(GLenum target, GLint level, GLint xoffset, GLint yoffset, GLsizei width, GLsizei height,GLenum format, GLenum type, const void * pixels);
using PFNGLBINDTEXTURE                      = void (GLAPIENTRY*)(GLenum target, GLuint texture);
using PFNGLTEXPARAMETERI                    = void (GLAPIENTRY*)(GLenum target, GLenum pname, GLint param);
using PFNGLGENTEXTURES                      = void (GLAPIENTRY*)(GLsizei n, GLuint *textures);
using PFNGLDELETETEXTURES                   = void (GLAPIENTRY*)(GLsizei n, const GLuint *textures);
using PFNGLMATRIXMODE                       = void (GLAPIENTRY*)(GLenum mode);
using PFNGLLOADIDENTITY                     = void (GLAPIENTRY*)();
using PFNGLORTHO                            = void (GLAPIENTRY*)(GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble near_val, GLdouble far_val);
using PFNGLFRUSTUM                          = void (GLAPIENTRY*)(GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble near_val, GLdouble far_val);
using PFNGLENABLE                           = void (GLAPIENTRY*)(GLenum cap);
using PFNGLDISABLE                          = void (GLAPIENTRY*)(GLenum cap);
using PFNGLBLENDFUNC                        = void (GLAPIENTRY*)(GLenum sfactor, GLenum dfactor);
using PFNGLCOLOR3UB                         = void (GLAPIENTRY*)(GLubyte red, GLubyte green, GLubyte blue);
using PFNGLCOLOR4UB                         = void (GLAPIENTRY*)(GLubyte red, GLubyte green, GLubyte blue, GLubyte alpha);
using PFNGLBEGIN                            = void (GLAPIENTRY*)(GLenum mode);
using PFNGLEND                              = void (GLAPIENTRY*)();
using PFNGLTEXCOORD2F                       = void (GLAPIENTRY*)(GLfloat s, GLfloat t);
using PFNGLTEXCOORD3F                       = void (GLAPIENTRY*)(GLfloat s, GLfloat t, GLfloat r);
using PFNGLVERTEX2I                         = void (GLAPIENTRY*)(GLint x, GLint y);
using PFNGLVERTEX2F                         = void (GLAPIENTRY*)(GLfloat x, GLfloat y);
using PFNGLVERTEX3F                         = void (GLAPIENTRY*)(GLfloat x, GLfloat y, GLfloat z);
using PFNGLCLEAR                            = void (GLAPIENTRY*)(GLbitfield mask);
using PFNGLCLEARCOLOR                       = void (GLAPIENTRY*)(GLclampf red, GLclampf green, GLclampf blue, GLclampf alpha);
using PFNGLGETPROGRAMINFOLOG                = void (GLAPIENTRY*)(GLuint program, GLsizei maxLength, GLsizei *length, GLchar *infoLog);

using PFNGLCREATESHADERPROC                 = GLuint (GLAPIENTRY*)(GLenum type);
using PFNGLCREATEPROGRAMPROC                = GLuint (GLAPIENTRY*)();
using PFNGLSHADERSOURCEPROC                 = void (GLAPIENTRY*)(GLuint shader, GLsizei count, const GLchar* const* strings, const GLint *lengths);
using PFNGLCOMPILESHADERPROC                = void (GLAPIENTRY*)(GLuint shader);
using PFNGLATTACHSHADERPROC                 = void (GLAPIENTRY*)(GLuint program, GLuint shader);
using PFNGLLINKPROGRAMPROC                  = void (GLAPIENTRY*)(GLuint program);
using PFNGLUSEPROGRAMPROC                   = void (GLAPIENTRY*)(GLuint program);
using PFNGLGETPROGRAMIVPROC                 = void (GLAPIENTRY*)(GLuint program, GLenum pname, GLint *param);
using PFNGLGETSHADERIVPROC                  = void (GLAPIENTRY*)(GLuint shader, GLenum pname, GLint *param);
using PFNGLGETINFOLOGARBPROC                = void (GLAPIENTRY*)(GLhandleARB obj, GLsizei maxLength, GLsizei *length, GLcharARB *infoLog);
using PFNGLDETACHSHADERPROC                 = void (GLAPIENTRY*)(GLuint program, GLuint shader);
using PFNGLDELETEPROGRAMPROC                = void (GLAPIENTRY*)(GLuint program);
using PFNGLDELETESHADERPROC                 = void (GLAPIENTRY*)(GLuint shader);
using PFNGLGETUNIFORMLOCATIONPROC           = GLint (GLAPIENTRY*)(GLuint program, const GLchar *name);
using PFNGLUNIFORM1FPROC                    = void (GLAPIENTRY*)(GLint location, GLfloat v0);
using PFNGLUNIFORM2FPROC                    = void (GLAPIENTRY*)(GLint location, GLfloat v0, GLfloat v1);
using PFNGLUNIFORM3FPROC                    = void (GLAPIENTRY*)(GLint location, GLfloat v0, GLfloat v1, GLfloat v2);
using PFNGLUNIFORM4FPROC                    = void (GLAPIENTRY*)(GLint location, GLfloat v0, GLfloat v1, GLfloat v2, GLfloat v3);
using PFNGLUNIFORM1IPROC                    = void (GLAPIENTRY*)(GLint location, GLint v0);
using PFNGLUNIFORM2IPROC                    = void (GLAPIENTRY*)(GLint location, GLint v0, GLint v1);
using PFNGLUNIFORM3IPROC                    = void (GLAPIENTRY*)(GLint location, GLint v0, GLint v1, GLint v2);
using PFNGLUNIFORM4IPROC                    = void (GLAPIENTRY*)(GLint location, GLint v0, GLint v1, GLint v2, GLint v3);
using PFNGLACTIVETEXTUREPROC                = void (GLAPIENTRY*)(GLuint texture);
using PFNGLTEXIMAGE3DPROC                   = void (GLAPIENTRY*)(GLenum, GLint, GLint, GLsizei, GLsizei, GLsizei, GLint, GLenum, GLenum, const GLvoid *);
using PFNGLTEXSUBIMAGE3DPROC                = void (GLAPIENTRY*)(GLenum, GLint, GLint, GLint, GLint, GLsizei, GLsizei, GLsizei, GLenum, GLenum, const GLvoid *);
using PFNGLGENPROGRAMSARBPROC               = void (GLAPIENTRY*)(GLsizei n, GLuint *programs);
using PFNGLBINDPROGRAMARBPROC               = void (GLAPIENTRY*)(GLenum target, GLuint program);
using PFNGLGETPROGRAMSTRINGARBPROC          = void (GLAPIENTRY*)(GLenum target, GLenum pname, GLvoid *string);
using PFNGLPROGRAMSTRINGARBPROC             = void (GLAPIENTRY*)(GLenum target, GLenum format, GLsizei len, const GLvoid *string);
using PFNGLPROGRAMLOCALPARAMETER4FARBPROC   = void (GLAPIENTRY*)(GLenum target, GLuint index, GLfloat x, GLfloat y, GLfloat z, GLfloat w);
using PFNGLPROGRAMENVPARAMETER4FARBPROC     = void (GLAPIENTRY*)(GLenum target, GLuint index, GLfloat x, GLfloat y, GLfloat z, GLfloat w);
using PFNGLDELETEPROGRAMSARBPROC            = void (GLAPIENTRY*)(GLsizei n, const GLuint *programs);
using PFNGLBINDBUFFER                       = void (GLAPIENTRY*)(GLenum, GLuint);
using PFNGLDELETEBUFFERS                    = void (GLAPIENTRY*)(GLsizei, const GLuint *);
using PFNGLGENBUFFERS                       = void (GLAPIENTRY*)(GLsizei, GLuint *);
using PFNGLBUFFERDATA                       = void (GLAPIENTRY*)(GLenum, GLsizeiptr, const GLvoid *, GLenum);
using PFNGLBUFFERSUBDATA                    = void (GLAPIENTRY*)(GLenum, GLintptr, GLsizeiptr, const GLvoid *);
using PFNGLMAPBUFFER                        = GLvoid* (GLAPIENTRY*)(GLenum, GLenum);
using PFNGLUNMAPBUFFER                      = GLboolean (GLAPIENTRY*)(GLenum);
using PFNGLBINDFRAMEBUFFEREXT               = void (GLAPIENTRY*)(GLenum, GLuint);
using PFNGLDELETEFRAMEBUFFERSEXT            = void (GLAPIENTRY*)(GLsizei, const GLuint *);
using PFNGLGENFRAMEBUFFERSEXT               = void (GLAPIENTRY*)(GLsizei, GLuint *);
using PFNGLFRAMEBUFFERTEXTURE2DEXT          = void (GLAPIENTRY*)(GLenum, GLenum, GLenum, GLuint, GLint);
using PFNGLTEXBUFFEREXT                     = void (GLAPIENTRY*)(GLenum target, GLenum internalformat, GLuint buffer);

using PFNWGLSWAPINTERVALEXTPROC             = GLint (GLAPIENTRY*)(GLint interval);

const static char *glnames[] =
{
    "glGetError",
    "glGetString",
    "glTexImage1D",
    "glTexImage2D",
    "glTexSubImage1D",
    "glTexSubImage2D",
    "glBindTexture",
    "glTexParameteri",
    "glGenTextures",
    "glDeleteTextures",
    "glMatrixMode",
    "glLoadIdentity",
    "glOrtho",
    "glFrustum",
    "glEnable",
    "glDisable",
    "glBlendFunc",
    "glColor3ub",
    "glColor4ub",
    "glBegin",
    "glEnd",
    "glTexCoord2f",
    "glTexCoord3f",
    "glVertex2i",
    "glVertex2f",
    "glVertex3f",
    "glClear",
    "glClearColor",
    "glGetProgramInfoLog",

    "glGenProgramsARB","glBindProgramARB",                                      //ARB ASM...
    "glGetProgramStringARB","glProgramStringARB",
    "glProgramLocalParameter4fARB",
    "glProgramEnvParameter4fARB",
    "glDeleteProgramsARB",

    "glActiveTexture","glTexImage3D","glTexSubImage3D",                         //multi/extended texture
    "wglSwapIntervalEXT",                                                       //limit refresh/sleep

    "glCreateShader","glCreateProgram",                                         //compile
    "glShaderSource","glCompileShader",                                         //
    "glAttachShader","glLinkProgram","glUseProgram",                            //link
    "glGetShaderiv","glGetProgramiv","glGetInfoLogARB",                         //get info
    "glDetachShader","glDeleteProgram","glDeleteShader",                        //decompile
    "glGetUniformLocation",                                                     //host->shader
    "glUniform1f" ,"glUniform2f" ,"glUniform3f" ,"glUniform4f",
    "glUniform1i" ,"glUniform2i" ,"glUniform3i" ,"glUniform4i",

    "glBindBuffer", "glDeleteBuffers", "glGenBuffers",                          //stuff for virtual mapping
    "glBufferData", "glBufferSubData",
    "glMapBuffer", "glUnmapBuffer",
    "glBindFramebufferEXT", "glDeleteFramebuffersEXT",
    "glGenFramebuffersEXT", "glFramebufferTexture2DEXT",

    "glTexBufferEXT",
};
const static char *glnames_old[] =
{
    "glGetError",
    "glGetString",
    "glTexImage1D",
    "glTexImage2D",
    "glTexSubImage1D",
    "glTexSubImage2D",
    "glBindTexture",
    "glTexParameteri",
    "glGenTextures",
    "glDeleteTextures",
    "glMatrixMode",
    "glLoadIdentity",
    "glOrtho",
    "glFrustum",
    "glEnable",
    "glDisable",
    "glBlendFunc",
    "glColor3ub",
    "glColor4ub",
    "glBegin",
    "glEnd",
    "glTexCoord2f",
    "glTexCoord3f",
    "glVertex2i",
    "glVertex2f",
    "glVertex3f",
    "glClear",
    "glClearColor",
    "glGetProgramInfoLog",

    "glGenProgramsARB","glBindProgramARB",                                      //ARB ASM...
    "glGetProgramStringARB","glProgramStringARB",
    "glProgramLocalParameter4fARB",
    "glProgramEnvParameter4fARB",
    "glDeleteProgramsARB",

    "glActiveTexture","glTexImage3D","glTexSubImage3D",                         //texture unit
    "wglSwapIntervalEXT",                                                       //limit refresh/sleep

    "glCreateShaderObjectARB","glCreateProgramObjectARB",                       //compile
    "glShaderSourceARB","glCompileShaderARB",                                   //
    "glAttachObjectARB","glLinkProgramARB","glUseProgramObjectARB",             //link
    "glGetObjectParameterivARB","glGetObjectParameterivARB","glGetInfoLogARB",  //get info
    "glDetachObjectARB","glDeleteObjectARB","glDeleteObjectARB",                //decompile
    "glGetUniformLocationARB",                                                  //host->shader
    "glUniform1fARB" ,"glUniform2fARB" ,"glUniform3fARB" ,"glUniform4fARB",
    "glUniform1iARB" ,"glUniform2iARB" ,"glUniform3iARB" ,"glUniform4iARB",

    "glBindBuffer", "glDeleteBuffers", "glGenBuffers",                          //stuff for virtual mapping
    "glBufferData", "glBufferSubData",
    "glMapBuffer", "glUnmapBuffer",
    "glBindFramebufferEXT", "glDeleteFramebuffersEXT",
    "glGenFramebuffersEXT", "glFramebufferTexture2DEXT",

    "glTexBufferEXT",
};
enum
{
    glGetError,
    glGetString,
    glTexImage1D,
    glTexImage2D,
    glTexSubImage1D,
    glTexSubImage2D,
    glBindTexture,
    glTexParameteri,
    glGenTextures,
    glDeleteTextures,
    glMatrixMode,
    glLoadIdentity,
    glOrtho,
    glFrustum,
    glEnable,
    glDisable,
    glBlendFunc,
    glColor3ub,
    glColor4ub,
    glBegin,
    glEnd,
    glTexCoord2f,
    glTexCoord3f,
    glVertex2i,
    glVertex2f,
    glVertex3f,
    glClear,
    glClearColor,
    glGetProgramInfoLog,

    glGenProgramsARB,glBindProgramARB,             //ARB ASM...
    glGetProgramStringARB,glProgramStringARB,
    glProgramLocalParameter4fARB,
    glProgramEnvParameter4fARB,
    glDeleteProgramsARB,

    glActiveTexture,glTexImage3D,glTexSubImage3D,  //multi/extended texture
    wglSwapIntervalEXT,                            //limit refesh/sleep

    glCreateShader,glCreateProgram,                //compile
    glShaderSource,glCompileShader,                //
    glAttachShader,glLinkProgram,glUseProgram,     //link
    glGetShaderiv,glGetProgramiv,glGetInfoLogARB,  //get info
    glDetachShader,glDeleteProgram,glDeleteShader, //decompile
    glGetUniformLocation,                          //host->shader
    glUniform1f, glUniform2f, glUniform3f, glUniform4f,
    glUniform1i, glUniform2i, glUniform3i, glUniform4i,

    glBindBuffer, glDeleteBuffers, glGenBuffers,   //stuff for virtual mapping
    glBufferData, glBufferSubData,
    glMapBuffer, glUnmapBuffer,
    glBindFramebufferEXT, glDeleteFramebuffersEXT,
    glGenFramebuffersEXT, glFramebufferTexture2DEXT,

    glTexBufferEXT,

    NUMGLFUNC
};

///////////////////////////////////////////////////////////////////////////////

#if defined(_WIN32)
    #define NOMINMAX
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
#endif

// TODO: REMOVE!
#include "macros.hpp"
#include <cstdlib>
// TODO: REMOVE!

///////////////////////////////////////////////////////////////////////////////

#if defined(_WIN32)
    //using proc = int (__stdcall*)();  // returns 32 or 64 bit maybe?
    //extern "C" void* __stdcall LoadLibraryA(const char*);
    //extern "C" proc __stdcall GetProcAddress(void*, const char*);
    //extern "C" proc __stdcall wglGetProcAddress(const char*);
    //extern "C" bool __stdcall FreeLibrary(void*);
#else
    using proc = void (*)();
    extern "C" proc glXGetProcAddressARB(const GLubyte *procName);
#endif

static void* GetAnyGLFuncAddress(const char* name)
{
#if defined(_WIN32)
    void *p = (void *)wglGetProcAddress(name);
    if(p == 0 || (p == (void*)0x1) || (p == (void*)0x2) || (p == (void*)0x3) || (p == (void*)-1))
    {
        HMODULE module = LoadLibraryA("opengl32.dll");
        p = (void *)GetProcAddress(module, name);
        FreeLibrary(module);
    }

    return p;
#else
    return (void*)glXGetProcAddressARB((const GLubyte*)name);
#endif
}

///////////////////////////////////////////////////////////////////////////////

using glfp_t = void (*)();
extern glfp_t glfp[NUMGLFUNC];

static bool LoadGL() {
    bool useoldglfuncs = false;
    for (int i = 0; i < NUMGLFUNC; i++) {
        if (!useoldglfuncs) {
            glfp[i] = (glfp_t)GetAnyGLFuncAddress(glnames[i]);
            if (glfp[i]) {
                continue;
            }
            useoldglfuncs = true;
        }
        glfp[i] = (glfp_t)GetAnyGLFuncAddress(glnames_old[i]);
        if (glfp[i]) {
            continue;
        }
        fprintf(stderr, "%s() / %s() not supported. :/\n", glnames[i], glnames_old[i]);
        return false;
    }
    return true;
}

///////////////////////////////////////////////////////////////////////////////

static void kgluninit (HWND hWnd, HDC hDC, HGLRC hRC)
{
#ifdef _WIN32
    wglMakeCurrent(0,0);
    wglDeleteContext(hRC);
    ReleaseDC(hWnd,hDC);
#endif
}

static int kglinit (HWND hwnd, int swapinterval, HDC *hDC, HGLRC *hRC)
{
#ifdef _WIN32
    PIXELFORMATDESCRIPTOR pfd;
    int format;

    *hDC = GetDC(hwnd);

    ZeroMemory(&pfd,sizeof(pfd));
    pfd.nSize = sizeof(pfd);
    pfd.nVersion = 1;
    pfd.dwFlags = PFD_DRAW_TO_WINDOW|PFD_SUPPORT_OPENGL|PFD_DOUBLEBUFFER;
    pfd.iPixelType = PFD_TYPE_RGBA;
    pfd.cColorBits = 24;
    pfd.cDepthBits = 32;
    pfd.iLayerType = PFD_MAIN_PLANE;
    format = ChoosePixelFormat(*hDC,&pfd);
    SetPixelFormat(*hDC,format,&pfd);

    *hRC = wglCreateContext(*hDC);
    wglMakeCurrent(*hDC,*hRC);

    if (!GetAnyGLFuncAddress("glCreateShaderObjectARB"))
    {
        //oct_useglsl = 0;
        std::abort();
    }
    else
    {
        glfp[glGetString] = (glfp_t)GetAnyGLFuncAddress("glGetString");
        const unsigned char* cptr = ((PFNGLGETSTRING)glfp[glGetString])(0x8B8C); // GL_SHADING_LANGUAGE_VERSION
        if (cptr[0] < '3') {
            kgluninit(hwnd, *hDC, *hRC);
            std::abort();
        }
    }

    if (!LoadGL()) {
        return -1;
    }

    if (glfp[wglSwapIntervalEXT]) {
        ((PFNWGLSWAPINTERVALEXTPROC)glfp[wglSwapIntervalEXT])(swapinterval);
    }
#else
    if (!GetAnyGLFuncAddress("glCreateShaderObjectARB"))
    {
        //oct_useglsl = 0;
        std::abort();
    }
    else
    {
        glfp[glGetString] = (glfp_t)GetAnyGLFuncAddress("glGetString");
        const unsigned char* cptr = ((PFNGLGETSTRING)glfp[glGetString])(0x8B8C); // GL_SHADING_LANGUAGE_VERSION
        if (cptr[0] < '3') {
            kgluninit(hwnd, *hDC, *hRC);
            std::abort();
        }
    }

    if (!LoadGL()) {
        return -1;
    }

    if (glfp[wglSwapIntervalEXT]) {
        ((PFNWGLSWAPINTERVALEXTPROC)glfp[wglSwapIntervalEXT])(swapinterval);
    }
#endif
    return (0);
}
///////////////////////////////////////////////////////////////////////////////

//static bool glcheckext(char* extnam) {
//    // Tiny strlen.
//    auto strlen = [](const char* string) -> int {
//        const char* start = string;
//        while (*string) {
//            ++string;
//        }
//        return static_cast<int>(string - start);
//    };
//
//    const char* st = (char*)((PFNGLGETSTRING)glfp[glGetString])(GL_EXTENSIONS);
//
//    st = strstr(st, extnam);
//    if (!st) {
//        return (0);
//    }
//    return (st[strlen(extnam)] <= 32);
//}

///////////////////////////////////////////////////////////////////////////////

static bool compileshaderGLSL(const char* source_vertex, const char* source_fragment, GLenum& shader_vertex, GLenum& shader_fragment, GLenum& shader_program)
{
    // Initialise outputs.
    shader_vertex = 0;
    shader_fragment = 0;
    shader_program = 0;

    // Ensure we are in GLSL mode.
    ((PFNGLDISABLE)glfp[glDisable])(GL_VERTEX_PROGRAM_ARB);
    ((PFNGLDISABLE)glfp[glDisable])(GL_FRAGMENT_PROGRAM_ARB);

    // Create the vertex program.
    shader_vertex = ((PFNGLCREATESHADERPROC)(glfp[glCreateShader]))(GL_VERTEX_SHADER);
    if (shader_vertex == 0) {
        fprintf(stderr, "Failed to create vertex shader.\n");
        return false;
    }

    // Compile the vertex program.
    ((PFNGLSHADERSOURCEPROC)glfp[glShaderSource])(shader_vertex, 1, (const char**)&source_vertex, 0);
    ((PFNGLCOMPILESHADERPROC)glfp[glCompileShader])(shader_vertex);

    // Check the compilation.
    GLint compiled_vertex = 0;
    ((PFNGLGETSHADERIVPROC)glfp[glGetShaderiv])(shader_vertex, GL_COMPILE_STATUS, &compiled_vertex);
    if (compiled_vertex == 0) {
        char error_vertex[4096];
        ((PFNGLGETINFOLOGARBPROC)glfp[glGetInfoLogARB])(shader_vertex, sizeof(error_vertex), 0, error_vertex);
        fprintf(stderr, "Failed to compile vertex shader:\n%s\n", error_vertex);
        // Cleanup.
        ((PFNGLDELETESHADERPROC)glfp[glDeleteShader])(shader_vertex);
        return false;
    }

    // Create the fragment program.
    shader_fragment = ((PFNGLCREATESHADERPROC)(glfp[glCreateShader]))(GL_FRAGMENT_SHADER);
    if (shader_fragment == 0) {
        fprintf(stderr, "Failed to create fragment shader.\n");
        // Cleanup.
        ((PFNGLDELETESHADERPROC)glfp[glDeleteShader])(shader_vertex);
        return false;
    }

    // Compile the fragment program.
    ((PFNGLSHADERSOURCEPROC)glfp[glShaderSource])(shader_fragment, 1, (const char**)&source_fragment, 0);
    ((PFNGLCOMPILESHADERPROC)glfp[glCompileShader])(shader_fragment);

    // Check the compilation.
    GLint compiled_fragment = 0;
    ((PFNGLGETSHADERIVPROC)glfp[glGetShaderiv])(shader_fragment, GL_COMPILE_STATUS, &compiled_fragment);
    if (compiled_fragment == 0) {
        char error_fragment[4096];
        ((PFNGLGETINFOLOGARBPROC)glfp[glGetInfoLogARB])(shader_fragment, sizeof(error_fragment), 0, error_fragment);
        fprintf(stderr, "Failed to compile fragment shader:\n%s\n", error_fragment);
        // Cleanup.
        ((PFNGLDELETESHADERPROC)glfp[glDeleteShader])(shader_vertex);
        ((PFNGLDELETESHADERPROC)glfp[glDeleteShader])(shader_fragment);
        return false;
    }

    // Create the shader program.
    shader_program = ((PFNGLCREATEPROGRAMPROC)glfp[glCreateProgram])();
    if (shader_program == 0) {
        fprintf(stderr, "Failed to create shader program.\n");
        // Cleanup.
        ((PFNGLDELETESHADERPROC)glfp[glDeleteShader])(shader_vertex);
        ((PFNGLDELETESHADERPROC)glfp[glDeleteShader])(shader_fragment);
        return false;
    }

    // Attach the shaders.
    ((PFNGLATTACHSHADERPROC)glfp[glAttachShader])(shader_program, shader_vertex);
    ((PFNGLATTACHSHADERPROC)glfp[glAttachShader])(shader_program, shader_fragment);

    // Link the program.
    ((PFNGLLINKPROGRAMPROC)glfp[glLinkProgram])(shader_program);

    // Check the linking.
    GLint linked_program = 0;
    ((PFNGLGETPROGRAMIVPROC)glfp[glGetProgramiv])(shader_program, GL_LINK_STATUS, &linked_program);
    if (linked_program == 0) {
        char error_program[4096];
        ((PFNGLGETPROGRAMINFOLOG)glfp[glGetProgramInfoLog])(shader_program, sizeof(error_program), 0, error_program);
        fprintf(stderr, "Failed to link shader program:\n%s\n", error_program);
        // Cleanup.
        ((PFNGLDELETESHADERPROC)glfp[glDeleteShader])(shader_vertex);
        ((PFNGLDELETESHADERPROC)glfp[glDeleteShader])(shader_fragment);
        ((PFNGLDELETEPROGRAMPROC)glfp[glDeleteProgram])(shader_program);
        return false;
    }

    return true;
}

static bool compileshaderASM(const char* source_vertex, const char* source_fragment, GLenum& shader_vertex, GLenum& shader_fragment)
{
    // Initialise outputs.
    shader_vertex = 0;
    shader_fragment = 0;

    // Tiny strlen.
    auto strlen = [](const char* string) -> int {
        const char* start = string;
        while (*string) {
            ++string;
        }
        return static_cast<int>(string - start);
    };

    // Ensure we are in ASM mode, aka not using a program.
    if (glfp[glUseProgram]) {
        ((PFNGLUSEPROGRAMPROC)glfp[glUseProgram])(0);
    }
    ((PFNGLENABLE)glfp[glEnable])(GL_VERTEX_PROGRAM_ARB);
    ((PFNGLENABLE)glfp[glEnable])(GL_FRAGMENT_PROGRAM_ARB);

    //Flush errors.
    ((PFNGLGETERROR)glfp[glGetError])();

    // Create vertex program.
    ((PFNGLGENPROGRAMSARBPROC)glfp[glGenProgramsARB])(1, &shader_vertex);

    // Load the program string.
    ((PFNGLBINDPROGRAMARBPROC)glfp[glBindProgramARB])(GL_VERTEX_PROGRAM_ARB, shader_vertex);
    ((PFNGLPROGRAMSTRINGARBPROC)glfp[glProgramStringARB])(GL_VERTEX_PROGRAM_ARB, GL_PROGRAM_FORMAT_ASCII_ARB, strlen(source_vertex), source_vertex);
    if (((PFNGLGETERROR)glfp[glGetError])() != GL_NO_ERROR) {
        fprintf(stderr, "Failed to compile vertex program:\n%s\n", reinterpret_cast<const char*>(((PFNGLGETSTRING)glfp[glGetString])(GL_PROGRAM_ERROR_STRING_ARB)));
        return false;
    }

    // Create the fragment program.
    ((PFNGLGENPROGRAMSARBPROC)glfp[glGenProgramsARB])(1, &shader_fragment);

    // Load the program string.
    ((PFNGLBINDPROGRAMARBPROC)glfp[glBindProgramARB])(GL_FRAGMENT_PROGRAM_ARB, shader_fragment);
    ((PFNGLPROGRAMSTRINGARBPROC)glfp[glProgramStringARB])(GL_FRAGMENT_PROGRAM_ARB, GL_PROGRAM_FORMAT_ASCII_ARB, strlen(source_fragment), source_fragment);
    if (((PFNGLGETERROR)glfp[glGetError])() != GL_NO_ERROR) {
        fprintf(stderr, "Failed to compile fragment program:\n%s\n", reinterpret_cast<const char*>(((PFNGLGETSTRING)glfp[glGetString])(GL_PROGRAM_ERROR_STRING_ARB)));
        return false;
    }

    return true;
}

///////////////////////////////////////////////////////////////////////////////

//Pixel Buffer Object (PBO); see: http://www.mathematik.uni-dortmund.de/~goeddeke/gpgpu/tutorial3.html
static void bo_uninit (unsigned int bufid) {
    ((PFNGLDELETEBUFFERS)glfp[glDeleteBuffers])(1, &bufid);
}

static unsigned int bo_init (int nbytes)
{
    unsigned int bufid;
    ((PFNGLGENBUFFERS)glfp[glGenBuffers])(1,&bufid);
    ((PFNGLBINDBUFFER)glfp[glBindBuffer])(GL_PIXEL_UNPACK_BUFFER,bufid);
    ((PFNGLBUFFERDATA)glfp[glBufferData])(GL_PIXEL_UNPACK_BUFFER,nbytes,NULL,GL_DYNAMIC_COPY);
    return(bufid);
}

static void* bo_begin (int bufid, int nbytes)
{
    ((PFNGLBINDBUFFER)glfp[glBindBuffer])(GL_PIXEL_UNPACK_BUFFER,bufid);
    if (nbytes)
    {     //NOTE:GL_STREAM_READ slightly faster than GL_STREAM_DRAW, but seems strange?
        ((PFNGLBUFFERDATA)glfp[glBufferData])(GL_PIXEL_UNPACK_BUFFER,nbytes,NULL,GL_STREAM_READ);
    }
    void* v = ((PFNGLMAPBUFFER)glfp[glMapBuffer])(GL_PIXEL_UNPACK_BUFFER,GL_WRITE_ONLY);
    if (!v) {
        fprintf(stderr, "((PFNGLMAPBUFFER)glfp[glMapBuffer])() error 0x%08x\n", ((PFNGLGETERROR)glfp[glGetError])());
    }
    return (v);
}

static void bo_end (int bufid, int xoffs, int yoffs, int xsiz, int ysiz, int k0, int k1, int offs)
{
    ((PFNGLBINDBUFFER)glfp[glBindBuffer])(GL_PIXEL_UNPACK_BUFFER, bufid);
    ((PFNGLUNMAPBUFFER)glfp[glUnmapBuffer])(GL_PIXEL_UNPACK_BUFFER);
    if (xsiz) {
        ((PFNGLTEXSUBIMAGE2D)glfp[glTexSubImage2D])(GL_TEXTURE_2D, 0, xoffs, yoffs, xsiz, ysiz, k0, k1, ((char*)NULL + (offs))); //fast - pure GPU copy
    }
    ((PFNGLBINDBUFFER)glfp[glBindBuffer])(GL_PIXEL_UNPACK_BUFFER, 0); //MUST unbind buffer object (failure to do so causes future unrelated ((PFNGLTEXSUBIMAGE2D)glfp[glTexSubImage2D])() to be corrupt!)
}

///////////////////////////////////////////////////////////////////////////////

static void kglActiveTexture(int texunit) {
    if (glfp[glActiveTexture]) {
        ((PFNGLACTIVETEXTUREPROC)glfp[glActiveTexture])((texunit & 3) + GL_TEXTURE0);
    }
}

static void kglProgramLocalParam(unsigned ind, float a, float b, float c, float d) {
    ((PFNGLPROGRAMLOCALPARAMETER4FARBPROC)glfp[glProgramLocalParameter4fARB])(GL_VERTEX_PROGRAM_ARB, ind, a, b, c, d);
    ((PFNGLPROGRAMLOCALPARAMETER4FARBPROC)glfp[glProgramLocalParameter4fARB])(GL_FRAGMENT_PROGRAM_ARB, ind, a, b, c, d);
}

static void kglProgramEnvParam(unsigned ind, float a, float b, float c, float d) {
    ((PFNGLPROGRAMENVPARAMETER4FARBPROC)glfp[glProgramEnvParameter4fARB])(GL_VERTEX_PROGRAM_ARB, ind, a, b, c, d);
    ((PFNGLPROGRAMENVPARAMETER4FARBPROC)glfp[glProgramEnvParameter4fARB])(GL_FRAGMENT_PROGRAM_ARB, ind, a, b, c, d);
}

static int kglGetUniformLoc(unsigned ind, char* shadvarnam) {
    return (((PFNGLGETUNIFORMLOCATIONPROC)glfp[glGetUniformLocation])(ind, shadvarnam));
}

static void kglUniform1f(unsigned sh, float v0) {
    ((PFNGLUNIFORM1FPROC)glfp[glUniform1f])(sh, v0);
}

static void kglUniform2f(unsigned sh, float v0, float v1) {
    ((PFNGLUNIFORM2FPROC)glfp[glUniform2f])(sh, v0, v1);
}

static void kglUniform3f(unsigned sh, float v0, float v1, float v2) {
    ((PFNGLUNIFORM3FPROC)glfp[glUniform3f])(sh, v0, v1, v2);
}

static void kglUniform4f(unsigned sh, float v0, float v1, float v2, float v3) {
    ((PFNGLUNIFORM4FPROC)glfp[glUniform4f])(sh, v0, v1, v2, v3);
}

static void kglUniform1i(unsigned sh, int v0) {
    ((PFNGLUNIFORM1IPROC)glfp[glUniform1i])(sh, v0);
}

static void kglUniform2i(unsigned sh, int v0, int v1) {
    ((PFNGLUNIFORM2IPROC)glfp[glUniform2i])(sh, v0, v1);
}

static void kglUniform3i(unsigned sh, int v0, int v1, int v2) {
    ((PFNGLUNIFORM3IPROC)glfp[glUniform3i])(sh, v0, v1, v2);
}

static void kglUniform4i(unsigned sh, int v0, int v1, int v2, int v3) {
    ((PFNGLUNIFORM4IPROC)glfp[glUniform4i])(sh, v0, v1, v2, v3);
}

///////////////////////////////////////////////////////////////////////////////

enum {
    KGL_BGRA32=0,
    KGL_RGBA32, // faster xfer than KGL_BGRA32!
    KGL_CHAR,
    KGL_SHORT,
    KGL_FLOAT,
    KGL_VEC4,
    KGL_NUM
};

enum {
    KGL_LINEAR = (0<<4),
    KGL_NEAREST = (1<<4),
    KGL_MIPMAP = (2<<4),
    KGL_MIPMAP3 = (2<<4),
    KGL_MIPMAP2 = (3<<4),
    KGL_MIPMAP1 = (4<<4),
    KGL_MIPMAP0 = (5<<4)
};

enum {
    KGL_REPEAT = (0<<8),
    KGL_MIRRORED_REPEAT = (1<<8),
    KGL_CLAMP = (2<<8),
    KGL_CLAMP_TO_EDGE = (3<<8)
};

static const int cubemapconst[6] =
{
    GL_TEXTURE_CUBE_MAP_POSITIVE_X, GL_TEXTURE_CUBE_MAP_NEGATIVE_X,
    GL_TEXTURE_CUBE_MAP_POSITIVE_Y, GL_TEXTURE_CUBE_MAP_NEGATIVE_Y,
    GL_TEXTURE_CUBE_MAP_POSITIVE_Z, GL_TEXTURE_CUBE_MAP_NEGATIVE_Z,
};
static const int cubemapindex[6] = {1,3,4,5,0,2};
static void kglalloctex (int itex, void *p, int xs, int ys, int zs, int icoltype)
{
    int i, targ, format, type, internalFormat = 0;

    targ = GL_TEXTURE_3D;
    if (zs == 1)
    {
        targ = GL_TEXTURE_2D;
        if (xs*6 == ys)
        {
            targ = GL_TEXTURE_CUBE_MAP;
            icoltype = (icoltype&~0xf00)|KGL_CLAMP_TO_EDGE;
            if ((icoltype&0xf0) >= KGL_MIPMAP) icoltype = (icoltype&~0xf0)|KGL_LINEAR;
        }
        else if (ys == 1) {
            targ = GL_TEXTURE_1D;
        }
    }

    if constexpr (oct_usegpubo) {
        if (!glfp[glBindBuffer]) {
            fprintf(stderr, "glBindBuffer() not supported .. enjoy the impending crash.\n");
            std::abort();
        }
        ((PFNGLBINDBUFFER)glfp[glBindBuffer])(GL_PIXEL_UNPACK_BUFFER, 0); //MUST unbind buffer object here, before usage of kglalloctex!
    }

    //glEnable(targ);
    ((PFNGLBINDTEXTURE)glfp[glBindTexture])(targ,itex);
    switch (icoltype&0xf0)
    {
        case KGL_LINEAR: default:
            ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
            ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
            break;
        case KGL_NEAREST:
            ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
            ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
            break;
        case KGL_MIPMAP0: case KGL_MIPMAP1: case KGL_MIPMAP2: case KGL_MIPMAP3:
            switch(icoltype&0xf0)
            {
                case KGL_MIPMAP0:
                    ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_MIN_FILTER,GL_NEAREST_MIPMAP_NEAREST);
                    break;
                case KGL_MIPMAP1:
                    ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_MIN_FILTER,GL_NEAREST_MIPMAP_LINEAR);
                    break;
                case KGL_MIPMAP2:
                    ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_MIN_FILTER,GL_LINEAR_MIPMAP_NEAREST);
                    break;
                case KGL_MIPMAP3:
                    ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_MIN_FILTER,GL_LINEAR_MIPMAP_LINEAR);
                    break;
            }
            ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
            break;
    }
    switch(icoltype&0xf00)
    {
        case KGL_REPEAT: default:
            ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_WRAP_S,GL_REPEAT);
            ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_WRAP_T,GL_REPEAT);
            ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_WRAP_R,GL_REPEAT);
            break;
        case KGL_MIRRORED_REPEAT:
            ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_WRAP_S,GL_MIRRORED_REPEAT);
            ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_WRAP_T,GL_MIRRORED_REPEAT);
            ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_WRAP_R,GL_MIRRORED_REPEAT);
            break;
        case KGL_CLAMP:
            ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_WRAP_S,GL_CLAMP);
            ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_WRAP_T,GL_CLAMP);
            ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_WRAP_R,GL_CLAMP);
            break;
        case KGL_CLAMP_TO_EDGE:
            ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
            ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
            ((PFNGLTEXPARAMETERI)glfp[glTexParameteri])(targ,GL_TEXTURE_WRAP_R,GL_CLAMP_TO_EDGE);
            break;
    }
    //glTexParameterf(targ,GL_TEXTURE_MAX_ANISOTROPY_EXT,1.f); //disable anisotropy for CPU raycast/GPU shader algo! (values > 1.f cause horrible 2x2 edge artifacts!)
    switch(icoltype&15)
    {
        case KGL_BGRA32: internalFormat =                   4; format =  GL_BGRA_EXT; type = GL_UNSIGNED_BYTE; break;
        case KGL_RGBA32: internalFormat =                   4; format =  GL_RGBA    ; type = GL_UNSIGNED_BYTE; break;
        case KGL_CHAR:   internalFormat =       GL_LUMINANCE8; format = GL_LUMINANCE; type = GL_UNSIGNED_BYTE; break;
        case KGL_SHORT:  internalFormat =      GL_LUMINANCE16; format = GL_LUMINANCE; type =GL_UNSIGNED_SHORT; break;
        case KGL_FLOAT:  internalFormat = GL_LUMINANCE32F_ARB; format = GL_LUMINANCE; type =         GL_FLOAT; break;
        case KGL_VEC4:   internalFormat =      GL_RGBA32F_ARB; format =      GL_RGBA; type =         GL_FLOAT; break;
    }
    switch(targ)
    {
        case GL_TEXTURE_1D: ((PFNGLTEXIMAGE1D)glfp[glTexImage1D])(targ,0,internalFormat,xs,   0,format,type,p); break;
        case GL_TEXTURE_2D: ((PFNGLTEXIMAGE2D)glfp[glTexImage2D])(targ,0,internalFormat,xs,ys,0,format,type,p); break;
        case GL_TEXTURE_3D: ((PFNGLTEXIMAGE3DPROC)glfp[glTexImage3D])(targ,0,internalFormat,xs,ys,zs,0,format,type,p); break;
        case GL_TEXTURE_CUBE_MAP:
            for(i=0;i<6;i++) {
                ((PFNGLTEXIMAGE2D)glfp[glTexImage2D])(cubemapconst[i],0,internalFormat,xs,xs,0,format,type,(void *)((INT_PTR)p + xs*xs*4*cubemapindex[i]));
            }
            //for(i=0;i<6;i++) {
            //    ((PFNGLTEXSUBIMAGE2D)glfp[glTexSubImage2D])(cubemapconst[i],0,0,0,xs,xs,format,type,(void *)(((INT_PTR)p)+xs*xs*4*cubemapindex[i]));
            //}
            break;
    }
}

#endif // GL_HPP
