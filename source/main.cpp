// This file has been modified from Ken Silverman's original release

#include "octree_physics.hpp"
#include "octree_modify.hpp"
#include "octree_renderer.hpp"

#include "brushes.hpp"
#include "types/sprite.hpp"
#include "utilities/gl.hpp"
#include "utilities/timestamp.hpp"
#include "utilities/window.hpp"
#include "utilities/macros.hpp"

#if defined(_WIN32)
    #define WIN32_LEAN_AND_MEAN
    #define NOMINMAX
    #include <windows.h>
#endif

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <memory>

static point3d ipos, irig, idow, ifor, ihaf;

#define OCTMAX 256
static oct_t oct[OCTMAX];
static int octnum = 0;

#define SPRMAX 256
static spr_t spr[SPRMAX];
static int sprnum = 0;

static INT_PTR htex = 0;

#define SURFPTR(oct,ind) (&((surf_t *)oct->sur.buf)[ind])

int WINAPI WinMain(HINSTANCE hinst, HINSTANCE hpinst, LPSTR cmdline, int ncmdshow);

int main(int argc, char* argv[]) {
    static_cast<void>(argc);
    static_cast<void>(argv);
    return WinMain(0, 0, 0, 0);
}

int WINAPI WinMain(HINSTANCE hinst, HINSTANCE hpinst, LPSTR cmdline, int ncmdshow) {
    static_cast<void>(hpinst);
    static_cast<void>(cmdline);
    static_cast<void>(ncmdshow);

    point3d fp, vec;
    point3i hit;
    double d, tim = 0, otim, dtim, avgdtim = 0.0;
    float f;
    int isurf, obstatus = 0;

    window_state window;

    window.xres = 1024;
    window.yres = 768;
    window.prognam = "HX";
    if (!window.initapp(hinst)) {
        return (-1);
    }

    oct_initonce(window.ghwnd, 256.0);

#if 0
    gskyid = oct_loadtex("kensky.jpg", KGL_BGRA32 + KGL_LINEAR);
#else
    int image_width = 256;
    int image_height = 1536;
    std::unique_ptr<unsigned int[]> image = std::unique_ptr<unsigned int[]>(new unsigned int[image_width * image_height]);
    for (int image_y = 0; image_y < image_height; ++image_y) {
        for (int image_x = 0; image_x < image_width; ++image_x) {
            image[image_y * image_width + image_x] = (0x00 << 24) | (image_y << 16) | (image_x << 8) | 0xFF;
        }
    }
    gskyid = oct_loadtex((INT_PTR)image.get(), image_width, image_height, KGL_RGBA32 + KGL_LINEAR);
#endif

#if 1
    oct_fogdist = 1e32;
    oct_fogcol = 0xff808080; //sky
#else
    oct_fogdist = 10.0;
    oct_fogcol = 0x80808080; //fog
#endif

    ipos.x = 0.0;
    ipos.y = 0.0;
    ipos.z = -1.0;
    irig.x = 1.0;
    irig.y = 0.0;
    irig.z = 0.0;
    idow.x = 0.0;
    idow.y = 1.0;
    idow.z = 0.0;
    ifor.x = 0.0;
    ifor.y = 0.0;
    ifor.z = 1.0;
    ihaf.x = window.xres * .5;
    ihaf.y = window.yres * .5;
    ihaf.z = ihaf.x;
    sprnum = 0;
    octnum = 0;
    //-----------------------------------------------------------------------------------------------

#if 0
    oct_load(window.xres, window.yres, &oct[octnum], "caco.kvo", &sp, &sr, &sd, &sf);
#else
    int oct_size = 7;// 11 works, 10 previous default, 7 fast.
    oct_new(window.xres, window.yres, &oct[octnum], oct_size, 0, 0, 0, 0);

    d = 1.0 / (double)oct[octnum].sid;
    spr[sprnum].br.x *= d;
    spr[sprnum].br.y *= d;
    spr[sprnum].br.z *= d;
    spr[sprnum].bd.x *= d;
    spr[sprnum].bd.y *= d;
    spr[sprnum].bd.z *= d;
    spr[sprnum].bf.x *= d;
    spr[sprnum].bf.y *= d;
    spr[sprnum].bf.z *= d;

    brush_sph_t sph;
    brush_box_t box;
    brush_cone_t cone;

    for (int i = (64 - 1); i >= 0; i--) {
        int x = rand() & (oct[octnum].sid - 1);
        int y = rand() & (oct[octnum].sid - 1);
        int z = rand() & (oct[octnum].sid - 1);
        int s = (rand() & ((oct[octnum].sid >> 4) - 1)) + (oct[octnum].sid >> 4);

        if (((float)rand() / (float)RAND_MAX) > 0.8) {
            brush_sph_init(&sph, x, y, z, s, 1);
            sph.col = rand() | 0xFF000000;
            oct_mod(&oct[octnum], (brush_t*)&sph, 1 + 2);
        }
        else if (((float)rand() / (float)RAND_MAX) > 0.8) {
            brush_box_init(&box, x, y, z, x+s, y+s, z+s, 1);
            box.col = rand() | 0xFF000000;
            oct_mod(&oct[octnum], (brush_t*)&box, 1 + 2);
        }
        else {
            brush_cone_init(&cone, x, y, z, 1, x, y+s, z, 5);
            cone.col = rand() | 0xFF000000;
            oct_mod(&oct[octnum], (brush_t*)&cone, 1 + 2);
        }
    }
#endif

    oct[octnum].tilid = 0;
    spr[sprnum].mixval = 0.5;
    spr[sprnum].imulcol = 0xff404040;
    spr[sprnum].oct = &oct[octnum];
    d = 1.0 / (double)oct[octnum].sid;
    spr[sprnum].p.x = -.5;
    spr[sprnum].p.y = -.5;
    spr[sprnum].p.z = -.5;
    spr[sprnum].r.x = d;
    spr[sprnum].r.y = 0;
    spr[sprnum].r.z = 0;
    spr[sprnum].d.x = 0;
    spr[sprnum].d.y = d;
    spr[sprnum].d.z = 0;
    spr[sprnum].f.x = 0;
    spr[sprnum].f.y = 0;
    spr[sprnum].f.z = d;
    sprnum++;
    octnum++;
    //-----------------------------------------------------------------------------------------------

    // Clip cursor to window
    RECT cursoldclip;
    POINT p;
    RECT r;
    GetClipCursor(&cursoldclip);
    p.x = (window.xres >> 1);
    p.y = (window.yres >> 1);
    ClientToScreen(window.ghwnd, &p);
    r.left = p.x;
    r.right = p.x + 1;
    r.top = p.y;
    r.bottom = p.y + 1;
    ClipCursor(&r);

    // Hide cursor
    ShowCursor(0);

    obstatus = 0;
    while (!window.breath())
    {
#if 0
        brush_sph_t sph;
        for (int i = 1; i >= 0; i--) {
            int x = rand() & (oct[0].sid - 1);
            int y = rand() & (oct[0].sid - 1);
            int z = rand() & (oct[0].sid - 1);
            int s = (rand() & ((oct[0].sid >> 4) - 1)) + (oct[0].sid >> 4);
            brush_sph_init(&sph, x, y, z, s, 1);
            oct_mod(&oct[0], (brush_t*)&sph, ((float)rand() / (float)RAND_MAX) > 0.1 ? 2 : 1 + 2);
        }
#elif 1
        brush_box_t box;
        for (int i = 1; i >= 0; i--) {
            int x = rand() & (oct[0].sid - 1);
            int y = rand() & (oct[0].sid - 1);
            int z = rand() & (oct[0].sid - 1);
            int s = (rand() & ((oct[0].sid >> 4) - 1)) + (oct[0].sid >> 4);
            brush_box_init(&box, x, y, z, x+s, y+s, z+s, 1);
            box.col = rand() | 0xFF000000;
            oct_mod(&oct[0], (brush_t*)&box, ((float)rand() / (float)RAND_MAX) > 0.1 ? 2 : 1 + 2);
        }
#endif

        //Read timer/mouse/keyboard
        otim = tim;
        tim = time_now_seconds();
        dtim = tim - otim;

        //Handle rotation
        if (!(window.bstatus & 2)) {
            rotvex(window.dmousx * 0.01, &ifor, &irig);
            f = irig.y * 0.1f;
        }
        else {
            f = (float)window.dmousx * -0.01f;
        }
        rotvex(window.dmousy * 0.01, &ifor, &idow);
        rotvex((double)f, &idow, &irig);
        window.dmousx = 0;
        window.dmousy = 0;
        window.dmousz = 0;

        //Handle movement (calculate movement vector from arrows&shifts)
        f = dtim;

        // Slow down.
        f *= 1.0/10.0;

        if (window.keystatus[0x2a]) f *= 1.0/16.0; //LShift
        if (window.keystatus[0x36]) f *= 16.0/1.0; //RShift
        fp.x = (window.keystatus[0xcd]-window.keystatus[0xcb])*f; //Right-Left
        fp.y = (window.keystatus[0x52]-window.keystatus[0x9d])*f; //KP0-RCtrl
        fp.z = (window.keystatus[0xc8]-window.keystatus[0xd0])*f; //Up-Down
        vec.x = irig.x*fp.x + idow.x*fp.y + ifor.x*fp.z;
        vec.y = irig.y*fp.x + idow.y*fp.y + ifor.y*fp.z;
        vec.z = irig.z*fp.x + idow.z*fp.y + ifor.z*fp.z;

        //Collision detection
        oct_world2voxpos(&spr[0],&ipos,&ipos);
        oct_world2voxdir(&spr[0],&vec,&vec);
        oct_slidemove(spr[0].oct,&ipos,&vec,0.5,&ipos);
        oct_vox2worldpos(&spr[0],&ipos,&ipos);

        surf_t cursurf;
        cursurf.b = rand()&255;
        cursurf.g = rand()&255;
        cursurf.r = rand()&255;
        cursurf.tex = rand()%108;

        while (int i = window.keyread())
        {
            switch((i>>8)&255)
            {
                case 0x01: window.quitloop(); break; //ESC
                case 0x0f: //Tab (grab surf)
                    oct_world2voxpos(&spr[0],&ipos,&fp);
                    oct_world2voxdir(&spr[0],&ifor,&vec); f = 4096.0; vec.x *= f; vec.y *= f; vec.z *= f;
                    isurf = oct_hitscan(spr[0].oct,&fp,&vec,&hit,&i,0);
                    if (isurf != -1) memcpy(&cursurf,SURFPTR(spr[0].oct,isurf),sizeof(surf_t));
                    break;
                case 0x39: //Space (paint surf)
                    oct_world2voxpos(&spr[0],&ipos,&fp);
                    oct_world2voxdir(&spr[0],&ifor,&vec); f = 4096.0; vec.x *= f; vec.y *= f; vec.z *= f;
                    isurf = oct_hitscan(spr[0].oct,&fp,&vec,&hit,&i,0);
                    if (isurf != -1) oct_writesurf(spr[0].oct,isurf,&cursurf);
                    break;
                case 0x13: //R (paint surf randomly)
                    oct_world2voxpos(&spr[0],&ipos,&fp);
                    oct_world2voxdir(&spr[0],&ifor,&vec); f = 4096.0; vec.x *= f; vec.y *= f; vec.z *= f;
                    isurf = oct_hitscan(spr[0].oct,&fp,&vec,&hit,&i,0);
                    if (isurf != -1)
                    {
                        cursurf.b = rand()&255;
                        cursurf.g = rand()&255;
                        cursurf.r = rand()&255;
                        cursurf.tex = rand()%108;
                        oct_writesurf(spr[0].oct,isurf,&cursurf);
                    }
                    break;
                case 0xd2: //Ins (insert voxels)
                    oct_world2voxpos(&spr[0],&ipos,&fp);
                    oct_world2voxdir(&spr[0],&ifor,&vec); f = 4096.0; vec.x *= f; vec.y *= f; vec.z *= f;
                    isurf = oct_hitscan(spr[0].oct,&fp,&vec,&hit,&i,0);
                    if (isurf != -1)
                    {
                        ((int *)&hit.x)[i>>1] += (i&1)*2-1;
                        int x, y, z;
                        if (oct_rebox(spr[0].oct,hit.x,hit.y,hit.z,hit.x+1,hit.y+1,hit.z+1,&x,&y,&z)) //grow octree if necessary
                        {
                            hit.x += x;
                            hit.y += y;
                            hit.z += z;
                            spr[0].p.x -= (spr[0].r.x*x + spr[0].d.x*y + spr[0].f.x*z);
                            spr[0].p.y -= (spr[0].r.y*x + spr[0].d.y*y + spr[0].f.y*z);
                            spr[0].p.z -= (spr[0].r.z*x + spr[0].d.z*y + spr[0].f.z*z);
                        }
                        oct_setvox(spr[0].oct,hit.x,hit.y,hit.z,SURFPTR(spr[0].oct,isurf),1+2);
                    }
                    break;
                case 0xd3: //Del (delete voxels)
                    oct_world2voxpos(&spr[0],&ipos,&fp);
                    oct_world2voxdir(&spr[0],&ifor,&vec);
                    f = 4096.0;
                    vec.x *= f;
                    vec.y *= f;
                    vec.z *= f;
                    isurf = oct_hitscan(spr[0].oct,&fp,&vec,&hit,&i,0);
                    if (isurf != -1) {
                        oct_setvox(spr[0].oct,hit.x,hit.y,hit.z,SURFPTR(spr[0].oct,isurf),2);
                    }
                    break;
            }
        }

        //Set window for drawing (handles both DirectDraw/OpenGL)
        if (oct_startdraw(window, &ipos,&irig,&idow,&ifor,ihaf.x,ihaf.y,ihaf.z) >= 0)
        {
                //Draw sprites
            for(int i=0;i<sprnum;i++) {
                oct_drawoct(spr[i].oct,&spr[i].p,&spr[i].r,&spr[i].d,&spr[i].f,spr[i].mixval,spr[i].imulcol);
            }

            oct_drawtext6x8(window.xres, window.yres, (window.xres>>1)-3,(window.yres>>1)-4,0xffffffff,0,"o");
            oct_drawline(window.xres, window.yres, window.xres-68, 0,window.xres-68,10,0xffc0c0c0);
            oct_drawline(window.xres, window.yres, window.xres-68,10,window.xres   ,10,0xffc0c0c0);
            avgdtim += (dtim-avgdtim)*.01;
            oct_drawtext6x8(window.xres, window.yres, window.xres-64,0,0xff00ffff,0xff0000ff,"%6.2f fps",1.0/avgdtim);

            oct_stopdraw(window);
        }
        obstatus = window.bstatus;
    }

    oct_freetex(gskyid);
    for(int i=octnum-1;i>=0;i--) {
        oct_free(&oct[i]);
    }
    oct_uninitonce(window.ghwnd);
    window.uninitapp();
    return(0);
}
