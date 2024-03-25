// This file has been modified from Ken Silverman's original release

#include "octree_physics.hpp"
#include "octree_modify.hpp"
#include "octree_renderer.hpp"

#include "brushes.hpp"
#include "types/parallelepiped.hpp"
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

static void orthofit (double *m, long c)
{
    double d, nm[9];
    int i;

        //(Cheap & simplified version of the 03/18/2006 algo)
        //Note: this version assumes input matrix has positive determinant
    for(;c>0;c--)
    {
        for(i=9-1;i>=0;i--) nm[i] = m[i];
        m[0] += nm[4]*nm[8] - nm[5]*nm[7];
        m[1] += nm[5]*nm[6] - nm[3]*nm[8];
        m[2] += nm[3]*nm[7] - nm[4]*nm[6];
        m[3] += nm[7]*nm[2] - nm[8]*nm[1];
        m[4] += nm[8]*nm[0] - nm[6]*nm[2];
        m[5] += nm[6]*nm[1] - nm[7]*nm[0];
        m[6] += nm[1]*nm[5] - nm[2]*nm[4];
        m[7] += nm[2]*nm[3] - nm[0]*nm[5];
        m[8] += nm[0]*nm[4] - nm[1]*nm[3];
        d = 1.0 / sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
        for(i=9-1;i>=0;i--) m[i] *= d;
    }
}

static void spr_update_rdf (spr_t *spr)
{
    spr->r.x = spr->br.x*spr->ori[0] + spr->br.y*spr->ori[1] + spr->br.z*spr->ori[2];
    spr->r.y = spr->br.x*spr->ori[3] + spr->br.y*spr->ori[4] + spr->br.z*spr->ori[5];
    spr->r.z = spr->br.x*spr->ori[6] + spr->br.y*spr->ori[7] + spr->br.z*spr->ori[8];
    spr->d.x = spr->bd.x*spr->ori[0] + spr->bd.y*spr->ori[1] + spr->bd.z*spr->ori[2];
    spr->d.y = spr->bd.x*spr->ori[3] + spr->bd.y*spr->ori[4] + spr->bd.z*spr->ori[5];
    spr->d.z = spr->bd.x*spr->ori[6] + spr->bd.y*spr->ori[7] + spr->bd.z*spr->ori[8];
    spr->f.x = spr->bf.x*spr->ori[0] + spr->bf.y*spr->ori[1] + spr->bf.z*spr->ori[2];
    spr->f.y = spr->bf.x*spr->ori[3] + spr->bf.y*spr->ori[4] + spr->bf.z*spr->ori[5];
    spr->f.z = spr->bf.x*spr->ori[6] + spr->bf.y*spr->ori[7] + spr->bf.z*spr->ori[8];
}

//10/26/2011:optimized algo :)
//WARNING:Assumes axis is unit length!
static void axisrotate (point3d *p, point3d *ax, double w)
{
    point3d op;
    double c, s, d;

    c = cos(w); s = sin(w);
    //P = cross(AX,P)*s + dot(AX,P)*(1-c)*AX + P*c;
    op = (*p);
    d = (op.x*ax->x + op.y*ax->y + op.z*ax->z)*(1.0-c);
    p->x = (ax->y*op.z - ax->z*op.y)*s + ax->x*d + op.x*c;
    p->y = (ax->z*op.x - ax->x*op.z)*s + ax->y*d + op.y*c;
    p->z = (ax->x*op.y - ax->y*op.x)*s + ax->z*d + op.z*c;
}

static void collide_binsearch (oct_t *o0, point3d *op, point3d *or_, point3d *od, point3d *of,
                                          point3d *np, point3d *nr, point3d *nd, point3d *nf,
                               oct_t *o1, point3d *p1, point3d *r1, point3d *d1, point3d *f1, oct_hit_t *ohit)
{
    point3d hitp, hitr, hitd, hitf;
    int k;

        //Binary search..
    for(k=8;k>0;k--)
    {
        hitp.x = (op->x + np->x)*.5;
        hitp.y = (op->y + np->y)*.5;
        hitp.z = (op->z + np->z)*.5;
        hitr = *or_;
        hitd = *od;
        hitf = *of; //FIX:hack!

        if (oct_touch_oct(o0,&hitp,&hitr,&hitd,&hitf,o1,p1,r1,d1,f1,ohit))
              { *np = hitp; *nr = hitr; *nd = hitd; *nf = hitf; }
        else { *op = hitp; *or_ = hitr; *od = hitd; *of = hitf; }
    }
}

//Binary search collision..
//  s0:     read: 1st object, starting posori, proposed delta movement
//  s1:     read: 2nd (static) object to test collision with
//  sfree: write: max posori not colliding
//  shit:  write: min posori     colliding (= {0} + 1 lsd binsrch unit)
//ohit: write: octree nodes, in order of {o0,o1}
//Returns ratio travelled {0.0-1.0}
static double collide_binsearch (spr_t *s0, spr_t *s1, spr_t *sfree, spr_t *shit, double dtim, oct_hit_t *ohit)
{
    spr_t tspr;
    point3d uax;
    double f, t, nt, tstep, mat[9], w;
    int j, k;

    w = s0->ax.x*s0->ax.x + s0->ax.y*s0->ax.y + s0->ax.z*s0->ax.z;
    if (w != 0.0)
    {
        w = sqrt(w);
        f = 1.0/w; uax.x = s0->ax.x*f; uax.y = s0->ax.y*f; uax.z = s0->ax.z*f;
    }

    t = 0.0; tstep = 1.0;
    for(k=8;k>0;k--,tstep*=.5)
    {
        nt = t+tstep;
        //---------------------------------------------------
        tspr.p.x = s0->p.x + s0->v.x*nt*dtim;
        tspr.p.y = s0->p.y + s0->v.y*nt*dtim;
        tspr.p.z = s0->p.z + s0->v.z*nt*dtim;

            //adjust pivot before rotation
        tspr.p.x += (s0->r.x*s0->cen.x + s0->d.x*s0->cen.y + s0->f.x*s0->cen.z);
        tspr.p.y += (s0->r.y*s0->cen.x + s0->d.y*s0->cen.y + s0->f.y*s0->cen.z);
        tspr.p.z += (s0->r.z*s0->cen.x + s0->d.z*s0->cen.y + s0->f.z*s0->cen.z);

        if (w != 0.0)
        {
            f = w*nt*dtim;
            for(j=0;j<9;j++) mat[j] = s0->ori[(j/3)+(j%3)*3];
            axisrotate((point3d *)&mat[0],&uax,f); //screw in
            axisrotate((point3d *)&mat[3],&uax,f);
            axisrotate((point3d *)&mat[6],&uax,f);
            for(j=0;j<9;j++) tspr.ori[j] = mat[(j/3)+(j%3)*3];
        }
        else
        {
            memcpy(tspr.ori,s0->ori,sizeof(tspr.ori));
        }
        tspr.br = s0->br; tspr.bd = s0->bd; tspr.bf = s0->bf;
        spr_update_rdf(&tspr);

            //adjust pivot after rotation
        tspr.p.x -= (tspr.r.x*s0->cen.x + tspr.d.x*s0->cen.y + tspr.f.x*s0->cen.z);
        tspr.p.y -= (tspr.r.y*s0->cen.x + tspr.d.y*s0->cen.y + tspr.f.y*s0->cen.z);
        tspr.p.z -= (tspr.r.z*s0->cen.x + tspr.d.z*s0->cen.y + tspr.f.z*s0->cen.z);
        //---------------------------------------------------
        j = oct_touch_oct(s0->oct,&tspr.p,&tspr.r,&tspr.d,&tspr.f,s1->oct,&s1->p,&s1->r,&s1->d,&s1->f,ohit);
        if (!j) { (*sfree) = tspr; t = nt; if (t == 1.0) break; }
            else { (*shit)  = tspr; }
    }

    if (t == 0.0)
    {
        (*sfree) = (*s0);
    }
    else
    {
        orthofit(sfree->ori,1); //not required, but doesn't hurt
        spr_update_rdf(sfree);
    }
    return(t);
}

//see rigid1.pdf&rigid2.pdf for derivation
    //---------------------------------------------------------
    //e      r   Coefficient of restitution: 0=plastic, 1=elastic
    //cp,cn  r   Contact point and unit vector nomal (provided by caller)
    //?pos   r   Centroid of body, world space
    //?ori   r   3x3 rotation matrix <rig,dow,for>, body to world space
    //?vel   rw  Translational velocity, world space
    //?rax   rw  Rotation axis scaled by angular velocity (rev/s), world space
    //?rmas  r   Reciprocal of mass
    //?rmoi  r   Inverse moment of inertia 3x3 matrix, ->BODY<- space.
static void doimpulse3d (double e, point3d *cp, point3d *cn,
    point3d *apos, double *aori, point3d *avel, point3d *arax, double armas, double *armoi,
    point3d *bpos, double *bori, point3d *bvel, point3d *brax, double brmas, double *brmoi)
{
    point3d fp, padot, pbdot, ra, rb, ta, tb;
    double fj, vrel, num, den, namoi[9], nbmoi[9];
        //FIXFIXFIXFIX:remove!
    printf("e:%g cp:<%g %g %g> cn:<%g %g %g>\n"
             "apos:<%g %g %g> aori:<%g %g %g %g %g %g %g %g %g> armas:%g armoi:<%g %g %g %g %g %g %g %g %g>\n"
             "bpos:<%g %g %g> bori:<%g %g %g %g %g %g %g %g %g> brmas:%g brmoi:<%g %g %g %g %g %g %g %g %g>\n",
        e,cp->x,cp->y,cp->z,cn->x,cn->y,cn->z,
        apos->x,apos->y,apos->z,aori[0],aori[1],aori[2],aori[3],aori[4],aori[5],aori[6],aori[7],aori[8],armas,armoi[0],armoi[1],armoi[2],armoi[3],armoi[4],armoi[5],armoi[6],armoi[7],armoi[8],
        bpos->x,bpos->y,bpos->z,bori[0],bori[1],bori[2],bori[3],bori[4],bori[5],bori[6],bori[7],bori[8],brmas,brmoi[0],brmoi[1],brmoi[2],brmoi[3],brmoi[4],brmoi[5],brmoi[6],brmoi[7],brmoi[8]);
    printf("avel:<%g %g %g> arax:<%g %g %g> bvel:<%g %g %g> brax:<%g %g %g>\n",avel->x,avel->y,avel->z,arax->x,arax->y,arax->z,bvel->x,bvel->y,bvel->z,brax->x,brax->y,brax->z);
        //Calc moi in world coords
        //I_world^-1 = R * I_body^-1 * R^T
    simxform(namoi,aori,armoi);
    simxform(nbmoi,bori,brmoi);
        //Calc impulse magnitude (j)
        //j = (-e-1) * (((Va + Wa x Rap) - (Vb + Wb x Rbp)) dot N) /
        //   ( 1/ma + N dot ((Ia^-1 * (Rap cross N)) cross Rap)
        //   + 1/mb + N dot ((Ib^-1 * (Rbp cross N)) cross Rbp))
        //
        //j = num/(armas + brmas + cn . ((a->Iinv * (ra x cn)) x ra)
        //                         cn . ((b->Iinv * (rb x cn)) x rb));
        //Calc relative velocity (vrel), parallel to contact normal (cn) at contact point (cp)
    vecsub(&ra,cp,apos); //ra = cp - a->x;
    vecsub(&rb,cp,bpos); //rb = cp - b->x;
    veccross(&fp,arax,&ra); vecadd(&padot,avel,&fp); //padot = a->v + (a->w x ra)
    veccross(&fp,brax,&rb); vecadd(&pbdot,bvel,&fp); //pbdot = b->v + (b->w x rb)
    vecsub(&fp,&padot,&pbdot); vrel = vecdot(cn,&fp); //vrel = cn . (padot - pbdot);
    num = (-e-1.0)*vrel;
    den = armas+brmas;
    veccross(&fp,&ra,cn); matvecmul(&ta,namoi,&fp); veccross(&fp,&ta,&ra); den += vecdot(cn,&fp);
    veccross(&fp,&rb,cn); matvecmul(&tb,nbmoi,&fp); veccross(&fp,&tb,&rb); den += vecdot(cn,&fp);
    fj = num/den;
        //Do impulse:translation
    vecscale(&fp,cn,fj*armas); vecadd(avel,avel,&fp); //nVa = Va + (j*N)/ma;
    vecscale(&fp,cn,fj*brmas); vecsub(bvel,bvel,&fp); //nVb = Vb - (j*N)/mb;
        //Do impulse:rotation
    vecscale(&fp,&ta,fj); vecadd(arax,arax,&fp);       //nWa = Wa + Ia^-1 * (Rap cross (j*N));
    vecscale(&fp,&tb,fj); vecsub(brax,brax,&fp);       //nWb = Wb - Ib^-1 * (Rbp cross (j*N));
        //FIXFIXFIXFIX:remove!
    printf("avel:<%g %g %g> arax:<%g %g %g> bvel:<%g %g %g> brax:<%g %g %g>\n",avel->x,avel->y,avel->z,arax->x,arax->y,arax->z,bvel->x,bvel->y,bvel->z,brax->x,brax->y,brax->z);
}

static void sprite_del (int k)
{
    int i, j;

    sprnum--; spr[k] = spr[sprnum];
    for(i=octnum-1;i>=0;i--)
    {
        for(j=sprnum-1;j>=0;j--) if (spr[j].oct == &oct[i]) break;
        if (j >= 0) continue;
        oct_free(&oct[i]);
        octnum--; oct[i] = oct[octnum];
        for(j=sprnum-1;j>=0;j--) if (spr[j].oct == &oct[octnum]) spr[j].oct = &oct[i];
    }
}

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
    static oct_hit_t ohit[2];

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

            for(int i=sprnum-1;i>0;i--)
            {
                point3d lvec;
                double dt, bdt, mat[9];
                spr_t sfree, shit, bsfree, bshit;
                oct_hit_t bohit[2];
                int bj;

                if (spr[i].cnt < 3) spr[i].v.y += dtim*.25; //gravity

                oct_drawtext6x8(window.xres, window.yres, window.xres>>1,30+i*10,0xffffffff,0,"spr[%d].v: %f %f %f",i,spr[i].v.x,spr[i].v.y,spr[i].v.z); //FIXFIXFIXFIX:remove

                bdt = 1.0;
                for(int j=sprnum-1;j>=0;j--) ////FIX:doesn't handle hitting multiple sprites simultaneously properly :/
                {
                    if (i == j) continue;
                    dt = collide_binsearch(&spr[i],&spr[j],&sfree,&shit,dtim,ohit);
                    if (dt <= bdt) { bdt = dt; bj = j; bsfree = sfree; bshit = shit; bohit[0] = ohit[0]; bohit[1] = ohit[1]; }
                }

                    //precession (AngMom = Iworld*AngVel); see FIZEX3D.KC for derivation
                simxform(mat,spr[i].ori,spr[i].moi); matvecmul(&lvec,mat,&spr[i].ax);

                spr[i].p = bsfree.p; spr[i].r = bsfree.r; spr[i].d = bsfree.d; spr[i].f = bsfree.f;
                memcpy(&spr[i].ori,&bsfree.ori,sizeof(spr[i].ori));
                spr_update_rdf(&spr[i]);

                    //precession (AngVel = Iworld^-1*AngMom); see FIZEX3D.KC for derivation
                simxform(mat,spr[i].ori,spr[i].rmoi); matvecmul(&spr[i].ax,mat,&lvec);

                if (bdt < 1.0)
                {
                    pgram3d_t bi;
                    point3d dp0, dr0, dd0, df0, dp1, dr1, dd1, df1, dpos, hit, norm;

                    //if (bdt == 0.0) { spr[i].v.x = 0; spr[i].v.y = 0.0; spr[i].v.z = 0.0; } //FIXFIXFIXFIX:prevents daffy ducking, but makes more stucks :/
                    if (bdt == 0.0) spr[i].cnt++; else spr[i].cnt = 0;

                    f = (1<<bohit[0].ls)*.5;
                    dr0.x = bshit.r.x*f; dd0.x = bshit.d.x*f; df0.x = bshit.f.x*f;
                    dr0.y = bshit.r.y*f; dd0.y = bshit.d.y*f; df0.y = bshit.f.y*f;
                    dr0.z = bshit.r.z*f; dd0.z = bshit.d.z*f; df0.z = bshit.f.z*f;
                    fp.x = bohit[0].x+f; fp.y = bohit[0].y+f; fp.z = bohit[0].z+f; oct_vox2worldpos(&bshit,&fp,&dp0);
                    f = (1<<bohit[1].ls)*.5;
                    dr1.x = spr[bj].r.x*f; dd1.x = spr[bj].d.x*f; df1.x = spr[bj].f.x*f;
                    dr1.y = spr[bj].r.y*f; dd1.y = spr[bj].d.y*f; df1.y = spr[bj].f.y*f;
                    dr1.z = spr[bj].r.z*f; dd1.z = spr[bj].d.z*f; df1.z = spr[bj].f.z*f;
                    fp.x = bohit[1].x+f; fp.y = bohit[1].y+f; fp.z = bohit[1].z+f; oct_vox2worldpos(&spr[bj],&fp,&dp1);

                    pgram3d_init(&bi,&dr0,&dd0,&df0,&dr1,&dd1,&df1);
                    if (pgram3d_gethit(&bi,&dp0,&dp1,&dpos,&hit,&norm))
                    {
                        if (!bj) { spr[0].rmas = 0.0; for(int k=0;k<9;k++) spr[0].rmoi[k] = 0.0; } //Hack to prevent spr[0] from moving
                        //if (keystatus[0x38]) { gravity = 0; fixhitit = 1; curspr = 1; hitp = bshit.p; hitr = bshit.r; hitd = bshit.d; hitf = bshit.f; break; } //FIX:debug only!

                        oct_vox2worldpos(&spr[ i],&spr[ i].cen,&dp0);
                        oct_vox2worldpos(&spr[bj],&spr[bj].cen,&dp1);
                        doimpulse3d(0.5,&hit,&norm,&dp0,spr[ i].ori,&spr[ i].v,&spr[ i].ax,spr[ i].rmas,spr[ i].rmoi,
                                                            &dp1,spr[bj].ori,&spr[bj].v,&spr[bj].ax,spr[bj].rmas,spr[bj].rmoi);
                    }

                    spr[i].tim += dtim;
                    if (spr[i].tim < 1.0)
                        { if (spr[i].tim >= 0.0) spr[i].imulcol = (spr[i].imulcol&0xffffff) + (std::min(std::max((int)(255.0-spr[i].tim*256.0),0),255)<<24); }
                    else {
                        sprite_del(i);
                    }
                }
                else
                {
                    spr[i].cnt = 0;
                    spr[i].tim = std::max(spr[i].tim-dtim*.5,-2.0);
                }
                if (spr[i].p.y > 4.0) { sprite_del(i); }
            }

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
