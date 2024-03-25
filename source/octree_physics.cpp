// This file has been modified from Ken Silverman's original release

#include "octree_physics.hpp"

// TODO: REMOVE!
#include "octree_renderer.hpp"
// TODO: REMOVE!

#include "types/parallelepiped.hpp"

#include "utilities/gl.hpp"
#include "utilities/macros.hpp"
#include "utilities/maths.hpp"
#include "utilities/popcount.hpp"
#include "utilities/thread_pool.hpp"

#include <algorithm>
#include <cstring>

//--------------------------------------------------------------------------------------------------
//              CPU:  GPU:
//pnd3d /thr=1: 290ms 312ms
//voxed /thr=1: 192ms
//pnd3d /thr=2: 145ms 158ms
//pnd3d /thr=3: 133ms 143ms
//pnd3d /thr=4:  91ms 109ms
//pnd3d /thr=5:  87ms 105ms
//pnd3d /thr=6:  67ms  85ms
//pnd3d /thr=7:  65ms  xxms
//pnd3d /thr=8:  62ms  xxms

//new algo (surf_t following)

void oct_bitcac_reset (oct_bitcac_t *bitcac) {
    #if (USENEWCAC == 0)
        bitcac->pt.x = 0x80000000;
    #else
        bitcac->n = 0;
        memset(bitcac->hashead,-1,CACHASHN*sizeof(int));
    #endif
}

void oct_estnorm (oct_t *loct, int nx, int ny, int nz, int rad, point3f *norm, oct_bitcac_t *bitcac)
{
    float f;
    int i, j, x, y, z, dx, dy, dz;

#if 1
    //Voxlap algo (faster)
    static const int bitsnum[32] = {
        0        ,1-(2<<16),1-(1<<16),2-(3<<16),
        1        ,2-(2<<16),2-(1<<16),3-(3<<16),
        1+(1<<16),2-(1<<16),2        ,3-(2<<16),
        2+(1<<16),3-(1<<16),3        ,4-(2<<16),
        1+(2<<16),2        ,2+(1<<16),3-(1<<16),
        2+(2<<16),3        ,3+(1<<16),4-(1<<16),
        2+(3<<16),3+(1<<16),3+(2<<16),4,
        3+(3<<16),4+(1<<16),4+(2<<16),5
    };
    unsigned int bitsol[32*32], *bitsolp, *bitsolp2;
    int xx, yy, zz, dia, b[5];

    if (bitcac)
    {
        dia = 32; i = 32-rad*2;

        #if (USENEWCAC == 0)
            if (((unsigned)(nx-rad-bitcac->pt.x) >= (unsigned)i) ||
                 ((unsigned)(ny-rad-bitcac->pt.y) >= (unsigned)i) ||
                 ((unsigned)(nz-rad-bitcac->pt.z) >= (unsigned)i))
            {
                bitcac->pt.x = (nx-rad-6)&-8;
                bitcac->pt.y = (ny-rad-6)&-8;
                bitcac->pt.z = (nz-rad-6)&-8;
                oct_sol2bit(loct,bitcac->buf,bitcac->pt.x,bitcac->pt.y,bitcac->pt.z,32,32,32,1);
            }
            x = nx-bitcac->pt.x; y = ny-bitcac->pt.y; z = nz-bitcac->pt.z; bitsolp = bitcac->buf;
        #else
            if ((!bitcac->n) ||
                 ((unsigned)(nx-rad-bitcac->pt[bitcac->cur].x) >= (unsigned)i) ||
                 ((unsigned)(ny-rad-bitcac->pt[bitcac->cur].y) >= (unsigned)i) ||
                 ((unsigned)(nz-rad-bitcac->pt[bitcac->cur].z) >= (unsigned)i))
            {
                x = (nx-rad-6)&-8;
                y = (ny-rad-6)&-8;
                z = (nz-rad-6)&-8;

                hash = ((x>>3) + (y>>3)*3 + (z>>3)*5)&(CACHASHN-1);
                bitcac->cur = -1;
                for(i=bitcac->hashead[hash];i>=std::max(bitcac->n-CACN,0);i=bitcac->ptn[i&(CACN-1)])
                    if ((bitcac->pt[i].x == x) && (bitcac->pt[i].y == y) && (bitcac->pt[i].z == z)) { bitcac->cur = (i&(CACN-1)); break; }
                if (bitcac->cur < 0)
                {
                    bitcac->cur = (bitcac->n&(CACN-1));
                    bitcac->ptn[bitcac->cur] = bitcac->hashead[hash]; bitcac->hashead[hash] = bitcac->n; bitcac->n++;
                    bitcac->pt[bitcac->cur].x = x;
                    bitcac->pt[bitcac->cur].y = y;
                    bitcac->pt[bitcac->cur].z = z;
                    oct_sol2bit(loct,bitcac->buf[bitcac->cur],x,y,z,32,32,32,1);
                }
            }
            x = nx-bitcac->pt[bitcac->cur].x;
            y = ny-bitcac->pt[bitcac->cur].y;
            z = nz-bitcac->pt[bitcac->cur].z; bitsolp = bitcac->buf[bitcac->cur];
        #endif
    }
    else {
        dia = rad*2+1;
        oct_sol2bit(loct,(unsigned int *)bitsol,nx-rad,ny-rad,nz-rad,dia,dia,dia,1);
        x = rad; y = rad; z = rad; bitsolp = bitsol;
    }

    dx = 0; dy = 0; dz = 0;
    if (rad == 2)
    {
        bitsolp2 = &bitsolp[(z-2)*dia + y-2]; i = x-2;
        for(yy=-2;yy<=2;yy++)
        {
            b[0] = ((bitsolp2[0]>>i)&31);
            b[1] = ((bitsolp2[1]>>i)&31);
            b[2] = ((bitsolp2[2]>>i)&31);
            b[3] = ((bitsolp2[3]>>i)&31);
            b[4] = ((bitsolp2[4]>>i)&31); bitsolp2 += dia;

            dy += (popcount8(b[4])-popcount8(b[0]))*2 + popcount8(b[3])-popcount8(b[1]);
            j = bitsnum[b[0]] + bitsnum[b[1]] + bitsnum[b[2]] + bitsnum[b[3]] + bitsnum[b[4]];
            dx += j; dz += (*(signed short *)&j)*yy;
        }
        dx >>= 16;
    }
    else
    {
        for(xx=-2;xx<=2;xx++)
            for(yy=-2;yy<=2;yy++)
                for(zz=-2;zz<=2;zz++)
                    if (bitsolp[(z+zz)*dia + (y+yy)]&(1<<(x+xx))) { dx += xx; dy += yy; dz += zz; }
    }
#else
    //New algo: search surfaces (bfs) - not faster, but better quality!
    //be careful: rad must be <= 8!
    //12 edges + 6 faces = 18
    static const int dirx[18] = { 0, 0,-1,+1, 0, 0,   0,-1,+1, 0, -1,+1,-1,+1,  0,-1,+1, 0};
    static const int diry[18] = { 0,-1, 0, 0,+1, 0,  -1, 0, 0,+1, -1,-1,+1,+1, -1, 0, 0,+1};
    static const int dirz[18] = {-1, 0, 0, 0, 0,+1,  -1,-1,-1,-1,  0, 0, 0, 0, +1,+1,+1,+1};

    #define ESFIFMAX 4096
    typedef struct { int x, y, z, d; } fif_t;
    fif_t fif[ESFIFMAX];
    int xx, yy, zz, dia, fifw, fifr;
    unsigned int bitsol[32*32], *bitsolp, *bitsolp2, bitgot[32*32];
    static const int distweight[] = //must have >= rad+1 entries defined
    {
        (int)(exp(0.0*0.0*-.1)*65536.0), //65536
        (int)(exp(1.0*1.0*-.1)*65536.0), //59299
        (int)(exp(2.0*2.0*-.1)*65536.0), //43930
        (int)(exp(3.0*3.0*-.1)*65536.0), //26644
        (int)(exp(4.0*4.0*-.1)*65536.0), //13231
        (int)(exp(5.0*5.0*-.1)*65536.0), // 5379
        (int)(exp(6.0*6.0*-.1)*65536.0), // 1790
        (int)(exp(7.0*7.0*-.1)*65536.0), //  488
        (int)(exp(8.0*8.0*-.1)*65536.0), //  108
    };

        //       face:       |       edge:      |  edge#2:
        //  pass:     fail:  |  pass:    fail:  |   fail:
        //                   |     X        .   |      X
        //   . .       . X   |   . 1 X    X 1 . |    . 1 X
        // X 0 1 X   X 0 1 . | X 0 X    . 0 X   |  X 0 .
        //   X X       . X   |   X        .     |    X

        //  6    10  1 11     14
        //7 0 8   2     3  15  5 16
        //  9    12  4 13     17
        //
        //0 1 2   9 10 11  18 19 20
        //3 4 5  12 13 14  21 22 23
        //6 7 8  15 16 17  24 25 26
    static const int neighmsk[18][3] =
    {
            //6 faces
        (1<< 4),9,(1<< 1)+(1<< 3)+(1<< 5)+(1<< 7),
        (1<<10),3,(1<< 1)+(1<< 9)+(1<<11)+(1<<19),
        (1<<12),1,(1<< 3)+(1<< 9)+(1<<15)+(1<<21),
        (1<<14),1,(1<< 4)+(1<<10)+(1<<16)+(1<<22),
        (1<<16),3,(1<< 4)+(1<<12)+(1<<14)+(1<<22),
        (1<<22),9,(1<<10)+(1<<12)+(1<<14)+(1<<16),

            //12 edges
        (1<< 1),(1<< 4),(1<<10),
        (1<< 3),(1<< 4),(1<<12),
        (1<< 5),(1<< 4),(1<<14),
        (1<< 7),(1<< 4),(1<<16),
        (1<< 9),(1<<10),(1<<12),
        (1<<11),(1<<10),(1<<14),
        (1<<15),(1<<12),(1<<16),
        (1<<17),(1<<14),(1<<16),
        (1<<19),(1<<10),(1<<22),
        (1<<21),(1<<12),(1<<22),
        (1<<23),(1<<14),(1<<22),
        (1<<25),(1<<16),(1<<22),
    };

    if (bitcac)
    {
        dia = 32;
        if ((nx-rad-1 < bitcac->pt.x) || (nx+rad+1 >= bitcac->pt.x+32) ||
             (ny-rad-1 < bitcac->pt.y) || (ny+rad+1 >= bitcac->pt.y+32) ||
             (nz-rad-1 < bitcac->pt.z) || (nz+rad+1 >= bitcac->pt.z+32))
        {
            bitcac->pt.x = ((nx-rad-1)-5)&-8;
            bitcac->pt.y = ((ny-rad-1)-5)&-8;
            bitcac->pt.z = ((nz-rad-1)-5)&-8;
            oct_sol2bit(loct,bitcac->buf,bitcac->pt.x,bitcac->pt.y,bitcac->pt.z,32,32,32,1);
        }

        x = (nx-rad-1)-bitcac->pt.x;
        y = (ny-rad-1)-bitcac->pt.y;
        z = (nz-rad-1)-bitcac->pt.z;

        memset(&bitgot[z*32],0,(rad*2+3)*32*sizeof(bitgot[0]));
        x += rad+1; y += rad+1; z += rad+1; bitsolp = bitcac->buf;
    }
    else
    {
        dia = rad*2+3;
        oct_sol2bit(loct,(unsigned int *)bitsol,nx-rad-1,ny-rad-1,nz-rad-1,dia,dia,dia,1);
        memset(bitgot,0,dia*dia*sizeof(bitgot[0]));
        x = rad+1; y = rad+1; z = rad+1; bitsolp = bitsol;
    }

    dx = 0; dy = 0; dz = 0;
    bitgot[z*dia+y] |= (1<<x);
    fif[0].x = x; fif[0].y = y; fif[0].z = z; fif[0].d = 0; fifw = 1;
    for(fifr=0;fifr<fifw;fifr++)
    {
        x = fif[fifr].x; y = fif[fifr].y; z = fif[fifr].z;

            //0 1 2   9 10 11  18 19 20
            //3 4 5  12 13 14  21 22 23
            //6 7 8  15 16 17  24 25 26
        bitsolp2 = &bitsolp[z*dia + y]; i = x-1;
        j = (((bitsolp2[-dia-1]>>i)&7)    ) + (((bitsolp2[-dia+0]>>i)&7)<< 3) + (((bitsolp2[-dia+1]>>i)&7)<< 6) +
             (((bitsolp2[    -1]>>i)&7)<< 9) + (((bitsolp2[    +0]>>i)&7)<<12) + (((bitsolp2[    +1]>>i)&7)<<15) +
             (((bitsolp2[+dia-1]>>i)&7)<<18) + (((bitsolp2[+dia+0]>>i)&7)<<21) + (((bitsolp2[+dia+1]>>i)&7)<<24);

        i = distweight[fif[fifr].d]; //give surfaces near nucleus higher weight (solves some issues in addition to looking nice)
        if (j&(1<<12)) dx -= i;
        if (j&(1<<14)) dx += i;
        if (j&(1<<10)) dy -= i;
        if (j&(1<<16)) dy += i;
        if (j&(1<< 4)) dz -= i;
        if (j&(1<<22)) dz += i;

        if (fif[fifr].d >= rad) continue;

        for(i=18-1;i>=0;i--)
        {
            if (!(j&neighmsk[i][0])) continue;
            if (i < 6) { if ((((j>>neighmsk[i][1])|j)&neighmsk[i][2]) == neighmsk[i][2]) continue; }
                    else { if (!(j&neighmsk[i][1]) == !(j&neighmsk[i][2])) continue; }

            xx = dirx[i]+x; yy = diry[i]+y; zz = dirz[i]+z; if (bitgot[zz*dia+yy]&(1<<xx)) continue;
            bitgot[zz*dia+yy] |= (1<<xx);

            fif[fifw].x = xx; fif[fifw].y = yy; fif[fifw].z = zz; fif[fifw].d = fif[fifr].d+1; fifw++;
        }
    }
#endif
    if (!(dx|dy|dz)) { norm->x = 0.0; norm->y = 0.0; norm->z = 0.0; return; }

    f = (float)dx*dx + (float)dy*dy + (float)dz*dz;

    #if !defined(_WIN32)
        f = 1.0/sqrt(f);
    #else
        _asm
        {
            rsqrtss xmm0, f
            movss f, xmm0
        }
    #endif

    norm->x = (float)dx*f;
    norm->y = (float)dy*f;
    norm->z = (float)dz*f;
}

void oct_estnorm (oct_t *loct, int nx, int ny, int nz, int rad, point3d *dnorm, oct_bitcac_t *bitcac)
{
    point3f norm;
    oct_estnorm(loct,nx,ny,nz,rad,&norm,bitcac);
    dnorm->x = norm.x; dnorm->y = norm.y; dnorm->z = norm.z;
}

typedef struct { oct_t *loct; int dist, x0, y0, z0, x1, y1, z1; } oct_rn_t;
static void oct_refreshnorms_mt (int i, void *_)
{
    oct_rn_t *rn;
    typedef struct { octv_t *ptr; int x, y, z, j; } stk_t;
    stk_t stk[OCT_MAXLS];
    point3f norm;
    surf_t *psurf;
    octv_t *ptr;
    int j, x, y, z, nx, ny, nz, ls, s;
    oct_bitcac_t bitcac;

    rn = &((oct_rn_t *)_)[i];

    oct_bitcac_reset(&bitcac);

    ls = rn->loct->lsid-1; s = (1<<ls); ptr = &((octv_t *)rn->loct->nod.buf)[rn->loct->head]; x = 0; y = 0; z = 0; j = 0/*Don't reverse direction - would mess up bitcac inside oct_estnorm()*/;
    while (1)
    {
        i = (1<<j); if (!(ptr->chi&i)) goto tosibly;
        nx = (((j   )&1)<<ls)+x; if ((nx >= rn->x1) || (nx+s <= rn->x0)) goto tosibly;
        ny = (((j>>1)&1)<<ls)+y; if ((ny >= rn->y1) || (ny+s <= rn->y0)) goto tosibly;
        nz = (((j>>2)&1)<<ls)+z; if ((nz >= rn->z1) || (nz+s <= rn->z0)) goto tosibly;

        i = popcount8(ptr->chi&(i-1)) + ptr->ind;

        if (ls <= 0)
        {
            psurf = &((surf_t *)rn->loct->sur.buf)[i];
            oct_estnorm(rn->loct,nx,ny,nz,rn->dist,&norm,&bitcac);
            psurf->norm[0] = (signed char)(norm.x*127.0);
            psurf->norm[1] = (signed char)(norm.y*127.0);
            psurf->norm[2] = (signed char)(norm.z*127.0);

            if constexpr (oct_usegpu && oct_usegpubo) {
                memcpy(&rn->loct->gsurf[i],psurf,rn->loct->sur.siz);
            }
            goto tosibly;
        }

        stk[ls].ptr = ptr; stk[ls].x = x; stk[ls].y = y; stk[ls].z = z; stk[ls].j = j; ls--; s >>= 1; //2child
        ptr = &((octv_t *)rn->loct->nod.buf)[i]; x = nx; y = ny; z = nz; j = 0;
        continue;

tosibly:
        j++; if (j < 8) continue;
        do { ls++; s <<= 1; if (ls >= rn->loct->lsid) return; j = stk[ls].j+1; } while (j >= 8); //2parent
        ptr = stk[ls].ptr; x = stk[ls].x; y = stk[ls].y; z = stk[ls].z;
    }
}

void oct_refreshnorms (oct_t *loct, int dist, int x0, int y0, int z0, int x1, int y1, int z1)
{
    #define LXSIZ 3
    #define LYSIZ 3
    #define LZSIZ 3
    static oct_rn_t octrn[1<<(LXSIZ+LYSIZ+LZSIZ)];
    int i, x, y, z;

    if constexpr (oct_usegpu && oct_usegpubo) {
        if (!loct->gsurf) {
            loct->gsurf = (surf_t *)bo_begin(loct->bufid,0);
        }
    }

    i = 0;
    for(x=(1<<LXSIZ)-1;x>=0;x--) {
        for(y=(1<<LYSIZ)-1;y>=0;y--) {
            for(z=(1<<LZSIZ)-1;z>=0;z--,i++)
            {
                octrn[i].loct = loct; octrn[i].dist = dist;
                octrn[i].x0 = (((x1-x0)*(x+0))>>LXSIZ) + x0;
                octrn[i].y0 = (((y1-y0)*(y+0))>>LYSIZ) + y0;
                octrn[i].z0 = (((z1-z0)*(z+0))>>LZSIZ) + z0;
                octrn[i].x1 = (((x1-x0)*(x+1))>>LXSIZ) + x0;
                octrn[i].y1 = (((y1-y0)*(y+1))>>LYSIZ) + y0;
                octrn[i].z1 = (((z1-z0)*(z+1))>>LZSIZ) + z0;
            }
        }
    }
    htrun(oct_refreshnorms_mt,octrn,0,1<<(LXSIZ+LYSIZ+LZSIZ),true);

    if constexpr (oct_usegpu && oct_usegpubo) {
        ((PFNGLBINDTEXTURE)glfp[glBindTexture])(GL_TEXTURE_2D,loct->octid);
        ((PFNGLTEXSUBIMAGE2D)glfp[glTexSubImage2D])(GL_TEXTURE_2D,0,0,0,loct->gxsid,(loct->sur.mal*2)>>loct->glxsid,GL_RGBA,GL_UNSIGNED_BYTE,(void *)loct->sur.buf);
    }
}

//NOTE:x0,etc. must be filled with bounded box to check on input. For whole object, set x0=0; x1=loct->sid; etc..
//NOTE:ins&outs are inclusive on x0/y0/z0 and exclusive on x1/y1/z1
int oct_getsolbbox (oct_t *loct, int *rx0, int *ry0, int *rz0, int *rx1, int *ry1, int *rz1)
{
    typedef struct { octv_t *ptr; int j; } stk_t;
    stk_t stk[OCT_MAXLS];
    octv_t *ptr;
    int i, j, ls, s, x, y, z, x0, y0, z0, x1, y1, z1;

    x0 = loct->sid; x1 = -1;
    y0 = loct->sid; y1 = -1;
    z0 = loct->sid; z1 = -1;
    ls = loct->lsid-1; s = (1<<ls); ptr = &((octv_t *)loct->nod.buf)[loct->head]; x = 0; y = 0; z = 0; j = 8-1;
    while (1)
    {
        i = (1<<j); if (!((ptr->chi|ptr->sol)&i)) goto tosibly;

        x = ((-((j   )&1))&s) + (x&(-s-s));
        y = ((-((j>>1)&1))&s) + (y&(-s-s));
        z = ((-((j>>2)&1))&s) + (z&(-s-s));

        if ((x >= (*rx1)) || (x+s <= (*rx0)) || //skip if fully outside user input bbox
             (y >= (*ry1)) || (y+s <= (*ry0)) ||
             (z >= (*rz1)) || (z+s <= (*rz0))) goto tosibly;
        if ((x >= x0) && (x+s <= x1) && //skip if cell fully inside current extents
             (y >= y0) && (y+s <= y1) &&
             (z >= z0) && (z+s <= z1)) goto tosibly;

        if (ptr->sol&i)
        {
            x0 = std::min(x0,x); x1 = std::max(x1,x+s);
            y0 = std::min(y0,y); y1 = std::max(y1,y+s);
            z0 = std::min(z0,z); z1 = std::max(z1,z+s);
            goto tosibly;
        }

        stk[ls].ptr = ptr; stk[ls].j = j; s >>= 1; ls--; //2child
        ptr = &((octv_t *)loct->nod.buf)[popcount8(ptr->chi&(i-1)) + ptr->ind]; j = 8-1;
        continue;

tosibly:;
        j--; if (j >= 0) continue;
        do { ls++; s <<= 1; if (ls >= loct->lsid) goto endit; j = stk[ls].j-1; } while (j < 0); //2parent
        ptr = stk[ls].ptr;
    }
endit:
    if (x0 >= x1) return(0);
    (*rx0) = std::max(*rx0,x0); (*rx1) = std::min(*rx1,x1);
    (*ry0) = std::max(*ry0,y0); (*ry1) = std::min(*ry1,y1);
    (*rz0) = std::max(*rz0,z0); (*rz1) = std::min(*rz1,z1);
    return(1);
}

//--------------------------------------------------------------------------------------------------

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
void doimpulse3d (double e, point3d *cp, point3d *cn,
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

//+-------+-------+
//|SOL    |    CHI|
//|       |       |
//|       |/^\    |
//|      /|air \  |
//+----/^\+----/^\+
//|  / chi|\ /     \
//|  \    |/^\ chi|/
//|    \ /|sol \ /|
//|CHI   \|    AIR|
//+-------+\-/----+
struct oto_t
{
    oct_t *loct[2];
    pgram3d_t bi[OCT_MAXLS*2+1];
    point3f p, r, d, f;
    oct_hit_t hit[2];
};

static float gpow2[OCT_MAXLS*2+1] = {};

static int oct_touch_oct_recur (oto_t *oto, octv_t *inode0, int x0, int y0, int z0, int ls0, octv_t *inode1, int x1, int y1, int z1, int ls1, int splitmsk)
{
    double fx, fy, fz;
    int i, j, k, v, x, y, z, s0, s1, ind;

    i = (ls0 > 0) && (splitmsk&1);
    j = (ls1 > 0) && (splitmsk&2);

    if ((i) && (j) && (ls0 < ls1)) i = 0; //split larger box when there's choice (absolutely necessary optimization to reduce recursion count!)

    if (i) //split 0
    {
        ls0--; s0 = (1<<ls0); s1 = (1<<ls1); ind = inode0->ind;
        fx = s0*.5 - (x1+s1*.5)*oto->r.x - (y1+s1*.5)*oto->d.x - (z1+s1*.5)*oto->f.x - oto->p.x;
        fy = s0*.5 - (x1+s1*.5)*oto->r.y - (y1+s1*.5)*oto->d.y - (z1+s1*.5)*oto->f.y - oto->p.y;
        fz = s0*.5 - (x1+s1*.5)*oto->r.z - (y1+s1*.5)*oto->d.z - (z1+s1*.5)*oto->f.z - oto->p.z;
        for(v=(inode0->sol|inode0->chi);v;v^=i)
        {
            i = (-v)&v;
            x = x0; if (i&0xaa) x += s0;
            y = y0; if (i&0xcc) y += s0;
            z = z0; if (i&0xf0) z += s0;

            j = pgram3d_isint(&oto->bi[ls1-ls0+OCT_MAXLS],(x+fx)*gpow2[OCT_MAXLS-ls0]*4.0,
                                                                         (y+fy)*gpow2[OCT_MAXLS-ls0]*4.0,
                                                                         (z+fz)*gpow2[OCT_MAXLS-ls0]*4.0);
            if (!j) continue;

            k = splitmsk; if (inode0->sol&i) k &= ~1;
            if (oct_touch_oct_recur(oto,&((octv_t *)oto->loct[0]->nod.buf)[popcount8(inode0->chi&(i-1))+ind],x,y,z,ls0,inode1,x1,y1,z1,ls1,k)) return(1);
        }
        return(0);
    }

    if (j) //split 1
    {
        ls1--; s0 = (1<<ls0); s1 = (1<<ls1); ind = inode1->ind;
        fx = x0+s0*.5 - (s1*.5)*oto->r.x - (s1*.5)*oto->d.x - (s1*.5)*oto->f.x - oto->p.x;
        fy = y0+s0*.5 - (s1*.5)*oto->r.y - (s1*.5)*oto->d.y - (s1*.5)*oto->f.y - oto->p.y;
        fz = z0+s0*.5 - (s1*.5)*oto->r.z - (s1*.5)*oto->d.z - (s1*.5)*oto->f.z - oto->p.z;
        for(v=(inode1->sol|inode1->chi);v;v^=i)
        {
            i = (-v)&v;
            x = x1; if (i&0xaa) x += s1;
            y = y1; if (i&0xcc) y += s1;
            z = z1; if (i&0xf0) z += s1;

            j = pgram3d_isint(&oto->bi[ls1-ls0+OCT_MAXLS],(fx - x*oto->r.x - y*oto->d.x - z*oto->f.x)*gpow2[OCT_MAXLS-ls0]*4.0,
                                                                         (fy - x*oto->r.y - y*oto->d.y - z*oto->f.y)*gpow2[OCT_MAXLS-ls0]*4.0,
                                                                         (fz - x*oto->r.z - y*oto->d.z - z*oto->f.z)*gpow2[OCT_MAXLS-ls0]*4.0);
            if (!j) continue;

            k = splitmsk; if (inode1->sol&i) k &= ~2;
            if (oct_touch_oct_recur(oto,inode0,x0,y0,z0,ls0,&((octv_t *)oto->loct[1]->nod.buf)[popcount8(inode1->chi&(i-1))+ind],x,y,z,ls1,k)) return(1);
        }
        return(0);
    }

    oto->hit[0].x = x0; oto->hit[0].y = y0; oto->hit[0].z = z0; oto->hit[0].ls = ls0;
    oto->hit[1].x = x1; oto->hit[1].y = y1; oto->hit[1].z = z1; oto->hit[1].ls = ls1;
    return(1);
}

int oct_touch_oct (oct_t *bloct0, point3d *bp0, point3d *br0, point3d *bd0, point3d *bf0, oct_t *bloct1, point3d *bp1, point3d *br1, point3d *bd1, point3d *bf1, oct_hit_t *hit)
{
    oto_t oto;
    octv_t *inode;
    oct_t *loct0, *loct1;
    point3d p0, r0, d0, f0, p1, r1, d1, f1;
    point3d cr, cd, cf, dr, dd, df;
    double f, nx, ny, nz, mat[9];
    int i, s, v, x, y, z, ind, ls0, ls1;

    if (gpow2[0] == 0.0) //FIX:should init elsewhere
    {
        gpow2[OCT_MAXLS] = 1.0;
        for(i=OCT_MAXLS+1;i<=OCT_MAXLS*2;i++) gpow2[i] = gpow2[i-1]*2.0;
        for(i=OCT_MAXLS-1;i>=          0;i--) gpow2[i] = gpow2[i+1]*0.5;
    }

        //Make sure that oct_touch_oct_recur() chooses to split 0 on 1st call
    if (bloct0->lsid >= bloct1->lsid)
    {
        loct0 = bloct0; p0 = (*bp0); r0 = (*br0); d0 = (*bd0); f0 = (*bf0);
        loct1 = bloct1; p1 = (*bp1); r1 = (*br1); d1 = (*bd1); f1 = (*bf1);
    }
    else
    {
        loct0 = bloct1; p0 = (*bp1); r0 = (*br1); d0 = (*bd1); f0 = (*bf1);
        loct1 = bloct0; p1 = (*bp0); r1 = (*br0); d1 = (*bd0); f1 = (*bf0);
    }

    //FIXFIXFIXFIX:can remove this inverse now that pgram3d handles it!
    //transform loct1 to make loct0 identity (p0=<0,0,0>, r0=<1,0,0>, d0=<0,1,0>, f0=<0,0,1>)
    //oto = loct0^-1 * loct1
    invert3x3(&r0,&d0,&f0,mat);
    oto.r.x = r1.x*mat[0] + r1.y*mat[1] + r1.z*mat[2];
    oto.r.y = r1.x*mat[3] + r1.y*mat[4] + r1.z*mat[5];
    oto.r.z = r1.x*mat[6] + r1.y*mat[7] + r1.z*mat[8];
    oto.d.x = d1.x*mat[0] + d1.y*mat[1] + d1.z*mat[2];
    oto.d.y = d1.x*mat[3] + d1.y*mat[4] + d1.z*mat[5];
    oto.d.z = d1.x*mat[6] + d1.y*mat[7] + d1.z*mat[8];
    oto.f.x = f1.x*mat[0] + f1.y*mat[1] + f1.z*mat[2];
    oto.f.y = f1.x*mat[3] + f1.y*mat[4] + f1.z*mat[5];
    oto.f.z = f1.x*mat[6] + f1.y*mat[7] + f1.z*mat[8];
    nx = p1.x-p0.x; ny = p1.y-p0.y; nz = p1.z-p0.z;
    oto.p.x = nx*mat[0] + ny*mat[1] + nz*mat[2];
    oto.p.y = nx*mat[3] + ny*mat[4] + nz*mat[5];
    oto.p.z = nx*mat[6] + ny*mat[7] + nz*mat[8];

    oto.loct[0] = loct0;
    oto.loct[1] = loct1;

    //FIX:optimize by precalcing only for: -std::max(loct0->lsid,loct1->lsid) .. +std::max(loct0->lsid,loct1->lsid)
    ls0 = loct0->lsid; ls1 = loct1->lsid;
    cr.x = 2.0; cr.y = 0.0; cr.z = 0.0;
    cd.x = 0.0; cd.y = 2.0; cd.z = 0.0;
    cf.x = 0.0; cf.y = 0.0; cf.z = 2.0;
    for(i=OCT_MAXLS-ls0;i<=OCT_MAXLS+ls1;i++)
    {
        f = gpow2[i]*2.0;
        dr.x = oto.r.x*f; dr.y = oto.r.y*f; dr.z = oto.r.z*f;
        dd.x = oto.d.x*f; dd.y = oto.d.y*f; dd.z = oto.d.z*f;
        df.x = oto.f.x*f; df.y = oto.f.y*f; df.z = oto.f.z*f;
        pgram3d_init(&oto.bi[i],&cr,&cd,&cf,&dr,&dd,&df);
    }

    //full size early out check
    if (!pgram3d_isint(&oto.bi[ls1-ls0+OCT_MAXLS],((1<<ls0)*.5 - ((1<<ls1)*.5*oto.r.x + (1<<ls1)*.5*oto.d.x + (1<<ls1)*.5*oto.f.x + oto.p.x))*gpow2[OCT_MAXLS-ls0]*4.0,
                                                                 ((1<<ls0)*.5 - ((1<<ls1)*.5*oto.r.y + (1<<ls1)*.5*oto.d.y + (1<<ls1)*.5*oto.f.y + oto.p.y))*gpow2[OCT_MAXLS-ls0]*4.0,
                                                                 ((1<<ls0)*.5 - ((1<<ls1)*.5*oto.r.z + (1<<ls1)*.5*oto.d.z + (1<<ls1)*.5*oto.f.z + oto.p.z))*gpow2[OCT_MAXLS-ls0]*4.0)) return(0);

    s = ((1<<ls1)>>1); inode = &((octv_t *)loct1->nod.buf)[loct1->head]; ind = inode->ind;
    for(v=(inode->sol|inode->chi);v;v^=i)
    {
        i = (-v)&v;
        x = 0; if (i&0xaa) x = s;
        y = 0; if (i&0xcc) y = s;
        z = 0; if (i&0xf0) z = s;
        if (!oct_touch_oct_recur(&oto,&((octv_t *)loct0->nod.buf)[loct0->head],0,0,0,ls0,&((octv_t *)loct1->nod.buf)[popcount8(inode->chi&(i-1))+ind],x,y,z,ls1-1,3)) continue;
        if (hit)
        {
            if (bloct0->lsid >= bloct1->lsid) { hit[0] = oto.hit[0]; hit[1] = oto.hit[1]; } //must preserve order for return value!
                                                  else { hit[0] = oto.hit[1]; hit[1] = oto.hit[0]; }
        }
        return(1);
    }
    return(0);
}

int oct_touch_oct (oct_t *bloct0, point3f *bp0, point3f *br0, point3f *bd0, point3f *bf0, oct_t *bloct1, point3f *bp1, point3f *br1, point3f *bd1, point3f *bf1, oct_hit_t *hit)
{
    point3d dp0, dr0, dd0, df0, dp1, dr1, dd1, df1;
    dp0.x = (double)bp0->x; dp0.y = (double)bp0->y; dp0.z = (double)bp0->z;
    dr0.x = (double)br0->x; dr0.y = (double)br0->y; dr0.z = (double)br0->z;
    dd0.x = (double)bd0->x; dd0.y = (double)bd0->y; dd0.z = (double)bd0->z;
    df0.x = (double)bf0->x; df0.y = (double)bf0->y; df0.z = (double)bf0->z;
    dp1.x = (double)bp1->x; dp1.y = (double)bp1->y; dp0.z = (double)bp0->z;
    dr1.x = (double)br1->x; dr1.y = (double)br1->y; dr0.z = (double)br0->z;
    dd1.x = (double)bd1->x; dd1.y = (double)bd1->y; dd0.z = (double)bd0->z;
    df1.x = (double)bf1->x; df1.y = (double)bf1->y; df0.z = (double)bf0->z;
    return(oct_touch_oct(bloct0,&dp0,&dr0,&dd0,&df0,bloct1,&dp1,&dr1,&dd1,&df1,hit));
}

//--------------------------------------------------------------------------------------------------
    //See VOX_MOI.KC for derivation
void oct_getmoi (oct_t *loct, double *mas, double *cx, double *cy, double *cz, double *ixx, double *iyy, double *izz, double *ixy, double *ixz, double *iyz)
{
    typedef struct { octv_t *ptr; int j; } stk_t;
    stk_t stk[OCT_MAXLS];
    octv_t *ptr;
    double f;

    //swwlut[ls] = ((((2^ls*2 - 3)*2^ls) + 1)*2^(ls*3))/6;  NOTE:swwlut[14] >= 2^64
    static const int64_t swwlut[OCT_MAXLS] = {
        0,4,224,8960,317440,10665984,349569024,11319377920,364359188480,11693786660864,374750392090624,12000804344954880,384166442167173120,static_cast<int64_t>(12295577674285318000),0,0
    };
    //swlut[ls] = (2^(ls*4-1)) - (2^(ls*3-1));
    static const int64_t swlut[OCT_MAXLS] = {
        0,4,96,1792,30720,507904,8257536,133169152,2139095040,34292629504,549218942976,8791798054912,140703128616960,2251524935778304,36026597995708416,576443160117379070
    };

    int64_t sx = 0, sy = 0, sz = 0, sxx = 0, syy = 0, szz = 0, sxy = 0, sxz = 0, syz = 0;
    int64_t qx, qy, qz, sw = 0, sww = 0, swx = 0, swy = 0, swz = 0;
    int i, j, k, ls, s, x, y, z, xx, yy, zz, sn = 0;
    #define FASTMOI 1

    ls = loct->lsid-1; s = (1<<ls); ptr = &((octv_t *)loct->nod.buf)[loct->head]; x = 0; y = 0; z = 0; j = 8-1;
    while (1)
    {
        if (!ls) //optimization: calculates all 8 leaf nodes simultaneously
        {
            k = ptr->sol; i = popcount8(k); sn += i;
            j = popcount8(k&0xaa); xx = i*x + j; sx += xx; sxx += (xx+j)*x + j;    //Be careful with j when moving these!
            j = popcount8(k&0xcc); yy = i*y + j; sy += yy; syy += (yy+j)*y + j; sxy += xx*y + j*x + popcount8(k&0x88);
            j = popcount8(k&0xf0); zz = i*z + j; sz += zz; szz += (zz+j)*z + j; sxz += xx*z + j*x + popcount8(k&0xa0);
                                                                                                     syz += yy*z + j*y + popcount8(k&0xc0);
            j = 0; goto tosibly;
        }

        i = (1<<j); if (!((ptr->chi|ptr->sol)&i)) goto tosibly; //skip pure air
        k = -s-s;
        x &= k; if (j&1) x += s;
        y &= k; if (j&2) y += s;
        z &= k; if (j&4) z += s;

        if (!(ptr->sol&i)) //hybrid
        {
            stk[ls].ptr = ptr; stk[ls].j = j; ls--; s >>= 1; //2child
            ptr = &((octv_t *)loct->nod.buf)[popcount8(ptr->chi&(i-1)) + ptr->ind]; j = 8-1;
            continue;
        }

             //pure solid; process node (x,y,z,ls)
#if (FASTMOI != 0)
        if (ls == 1) //faster because using 32-bit temps
        {
            sw += swlut[1]; sww += swwlut[1]; sn += 8;
            xx = (x<<3); sx += xx; swx += xx; sxx += xx*x; sxy += xx*y;
            yy = (y<<3); sy += yy; swy += yy; syy += yy*y; sxz += xx*z;
            zz = (z<<3); sz += zz; swz += zz; szz += zz*z; syz += yy*z;
        }
        else //general solution for 64-bit precision
        {
            sw += swlut[ls]; sww += swwlut[ls]; k = ls*3; sn += (1<<k); i = s-1;
            qx = (((int64_t)x)<<k); sx += qx; swx += qx*i; sxx += qx*x; sxy += qx*y;
            qy = (((int64_t)y)<<k); sy += qy; swy += qy*i; syy += qy*y; sxz += qx*z;
            qz = (((int64_t)z)<<k); sz += qz; swz += qz*i; szz += qz*z; syz += qy*z;
        }
#else
        for(zz=z;zz<z+s;zz++)
            for(yy=y;yy<y+s;yy++)
                for(xx=x;xx<x+s;xx++)
                    { sn++; sx += xx; sy += yy; sz += zz; sxx += xx*xx; syy += yy*yy; szz += zz*zz; sxy += xx*yy; sxz += xx*zz; syz += yy*zz; }
#endif

tosibly:;
        j--; if (j >= 0) continue;
        do { ls++; s <<= 1; if (ls >= loct->lsid) goto break2; j = stk[ls].j-1; } while (j < 0); //2parent
        ptr = stk[ls].ptr;
    }
break2:;
#if (FASTMOI != 0)
    sx += sw; sxx += swx+sww;
    sy += sw; syy += swy+sww;
    sz += sw; szz += swz+sww;
    sww = sww*3 - sw;
    sxy += (((swx+swy)*2 + sww)>>2);
    sxz += (((swx+swz)*2 + sww)>>2);
    syz += (((swy+swz)*2 + sww)>>2);
#endif
    (*mas) = sn; if (sn > 0.0) f = 1.0/sn; else f = 0.0; (*cx) = sx*f; (*cy) = sy*f; (*cz) = sz*f;
    (*ixx) = syy+szz - ((*cy)*(*cy) + (*cz)*(*cz))*(*mas); (*ixy) = (*cx)*(*cy)*(*mas) - sxy;
    (*iyy) = sxx+szz - ((*cx)*(*cx) + (*cz)*(*cz))*(*mas); (*ixz) = (*cx)*(*cz)*(*mas) - sxz;
    (*izz) = sxx+syy - ((*cx)*(*cx) + (*cy)*(*cy))*(*mas); (*iyz) = (*cy)*(*cz)*(*mas) - syz;
    (*cx) += .5; (*cy) += .5; (*cz) += .5;

        //Correct values for default model: 3.42269e6, 122.157, 125.582, 123.841, 3.49237e10, 3.64484e10, 3.83971e10, -3.1043e9, -3.2001e9, -6.77705e8
    //{ char tbuf[2048]; sprintf(tbuf,"%g\n%g\n%g\n%g\n%g\n%g\n%g\n%g\n%g\n%g",*mas,*cx,*cy,*cz,*ixx,*iyy,*izz,*ixy,*ixz,*iyz); MessageBox(ghwnd,tbuf,prognam,MB_OK); }
}
void oct_getmoi (oct_t *loct, float *mas, float *cx, float *cy, float *cz, float *ixx, float *iyy, float *izz, float *ixy, float *ixz, float *iyz)
{
    double dmas, dcx, dcy, dcz, dixx, diyy, dizz, dixy, dixz, diyz;
    oct_getmoi(loct,&dmas,&dcx,&dcy,&dcz,&dixx,&diyy,&dizz,&dixy,&dixz,&diyz);
    (*mas) = (float)dmas;
    (*cx ) = (float)dcx;  (*cy ) = (float)dcy;  (*cz ) = (float)dcz;
    (*ixx) = (float)dixx; (*iyy) = (float)diyy; (*izz) = (float)dizz;
    (*ixy) = (float)dixy; (*ixz) = (float)dixz; (*iyz) = (float)diyz;
}

typedef struct { octv_t *octv; int bx, by, bz, bdx, bdy, bdz, stat; } ogbs_t;
static void oct_getboxstate_recur (ogbs_t *ogbs, int x0, int y0, int z0, int s, int ind)
{
    octv_t *ptr;
    int i, x, y, z, v;

    ptr = &ogbs->octv[ind]; v = (1<<8)-1;
    x = x0+s-ogbs->bx; if (x <= 0) v &= 0xaa; else if (x >= ogbs->bdx) v &= 0x55;
    y = y0+s-ogbs->by; if (y <= 0) v &= 0xcc; else if (y >= ogbs->bdy) v &= 0x33;
    z = z0+s-ogbs->bz; if (z <= 0) v &= 0xf0; else if (z >= ogbs->bdz) v &= 0x0f;
    if ( v&          ptr->sol       ) { ogbs->stat |= 2; if (ogbs->stat == 3) return; }
    if ((v&(ptr->chi|ptr->sol)) != v) { ogbs->stat |= 1; if (ogbs->stat == 3) return; }
    if (s == 1) return;
    v = (v&ptr->chi&(~ptr->sol));
    for(;v;v^=i)
    {
        i = (-v)&v;
        x = x0; if (i&0xaa) x += s;
        y = y0; if (i&0xcc) y += s;
        z = z0; if (i&0xf0) z += s;
        oct_getboxstate_recur(ogbs,x,y,z,s>>1,popcount8((i-1)&ptr->chi)+ptr->ind);
        if (ogbs->stat == 3) return;
    }
}
    //returns: 0:pure air, 1:some air&sol, 2:pure sol
int oct_getboxstate (oct_t *loct, int bx0, int by0, int bz0, int bx1, int by1, int bz1)
{
    ogbs_t ogbs;

    ogbs.octv = (octv_t *)loct->nod.buf;
    ogbs.bx = bx0; ogbs.bdx = bx1-bx0;
    ogbs.by = by0; ogbs.bdy = by1-by0;
    ogbs.bz = bz0; ogbs.bdz = bz1-bz0;
    ogbs.stat = 0;
    oct_getboxstate_recur(&ogbs,0,0,0,loct->sid>>1,loct->head);
    if (ogbs.stat == 3) return(1);
    return(ogbs.stat&2);
}

//Test intersection between infinite ray and box; returns 1 if intersect else 0
//See INBOX.KC for test program
static int inbox3d (point3d *p, point3d *d, double bx, double by, double bz, double bs)
{
    double hs, x, y, z;
        //x = d->x*t + p->x;
        //y = d->y*t + p->y;
        //z = d->z*t + p->z;
    hs = bs*.5; bx += hs-p->x; by += hs-p->y; bz += hs-p->z;

    if (fabs(d->x) > 0.0)
    {
        x = fabs(d->x)*hs; y = d->x*by; z = d->x*bz;
        if (std::max(fabs((bx-hs)*d->y - y),fabs((bx-hs)*d->z - z)) <= x) return(1);
        if (std::max(fabs((bx+hs)*d->y - y),fabs((bx+hs)*d->z - z)) <= x) return(1);
    }

    if (fabs(d->y) > 0.0)
    {
        y = fabs(d->y)*hs; x = d->y*bx; z = d->y*bz;
        if (std::max(fabs((by-hs)*d->x - x),fabs((by-hs)*d->z - z)) <= y) return(1);
        if (std::max(fabs((by+hs)*d->x - x),fabs((by+hs)*d->z - z)) <= y) return(1);
    }

    if (fabs(d->z) > 0.0)
    {
        z = fabs(d->z)*hs; x = d->z*bx; y = d->z*by;
        if (std::max(fabs((bz-hs)*d->x - x),fabs((bz-hs)*d->y - y)) <= z) return(1);
        if (std::max(fabs((bz+hs)*d->x - x),fabs((bz+hs)*d->y - y)) <= z) return(1);
    }

    return(0);
}

int oct_hitscan (oct_t *loct, point3d *p, point3d *padd, point3i *hit, int *rhitdir, double *fracwent)
{
    typedef struct { octv_t *ptr; int x, y, z, j; } stk_t;
    stk_t stk[OCT_MAXLS];
    octv_t *ptr;
    double d, f, fx, fy, fz;
    int i, j, k, ls, x, y, z, nx, ny, nz, ix, iy, iz, ox, oy, oz, ox2, oy2, oz2, ord;

    if (fabs(padd->x) < 1e-8) { fx = 0.0; ix = 0; } else { fx = padd->x; ix = (fx>0)*2-1; }
    if (fabs(padd->y) < 1e-8) { fy = 0.0; iy = 0; } else { fy = padd->y; iy = (fy>0)*2-1; }
    if (fabs(padd->z) < 1e-8) { fz = 0.0; iz = 0; } else { fz = padd->z; iz = (fz>0)*2-1; }

    ox = (int)floor(p->x); ox2 = (int)floor(p->x+fx);
    oy = (int)floor(p->y); oy2 = (int)floor(p->y+fy);
    oz = (int)floor(p->z); oz2 = (int)floor(p->z+fz);
    ord = (fz >= 0.0)*4 + (fy >= 0.0)*2 + (fx >= 0.0);

    ls = loct->lsid-1; ptr = &((octv_t *)loct->nod.buf)[loct->head]; x = 0; y = 0; z = 0; j = 8-1;
    while (1)
    {
        k = j^ord;

        i = (1<<k); if (!(ptr->chi&i)) goto tosibly;

        nx = (( k    &1)<<ls)+x;
        ny = (((k>>1)&1)<<ls)+y;
        nz = (((k>>2)&1)<<ls)+z;

        if (!inbox3d(p,padd,nx,ny,nz,1<<ls)) goto tosibly;
              if (ix > 0) { if ((nx+(1<<ls) <= ox) || (nx         >  ox2)) goto tosibly; }
        else if (ix < 0) { if ((nx         >  ox) || (nx+(1<<ls) <= ox2)) goto tosibly; }
              if (iy > 0) { if ((ny+(1<<ls) <= oy) || (ny         >  oy2)) goto tosibly; }
        else if (iy < 0) { if ((ny         >  oy) || (ny+(1<<ls) <= oy2)) goto tosibly; }
              if (iz > 0) { if ((nz+(1<<ls) <= oz) || (nz         >  oz2)) goto tosibly; }
        else if (iz < 0) { if ((nz         >  oz) || (nz+(1<<ls) <= oz2)) goto tosibly; }

        if (ls <= 0)
        {
            if (hit) { hit->x = nx; hit->y = ny; hit->z = nz; }
            if ((rhitdir) || (fracwent))
            {
                d = 0.0;
                if (rhitdir) (*rhitdir) = 0;
                if (fx != 0.0) { f = (nx+(fx<0.0)-p->x)/fx; if (f > d) { d = f; if (rhitdir) (*rhitdir) = (fx<0.0)+0; } }
                if (fy != 0.0) { f = (ny+(fy<0.0)-p->y)/fy; if (f > d) { d = f; if (rhitdir) (*rhitdir) = (fy<0.0)+2; } }
                if (fz != 0.0) { f = (nz+(fz<0.0)-p->z)/fz; if (f > d) { d = f; if (rhitdir) (*rhitdir) = (fz<0.0)+4; } }
                if (fracwent) (*fracwent) = d;
            }
            return(popcount8(ptr->chi&(i-1)) + ptr->ind);
        }

        stk[ls].ptr = ptr; stk[ls].x = x; stk[ls].y = y; stk[ls].z = z; stk[ls].j = j; ls--; //2child
        ptr = &((octv_t *)loct->nod.buf)[popcount8(ptr->chi&(i-1)) + ptr->ind]; x = nx; y = ny; z = nz; j = 8-1;
        continue;

tosibly:;
        j--; if (j >= 0) continue;
        do { ls++; if (ls >= loct->lsid) { if (fracwent) (*fracwent) = 1.0; return(-1); } j = stk[ls].j-1; } while (j < 0); //2parent
        ptr = stk[ls].ptr; x = stk[ls].x; y = stk[ls].y; z = stk[ls].z;
    }
}
int oct_hitscan (oct_t *loct, point3f *fp, point3f *fpadd, point3i *hit, int *rhitdir, float *fracwent)
{
    point3d dp, dpadd;
    double dfracwent;
    int hitind;

    dp.x = fp->x; dp.y = fp->y; dp.z = fp->z;
    dpadd.x = fpadd->x; dpadd.y = fpadd->y; dpadd.z = fpadd->z;
    hitind = oct_hitscan(loct,&dp,&dpadd,hit,rhitdir,&dfracwent);
    (*fracwent) = dfracwent;
    return(hitind);
}

//--------------------------------------------------------------------------------------------------
    //06/06/2011:Ported from CYL_INTERSECT.KC (see .KC for derivation)
struct int_cube_cyl_3d_t
{
    double x0, y0, z0, dx, dy, dz, rad, rad2, rdx, rdy, rdz;
    double dyz2, dxz2, dxy2, dxyz2, rdyz2, rdxz2, rdxy2, rdxyz2;
    int sdx, sdy, sdz, dumalign;
};
static void int_cube_cyl_3d_init (int_cube_cyl_3d_t *ic, double x0, double y0, double z0, double x1, double y1, double z1, double rad)
{
    double dx2, dy2, dz2;
    ic->x0 = x0; ic->dx = x1-x0; if (ic->dx) ic->rdx = 1.0/ic->dx;
    ic->y0 = y0; ic->dy = y1-y0; if (ic->dy) ic->rdy = 1.0/ic->dy;
    ic->z0 = z0; ic->dz = z1-z0; if (ic->dz) ic->rdz = 1.0/ic->dz;
    if (ic->dx > 0.0) ic->sdx = 1; else if (ic->dx < 0.0) ic->sdx = -1; else ic->sdx = 0;
    if (ic->dy > 0.0) ic->sdy = 1; else if (ic->dy < 0.0) ic->sdy = -1; else ic->sdy = 0;
    if (ic->dz > 0.0) ic->sdz = 1; else if (ic->dz < 0.0) ic->sdz = -1; else ic->sdz = 0;
    dx2 = ic->dx*ic->dx; dy2 = ic->dy*ic->dy; dz2 = ic->dz*ic->dz;
    ic->dyz2 = dy2 + dz2; if (ic->dyz2) ic->rdyz2 = 1.0/ic->dyz2;
    ic->dxz2 = dx2 + dz2; if (ic->dxz2) ic->rdxz2 = 1.0/ic->dxz2;
    ic->dxy2 = dx2 + dy2; if (ic->dxy2) ic->rdxy2 = 1.0/ic->dxy2;
    ic->dxyz2 = dx2 + ic->dyz2; if (ic->dxyz2) ic->rdxyz2 = 1.0/ic->dxyz2;
    ic->rad = rad; ic->rad2 = rad*rad;
}
    //returns t, where:
    //   goalx = ic->dx*t + ic->x0
    //   goaly = ic->dy*t + ic->y0
    //   goalz = ic->dz*t + ic->z0
    //and hitx/y/z is nearest point on cube (bx,by,bz,bs) to goalx/y/z
static double int_cube_cyl_3d_gett (int_cube_cyl_3d_t *ic, double bx, double by, double bz, double bs, double tmin)
{
    double f, t, xx, yy, zz, Zb, Zc, insqr, hbs;
    int x, y, z, u, v;

    hbs = bs*.5;
    bx += hbs-ic->x0;
    by += hbs-ic->y0;
    bz += hbs-ic->z0;

        //early exit bounding box
    f = hbs+ic->rad;
    t = ic->dx*tmin; if ((bx > std::max(t,0.0)+f) || (bx < std::min(t,0.0)-f)) return(tmin);
    t = ic->dy*tmin; if ((by > std::max(t,0.0)+f) || (by < std::min(t,0.0)-f)) return(tmin);
    t = ic->dz*tmin; if ((bz > std::max(t,0.0)+f) || (bz < std::min(t,0.0)-f)) return(tmin);

        //tests rounded cube at start in case it already intersects..
    xx = fabs(bx); xx -= std::min(xx,hbs);
    yy = fabs(by); yy -= std::min(yy,hbs);
    zz = fabs(bz); zz -= std::min(zz,hbs);
    if (xx*xx + yy*yy + zz*zz < ic->rad2) return(0.0);

        //test 3/6 faces
    f = hbs+ic->rad;
    if (ic->sdx) { t = (bx-(double)ic->sdx*f)*ic->rdx; if ((t >= 0.0) && (t < tmin) && (std::max(fabs(ic->dy*t-by),fabs(ic->dz*t-bz)) < hbs)) tmin = t; }
    if (ic->sdy) { t = (by-(double)ic->sdy*f)*ic->rdy; if ((t >= 0.0) && (t < tmin) && (std::max(fabs(ic->dx*t-bx),fabs(ic->dz*t-bz)) < hbs)) tmin = t; }
    if (ic->sdz) { t = (bz-(double)ic->sdz*f)*ic->rdz; if ((t >= 0.0) && (t < tmin) && (std::max(fabs(ic->dx*t-bx),fabs(ic->dy*t-by)) < hbs)) tmin = t; }

        //test 7/8 verts
    if (ic->dxyz2)
        for(x=-1;x<=1;x+=2)
            for(y=-1;y<=1;y+=2)
                for(z=-1;z<=1;z+=2)
                {
                    if ((x == ic->sdx) && (y == ic->sdy) && (z == ic->sdz)) continue;
                    xx = (double)x*hbs+bx; yy = (double)y*hbs+by; zz = (double)z*hbs+bz; //(xx - ic->dx*t)^2+ "y + "z = ic->rad2)
                    Zb = xx*ic->dx + yy*ic->dy + zz*ic->dz;
                    Zc = xx*xx + yy*yy + zz*zz - ic->rad2;
                    insqr = Zb*Zb - ic->dxyz2*Zc; if (insqr < 0.0) continue;
                    t = (Zb-sqrt(insqr))*ic->rdxyz2; if ((t >= 0.0) && (t < tmin)) tmin = t;
                }

        //test 9/12 edges
    for(v=-1;v<=1;v+=2)
        for(u=-1;u<=1;u+=2)
        {
            if ((ic->dyz2) && ((u != ic->sdy) || (v != ic->sdz)))
            {
                yy = hbs*(double)u+by; zz = hbs*(double)v+bz; //(yy - ic->dy*t)^2 + "z = rad2
                Zb = yy*ic->dy + zz*ic->dz;
                Zc = yy*yy + zz*zz - ic->rad2;
                insqr = Zb*Zb - ic->dyz2*Zc;
                if (insqr >= 0.0) { t = (Zb-sqrt(insqr))*ic->rdyz2; if ((t >= 0.0) && (t < tmin) && (fabs(ic->dx*t-bx) < hbs)) tmin = t; }
            }
            if ((ic->dxz2) && ((u != ic->sdx) || (v != ic->sdz)))
            {
                xx = hbs*(double)u+bx; zz = hbs*(double)v+bz; //(xx - ic->dx*t)^2 + "z = rad2
                Zb = xx*ic->dx + zz*ic->dz;
                Zc = xx*xx + zz*zz - ic->rad2;
                insqr = Zb*Zb - ic->dxz2*Zc;
                if (insqr >= 0.0) { t = (Zb-sqrt(insqr))*ic->rdxz2; if ((t >= 0.0) && (t < tmin) && (fabs(ic->dy*t-by) < hbs)) tmin = t; }
            }
            if ((ic->dxy2) && ((u != ic->sdx) || (v != ic->sdy)))
            {
                xx = hbs*(double)u+bx; yy = hbs*(double)v+by; //(xx - ic->dx*t)^2 + "y = rad2
                Zb = xx*ic->dx + yy*ic->dy;
                Zc = xx*xx + yy*yy - ic->rad2;
                insqr = Zb*Zb - ic->dxy2*Zc;
                if (insqr >= 0.0) { t = (Zb-sqrt(insqr))*ic->rdxy2; if ((t >= 0.0) && (t < tmin) && (fabs(ic->dz*t-bz) < hbs)) tmin = t; }
            }
        }

    return(tmin);
}

static double int_tri_cyl_3d_gett (const point3d *c0, double cr, const point3d *cv, const point3d *tri, double tmin, point3d *hit)
{
    point3d ntri[3];
    double f, t, u, v, cr2, ux, uy, uz, vx, vy, vz, k0, k1, k2, k3, k4, k5, k6, k7;
    double k8, k9, ka, kb, kc, kd, ke, kf, kg, kh, Za, Zb, Zc, insqr, rZa, px, py, pz, dx, dy, dz;
    int i, j, s;

    cr2 = cr*cr;
    for(i=3-1;i>=0;i--)
    {
        ntri[i].x = tri[i].x-c0->x;
        ntri[i].y = tri[i].y-c0->y;
        ntri[i].z = tri[i].z-c0->z;
    }

        //Check plane..
    ux = ntri[1].x-ntri[0].x; vx = ntri[2].x-ntri[0].x;
    uy = ntri[1].y-ntri[0].y; vy = ntri[2].y-ntri[0].y;
    uz = ntri[1].z-ntri[0].z; vz = ntri[2].z-ntri[0].z;
        //hx = ux*u + vx*v + ntri[0].x, ix = cv->x*t
        //hy = uy*u + vy*v + ntri[0].y, iy = cv->y*t
        //hz = uz*u + vz*v + ntri[0].z, iz = cv->z*t
        //(hx-ix)*ux + (hy-iy)*uy + (hz-iz)*uz = 0
        //(hx-ix)*vx + (hy-iy)*vy + (hz-iz)*vz = 0
        //(hx-ix)^2  + (hy-iy)^2  + (hz-iz)^2 = cr2
        //-----------------------------------------
        //(ux*u + vx*v - cv->x*t + ntri[0].x)*ux +
        //(uy*u + vy*v - cv->y*t + ntri[0].y)*uy +
        //(uz*u + vz*v - cv->z*t + ntri[0].z)*uz = 0
        //
        //(ux*u + vx*v - cv->x*t + ntri[0].x)*vx +
        //(uy*u + vy*v - cv->y*t + ntri[0].y)*vy +
        //(uz*u + vz*v - cv->z*t + ntri[0].z)*vz = 0
        //
        //(ux*u + vx*v - cv->x*t + ntri[0].x)^2 +
        //(uy*u + vy*v - cv->y*t + ntri[0].y)^2 +
        //(uz*u + vz*v - cv->z*t + ntri[0].z)^2 = cr2
    k0 = ux*ux + uy*uy + uz*uz;
    k1 = ux*vx + uy*vy + uz*vz;
    k5 = vx*vx + vy*vy + vz*vz;
    k2 = cv->x*ux + cv->y*uy + cv->z*uz; k3 = ntri[0].x*ux + ntri[0].y*uy + ntri[0].z*uz;
    k6 = cv->x*vx + cv->y*vy + cv->z*vz; k7 = ntri[0].x*vx + ntri[0].y*vy + ntri[0].z*vz;
        //k0*u + k1*v + -k2*t = -k3
        //k1*u + k5*v + -k6*t = -k7
    f = k0*k5 - k1*k1;
    k8 = k2*k5 - k1*k6; k9 = k1*k7 - k3*k5;
    ka = k0*k6 - k1*k2; kb = k1*k3 - k0*k7;
        //u = (t*k8 + k9)/f
        //v = (t*ka + kb)/f
        //(ux*(t*k8 + k9) + vx*(t*ka + kb) - cv->x*t*f + ntri[0].x*f)^2 +
        //(uy*(t*k8 + k9) + vy*(t*ka + kb) - cv->y*t*f + ntri[0].y*f)^2 +
        //(uz*(t*k8 + k9) + vz*(t*ka + kb) - cv->z*t*f + ntri[0].z*f)^2 = cr2*f*f
    kc = ux*k8 + vx*ka - cv->x*f; kd = ux*k9 + vx*kb + ntri[0].x*f;
    ke = uy*k8 + vy*ka - cv->y*f; kf = uy*k9 + vy*kb + ntri[0].y*f;
    kg = uz*k8 + vz*ka - cv->z*f; kh = uz*k9 + vz*kb + ntri[0].z*f;
        //(kc*t + kd)^2 +
        //(ke*t + kf)^2 +
        //(kg*t + kh)^2 = cr2*f*f
    Za = kc*kc + ke*ke + kg*kg;
    Zb = kc*kd + ke*kf + kg*kh;
    Zc = kd*kd + kf*kf + kh*kh - cr2*f*f;
    insqr = Zb*Zb - Za*Zc;
    if ((insqr >= 0.0) && (Za != 0.0))
    {
        rZa = 1.0/Za;
        for(s=-1;s<=1;s+=2)
        {
            t = (sqrt(insqr)*(double)s - Zb)*rZa; if ((t < 0.0) || (t >= tmin)) continue;
            u = t*k8 + k9; if (u < 0.0) continue;
            v = t*ka + kb; if ((v < 0.0) || (u+v > f)) continue;
            tmin = t; f = 1.0/f; u *= f; v *= f;
            hit->x = ux*u + vx*v + tri[0].x;
            hit->y = uy*u + vy*v + tri[0].y;
            hit->z = uz*u + vz*v + tri[0].z;
        }
    }

    for(j=2,i=0;i<3;j=i,i++)
    {
        px = ntri[i].x; dx = ntri[j].x-px;
        py = ntri[i].y; dy = ntri[j].y-py;
        pz = ntri[i].z; dz = ntri[j].z-pz;

            //Check 3 line segments..
            //ix = dx*u + px
            //iy = dy*u + py
            //iz = dz*u + pz
            //
            //   //if sphere hits inside line segment,
            //   //(ix,iy,iz)-(cv->x*t,cv->y*t,cv->z*t) must be
            //   //perpendicular to line ntri[i]-ntri[j]
            //(ix-cv->x*t)*dx + "y + "z = 0
            //(ix-cv->x*t)^2 + "y + "z = cr^2, 0<t<tmin
            //
            //   2eq/2unk:
            //(dx*u + px - cv->x*t)*dx + "y + "z = 0
            //(dx*u + px - cv->x*t)^2 + "y + "z = cr^2
        f = cv->x*dx + cv->y*dy + cv->z*dz;
        k0 = dx*dx + dy*dy + dz*dz;
        k1 = px*dx + py*dy + pz*dz;
            //t = (k0*u + k1)/f;
            //
            //(dx*u*f + px*f - cv->x*(k0*u + k1))^2 + "y + "z = cr^2*f*f
            //((dx*f - cv->x*k0)*u + (px*f - cv->x*k1))^2 + "y + "z = cr^2*f*f
        k2 = dx*f - cv->x*k0; k5 = px*f - cv->x*k1;
        k3 = dy*f - cv->y*k0; k6 = py*f - cv->y*k1;
        k4 = dz*f - cv->z*k0; k7 = pz*f - cv->z*k1;
        Za = k2*k2 + k3*k3 + k4*k4;
        Zb = k2*k5 + k3*k6 + k4*k7;
        Zc = k5*k5 + k6*k6 + k7*k7 - cr2*f*f;
        insqr = Zb*Zb - Za*Zc;
        if ((insqr >= 0.0) && (Za != 0.0))
        {
            rZa = 1.0/Za;
            for(s=-1;s<=1;s+=2)
            {
                u = (sqrt(insqr)*(double)s - Zb)*rZa; if ((u < 0.0) || (u > 1.0)) continue;
                t = (k0*u + k1)/f; if ((t < 0.0) || (t >= tmin)) continue;
                tmin = t;
                hit->x = dx*u + tri[i].x;
                hit->y = dy*u + tri[i].y;
                hit->z = dz*u + tri[i].z;
            }
        }

            //Check 3 endpoints..
            //(cv->x*t - px)^2 + "y + "z = cr^2
        Za = cv->x*cv->x + cv->y*cv->y + cv->z*cv->z;
        Zb = cv->x*px + cv->y*py + cv->z*pz;
        Zc = px*px + py*py + pz*pz - cr2;
        insqr = Zb*Zb - Za*Zc;
        if ((insqr >= 0.0) && (Za != 0.0))
        {
            t = (Zb - sqrt(insqr))/Za;
            if ((t >= 0.0) && (t < tmin)) { tmin = t; (*hit) = tri[i]; }
        }
    }

    return(tmin);
}

surf_t *oct_sphtrace (oct_t *loct, point3d *p, point3d *padd, double rad, point3i *hit, point3d *pgoal, point3d *hitnorm)
{
    typedef struct { octv_t *ptr; int x, y, z, j; } stk_t;
    stk_t stk[OCT_MAXLS];
    octv_t *ptr;
    surf_t *lsurf = 0;
    int_cube_cyl_3d_t ic;
    double d, dx, dy, dz, tmin, ntmin;
    int i, j, k, ls, x, y, z, nx, ny, nz, hnx, hny, hnz, ord;

    if (pgoal)
    {
        pgoal->x = p->x+padd->x;
        pgoal->y = p->y+padd->y;
        pgoal->z = p->z+padd->z;
    }

    int_cube_cyl_3d_init(&ic,p->x,p->y,p->z,p->x+padd->x,p->y+padd->y,p->z+padd->z,rad);
    tmin = 1.0; ord = (padd->z >= 0.0)*4 + (padd->y >= 0.0)*2 + (padd->x >= 0.0);
    ls = loct->lsid-1; ptr = &((octv_t *)loct->nod.buf)[loct->head]; x = 0; y = 0; z = 0; j = 8-1;
    while (1)
    {
        k = j^ord;

        i = (1<<k); if (!(ptr->chi&i)) goto tosibly;

        nx = (( k    &1)<<ls)+x;
        ny = (((k>>1)&1)<<ls)+y;
        nz = (((k>>2)&1)<<ls)+z;

        ntmin = int_cube_cyl_3d_gett(&ic,nx,ny,nz,1<<ls,tmin);
        if (ntmin >= tmin) goto tosibly;

        if (ls <= 0)
        {
            hnx = nx; hny = ny; hnz = nz;
            tmin = ntmin;
            lsurf = &((surf_t *)loct->sur.buf)[popcount8(ptr->chi&(i-1)) + ptr->ind];
            goto tosibly;
        }

        stk[ls].ptr = ptr; stk[ls].x = x; stk[ls].y = y; stk[ls].z = z; stk[ls].j = j; ls--; //2child
        ptr = &((octv_t *)loct->nod.buf)[popcount8(ptr->chi&(i-1)) + ptr->ind]; x = nx; y = ny; z = nz; j = 8-1;
        continue;

tosibly:;
        j--; if (j >= 0) continue;
        do
        {
            ls++;
            if (ls >= loct->lsid)
            {
                dx = padd->x*tmin + p->x;
                dy = padd->y*tmin + p->y;
                dz = padd->z*tmin + p->z;
                if (pgoal)
                {
                    pgoal->x = dx;
                    pgoal->y = dy;
                    pgoal->z = dz;
                }
                if (!lsurf) return(0);
                if (hit) { hit->x = hnx; hit->y = hny; hit->z = hnz; }
                if (hitnorm)
                {
                        //hitx/y/z is nearest point on cube (hit->x,hit->y,hit->z,1) to pgoal
                    d = 1.0/rad;
                    hitnorm->x = (std::min(std::max(dx,(double)hnx),(double)(hnx+1)) - dx)*d;
                    hitnorm->y = (std::min(std::max(dy,(double)hny),(double)(hny+1)) - dy)*d;
                    hitnorm->z = (std::min(std::max(dz,(double)hnz),(double)(hnz+1)) - dz)*d;
                }
                return(lsurf);
            }
            j = stk[ls].j-1;
        } while (j < 0); //2parent
        ptr = stk[ls].ptr; x = stk[ls].x; y = stk[ls].y; z = stk[ls].z;
    }
}
surf_t *oct_sphtrace (oct_t *loct, point3f *fp, point3f *fpadd, float rad, point3i *hit, point3f *fpgoal, point3f *fhitnorm)
{
    point3d dp, dpadd, dpgoal, dhitnorm;
    surf_t *retval;
    dp.x = fp->x; dp.y = fp->y; dp.z = fp->z;
    dpadd.x = fpadd->x; dpadd.y = fpadd->y; dpadd.z = fpadd->z;
    retval = oct_sphtrace(loct,&dp,&dpadd,rad,hit,&dpgoal,&dhitnorm);
    fpgoal->x = dpgoal.x; fpgoal->y = dpgoal.y; fpgoal->z = dpgoal.z;
    fhitnorm->x = dhitnorm.x; fhitnorm->y = dhitnorm.y; fhitnorm->z = dhitnorm.z;
    return(retval);
}
//--------------------------------------------------------------------------------------------------

    //   x,y,z: vector difference of test point to cube's top-left-etc corner
    //     dia: cube's diameter
    // returns: distance squared from point to nearest point on cube (0.0 if inside)
static double dist2cube2_nearest (double x, double y, double z, double dia)
{
    double rad = dia*.5;
    x = fabs(x+rad); x -= std::min(x,rad);
    y = fabs(y+rad); y -= std::min(y,rad);
    z = fabs(z+rad); z -= std::min(z,rad);
    return(x*x + y*y + z*z);
}
static double dist2cube2_farthest (double x, double y, double z, double dia)
{
    double rad = dia*.5;
    x = fabs(x+rad)+dia;
    y = fabs(y+rad)+dia;
    z = fabs(z+rad)+dia;
    return(x*x + y*y + z*z);
}
static double dist2tri2_nearest_2d (point2d *p, point2d *p0, point2d *p1, point2d *p2)
{
    double x0, y0, x1, y1, x2, y2, x10, y10, x21, y21, x02, y02, t0, t1, t2;
    int a0, a1, a2, a3, a4, a5, sid;

    x0 = p->x-p0->x; y0 = p->y-p0->y; x10 = p1->x-p0->x; y10 = p1->y-p0->y;
    x1 = p->x-p1->x; y1 = p->y-p1->y; x21 = p2->x-p1->x; y21 = p2->y-p1->y;
    x2 = p->x-p2->x; y2 = p->y-p2->y; x02 = p0->x-p2->x; y02 = p0->y-p2->y;
    a5 = (x0*x02 + y0*y02 > 0.0);
    a0 = (x0*x10 + y0*y10 > 0.0); if (a5 > a0) return(x0*x0 + y0*y0);
    a1 = (x1*x10 + y1*y10 > 0.0);
    a2 = (x1*x21 + y1*y21 > 0.0); if (a1 > a2) return(x1*x1 + y1*y1);
    a3 = (x2*x21 + y2*y21 > 0.0);
    a4 = (x2*x02 + y2*y02 > 0.0); if (a3 > a4) return(x2*x2 + y2*y2);
    sid = (x10*y21 < y10*x21);
    t0 = x0*y10 - y0*x10;
    t1 = x1*y21 - y1*x21;
    t2 = x2*y02 - y2*x02;
    if (((t0<0.0) == sid) && (a0 > a1)) return(t0*t0/(x10*x10 + y10*y10));
    if (((t1<0.0) == sid) && (a2 > a3)) return(t1*t1/(x21*x21 + y21*y21));
    if (((t2<0.0) == sid) && (a4 > a5)) return(t2*t2/(x02*x02 + y02*y02));
    return(0.0);
}
static double dist2tri2_nearest_3d (point3d *p, point3d *p0, point3d *p1, point3d *p2)
{
    point3d a, b, n;
    point2d q, q0, q1, q2;
    double f;

    n.x = (p1->y-p0->y)*(p2->z-p0->z) - (p1->z-p0->z)*(p2->y-p0->y);
    n.y = (p1->z-p0->z)*(p2->x-p0->x) - (p1->x-p0->x)*(p2->z-p0->z);
    n.z = (p1->x-p0->x)*(p2->y-p0->y) - (p1->y-p0->y)*(p2->x-p0->x);
    f = 1.0/sqrt(n.x*n.x + n.y*n.y + n.z*n.z); n.x *= f; n.y *= f; n.z *= f;
    a.x = (p1->y-p0->y)*n.z - (p1->z-p0->z)*n.y;
    a.y = (p1->z-p0->z)*n.x - (p1->x-p0->x)*n.z;
    a.z = (p1->x-p0->x)*n.y - (p1->y-p0->y)*n.x;
    b.x = n.y*a.z - n.z*a.y;
    b.y = n.z*a.x - n.x*a.z;
    b.z = n.x*a.y - n.y*a.x;

    f   = (p->x-p0->x)*n.x + (p->y-p0->y)*n.y + (p->z-p0->z)*n.z;
    q.x = (p->x-p0->x)*a.x + (p->y-p0->y)*a.y + (p->z-p0->z)*a.z;
    q.y = (p->x-p0->x)*b.x + (p->y-p0->y)*b.y + (p->z-p0->z)*b.z;
    q0.x = 0.0;
    q0.y = 0.0;
    q1.x = 0.0;
    q1.y = (p1->x-p0->x)*b.x + (p1->y-p0->y)*b.y + (p1->z-p0->z)*b.z;
    q2.x = (p2->x-p0->x)*a.x + (p2->y-p0->y)*a.y + (p2->z-p0->z)*a.z;
    q2.y = (p2->x-p0->x)*b.x + (p2->y-p0->y)*b.y + (p2->z-p0->z)*b.z;
    return(dist2tri2_nearest_2d(&q,&q0,&q1,&q2)/(a.x*a.x + a.y*a.y + a.z*a.z) + f*f);
}

    //find closest distance from point to any voxel (renamed from findmaxcr)
double oct_balloonrad (oct_t *loct, point3d *p, double cr, point3i *hit, surf_t **hitsurf)
{
    typedef struct { octv_t *ptr; int x, y, z, j, ord; } stk_t;
    stk_t stk[OCT_MAXLS];
    octv_t *ptr;
    double cr2, ncr2;
    int i, j, k, s, ls, x, y, z, nx, ny, nz, ord, ipx, ipy, ipz;

    if (cr <= 0.0) return(-1.0);

    ipx = (int)floor(p->x); ipy = (int)floor(p->y); ipz = (int)floor(p->z);

    cr2 = cr*cr; hit->x = -1;
    ls = loct->lsid-1; s = (1<<ls); ptr = &((octv_t *)loct->nod.buf)[loct->head];
    x = 0; y = 0; z = 0; j = 8-1; ord = (ipz < z+s)*4 + (ipy < y+s)*2 + (ipx < x+s);
    while (1)
    {
        k = (j^ord);

        i = (1<<k); if (!(ptr->chi&i)) goto tosibly;

        nx = (( k    &1)<<ls)+x;
        ny = (((k>>1)&1)<<ls)+y;
        nz = (((k>>2)&1)<<ls)+z;

        ncr2 = dist2cube2_nearest((double)nx-p->x,(double)ny-p->y,(double)nz-p->z,(double)s);
        if (ncr2 >= cr2) goto tosibly;

        if (ls <= 0)
        {
            cr2 = ncr2;
            if (hit) { hit->x = nx; hit->y = ny; hit->z = nz; }
            if (hitsurf) (*hitsurf) = &((surf_t *)loct->sur.buf)[popcount8(ptr->chi&(i-1)) + ptr->ind];
            goto tosibly;
        }

        if (!(ptr->sol&i)) cr2 = std::min(cr2,dist2cube2_farthest((double)nx-p->x,(double)ny-p->y,(double)nz-p->z,(double)s));

        stk[ls].ptr = ptr; stk[ls].x = x; stk[ls].y = y; stk[ls].z = z; stk[ls].j = j; stk[ls].ord = ord; ls--; s >>= 1; //2child
        ptr = &((octv_t *)loct->nod.buf)[popcount8(ptr->chi&(i-1)) + ptr->ind]; x = nx; y = ny; z = nz; j = 8-1; ord = (ipz < z+s)*4 + (ipy < y+s)*2 + (ipx < x+s);
        continue;

tosibly:;
        j--; if (j >= 0) continue;
        do { ls++; s <<= 1; if (ls >= loct->lsid) return(sqrt(cr2)); j = stk[ls].j-1; } while (j < 0); //2parent
        ptr = stk[ls].ptr; x = stk[ls].x; y = stk[ls].y; z = stk[ls].z; ord = stk[ls].ord;
    }
}
double oct_balloonrad (oct_t *loct, point3f *fp, double cr, point3i *hit, surf_t **hitsurf)
{
    point3d dp;
    dp.x = fp->x; dp.y = fp->y; dp.z = fp->z;
    return(oct_balloonrad(loct,&dp,cr,hit,hitsurf));
}

//--------------------------------------------------------------------------------------------------

int oct_slidemove (oct_t *loct, point3d *p, point3d *padd, double fat, point3d *pgoal)
{
    point3d fp, fa, hit[3], hnorm[3];
    point3i lp;
    surf_t *psurf;
    double d;
    int i;

    fp = *p; fa = *padd; fat = oct_balloonrad(loct,&fp,fat,&lp,&psurf);
    for(i=0;1;i++)
    {
        fat -= 1e-7; if (!oct_sphtrace(loct,&fp,&fa,fat,0,&hit[i],&hnorm[i])) { *pgoal = hit[i]; return(i); }
        if (i == 2) break;

        fa.x += fp.x-hit[i].x;
        fa.y += fp.y-hit[i].y;
        fa.z += fp.z-hit[i].z;
        fp = hit[i];

        if (!i)
        {
            d = fa.x*hnorm[0].x + fa.y*hnorm[0].y + fa.z*hnorm[0].z;
            fa.x -= hnorm[0].x*d;
            fa.y -= hnorm[0].y*d;
            fa.z -= hnorm[0].z*d;
        }
        else
        {
            hnorm[2].x = hnorm[0].y*hnorm[1].z - hnorm[0].z*hnorm[1].y;
            hnorm[2].y = hnorm[0].z*hnorm[1].x - hnorm[0].x*hnorm[1].z;
            hnorm[2].z = hnorm[0].x*hnorm[1].y - hnorm[0].y*hnorm[1].x;
            d = hnorm[2].x*hnorm[2].x + hnorm[2].y*hnorm[2].y + hnorm[2].z*hnorm[2].z; if (d == 0.0) return(1);
            d = (fa.x*hnorm[2].x + fa.y*hnorm[2].y + fa.z*hnorm[2].z)/d;
            fa.x = hnorm[2].x*d;
            fa.y = hnorm[2].y*d;
            fa.z = hnorm[2].z*d;
        }
    }
    *pgoal = hit[2]; return(3);
}
int oct_slidemove (oct_t *loct, point3f *fp, point3f *fpadd, double fat, point3f *fpgoal)
{
    point3d dp, dpadd, dpgoal;
    int retval;
    dp.x = fp->x; dp.y = fp->y; dp.z = fp->z;
    dpadd.x = fpadd->x; dpadd.y = fpadd->y; dpadd.z = fpadd->z;
    retval = oct_slidemove(loct,&dp,&dpadd,fat,&dpgoal);
    fpgoal->x = dpgoal.x; fpgoal->y = dpgoal.y; fpgoal->z = dpgoal.z;
    return(retval);
}
