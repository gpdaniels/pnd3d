# Major Changes of this Refactor

- Replaced window code with a platform generic version.
- Changed the default application to a stress test that adds/removes loads of voxels.
- Fixed loads of compiler warnings.
- Fixed a fair few issues spotted by the many compiler sanitisers.
- Removed message box logging, replaced with printfs.
- Got the code to build for 64 bit platforms.
- Removed bilinear filtering / polygon drawing / image loading code.
- Added C++ versions all needed ASM functions.
- Split apart files as much as possible.
- Fixed some issues with the GLSL code for non-textured rendering.
- Removed texturing support, only single colour per voxel supported.
- Fall back to single threaded on platforms other than windows.
- Removed all of the marching cubes sections.
- Removed loads of unneeded / not-running code.
- Got code building on a modern version of visual studio.

---

# Archived Original Release Information

>  Amalgamated from the original pnd3d readme and Ken Silverman's website: http://advsys.net/ken/voxlap/pnd3d.htm

##  PND3D Demo and Source Code

```
06/24/2014: First public release of PND3D engine demo
03/07/2018: First public release of PND3D source code
09/10/2018: Remove WINMAIN.CPP/SYSMAIN.H; fix mouse code for Windows 10
```

The PND3D engine is something I started, really in response to the hype of Unlimited Detail. Initially my goal was to figure out their algorithm. I started my first prototype in August 2008. It wasn't until December 2009 that I found a rendering algorithm I was happy with (and by that I mean a fast enough worst case). Then in 2011, Ace of Spades happened, and so I shifted my goal to making a successor to the Voxlap engine. Voxlap's main feature is the ability to support unlimited modification to voxels on the grid. Note that this design goal is in conflict to tricks like instancing, which can save greatly on memory.

I did a lot of stuff over the next year, then AoS changed hands and I lost interest in the project. I got frustrated over silly things like physics code and shifted into other projects. I figured today was a good day to end the silence.

### Some major improvements over Voxlap are:

- Texture mapping
- Support up to 4096^3 grids (higher than that works but with artifacts due to an optimization of structure size)
- Support multiple voxel grids; world and sprites are handled the same
- Fast, voxel-accurate collision detection
- More compact file format: *.KVO (Ken's Voxel Octree)

### Some less noticable improvements:

- Ability to quickly double or halve voxel grid size if necessary
- Ability to safely travel inside solid (Voxlap would crash if you tried!)
- Ability to render from any viewpoint (Voxlap required you to be inside or above map!)
- Fewer visual artifacts (no jagged edges; no artifacts when looking up or down; no mip-mapping required)
- Multithread support
- Some GPU acceleration
- Effects such as transparency and entire sprite shading offset
- Built-in drawsph() and drawcone() functions - handy for debugging
- Faster / more accurate normal estimation

### Some of my disappointments with the new engine are:

- Mod speed. This is mostly the fault of copying from CPU to GPU memory. For games that require less frequent, small changes, this is not an issue. Also note that it is not an issue when using /cpu mode. I suppose this will no longer be an issue once unified CPU&GPU memory becomes common.
- Physics. I was greedy - I wanted it all. I never figured out how to handle multiple simultaneous contacts correctly - not even a hack. Currently, objects get stuck and disappear when this happens. I found a paper about it and I admit the math is beyond me. The good news is things are finished on the voxel engine side - meaning I have a fully working function that detects exact cube voxel to voxel collision, which returns contact and normal vector.

### Rendering tricks:

I often get asked what makes PND3D's rendering so fast. You think you can do as good as me? Well, you better have a thorough understanding of assembly language, multithread techniques, and GPU shaders. I will do you a favor and summarize all the tricks used in PND3D:

- While not necessarily a trick, it should first be stated that PND3D uses a standard 8-way octree.
- Each non-leaf node is a 64-bit structure containing:
  - `char chi`: child byte mask: 1=has children - must recurse, 0=pure air or pure solid - don't recurse), and a solid byte
  - `char sol`: solid byte mask: 1=pure solid inside, 0=not pure solid inside
  - `char mrk, mrk2`: dummy bytes for structure alignment
  - `int ind`: an index pointer to the first child node. There can be anywhere from 0 to 8 children depending on the child mask.
 - During rendering, the octree is visited in front-to-back order, depth-first. This 'sorting' helps speed up rendering by allowing occlusion queries to skip large sections. Again, this is probably nothing new here, but it needs to be stated.
  - I calculate bounding rectangles on screen using a single 'divps' instruction during projection. Each cube's bounding box is determined by only 4 of the 8 vertices of a cube. These 4 vertices can be known in advance based on which rectangular tile of the screen it is inside.
  - The screen is divided into large rectangular tiles (typically around 16-48 tiles per screen). This is done for 2 reasons:

    1. To calculate bounding boxes of cubes more quickly. The brute force method to calculate a bounding box of a cube is to transform and project all 8 points and then find the minimum and maximum x and y coordinates. This can be sped up considerably if you happen to know which 4 of the 8 vertices actually contribute to an edge of the bounding box in advance. In fact, it turns out large sections of the screen have the same 4 vertices contributing to the same edges. Even better, these sections are axis-aligned rectangular tiles on the screen.

        These tiles are determined by the 3 'vanishing' points of the voxel grid (and on the front half of the screen space). To calculate a vanishing point, imagine trying to render a cube really far away on each axis, such as at voxel grid location: (0,0,-1e9), (0,+1e9,0), etc.. If you projected these to screen coordinates, you would get 6 sets of x/y coordinates. Simply cut the screen horizontally and vertically (whenever they happen to be within the bounds of the screen - most aren't), and you will have your screen nicely divided into areas where a cube has the same 4 vertices generating the same 4 edges of the bounding box.

    2. For multithread optimization. I make a few additional horizontal cuts (about 16). Each tile is rendered as an independent job, clipped to the tile like it was a viewport, and then rendered from the root of the octree. So because of these cuts / tiles, some voxels along the edge can actually be rendered multiple times, but it's well worth it.

- I do an occlusion query for every non-leaf node of the tree that contains some solid. If the projected bounding box is larger than 32 pixels of width, I exit early and simply assume it's visible and visit its children anyway. Note that there are not very many large or nearby octree nodes, but this does skip a lot of needless screen area when processing the query.

  For a bounding box of 32 pixels of width or less, I have a bit array (cover buffer) - 1 bit per pixel, where each bit determines whether that pixel has been rendered or not. Then using an 'and' instruction, I can test 32 pixels simultaneously. If the region is found to be fully occluded, I skip it and all its children of the octree.

- For leaf nodes only (i.e. solid surface voxels being rendered), I visit each pixel of the bounding rectangle and perform 6 dot products in the inner loop using 2 'addps' and 2 'movmskps' instructions (25% waste) to determine whether the pixel actually falls inside the projected cube. Note that a cube, when projected to screen, results in a 4- or 6- sided convex polygon.

- The above steps are all done on the CPU. A rendered pixel on the screen is not your typical 32-bit ARGB color. Instead, texture mapping is done during a second pass. I write a 64-bit pixel, which contains 12 bits x, 12 bits y, 12 bits z, and a 28-bit index to the leaf node which holds information about the surface voxel. The GPU is given the voxel position and a pointer to its surface structure. The shader basically raytraces to a cube which it already knows it will hit. The GPU determines what face of the cube with a few simple cases. The GPU then does the work of texture mapping with fancy filters like bilinear and mip mapping.
  
### Original Demo Download

[PND3D.ZIP](http://advsys.net/ken/voxlap/pnd3d.zip) (669,942 bytes, 09/10/2018) Includes demo and source code (mouse code now works on Windows 10!).

Note: The GLSL shaders have not been tested extensively on all video cards! I have an NVidia card so that tends to work best. If you have an old ATI or Intel based video card, you may need to override the default rendering mode with one of the following:
- `pnd3d /arbasm` is still full GPU acceleration and quality, but older style shader code. It may work in some cases.
- `pnd3d /cpu` is CPU-based emulation code. It should work on all machines. The CPU code is actually not too slow, but it does lack the highest quality texture-mapping (i.e. nearest+mipmap instead of trilinear filtering).

### Example Rendered Screenshots

|![PND3D ~tomland.png](https://raw.githubusercontent.com/gpdaniels/pnd3d/master/media/pnd3d00.jpg "PND3D ~tomland.png")|![PND3D /ls=12  //4096 voxels across](https://raw.githubusercontent.com/gpdaniels/pnd3d/master/media/pnd3d01.jpg "PND3D /ls=12  //4096 voxels across")|
|:--:|:--:|
|![PND3D /ils=4  //16^3 pig tanks (model from Duke 3D HRP)](https://raw.githubusercontent.com/gpdaniels/pnd3d/master/media/pnd3d02.jpg "PND3D /ils=4  //16^3 pig tanks (model from Duke 3D HRP)")|![PND3D untitled.vxl /ils=2  //4^3 1024x1024x256 voxel maps](https://raw.githubusercontent.com/gpdaniels/pnd3d/master/media/pnd3d03.jpg "PND3D untitled.vxl /ils=2  //4^3 1024x1024x256 voxel maps")|

## PND3D Engine non-commercial license:

```
[1] Any derivative works based on PND3D may be distributed as long as it is
    free of charge and through noncommercial means.

[2] You must give me proper credit. This line of text is sufficient:

       PND3D engine by Ken Silverman (http://advsys.net/ken)

    Make sure it is clearly visible somewhere in your archive.

[3] If you wish to release modified source code to your game, please add the
    following line to each source file changed:

   // This file has been modified from Ken Silverman's original release

[4] I am open to commercial applications based on PND3D, however you must
    consult with me first to acquire a commercial license. Using PND3D as a
    test platform or as an advertisement to another commercial game is
    commercial exploitation and prohibited without a commercial license.
```
