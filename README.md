# GSoC-18
Detailed Documentation of my work for BRL-CAD organisation under Google Summer of Code 2018.

# OpenCL GPGPU Raytracing
## Abstract
BRL-CAD objects are made from about two dozen primitives, ray-tracing of which can be quickly realised in a GPU owing to its parallel architecture. There are lot of function calls to see if the shooted ray is intersecting with the primitive and to check on which direction it should be reflected. These calls on different data can be effectively launched on GPU/multi-core CPU if ported into OpenCL to achieve higher speedups.

## Methodology
The project has two parts divided as follows:
* Accelerating primitives by porting the ray-primitive intersection (shot) code from C to OpenCL and provide the glue code which packs and unpacks that primitive from the C side to the OpenCL side.
* Developing BREP support involves serializing the ‘ON_BREP *’ of the primitives and passing it to the GPU and then porting the generic ray-BREP intersection code.

**Project Plan:**       https://brlcad.org/wiki/User:Sreyanshjainrkl/GSoC18/Project              
**Development Logs:**   https://https://brlcad.org/wiki/User:Sreyanshjainrkl/GSoC18/Log

## Patches
* SUPERELL:   https://sourceforge.net/p/brlcad/patches/477/
* HYP:        https://sourceforge.net/p/brlcad/patches/489/
* ARBN:       https://sourceforge.net/p/brlcad/patches/490/
* DSP:        https://sourceforge.net/p/brlcad/patches/492/
* VOL:        https://sourceforge.net/p/brlcad/patches/513/
* METABALL:
* PIPE:
* BREP:

