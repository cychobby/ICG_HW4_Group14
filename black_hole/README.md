# ICG 2025 HW4 - Black Hole (Group14)

Black hole simulation project for ICG 2025 HW4. This repo includes a 3D scene
with GPU geodesic ray tracing and a 2D lensing demo.

## Features
- 3D geodesic ray tracing in `geodesic.comp` with an accretion disk and object hits.
- Warped spacetime grid overlay.
- Geometry-shader orbital trails (`orbit_particles.vert/.geom/.frag`).
- Meteor event with tidal stretching, debris, and a separate trail.
- 2D lensing demo (`2D_lensing.cpp`).

## Controls (BlackHole3D)
- LMB or MMB drag: orbit camera around the black hole.
- Mouse wheel: zoom in/out.
- RMB hold: enable gravity while held.
- G: toggle gravity on/off.
- M: spawn meteor.
- Shift+M: spawn fast meteor.

## Build Requirements
- C++ compiler supporting C++17 or newer
- CMake
- vcpkg
- Git
- OpenGL 4.3 capable GPU

## Build Instructions
1. Clone the repository:
   - `git clone https://github.com/kavan010/black_hole.git`
2. CD into the newly cloned directory:
   - `cd ./black_hole`
3. Install dependencies with vcpkg:
   - `vcpkg install`
4. Get the vcpkg CMake toolchain file path:
   - `vcpkg integrate install`
5. Configure the project:
   - `cmake -B build -S . -DCMAKE_TOOLCHAIN_FILE=/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake`
6. Build the project:
   - `cmake --build build`

## Run
- Windows: `build\\BlackHole3D.exe`, `build\\BlackHole2D.exe`
- Linux/macOS: `./build/BlackHole3D`, `./build/BlackHole2D`

Note: Shader files are copied next to the `BlackHole3D` executable after build.

### Alternative: Debian/Ubuntu apt workaround
If you do not want to use vcpkg, install these packages and run the same CMake
steps above:

```bash
sudo apt update
sudo apt install build-essential cmake \
  libglew-dev libglfw3-dev libglm-dev libgl1-mesa-dev
```