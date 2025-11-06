Assignment 4 - Ray Tracing - Dianne Cao
=====================================
To compile and run the code:
1. Navigate to the build directory: cd build
2. Run: cmake .. -DCMAKE_POLICY_VERSION_MINIMUM=3.5
3. Run: make
4. Run: ./assignment4 -input <scene_file> -output <output_file> [options]
   - Example: ./assignment4 -input scene01_plane.txt -output 01.png -size 200 200
   - Options:
     * -size <width> <height>: Set image dimensions
     * -bounces <N>: Maximum recursion depth for reflections (default: 0)
     * -shadows: Enable shadow ray visibility testing
     * -aa <N> or -spp <N>: Supersampling samples per dimension (N×N samples per pixel)
     * -aafilter <box|tent>: Anti-aliasing filter type (default: box)
     * -camera <perspective|spiral>: Camera type (default: perspective)
     * -twist <float>: Twist strength for spiral camera in radians (default: 0.5)

COLLABORATION
-------------
I discussed the assignment with people around me, but there was
no official collaboration. Quinne and Yifan both helped me answer questions when I
encountered bugs.


RESOURCES
----------
- Cursor AI: Used for code structure explanation, debugging compilation issues,
  implementing ray tracing primitives, shading models, and extra features
- ChatGPT: Used for debugging parts of my code.
- Course materials and lecture notes on ray tracing, Phong shading, and
  Monte Carlo integration
- GLOO framework documentation for rendering and scene management
- GLM library documentation for vector and matrix operations

KNOWN ISSUES
-------------
None known at this time.

EXTRA FEATURES
---------------
1. Super Sampling Anti-Aliasing
   - Implements pixel-level supersampling with configurable samples per dimension
   - Supports two filter types:
     * Box filter: Uniform weighting (equal weights for all samples)
     * Tent filter: Weighted average emphasizing central samples
   - Command-line control via -aa/-spp and -aafilter flags
   - Reduces aliasing artifacts and produces smoother edges
Examples created for all 7 image.

2. Glossy Reflection
   - Implements Monte Carlo integration for glossy reflections
   - Adds roughness parameter to Material class for controlling reflection spread
   - Uses per-pixel seeded random sampling for deterministic noise
   - Samples random directions within a cone around the perfect reflection direction
   - Supports smooth transitions from perfect mirror (roughness=0) to glossy materials
Effect applied to 06_glossy.png.

3. Normal Mapping
   - Implements tangent-space normal mapping for surface detail
   - Extends Material class with normal map support
   - Computes tangent frame (T, B, N) at triangle intersections
   - Interpolates texture coordinates using barycentric coordinates
   - Transforms tangent-space normals to world space using TBN matrix
   - Includes robust handling of degenerate cases (NaN, zero vectors)
Applied to cube in scene02_cube.txt using normal_mapping_normal_map.png


4. Spiral Camera Effect
   - Implements a new SpiralCamera class that rotates rays around the optical axis
   - Creates a "vortex" effect where rays are twisted based on horizontal position
   - Twist strength is proportional to pixel distance from center
   - Command-line control via -camera spiral and -twist <float> flags
   - Maintains proper FOV scaling and camera basis calculations
Applied to assignment4/01_spiral.png.
