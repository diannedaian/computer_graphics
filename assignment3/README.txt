Assignment 3 - Physically-Based Simulation - Dianne Cao
=====================================
To compile and run the code:
1. Navigate to the build directory: cd build
2. Run: cmake .. -DCMAKE_POLICY_VERSION_MINIMUM=3.5
3. Run: make
4. Run: ./assignment3 <integrator> <timestep>
   - Integrators: 'e' (Euler), 't' (Trapezoidal), 'r' (RK4)
   - Example: ./assignment3 r 0.005 (RK4 with 5ms timestep)

COLLABORATION
-------------
I discussed the assignment with people around me, but there was
no official collaboration.


CONTROLS
--------
- 'R' key: Reset simulation to initial state
- 'W' key: Toggle wind on/off (cloth simulation only)
- Mouse: ArcBall camera controls for viewing

RESOURCES
----------
- Cursor AI: Used for code structure explanation, debugging compilation issues,
  implementing integrators, and building the cloth simulation system
- Course materials and lecture notes on numerical integration
- GLOO framework documentation for rendering and scene management

KNOWN ISSUES
-------------
There isn't a known issue other than the wind doesn't effect the cloth that much until after it bounce.

EXTRA FEATURES
---------------
- Wind system with realistic varying forces
- Wind visualization with directional arrow

FEEDBACK
--------
It was very fun and sad I dont have time to implement the cloth render.
