Assignment 2 - Hierarchical Modeling and SSD - Dianne Cao
=====================================
To compile and run the code:
1. Navigate to the build directory: cd build
2. Run: cmake ..
3. Run: make
4. Run: ./assignment2

Initially, I didn't realize that I needed to recompile the code every time I made changes
I spent about 1 hour debugging this issue during office hours before realizing I was missing the compile step.

COLLABORATION
-------------
I discussed the assignment with people around me during review sessions, but there was
no official collaboration.

CONTROLS
--------
- 'S' key: Toggle between Skeleton mode and SSD (Skinned) mode
- '1' key: Switch control to Character 1 (left character at X=-1.0)
- '2' key: Switch control to Character 2 (right character at X=1.0)
- UI Sliders: Adjust joint rotations for the currently active character

RESOURCE
----------
- Cursor AI: Used to explain code structure, debug issues, and implement
  the dual-character system for extra credit
- Course materials and lecture notes
- GLFW documentation for key codes
- GLOO framework documentation

KNOWN PROBLEMS
--------------
1. upon first compile there is a weird cylinder at the head of the skeleton, but
the once you refresh it with the ui toggle the weird display is gone and everything works fine
but I couldn't debug the diaplay


EXTRA CREDIT
------------------------------------
I implemented a dual-character system with two independent skeletal characters displayed
side-by-side

HOW IT WORKS:
- Press '1' to control the left character (Character 1)
- Press '2' to control the right character (Character 2)
- When you switch characters, the UI sliders automatically update to show the selected
  character's current pose
- Each character remembers its own pose independently, so you can pose Character 1,
  switch to Character 2 to pose it differently, then switch back to Character 1 and
  it will still have its original pose
- The 'S' key toggles both characters between Skeleton and SSD mode simultaneously

FEEDBACK
-------------------
The assignment was super fun!
