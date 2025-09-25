Assignment 1 - Spline Viewer - Dianne Cao
==========================================
To compile and run the code once I was in the folder:
1. Run: mkdir build && cd build
2. Run: cmake ..
3. Run: make
4.
Run: ./assignment1 ../assets/assignment1/curve.spline for the curve render
Run: ./assignment1 ../assets/assignment1/teapot.spline for the teapot render
Run: ./assignment1 ../assets/assignment1/teaspoon.spline for the spoon render
Run: ./assignment1 ../assets/assignment1/teacut.spline for the cup render


COLLABORATION
-------------
I asked for help from my classmate Yifan and he helped me a lot on the curve section of the code.
I was having a bug that I just could not fix where the curve was not updating when I pressed B and Z and
it was one step behind and Yifan helped me debug it.
I also used cursor to help me debug the curve section and worked with it to implement the material extra credit.

RESOURCE
----------
- Course materials and lecture notes
- http://holmes3d.net/graphics/teapot/ where I found the Bezier patches for the spoon and the cup
- Cursor AI


KNOWN PROBLEMS
--------------
- I tried to implement gradient material according to calculating curvature but failed, so I pivoted to calculating
the Z axis to color the models with different color matierials. Thats why in my code there is a lot of curvature calculation
that I didnt end up using. If I have more time I think I would have gone to office hour because I think trying to implement this
was a bit ambicious and I didn't test every step.
- I tried to implement a torus (donut shaped) spline model but also failed, but since I already have the spoon and the
cup so I left it there as well. I think if I had more time I would model the torus in blender and calculate the control points/
write code to calculate the patches to understand fully how the points are being calculated. 

EXTRA CREDIT
------------
I implemented an additional feature: multi-color teapot rendering where different parts of the teapot (lid, body, base, etc.) are colored based on their Z-coordinate position.
You can just run the teapot normally and the new materials will show up, I have also attached a screenshot!

FEEDBACK
-------------------
The end results of the pset was very fun, however it was very hard as someone that only has 6.1010 as background to deal with so much C++ right of the bat.
I struggled a lot in the curve part of the code. I think the assignment could have been more detailed in walking us through what to implement and then give us
creative freedom.
