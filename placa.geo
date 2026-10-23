// One quarter cylinder under internal pressure
ref=0.1;
Point(1) = {0, 0, 0, ref};
Point(2) = {1, 0, 0, ref};
Point(3) = {1.5, 0, 0, ref};
Point(4) = {0, 1.5, 0, ref};
Point(5) = {0, 1, 0, ref};
Line(1) = {2, 3};
Circle(2) = {3, 1, 4};
Line(3) = {4, 5};
Circle(4) = {5, 1, 2};
Line loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Physical Curve("fix1") = {1};
Physical Curve("fix2") = {3};
Physical Curve("loaded") = {4};
Physical Curve("free") = {2};
Physical Surface("surf1") = {1};