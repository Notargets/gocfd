Point(1) = {-0.5, 0, 0, 0.05};
Point(2) = {0, 0, 0, 0.05};
Point(3) = {1.0, 0.2679491924311227, 0, 0.05};
Point(4) = {1.0, 1.5, 0, 0.05};
Point(5) = {-0.5, 1.5, 0, 0.2};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};
Curve Loop(1) = {1, 2, 3, 4, 5};
Plane Surface(1) = {1};
Physical Curve("wall") = {1, 2};
//+
Physical Surface("flowfield", 6) = {1};
//+
Show "*";
//+
Show "*";
//+
Show "*";
//+
Physical Curve("out", 7) = {3};
//+
Physical Curve("far", 8) = {4};
//+
Physical Curve("in", 9) = {5};
//+
Show "*";
