Point(1) = {0.0, 0.0, 0};
Point(2) = {4.2, 0.0, 0};
Point(3) = {4.2, 1.3, 0};
Point(4) = {0.0, 1.3, 0};
Line(1) = {3, 4};
Line(2) = {4, 1};
Line(3) = {1, 2};
Line(4) = {2, 3};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
