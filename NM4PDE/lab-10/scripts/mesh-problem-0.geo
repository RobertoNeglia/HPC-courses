N = 20;
h = 1.0 / N;

Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 2, 0, h};
Point(4) = {0, 2, 0, h};

Line(1)   = {1, 2};
Line(2)   = {2, 3};
Line(3)   = {3, 4};
Line(4)   = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Set tags to the boundaries. //////////////////////////////////////////////////

Physical Curve(0) = {4};
Physical Curve(1) = {2};
Physical Curve(2) = {1};
Physical Curve(3) = {3};

Physical Surface(10) = {1};

// Generate mesh ////////////////////////////////////////////////////////////////

Mesh 2;
Save "../mesh/mesh-problem-0.msh";