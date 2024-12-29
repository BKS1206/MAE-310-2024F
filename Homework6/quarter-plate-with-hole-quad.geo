R = 0.3;
L = 1.0;

Point(1) = {L, -L, 0};
Point(2) = {L, L, 0};
Point(3) = {-L, L, 0};
Point(4) = {-L, -L, 0};
Point(5) = {-L + R, -L, 0};
Point(6) = {-L, -L + R, 0};
Point(7) = {-L + Cos(Pi/4) * R, -L + Sin(Pi/4) * R, 0};

// Draw the points of the geometry. Points 5,6,7 corresponds to the bottom,top and the middle points of the small hole, while point 4 is the center.

Circle(1) = {5, 4, 7};
Circle(2) = {7, 4, 6};

//Draw two arcs to form the 1/4 of the small hole, based on point 4(the center), point 5,7(the bottom and top of the arc)

Line(3) = {6, 3};
Line(4) = {3, 2};
Line(5) = {2, 1};
Line(6) = {1, 5};
Line(7) = {2, 7};

// Draw the lines based on points, line 7 is the diagonal line

Curve Loop(1) = {4, 7, 2, 3};
Plane Surface(1) = {1};

// First create a curve loop to form a boundary of surface, the numbers in {} is the line numbers
// Then create surface based on the loop created here. Here is the upper part of the two surfaces

Curve Loop(2) = {7, -1, -6, -5};
Plane Surface(2) = {2};

// Same as previous codes. The negetive sigh implies the direction of loop is opposite the line's direction(line's direction is defined when connect the points)
// Here is the lower parts

Transfinite Line{1, 2, 3, 4, 5, 6, 7} = 3;

//Transfinite Line specifies sides which are applied with structured mesh, and parameters 3 specifies 3 grid points will be used when drawing mesh lines between the two defined faces

Transfinite Surface{1};
Transfinite Surface{2};

//Surface 1 and 2 are specified to be mesh with a structured mesh by the syntax "Transfinite Surface {i}"

Recombine Surface{1};
Recombine Surface{2};

// "Recombine surface" specify that the mesh will made up of quadrilaterials. Otherwise it will use structured triangles mesh.

Mesh.ElementOrder = 1;

// Set the element order of mesh is 1(linear elements). If we want to refine the mesh, we can use quadrilateral order (element order = 2) or cubic.

Mesh.Algorithm = 8;

// Define the algorithm to use to generate the mesh. 8 = DelQuad algorithm

// EOF
//+
Coherence;
