//////////////////////////////////////////////////////////////////
//
//                             D--------C
//                            /|       /|
//                           H--------G |
//                           | A------|-B
//                           |/       |/
//                           E--------F
//
//
//      y
//      |
//      |
//      o -----x
//     /
//    /
//   z
//
//////////////////////////////////////////////////////////////////

dx = 0.1;                                                       // Setting side discretization length...
dy = 0.1;                                                       // Setting side discretization length...
dz = 0.1;                                                       // Setting side discretization length...

x_min = -1.0 - 0.5*dx;                                          // Setting "x_min"...
x_max = +1.0 + 0.5*dx;                                          // Setting "x_max"...
y_min = -1.0 - 0.5*dy;                                          // Setting "y_min"...
y_max = +1.0 + 0.5*dy;                                          // Setting "y_max"...
z_min = -1.0 - 0.5*dz;                                          // Setting "z_min"...
z_max = +1.0 + 0.5*dz;                                          // Setting "z_max"...

A = 1;                                                          // Point "A".
B = 2;                                                          // Point "B".
C = 3;                                                          // Point "C".
D = 4;                                                          // Point "D".
E = 5;                                                          // Point "E".
F = 6;                                                          // Point "F".
G = 7;                                                          // Point "G".
H = 8;                                                          // Point "H".

AB = 1;                                                         // Line "AB".
BC = 2;                                                         // Line "BC".
CD = 3;                                                         // Line "CD".
DA = 4;                                                         // Line "DA".
EF = 5;                                                         // Line "EF".
FG = 6;                                                         // Line "FG".
GH = 7;                                                         // Line "GH".
HE = 8;                                                         // Line "HE".
EA = 9;                                                         // Line "EA".
FB = 10;                                                        // Line "FB".
CG = 11;                                                        // Line "CG".
DH = 12;                                                        // Line "DH".

ABCD = 13;                                                      // Loop "ABCD".
EFGH = 14;                                                      // Loop "EFGH".
ADHE = 15;                                                      // Loop "ADHE".
BCGF = 16;                                                      // Loop "BCGF".
ABFE = 17;                                                      // Loop "ABFE".
DCGH = 18;                                                      // Loop "DCGH".

VOLUME = 1;                                                     // Entire volume.

Point(A) = {x_min, y_min, z_min};                           // Setting point "A"...
Point(B) = {x_max, y_min, z_min};                           // Setting point "B"...
Point(C) = {x_max, y_max, z_min};                           // Setting point "C"...
Point(D) = {x_min, y_max, z_min};                           // Setting point "D"...
Point(E) = {x_min, y_min, z_max};                           // Setting point "E"...
Point(F) = {x_max, y_min, z_max};                           // Setting point "F"...
Point(G) = {x_max, y_max, z_max};                           // Setting point "G"...
Point(H) = {x_min, y_max, z_max};                           // Setting point "H"...

Line(AB) = {A, B};                                              // Setting side "AB"...
Line(BC) = {B, C};                                              // Setting side "BC"...
Line(CD) = {C, D};                                              // Setting side "CD"...
Line(DA) = {D, A};                                              // Setting side "DA"...
Line(EF) = {E, F};                                              // Setting side "EF"...
Line(FG) = {F, G};                                              // Setting side "FG"...
Line(GH) = {G, H};                                              // Setting side "GH"...
Line(HE) = {H, E};                                              // Setting side "HE"...
Line(EA) = {E, A};                                              // Setting side "EA"...
Line(FB) = {F, B};                                              // Setting side "FB"...
Line(CG) = {C, G};                                              // Setting side "CG"...
Line(DH) = {D, H};                                              // Setting side "DH"...

Transfinite Curve {AB} = (x_max - x_min)/dx + 1;                // Applying transfinite algorithm...
Transfinite Curve {EF} = (x_max - x_min)/dx + 1;                // Applying transfinite algorithm...
Transfinite Curve {CD} = (x_max - x_min)/dx + 1;                // Applying transfinite algorithm...
Transfinite Curve {GH} = (x_max - x_min)/dx + 1;                // Applying transfinite algorithm...
Transfinite Curve {BC} = (y_max - y_min)/dy + 1;                // Applying transfinite algorithm...
Transfinite Curve {DA} = (y_max - y_min)/dy + 1;                // Applying transfinite algorithm...
Transfinite Curve {FG} = (y_max - y_min)/dy + 1;                // Applying transfinite algorithm...
Transfinite Curve {HE} = (y_max - y_min)/dy + 1;                // Applying transfinite algorithm...
Transfinite Curve {EA} = (z_max - z_min)/dz + 1;                // Applying transfinite algorithm...
Transfinite Curve {FB} = (z_max - z_min)/dz + 1;                // Applying transfinite algorithm...
Transfinite Curve {CG} = (z_max - z_min)/dz + 1;                // Applying transfinite algorithm...
Transfinite Curve {DH} = (z_max - z_min)/dz + 1;                // Applying transfinite algorithm...

Curve Loop(ABCD) = {AB, BC, CD, DA};                            // Setting perimeter "ABCD"...
Curve Loop(EFGH) = {EF, FG, GH, HE};                            // Setting perimeter "EFGH"...
Curve Loop(ADHE) = {-DA, DH, HE, EA};                           // Setting perimeter "ADHE"...
Curve Loop(BCGF) = {BC, CG, -FG, FB};                           // Setting perimeter "BCGF"...
Curve Loop(ABFE) = {AB, -FB, -EF, EA};                          // Setting perimeter "ABFE"...
Curve Loop(DCGH) = {-CD, CG, GH, -DH};                          // Setting perimeter "DCGH"...

Plane Surface(ABCD) = {ABCD};                                   // Setting surface "ABCD"...
Plane Surface(EFGH) = {EFGH};                                   // Setting surface "EFGH"...
Plane Surface(ADHE) = {ADHE};                                   // Setting surface "ADHE"...
Plane Surface(BCGF) = {BCGF};                                   // Setting surface "BCGF"...
Plane Surface(ABFE) = {ABFE};                                   // Setting surface "ABFE"...
Plane Surface(DCGH) = {DCGH};                                   // Setting surface "DCGH"...

Transfinite Surface{ABCD} = {A, B, C, D};                       // Applying transfinite algorithm...
Recombine Surface {ABCD};                                       // Recombining triangles into quadrangles...

Transfinite Surface {EFGH} = {E, F, G, H};                      // Applying transfinite algorithm...
Recombine Surface {EFGH};                                       // Recombining triangles into quadrangles...

Transfinite Surface {ADHE} = {A, D, H, E};                      // Applying transfinite algorithm...
Recombine Surface {ADHE};                                       // Recombining triangles into quadrangles...

Transfinite Surface {BCGF} = {B, C, G, F};                      // Applying transfinite algorithm...
Recombine Surface {BCGF};                                       // Recombining triangles into quadrangles...

Transfinite Surface {ABFE} = {A, B, F, E};                      // Applying transfinite algorithm...
Recombine Surface {ABFE};                                       // Recombining triangles into quadrangles...

Transfinite Surface {DCGH} = {D, C, G, H};                      // Applying transfinite algorithm...
Recombine Surface {DCGH};                                       // Recombining triangles into quadrangles...

out[] = Extrude {0.0 , 0.0, (z_max - z_min)}                    // Creating extrusion along z-axis...
{
  Surface{ABCD};                                                // Setting surface to be extruded...
  Layers{(z_max - z_min)/dz};
  //Layers{2};
  Recombine;
};

Physical Volume(VOLUME) = out[1];                               // Setting physical group...

Physical Surface(ABCD) = {ABCD};                                // Setting group: perimeter "ABCD"...
Physical Surface(EFGH) = {EFGH};                                // Setting group: perimeter "EFGH"...
Physical Surface(ADHE) = {ADHE};                                // Setting group: perimeter "ADHE"...
Physical Surface(BCGF) = {BCGF};                                // Setting group: perimeter "BCGF"...
Physical Surface(ABFE) = {ABFE};                                // Setting group: perimeter "ABFE"...
Physical Surface(DCGH) = {DCGH};                                // Setting group: perimeter "DCGH"...

Mesh 5;                                                         // Setting mesh type: quadrangles...

Mesh.SaveAll = 1;                                               // Saving all mesh nodes...