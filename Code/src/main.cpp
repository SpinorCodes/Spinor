/// @file     main.cpp
/// @author   Erik ZORZIN
/// @date     12JAN2021
/// @brief    Single 1/2 spin spinor, simulated as a tangle of a 3D continuum body.

#define INTEROP       true                                                                           // "true" = use OpenGL-OpenCL interoperability.
#define SX            800                                                                            // Window x-size [px].
#define SY            600                                                                            // Window y-size [px].
#define NAME          "Neutrino - Spinor"                                                            // Window name.
#define ORB_X         0.0f                                                                           // x-axis orbit initial rotation.
#define ORB_Y         0.0f                                                                           // y-axis orbit initial rotation.
#define PAN_X         0.0f                                                                           // x-axis pan initial translation.
#define PAN_Y         0.0f                                                                           // y-axis pan initial translation.
#define PAN_Z         -2.0f                                                                          // z-axis pan initial translation.

#ifdef __linux__
  #define SHADER_HOME "../Spinor/Code/shader/"                                                       // Linux OpenGL shaders directory.
  #define KERNEL_HOME "../Spinor/Code/kernel/"                                                       // Linux OpenCL kernels directory.
  #define GMSH_HOME   "../Spinor/Code/mesh/"                                                         // Linux GMSH mesh directory.
#endif

#ifdef __APPLE__
  #define SHADER_HOME "../Spinor/Code/shader/"                                                       // Mac OpenGL shaders directory.
  #define KERNEL_HOME "../Spinor/Code/kernel/"                                                       // Mac OpenCL kernels directory.
  #define GMSH_HOME   "../Spinor/Code/mesh/"                                                         // Linux GMSH mesh directory.
#endif

#ifdef WIN32
  #define SHADER_HOME "..\\..\\..\\Spinor\\Code\\shader\\"                                           // Windows OpenGL shaders directory.
  #define KERNEL_HOME "..\\..\\..\\Spinor\\Code\\kernel\\"                                           // Windows OpenCL kernels directory.
  #define GMSH_HOME   "..\\..\\..\\Spinor\\Code\\mesh\\"                                             // Linux GMSH mesh directory.
#endif

#define SHADER_VERT   "voxel_vertex.vert"                                                            // OpenGL vertex shader.
#define SHADER_GEOM_1 "voxel_geometry_1.geom"                                                        // OpenGL geometry shader.
#define SHADER_GEOM_2 "voxel_geometry_2.geom"                                                        // OpenGL geometry shader.
#define SHADER_FRAG   "voxel_fragment.frag"                                                          // OpenGL fragment shader.
#define KERNEL_1      "spinor_kernel_1.cl"                                                           // OpenCL kernel source.
#define KERNEL_2      "spinor_kernel_2.cl"                                                           // OpenCL kernel source.
#define GMSH_MESH     "spinor.msh"                                                                   // GMSH mesh.

// INCLUDES:
#include "nu.hpp"                                                                                    // Neutrino header file.

int main ()
{
  // INDEXES:
  size_t              i;                                                                             // Index [#].
  size_t              j;                                                                             // Index [#].
  size_t              j_min;                                                                         // Index [#].
  size_t              j_max;                                                                         // Index [#].

  // MOUSE PARAMETERS:
  float               ms_orbit_rate  = 1.0f;                                                         // Orbit rotation rate [rev/s].
  float               ms_pan_rate    = 5.0f;                                                         // Pan translation rate [m/s].
  float               ms_decaytime   = 1.25f;                                                        // Pan LP filter decay time [s].

  // GAMEPAD PARAMETERS:
  float               gmp_orbit_rate = 1.0f;                                                         // Orbit angular rate coefficient [rev/s].
  float               gmp_pan_rate   = 1.0f;                                                         // Pan translation rate [m/s].
  float               gmp_decaytime  = 1.25f;                                                        // Low pass filter decay time [s].
  float               gmp_deadzone   = 0.1f;                                                         // Gamepad joystick deadzone [0...1].

  // NEUTRINO:
  opengl*             gl             = new opengl (NAME, SX, SY, ORB_X, ORB_Y, PAN_X, PAN_Y, PAN_Z); // OpenGL context.
  opencl*             cl             = new opencl (NU_GPU);                                          // OpenCL context.
  shader*             S              = new shader ();                                                // OpenGL shader program.
  kernel*             K1             = new kernel ();                                                // OpenCL kernel array.
  kernel*             K2             = new kernel ();                                                // OpenCL kernel array.

  // KERNEL VARIABLES:
  nu_float4*          color          = new nu_float4 (0);                                            // Color [].
  nu_float4*          position       = new nu_float4 (1);                                            // Position [m].
  nu_float4*          velocity       = new nu_float4 (2);                                            // Velocity [m/s].
  nu_float4*          acceleration   = new nu_float4 (3);                                            // Acceleration [m/s^2].
  nu_float4*          position_int   = new nu_float4 (4);                                            // Position (intermediate) [m].
  nu_float4*          velocity_int   = new nu_float4 (5);                                            // Velocity (intermediate) [m/s].
  nu_float*           stiffness      = new nu_float (6);                                             // Stiffness.
  nu_float*           resting        = new nu_float (7);                                             // Resting.
  nu_float*           friction       = new nu_float (8);                                             // Friction.
  nu_int*             central        = new nu_int (9);                                               // Central nodes.
  nu_int*             neighbour      = new nu_int (10);                                              // Neighbour.
  nu_int*             offset         = new nu_int (11);                                              // Offset.
  nu_int*             freedom        = new nu_int (12);                                              // Freedom.
  nu_float*           dt             = new nu_float (13);                                            // Time step [s].

  // MESH:
  mesh*               spinor         = new mesh (std::string (GMSH_HOME) + std::string (GMSH_MESH)); // Mesh cloth.
  size_t              nodes;                                                                         // Number of nodes.
  size_t              elements;                                                                      // Number of elements.
  size_t              groups;                                                                        // Number of groups.
  size_t              neighbours;                                                                    // Number of neighbours.
  std::vector<size_t> frame;                                                                         // Nodes on frame.
  size_t              frame_nodes;                                                                   // Number of frame nodes.
  float               x_min          = -0.9f;                                                        // "x_min" spatial boundary [m].
  float               x_max          = +0.9f;                                                        // "x_max" spatial boundary [m].
  float               y_min          = -0.9f;                                                        // "y_min" spatial boundary [m].
  float               y_max          = +0.9f;                                                        // "y_max" spatial boundary [m].
  float               z_min          = -0.9f;                                                        // "z_min" spatial boundary [m].
  float               z_max          = +0.9f;                                                        // "z_max" spatial boundary [m].
  float               ds             = 0.9f;                                                         // Cell size [m].
  float               pos_x;                                                                         // Position "x" component...
  float               pos_y;                                                                         // Position "y" component...
  float               pos_z;                                                                         // Position "z" component...
  float               link_x;                                                                        // Link "x" component...
  float               link_y;                                                                        // Link "y" component...
  float               link_z;                                                                        // Link "z" component...
  float               link_length;                                                                   // Link length...

  // SIMULATION VARIABLES:
  float               m              = 0.001f;                                                       // Node mass [kg].
  float               K              = 100.0f;                                                       // Link elastic constant [kg/s^2].
  float               B              = 1.0f;                                                         // Damping [kg*s*m].
  float               dt_critical;                                                                   // Critical time step [s].
  float               dt_simulation;                                                                 // Simulation time step [s].

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////// DATA INITIALIZATION ///////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // MESH:
  spinor->process (1, 3, NU_MSH_HEX_8);                                                              // Processing mesh...

  position->data  = spinor->node_coordinates;                                                        // Setting all node coordinates...
  neighbour->data = spinor->neighbour;                                                               // Setting neighbour indices...
  offset->data    = spinor->neighbour_offset;                                                        // Setting neighbour offsets...
  resting->data   = spinor->neighbour_length;                                                        // Setting resting distances...

  nodes           = spinor->node.size ();                                                            // Getting the number of nodes...
  elements        = spinor->element.size ();                                                         // Getting the number of elements...
  groups          = spinor->group.size ();                                                           // Getting the number of groups...
  neighbours      = spinor->neighbour.size ();                                                       // Getting the number of neighbours...

  std::cout << "nodes = " << nodes << std::endl;
  std::cout << "elements = " << elements << std::endl;
  std::cout << "groups = " << groups << std::endl;
  std::cout << "neighbours = " << neighbours << std::endl;
  std::cout << "offsets = " << spinor->neighbour_offset.size () << std::endl;
  std::cout << "lenghts = " << spinor->neighbour_length.size () << std::endl;
  std::cout << "links = " << spinor->neighbour_link.size () << std::endl;

  // Building central node vector:


  dt_critical   = sqrt (m/K);                                                                        // Critical time step [s].
  dt_simulation = 0.5f*dt_critical;                                                                  // Simulation time step [s].

  // SETTING NEUTRINO ARRAYS (parameters):
  friction->data.push_back (B);                                                                      // Setting friction...
  dt->data.push_back (dt_simulation);                                                                // Setting time step...

  // SETTING NEUTRINO ARRAYS ("nodes" depending):
  for(i = 0; i < nodes; i++)
  {
    position_int->data.push_back ({0.0f, 0.0f, 0.0f, 1.0f});                                         // Setting intermediate position...
    velocity->data.push_back ({0.0f, 0.0f, 0.0f, 1.0f});                                             // Setting velocity...
    velocity_int->data.push_back ({0.0f, 0.0f, 0.0f, 1.0f});                                         // Setting intermediate velocity...
    acceleration->data.push_back ({0.0f, 0.0f, 0.0f, 1.0f});                                         // Setting acceleration...
    freedom->data.push_back (1);                                                                     // Setting freedom flag...

    std::cout << "i = " << i << " offset = " << offset->data[i] <<
      " central = " << central->data[i] <<
      " x = " << position->data[central->data[i]].x <<
      " y = " << position->data[central->data[i]].y <<
      " z = " << position->data[central->data[i]].z <<
      " nodes: ";

    // Computing minimum element offset index:
    if(i == 0)
    {
      j_min = 0;                                                                                     // Setting minimum element offset index...
    }
    else
    {
      j_min = offset->data[i - 1];                                                                   // Setting minimum element offset index...
    }

    j_max = offset->data[i];                                                                         // Setting maximum element offset index...

    for(j = j_min; j < j_max; j++)
    {
      std::cout << neighbour->data[j] << " ";

      central->data.push_back (spinor->node[i]);                                                     // Building central node vector...
    }

    std::cout << std::endl;

    stiffness->data.push_back (K);                                                                   // Setting stiffness...

    if(resting->data[j] > 0.2)
    {
      color->data.push_back ({1.0f, 0.0f, 0.0f, 1.0f});                                              // Setting color...
    }
    else
    {
      color->data.push_back ({0.0f, 1.0f, 0.0f, 1.0f});                                              // Setting color...
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////// OPENCL KERNELS INITIALIZATION /////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  K1->addsource (std::string (KERNEL_HOME) + std::string (KERNEL_1));                                // Setting kernel source file...
  K1->build (nodes, 0, 0);                                                                           // Building kernel program...

  K2->addsource (std::string (KERNEL_HOME) + std::string (KERNEL_2));                                // Setting kernel source file...
  K2->build (nodes, 0, 0);                                                                           // Building kernel program...

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////// OPENGL SHADERS INITIALIZATION /////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  S->addsource (std::string (SHADER_HOME) + std::string (SHADER_VERT), NU_VERTEX);                   // Setting shader source file...
  S->addsource (std::string (SHADER_HOME) + std::string (SHADER_GEOM_1), NU_GEOMETRY);               // Setting shader source file...
  S->addsource (std::string (SHADER_HOME) + std::string (SHADER_GEOM_2), NU_GEOMETRY);               // Setting shader source file...
  S->addsource (std::string (SHADER_HOME) + std::string (SHADER_FRAG), NU_FRAGMENT);                 // Setting shader source file...
  S->build (neighbours);                                                                             // Building shader program...

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// SETTING OPENCL KERNEL ARGUMENTS /////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  cl->write ();

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// APPLICATION LOOP ////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  while(!gl->closed ())                                                                              // Opening window...
  {
    cl->get_tic ();                                                                                  // Getting "tic" [us]...
    cl->acquire ();
    cl->execute (K1, NU_WAIT);                                                                       // Executing OpenCL kernel...
    //cl->execute (K2, NU_WAIT);                                                                       // Executing OpenCL kernel...
    cl->release ();

    gl->clear ();                                                                                    // Clearing gl...
    gl->poll_events ();                                                                              // Polling gl events...
    gl->mouse_navigation (ms_orbit_rate, ms_pan_rate, ms_decaytime);
    gl->gamepad_navigation (gmp_orbit_rate, gmp_pan_rate, gmp_decaytime, gmp_deadzone);
    gl->plot (S);                                                                                    // Plotting shared arguments...

    gl->refresh ();                                                                                  // Refreshing gl...

    if(gl->button_CROSS)
    {
      gl->close ();                                                                                  // Closing gl...
    }

    cl->get_toc ();                                                                                  // Getting "toc" [us]...
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////// CLEANUP ////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  delete cl;                                                                                         // Deleting OpenCL context...
  delete color;                                                                                      // Deleting color data...
  delete position;                                                                                   // Deleting position data...
  delete position_int;                                                                               // Deleting intermediate position data...
  delete velocity;                                                                                   // Deleting velocity data...
  delete velocity_int;                                                                               // Deleting intermediate velocity data...
  delete acceleration;                                                                               // Deleting acceleration data...
  delete stiffness;                                                                                  // Deleting stiffness data...
  delete resting;                                                                                    // Deleting resting data...
  delete friction;                                                                                   // Deleting friction data...
  delete neighbour;                                                                                  // Deleting neighbours...
  delete offset;                                                                                     // Deleting offset...
  delete freedom;                                                                                    // Deleting freedom flag data...
  delete dt;                                                                                         // Deleting time step data...
  delete K1;                                                                                         // Deleting OpenCL kernel...
  delete K2;                                                                                         // Deleting OpenCL kernel...

  return 0;
}