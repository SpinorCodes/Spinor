/// @file     main.cpp
/// @author   Erik ZORZIN
/// @date     12JAN2021
/// @brief    Single 1/2 spinor.

#define INTEROP       true                                                                          // "true" = use OpenGL-OpenCL interoperability.
#define SX            800                                                                           // Window x-size [px].
#define SY            600                                                                           // Window y-size [px].
#define NAME          "Neutrino - Spinor"                                                           // Window name.
#define ORBX          0.0f                                                                          // x-axis orbit initial rotation.
#define ORBY          0.0f                                                                          // y-axis orbit initial rotation.
#define PANX          0.0f                                                                          // x-axis pan initial translation.
#define PANY          0.0f                                                                          // y-axis pan initial translation.
#define PANZ          -2.0f                                                                         // z-axis pan initial translation.

#ifdef __linux__
  #define SHADER_HOME "../../Code/shader/"                                                          // Linux OpenGL shaders directory.
  #define KERNEL_HOME "../../Code/kernel/"                                                          // Linux OpenCL kernels directory.
  #define GMSH_HOME   "../../Code/mesh/"                                                            // Linux GMSH mesh directory.
#endif

#ifdef WIN32
  #define SHADER_HOME "..\\..\\Code\\shader\\"                                                      // Windows OpenGL shaders directory.
  #define KERNEL_HOME "..\\..\\Code\\kernel\\"                                                      // Windows OpenCL kernels directory.
  #define GMSH_HOME   "..\\..\\Code\\mesh\\"                                                        // Linux GMSH mesh directory.
#endif

#define SHADER_VERT   "voxel_vertex.vert"                                                           // OpenGL vertex shader.
#define SHADER_GEOM   "voxel_geometry.geom"                                                         // OpenGL geometry shader.
#define SHADER_FRAG   "voxel_fragment.frag"                                                         // OpenGL fragment shader.
#define OVERLAY_VERT  "overlay_vertex.vert"                                                         // OpenGL vertex shader.
#define OVERLAY_GEOM  "overlay_geometry.geom"                                                       // OpenGL geometry shader.
#define OVERLAY_FRAG  "overlay_fragment.frag"                                                       // OpenGL fragment shader.
#define KERNEL_1      "spinor_kernel_1.cl"                                                          // OpenCL kernel source.
#define KERNEL_2      "spinor_kernel_2.cl"                                                          // OpenCL kernel source.
#define KERNEL_3      "spinor_kernel_3.cl"                                                          // OpenCL kernel source.
#define UTILITIES     "utilities.cl"                                                                // OpenCL utilities source.
#define MESH          "spacetime.msh"                                                               // GMSH mesh.

// INCLUDES:
#include "nu.hpp"                                                                                   // Neutrino header file.

int main ()
{
  // INDEXES:
  GLuint             i;                                                                             // Index [#].
  GLuint             j;                                                                             // Index [#].

  // MOUSE PARAMETERS:
  float              ms_orbit_rate  = 1.0f;                                                         // Orbit rotation rate [rev/s].
  float              ms_pan_rate    = 5.0f;                                                         // Pan translation rate [m/s].
  float              ms_decaytime   = 1.25f;                                                        // Pan LP filter decay time [s].

  // GAMEPAD PARAMETERS:
  float              gmp_orbit_rate = 1.0f;                                                         // Orbit angular rate coefficient [rev/s].
  float              gmp_pan_rate   = 1.0f;                                                         // Pan translation rate [m/s].
  float              gmp_decaytime  = 1.25f;                                                        // Low pass filter decay time [s].
  float              gmp_deadzone   = 0.30f;                                                        // Gamepad joystick deadzone [0...1].

  // OPENGL:
  nu::opengl*        gl             = new nu::opengl (NAME, SX, SY, ORBX, ORBY, PANX, PANY, PANZ);  // OpenGL context.
  nu::shader*        shader_1       = new nu::shader ();                                            // OpenGL shader program.
  nu::shader*        overlay        = new nu::shader ();                                            // OpenGL shader program.

  // OPENCL:
  nu::opencl*        cl             = new nu::opencl (NU_GPU);                                      // OpenCL context.
  nu::kernel*        kernel_1       = new nu::kernel ();                                            // OpenCL kernel array.
  nu::kernel*        kernel_2       = new nu::kernel ();                                            // OpenCL kernel array.
  nu::kernel*        kernel_3       = new nu::kernel ();                                            // OpenCL kernel array.

  nu::float4*        position       = new nu::float4 (0);                                           // vec4(position.xyz [m], freedom []).
  nu::float4*        position_int   = new nu::float4 (1);                                           // vec4(position (intermediate) [m], radiative energy [J]).
  nu::float4*        velocity       = new nu::float4 (2);                                           // vec4(velocity.xyz [m/s], friction [N*s/m]).
  nu::float4*        velocity_int   = new nu::float4 (3);                                           // Velocity (intermediate) [m/s].
  nu::float4*        acceleration   = new nu::float4 (4);                                           // vec4(acceleration.xyz [m/s^2], mass [kg]).

  nu::float4*        color          = new nu::float4 (5);                                           // vec4(color.xyz [], alpha []).
  nu::float1*        stiffness      = new nu::float1 (6);                                           // Stiffness.
  nu::float1*        resting        = new nu::float1 (7);                                           // Resting.
  nu::int1*          central        = new nu::int1 (8);                                             // Central nodes.
  nu::int1*          neighbour      = new nu::int1 (9);                                             // Neighbour.
  nu::int1*          offset         = new nu::int1 (10);                                            // Offset.

  nu::int1*          spinor         = new nu::int1 (11);                                            // Spinor.
  nu::int1*          spinor_num     = new nu::int1 (12);                                            // Spinor cells number.
  nu::float4*        spinor_pos     = new nu::float4 (13);                                          // Spinor cells position.
  nu::int1*          wall           = new nu::int1 (14);                                            // Wall.
  nu::int1*          wall_num       = new nu::int1 (15);                                            // Wall nodes number.
  nu::float4*        wall_pos       = new nu::float4 (16);                                          // Wall nodes position.

  nu::float1*        dispersion     = new nu::float1 (17);                                          // Dispersion fraction [-0.5...1.0].
  nu::float1*        dt             = new nu::float1 (18);                                          // Time step [s].

  // MESH:
  nu::mesh*          spacetime      = new nu::mesh (std::string (GMSH_HOME) + std::string (MESH));  // Spacetime mesh.
  size_t             nodes;                                                                         // Number of nodes.
  size_t             elements;                                                                      // Number of elements.
  size_t             groups;                                                                        // Number of groups.
  size_t             neighbours;                                                                    // Number of neighbours.
  std::vector<GLint> point;                                                                         // Point on frame.
  size_t             point_nodes;                                                                   // Number of point nodes.
  size_t             wall_nodes;                                                                    // Number of wall nodes.
  int                ABCD           = 13;                                                           // "ABCD" surface tag.
  int                EFGH           = 14;                                                           // "EFGH" surface tag.
  int                ADHE           = 15;                                                           // "ADHE" surface tag.
  int                BCGF           = 16;                                                           // "BCGF" surface tag.
  int                ABFE           = 17;                                                           // "ABFE" surface tag.
  int                DCGH           = 18;                                                           // "DCGH" surface tag.
  int                VOLUME         = 1;                                                            // Entire volume tag.
  std::vector<int>   boundary;                                                                      // Boundary array.
  float              px;
  float              py;
  float              pz;
  float              px_new;
  float              py_new;
  float              pz_new;

  // PLOT:
  bool               plot_overlay   = false;                                                        // Plot overlay flag.

  // SIMULATION VARIABLES:
  float              safety_CFL     = 0.2f;                                                         // Courant-Friedrichs-Lewy safety coefficient [].
  int                N              = 3;                                                            // Number of spatial dimensions of the MSM [].
  float              rho            = 200.0f;                                                       // Mass density [kg/m^3].
  float              E              = 100000.0f;                                                    // Young's modulus [Pa];
  float              nu             = -0.35f;                                                       // Poisson's ratio [];
  float              beta           = 10.0f;                                                        // Damping [kg*s*m].
  float              R              = 3;                                                            // Particle's radius [#cells].

  float              ds;                                                                            // Cell size [m].
  float              dV;                                                                            // Cell volume [m^3].
  float              k;                                                                             // Spring constant [N/m].
  float              K;                                                                             // Bulk modulus [Pa].
  float              c;                                                                             // Speed of pressure waves [m/s].
  float              dm;                                                                            // Node mass [kg].
  float              lambda;                                                                        // Lamé 1st parameter [Pa];
  float              mu;                                                                            // Lamé 2nd parameter [Pa];
  float              B;                                                                             // Dispersive pressure [Pa].
  float              Q;                                                                             // Dispersive to direct momentum flow ratio [].
  float              C;                                                                             // Interaction momentum carriers pressure [Pa].
  float              D;                                                                             // Dispersion fraction [-0.5...1.0].
  float              dt_CFL;                                                                        // Courant-Friedrichs-Lewy critical time step [s].
  float              dt_SIM;                                                                        // Simulation time step [s].

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////// DATA INITIALIZATION ///////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // MESH:
  spacetime->process (VOLUME, N, NU_MSH_HEX_8);                                                     // Processing mesh...
  position->data  = spacetime->node_coordinates;                                                    // Setting all node coordinates...
  neighbour->data = spacetime->neighbour;                                                           // Setting neighbour indices...
  central->data   = spacetime->neighbour_center;                                                    // Setting neighbour centers...
  offset->data    = spacetime->neighbour_offset;                                                    // Setting neighbour offsets...
  resting->data   = spacetime->neighbour_length;                                                    // Setting resting distances...
  nodes           = spacetime->node.size ();                                                        // Getting the number of nodes...
  elements        = spacetime->element.size ();                                                     // Getting the number of elements...
  groups          = spacetime->group.size ();                                                       // Getting the number of groups...
  neighbours      = spacetime->neighbour.size ();                                                   // Getting the number of neighbours...
  ds              = *std::min_element (std::begin (resting->data), std::end (resting->data));       // Getting cell size...

  dV              = pow (ds, N);                                                                    // Computing cell volume...
  dm              = rho*dV;                                                                         // Computing node mass...
  lambda          = (E*nu)/((nu + 1.0f)*(nu - N*nu + 1.0f));                                        // Computing 1st Lamé parameter...
  mu              = E/(2.0f*nu + 2.0f);                                                             // Computing 2nd Lamé parameter...
  B               = lambda - mu;                                                                    // Computing dispersive pressure...
  Q               = B/(mu*(1.0f + 2.0f/N));                                                         // Computing dispersive to direct momentum flow ratio...
  C               = mu + mu*abs (Q);                                                                // Computing interaction momentum carriers pressure...
  D               = Q/(1.0f + abs (Q));                                                             // Computing dispersion fraction...
  k               = 5.0f/(2.0f + 4.0f*sqrt (2.0f))*mu*dV/pow (ds, 2);                               // Computing spring constant...
  K               = E/(N + N*nu - (float)pow (N, 2)*nu);                                            // Computing bulk modulus...
  c               = sqrt (K/rho);                                                                   // Computing speed of pressure waves...
  dt_CFL          = 1.0f/(N*(c/ds));                                                                // Computing Courant-Friedrichs-Lewy critical time step [s]...
  dt_SIM          = safety_CFL*dt_CFL;                                                              // Setting simulation time step [s]...

  // SETTING NEUTRINO ARRAYS (parameters):
  dispersion->data.push_back (D);                                                                   // Setting dispersion fraction...
  dt->data.push_back (dt_SIM);                                                                      // Setting time step...

  // SETTING NEUTRINO ARRAYS ("nodes" depending):
  for(i = 0; i < nodes; i++)
  {
    position->data[i].w = 1.0f;                                                                     // Setting freedom flag...
    position_int->data.push_back (position->data[i]);                                               // Setting intermediate position...
    velocity->data.push_back ({0.0f, 0.0f, 0.0f, beta});                                            // Setting velocity...
    velocity_int->data.push_back ({0.0f, 0.0f, 0.0f, 1.0f});                                        // Setting intermediate velocity...
    acceleration->data.push_back ({0.0f, 0.0f, 0.0f, dm});                                          // Setting acceleration...

    // Finding spinor:
    if(
       (0 < sqrt (
                  pow (position->data[i].x, 2) +
                  pow (position->data[i].y, 2) +
                  pow (position->data[i].z, 2)
                 )) &&
       (sqrt (
              pow (position->data[i].x, 2) +
              pow (position->data[i].y, 2) +
              pow (position->data[i].z, 2)
             ) < (sqrt (3.0f)*ds*R + FLT_EPSILON)
       )
      )
    {
      spinor->data.push_back (i);                                                                   // Setting spinor index...
      spinor_pos->data.push_back (position->data[i]);                                               // Setting initial spinor's position...
      position->data[i].w = 1.0f;                                                                   // Resetting freedom flag... (EZOR 25APR2021: temporary set to 1)
    }
  }

  spinor_num->data.push_back (spinor->data.size ());

  // SETTING NEUTRINO ARRAYS ("neighbours" depending):
  for(i = 0; i < neighbours; i++)
  {
    // Building 3D isotropic 18-node cubic MSM:
    if(resting->data[i] < (sqrt (2.0f)*ds + FLT_EPSILON))
    {
      stiffness->data.push_back (k);                                                                // Setting 1st and 2nd nearest neighbour link stiffness...
    }
    else
    {
      stiffness->data.push_back (0.0f);                                                             // Setting 3rd nearest neighbour link stiffness...
    }

    // Showing only 1st neighbours:
    if(resting->data[i] < (ds + FLT_EPSILON))
    {
      color->data.push_back ({0.0f, 1.0f, 0.0f, 0.3f});                                             // Setting color...
    }
    else
    {
      color->data.push_back ({0.0f, 0.0f, 0.0f, 0.0f});                                             // Setting color...
    }
  }

  // SETTING MESH PHYSICAL CONSTRAINTS:
  boundary.push_back (ABCD);                                                                        // Setting boundary surface...
  boundary.push_back (EFGH);                                                                        // Setting boundary surface...
  //boundary.push_back (ADHE);                                                                        // Setting boundary surface...
  //boundary.push_back (BCGF);                                                                        // Setting boundary surface...
  //boundary.push_back (ABFE);                                                                        // Setting boundary surface...
  //boundary.push_back (DCGH);                                                                        // Setting boundary surface...

  for(i = 0; i < boundary.size (); i++)
  {
    spacetime->process (boundary[i], 2, NU_MSH_PNT);                                                // Processing mesh...
    point       = spacetime->node;                                                                  // Getting nodes on border...
    point_nodes = point.size ();                                                                    // Getting the number of nodes on border...

    for(j = 0; j < point_nodes; j++)
    {
      position->data[point[j]].w = 0.0f;                                                            // Resetting freedom flag...
    }

    point.clear ();
  }

  spacetime->process (boundary[0], 2, NU_MSH_PNT);                                                  // Processing mesh...
  wall->data = spacetime->node;                                                                     // Getting nodes on border...
  wall_nodes = wall->data.size ();                                                                  // Getting the number of nodes on border...

  for(j = 0; j < wall_nodes; j++)
  {
    wall->data.push_back (j);
    wall_pos->data.push_back (position->data[wall->data[j]]);
  }

  wall_num->data.push_back (wall_nodes);

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////// OPENCL KERNELS INITIALIZATION /////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  kernel_1->addsource (std::string (KERNEL_HOME) + std::string (UTILITIES));                        // Setting kernel source file...
  kernel_1->addsource (std::string (KERNEL_HOME) + std::string (KERNEL_1));                         // Setting kernel source file...
  kernel_1->build (nodes, 0, 0);                                                                    // Building kernel program...
  kernel_2->addsource (std::string (KERNEL_HOME) + std::string (UTILITIES));                        // Setting kernel source file...
  kernel_2->addsource (std::string (KERNEL_HOME) + std::string (KERNEL_2));                         // Setting kernel source file...
  kernel_2->build (nodes, 0, 0);                                                                    // Building kernel program...
  kernel_3->addsource (std::string (KERNEL_HOME) + std::string (UTILITIES));                        // Setting kernel source file...
  kernel_3->addsource (std::string (KERNEL_HOME) + std::string (KERNEL_3));                         // Setting kernel source file...
  kernel_3->build (nodes, 0, 0);                                                                    // Building kernel program...

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////// OPENGL SHADERS INITIALIZATION /////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  shader_1->addsource (std::string (SHADER_HOME) + std::string (SHADER_VERT), NU_VERTEX);           // Setting shader source file...
  shader_1->addsource (std::string (SHADER_HOME) + std::string (SHADER_GEOM), NU_GEOMETRY);         // Setting shader source file...
  shader_1->addsource (std::string (SHADER_HOME) + std::string (SHADER_FRAG), NU_FRAGMENT);         // Setting shader source file...
  shader_1->build (neighbours);                                                                     // Building shader program...
  overlay->addsource (std::string (SHADER_HOME) + std::string (OVERLAY_VERT), NU_VERTEX);           // Setting shader source file...
  overlay->addsource (std::string (SHADER_HOME) + std::string (OVERLAY_GEOM), NU_GEOMETRY);         // Setting shader source file...
  overlay->addsource (std::string (SHADER_HOME) + std::string (OVERLAY_FRAG), NU_FRAGMENT);         // Setting shader source file...
  overlay->build (nodes);                                                                           // Building shader program...

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// SETTING OPENCL KERNEL ARGUMENTS /////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  cl->write ();

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////// PRINTING SIMULATION PARAMETERS /////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << "============================== LATTICE PARAMETERS ==============================="; // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << "Lattice spatial dimension:                  N          = " << N << " []";           // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << "Cell size:                                  ds         = " << ds << " [m]";         // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << "Cell volume:                                dV         = " << dV << " [m^N]";       // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << "Node mass:                                  dm         = " << dm << " [kg]";        // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << "Spring constant:                            k          = " << k << " [N/m]";        // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << std::endl;                                                                           // Printing message...

  std::cout << "============================= MECHANICAL PARAMETERS ============================="; // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << "Young's modulus:                            E          = " << E << " [Pa]";         // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << "Poisson's ratio:                            nu         = " << nu << " []";          // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << "Lame's 1st parameter:                       lambda     = " << lambda << " [Pa]";    // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << "Lame's 2nd parameter:                       mu         = " << mu << " [Pa]";        // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << "Dispersive pressure:                        B          = " << B << " [Pa]";         // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << "Dispersive-to-direct momentum flow ratio:   Q          = " << Q << " []";           // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << "Interaction momentum carriers pressure:     C          = " << C << " [Pa]";         // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << "Dispersion fraction:                        D          = " << D << " []";           // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << std::endl;                                                                           // Printing message...

  std::cout << "============================== DYNAMICS PARAMETERS =============================="; // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << "Mass density:                               rho        = " << rho << " [kg/m^N]";   // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << "Bulk modulus:                               K          = " << K << " [Pa]";         // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << "Speed of pressure waves:                    c          = " << c << " [m/s]";        // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << "Courant-Friedrichs-Lewy critical time step: dt_CFL     = " << dt_CFL << " [s]";     // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << "Courant-Friedrichs-Lewy safety coefficient: safety_CFL = " << safety_CFL << " []";  // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << "Simulation time step:                       dt_sim     = " << dt_SIM << " [s]";     // Printing message...
  std::cout << std::endl;                                                                           // Printing message...
  std::cout << std::endl;                                                                           // Printing message...

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// APPLICATION LOOP ////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  while(!gl->closed ())                                                                             // Opening window...
  {
    cl->get_tic ();                                                                                 // Getting "tic" [us]...
    //cl->write (17);                                                                                 // Writing spinor's position...
    cl->write (16);                                                                                 // Writing wall position...
    cl->acquire ();                                                                                 // Acquiring variables...
    cl->execute (kernel_1, NU_WAIT);                                                                // Executing OpenCL kernel...
    cl->execute (kernel_2, NU_WAIT);                                                                // Executing OpenCL kernel...
    cl->execute (kernel_3, NU_WAIT);                                                                // Executing OpenCL kernel...
    cl->release ();                                                                                 // Releasing variables...

    gl->clear ();                                                                                   // Clearing gl...
    gl->poll_events ();                                                                             // Polling gl events...
    gl->mouse_navigation (ms_orbit_rate, ms_pan_rate, ms_decaytime);                                // Mouse navigation...
    gl->gamepad_navigation (gmp_orbit_rate, gmp_pan_rate, gmp_decaytime, gmp_deadzone);             // Gamepad navigation...
    gl->plot (shader_1);                                                                            // Plotting shared arguments...

    if(plot_overlay)
    {
      gl->plot (overlay);                                                                           // Plotting shared arguments...
    }

    gl->refresh ();                                                                                 // Refreshing gl...

    if(gl->button_DPAD_LEFT)
    {
      for(i = 0; i < spinor_num->data[0]; i++)
      {
        px                    = spinor_pos->data[i].x;
        py                    = spinor_pos->data[i].y;

        px_new                = +cos (0.01f)*px - sin (0.01f)*py;
        py_new                = +sin (0.01f)*px + cos (0.01f)*py;

        spinor_pos->data[i].x = px_new;
        spinor_pos->data[i].y = py_new;
      }
    }

    if(gl->button_DPAD_RIGHT)
    {
      for(i = 0; i < spinor_num->data[0]; i++)
      {
        px                    = spinor_pos->data[i].x;
        py                    = spinor_pos->data[i].y;

        px_new                = +cos (0.01f)*px + sin (0.01f)*py;
        py_new                = -sin (0.01f)*px + cos (0.01f)*py;

        spinor_pos->data[i].x = px_new;
        spinor_pos->data[i].y = py_new;
      }
    }

    if(gl->button_DPAD_DOWN)
    {
      for(i = 0; i < spinor_num->data[0]; i++)
      {
        py                    = spinor_pos->data[i].y;
        pz                    = spinor_pos->data[i].z;

        py_new                = +cos (0.01f)*py - sin (0.01f)*pz;
        pz_new                = +sin (0.01f)*py + cos (0.01f)*pz;

        spinor_pos->data[i].y = py_new;
        spinor_pos->data[i].z = pz_new;
      }
    }

    if(gl->button_DPAD_UP)
    {
      for(i = 0; i < spinor_num->data[0]; i++)
      {
        py                    = spinor_pos->data[i].y;
        pz                    = spinor_pos->data[i].z;

        py_new                = +cos (0.01f)*py + sin (0.01f)*pz;
        pz_new                = -sin (0.01f)*py + cos (0.01f)*pz;

        spinor_pos->data[i].y = py_new;
        spinor_pos->data[i].z = pz_new;
      }
    }

    if(gl->button_LEFT_BUMPER)
    {
      for(i = 0; i < spinor_num->data[0]; i++)
      {
        px                    = spinor_pos->data[i].x;
        py                    = spinor_pos->data[i].y;
        pz                    = spinor_pos->data[i].z;

        px_new                = px*0.99f;
        py_new                = py*0.99f;
        pz_new                = pz*0.99f;

        spinor_pos->data[i].x = px_new;
        spinor_pos->data[i].y = py_new;
        spinor_pos->data[i].z = pz_new;
      }
    }

    if(gl->button_RIGHT_BUMPER)
    {
      for(i = 0; i < spinor_num->data[0]; i++)
      {
        px                    = spinor_pos->data[i].x;
        py                    = spinor_pos->data[i].y;
        pz                    = spinor_pos->data[i].z;

        px_new                = px/0.99f;
        py_new                = py/0.99f;
        pz_new                = pz/0.99f;

        spinor_pos->data[i].x = px_new;
        spinor_pos->data[i].y = py_new;
        spinor_pos->data[i].z = pz_new;
      }
    }

    if(gl->button_SQUARE)
    {
      for(i = 0; i < wall_num->data[0]; i++)
      {
        pz                  = wall_pos->data[i].z;
        pz_new              = pz + ds/10.0f;
        wall_pos->data[i].z = pz_new;
      }
    }

    if(gl->button_CIRCLE)
    {
      for(i = 0; i < wall_num->data[0]; i++)
      {
        pz                  = wall_pos->data[i].z;
        pz_new              = pz - ds/10.0f;
        wall_pos->data[i].z = pz_new;
      }
    }

    if(gl->button_TRIANGLE)
    {
      plot_overlay = !plot_overlay;
    }

    if(gl->button_CROSS)
    {
      gl->close ();                                                                                 // Closing gl...
    }

    cl->get_toc ();                                                                                 // Getting "toc" [us]...
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////// CLEANUP ////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  delete cl;                                                                                        // Deleting OpenCL context...
  delete position;                                                                                  // Deleting position data...
  delete position_int;                                                                              // Deleting intermediate position data...
  delete velocity;                                                                                  // Deleting velocity data...
  delete velocity_int;                                                                              // Deleting intermediate velocity data...
  delete acceleration;                                                                              // Deleting acceleration data...
  delete color;                                                                                     // Deleting color data...
  delete stiffness;                                                                                 // Deleting stiffness data...
  delete resting;                                                                                   // Deleting resting data...
  delete central;                                                                                   // Deleting central...
  delete neighbour;                                                                                 // Deleting neighbours...
  delete offset;                                                                                    // Deleting offset...
  delete spinor;                                                                                    // Deleting spinor...
  delete spinor_num;                                                                                // Deleting spinor_num...
  delete spinor_pos;                                                                                // Deleting spinor_pos...
  delete wall;                                                                                      // Deleting wall...
  delete wall_num;                                                                                  // Deleting wall_num...
  delete wall_pos;                                                                                  // Deleting wall_pos...
  delete dt;                                                                                        // Deleting time step data...
  delete kernel_1;                                                                                  // Deleting OpenCL kernel...
  delete kernel_2;                                                                                  // Deleting OpenCL kernel...
  delete shader_1;                                                                                  // Deleting OpenGL shader...
  delete overlay;                                                                                   // Deleting OpenGL shader...
  delete spacetime;                                                                                 // Deleting spacetime mesh...
  delete cl;                                                                                        // Deleting OpenCL...
  delete gl;                                                                                        // Deleting OpenGL...

  return 0;
}