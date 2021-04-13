/// @file     main.cpp
/// @author   Erik ZORZIN
/// @date     12JAN2021
/// @brief    Single 1/2 spin spinor, simulated as a tangle of a 3D continuum body.

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
  #define SHADER_HOME "../Spinor/Code/shader/"                                                      // Linux OpenGL shaders directory.
  #define KERNEL_HOME "../Spinor/Code/kernel/"                                                      // Linux OpenCL kernels directory.
  #define GMSH_HOME   "../Spinor/Code/mesh/"                                                        // Linux GMSH mesh directory.
#endif

#ifdef WIN32
  #define SHADER_HOME "..\\..\\..\\Spinor\\Code\\shader\\"                                          // Windows OpenGL shaders directory.
  #define KERNEL_HOME "..\\..\\..\\Spinor\\Code\\kernel\\"                                          // Windows OpenCL kernels directory.
  #define GMSH_HOME   "..\\..\\..\\Spinor\\Code\\mesh\\"                                            // Linux GMSH mesh directory.
#endif

#define SHADER_VERT   "voxel_vertex.vert"                                                           // OpenGL vertex shader.
#define SHADER_GEOM   "voxel_geometry.geom"                                                         // OpenGL geometry shader.
#define SHADER_FRAG   "voxel_fragment.frag"                                                         // OpenGL fragment shader.
#define OVERLAY_VERT  "overlay_vertex.vert"                                                         // OpenGL vertex shader.
#define OVERLAY_GEOM  "overlay_geometry.geom"                                                       // OpenGL geometry shader.
#define OVERLAY_FRAG  "overlay_fragment.frag"                                                       // OpenGL fragment shader.
#define KERNEL_1      "spinor_kernel_1.cl"                                                          // OpenCL kernel source.
#define KERNEL_2      "spinor_kernel_2.cl"                                                          // OpenCL kernel source.
#define UTILITIES     "utilities.cl"                                                                // OpenCL utilities source.
#define MESH          "spinor.msh"                                                                  // GMSH mesh.

#define EPSILON       0.005f                                                                        // Float epsilon for mesh.

// INCLUDES:
#include "nu.hpp"                                                                                   // Neutrino header file.

int main ()
{
  // INDEXES:
  size_t             i;                                                                             // Index [#].
  size_t             j;                                                                             // Index [#].

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
  nu::shader*        S              = new nu::shader ();                                            // OpenGL shader program.
  nu::shader*        overlay        = new nu::shader ();                                            // OpenGL shader program.

  // OPENCL:
  nu::opencl*        cl             = new nu::opencl (NU_GPU);                                      // OpenCL context.
  nu::kernel*        K1             = new nu::kernel ();                                            // OpenCL kernel array.
  nu::kernel*        K2             = new nu::kernel ();                                            // OpenCL kernel array.
  nu::float4*        color          = new nu::float4 (0);                                           // Color [].
  nu::float4*        position       = new nu::float4 (1);                                           // Position [m].
  nu::float4*        velocity       = new nu::float4 (2);                                           // Velocity [m/s].
  nu::float4*        acceleration   = new nu::float4 (3);                                           // Acceleration [m/s^2].
  nu::float4*        position_int   = new nu::float4 (4);                                           // Position (intermediate) [m].
  nu::float4*        velocity_int   = new nu::float4 (5);                                           // Velocity (intermediate) [m/s].
  nu::float1*        stiffness      = new nu::float1 (6);                                           // Stiffness.
  nu::float1*        resting        = new nu::float1 (7);                                           // Resting.
  nu::float1*        friction       = new nu::float1 (8);                                           // Friction.
  nu::float1*        mass           = new nu::float1 (9);                                           // Mass.
  nu::int1*          central        = new nu::int1 (10);                                            // Central nodes.
  nu::int1*          neighbour      = new nu::int1 (11);                                            // Neighbour.
  nu::int1*          offset         = new nu::int1 (12);                                            // Offset.
  nu::int1*          freedom        = new nu::int1 (13);                                            // Freedom.
  nu::float1*        dt             = new nu::float1 (14);                                          // Time step [s].
  nu::int1*          particle       = new nu::int1 (15);                                            // Particle.
  nu::int1*          particle_num   = new nu::int1 (16);                                            // Particle number.
  nu::float4*        particle_pos   = new nu::float4 (17);                                          // Particle position.

  // MESH:
  nu::mesh*          spinor         = new nu::mesh (std::string (GMSH_HOME) + std::string (MESH));  // Mesh cloth.
  size_t             nodes;                                                                         // Number of nodes.
  size_t             elements;                                                                      // Number of elements.
  size_t             groups;                                                                        // Number of groups.
  size_t             neighbours;                                                                    // Number of neighbours.
  std::vector<GLint> point;                                                                         // Point on frame.
  size_t             point_nodes;                                                                   // Number of point nodes.
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
  float              ds             = 0.1f;                                                         // Cell size [m].
  float              m              = 0.01f;                                                        // Node mass [kg].
  float              K              = 100.0f;                                                       // Link elastic constant [kg/s^2].
  float              B              = 1.0f;                                                         // Damping [kg*s*m].
  float              R              = 10;                                                           // Particle's radius [#cells].
  float              dt_critical    = sqrt (m/K);                                                   // Critical time step [s].
  float              dt_simulation  = 0.2*dt_critical;                                              // Simulation time step [s].

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////// DATA INITIALIZATION ///////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // MESH:
  spinor->process (VOLUME, 3, NU_MSH_HEX_8);                                                        // Processing mesh...
  position->data  = spinor->node_coordinates;                                                       // Setting all node coordinates...
  neighbour->data = spinor->neighbour;                                                              // Setting neighbour indices...
  central->data   = spinor->neighbour_center;                                                       // Setting neighbour centers...
  offset->data    = spinor->neighbour_offset;                                                       // Setting neighbour offsets...
  resting->data   = spinor->neighbour_length;                                                       // Setting resting distances...
  nodes           = spinor->node.size ();                                                           // Getting the number of nodes...
  elements        = spinor->element.size ();                                                        // Getting the number of elements...
  groups          = spinor->group.size ();                                                          // Getting the number of groups...
  neighbours      = spinor->neighbour.size ();                                                      // Getting the number of neighbours...

  // SETTING NEUTRINO ARRAYS (parameters):
  friction->data.push_back (B);                                                                     // Setting friction...
  dt->data.push_back (dt_simulation);                                                               // Setting time step...

  // SETTING NEUTRINO ARRAYS ("nodes" depending):
  for(i = 0; i < nodes; i++)
  {
    position->data[i].x *= 1.0f;
    position->data[i].y *= 1.0f;
    position->data[i].z *= 1.0f;
    position_int->data.push_back (position->data[i]);                                               // Setting intermediate position...
    velocity->data.push_back ({0.0f, 0.0f, 0.0f, 1.0f});                                            // Setting velocity...
    velocity_int->data.push_back ({0.0f, 0.0f, 0.0f, 1.0f});                                        // Setting intermediate velocity...
    acceleration->data.push_back ({0.0f, 0.0f, 0.0f, 1.0f});                                        // Setting acceleration...
    mass->data.push_back (m);                                                                       // Setting mass...
    freedom->data.push_back (1);                                                                    // Setting freedom flag...

    // Finding particle:
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
             ) < (sqrt (3.0f)/2.0f)*R*ds + EPSILON)
      )
    {
      particle->data.push_back (i);                                                                 // Setting particle index...
      particle_pos->data.push_back (position->data[i]);                                             // Setting initial particle's position...
      freedom->data[i] = (0);                                                                       // Resetting freedom flag...
      std::cout << "Found particle: i = " << i << std::endl;
    }
  }

  particle_num->data.push_back (particle->data.size ());
  std::cout << "Particle num = " << particle_num->data[0] << std::endl;

  // SETTING NEUTRINO ARRAYS ("neighbours" depending):
  for(i = 0; i < neighbours; i++)
  {
    // Building 3D isotropic 18-node cubic MSM:
    if(resting->data[i] < ((sqrt (2.0f)/2.0f)*ds + EPSILON))
    {
      stiffness->data.push_back (K);                                                                // Setting link stiffness...
    }
    else
    {
      stiffness->data.push_back (K/200000.0f);                                                      // Setting link stiffness...
    }

    // Showing only [100] neighbours:
    if(resting->data[i] < (ds + EPSILON))
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
  boundary.push_back (ADHE);                                                                        // Setting boundary surface...
  boundary.push_back (BCGF);                                                                        // Setting boundary surface...
  boundary.push_back (ABFE);                                                                        // Setting boundary surface...
  boundary.push_back (DCGH);                                                                        // Setting boundary surface...

  for(i = 0; i < boundary.size (); i++)
  {
    spinor->process (boundary[i], 2, NU_MSH_PNT);                                                   // Processing mesh...
    point       = spinor->node;                                                                     // Getting nodes on border...
    point_nodes = point.size ();                                                                    // Getting the number of nodes on border...

    std::cout << "surface nodes = " << point_nodes << std::endl;

    for(j = 0; j < point_nodes; j++)
    {
      freedom->data[point[j]] = 0;                                                                  // Resetting freedom flag...
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////// OPENCL KERNELS INITIALIZATION /////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  K1->addsource (std::string (KERNEL_HOME) + std::string (UTILITIES));                              // Setting kernel source file...
  K1->addsource (std::string (KERNEL_HOME) + std::string (KERNEL_1));                               // Setting kernel source file...
  K1->build (nodes, 0, 0);                                                                          // Building kernel program...
  K2->addsource (std::string (KERNEL_HOME) + std::string (UTILITIES));                              // Setting kernel source file...
  K2->addsource (std::string (KERNEL_HOME) + std::string (KERNEL_2));                               // Setting kernel source file...
  K2->build (nodes, 0, 0);                                                                          // Building kernel program...

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////// OPENGL SHADERS INITIALIZATION /////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  S->addsource (std::string (SHADER_HOME) + std::string (SHADER_VERT), NU_VERTEX);                  // Setting shader source file...
  S->addsource (std::string (SHADER_HOME) + std::string (SHADER_GEOM), NU_GEOMETRY);                // Setting shader source file...
  S->addsource (std::string (SHADER_HOME) + std::string (SHADER_FRAG), NU_FRAGMENT);                // Setting shader source file...
  S->build (neighbours);                                                                            // Building shader program...
  overlay->addsource (std::string (SHADER_HOME) + std::string (OVERLAY_VERT), NU_VERTEX);           // Setting shader source file...
  overlay->addsource (std::string (SHADER_HOME) + std::string (OVERLAY_GEOM), NU_GEOMETRY);         // Setting shader source file...
  overlay->addsource (std::string (SHADER_HOME) + std::string (OVERLAY_FRAG), NU_FRAGMENT);         // Setting shader source file...
  overlay->build (nodes);                                                                           // Building shader program...

  std::cout << NU_NEAR_Z_CLIP << " " << NU_FAR_Z_CLIP << std::endl;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// SETTING OPENCL KERNEL ARGUMENTS /////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  cl->write ();

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////// APPLICATION LOOP ////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  while(!gl->closed ())                                                                             // Opening window...
  {
    cl->get_tic ();                                                                                 // Getting "tic" [us]...
    cl->write (17);                                                                                 // Writing particle's position...
    cl->acquire ();                                                                                 // Acquiring variables...
    cl->execute (K1, NU_WAIT);                                                                      // Executing OpenCL kernel...
    cl->execute (K2, NU_WAIT);                                                                      // Executing OpenCL kernel...
    cl->release ();                                                                                 // Releasing variables...

    gl->clear ();                                                                                   // Clearing gl...
    gl->poll_events ();                                                                             // Polling gl events...
    gl->mouse_navigation (ms_orbit_rate, ms_pan_rate, ms_decaytime);                                // Mouse navigation...
    gl->gamepad_navigation (gmp_orbit_rate, gmp_pan_rate, gmp_decaytime, gmp_deadzone);             // Gamepad navigation...
    gl->plot (S);                                                                                   // Plotting shared arguments...

    if(plot_overlay)
    {
      gl->plot (overlay);                                                                           // Plotting shared arguments...
    }

    gl->refresh ();                                                                                 // Refreshing gl...

    if(gl->button_DPAD_LEFT)
    {
      for(i = 0; i < particle_num->data[0]; i++)
      {
        px                      = particle_pos->data[i].x;
        py                      = particle_pos->data[i].y;

        px_new                  = +cos (0.01f)*px - sin (0.01f)*py;
        py_new                  = +sin (0.01f)*px + cos (0.01f)*py;

        particle_pos->data[i].x = px_new;
        particle_pos->data[i].y = py_new;
      }
    }

    if(gl->button_DPAD_RIGHT)
    {
      for(i = 0; i < particle_num->data[0]; i++)
      {
        px                      = particle_pos->data[i].x;
        py                      = particle_pos->data[i].y;

        px_new                  = +cos (0.01f)*px + sin (0.01f)*py;
        py_new                  = -sin (0.01f)*px + cos (0.01f)*py;

        particle_pos->data[i].x = px_new;
        particle_pos->data[i].y = py_new;
      }
    }

    if(gl->button_DPAD_DOWN)
    {
      for(i = 0; i < particle_num->data[0]; i++)
      {
        py                      = particle_pos->data[i].y;
        pz                      = particle_pos->data[i].z;

        py_new                  = +cos (0.01f)*py - sin (0.01f)*pz;
        pz_new                  = +sin (0.01f)*py + cos (0.01f)*pz;

        particle_pos->data[i].y = py_new;
        particle_pos->data[i].z = pz_new;
      }
    }

    if(gl->button_DPAD_UP)
    {
      for(i = 0; i < particle_num->data[0]; i++)
      {
        py                      = particle_pos->data[i].y;
        pz                      = particle_pos->data[i].z;

        py_new                  = +cos (0.01f)*py + sin (0.01f)*pz;
        pz_new                  = -sin (0.01f)*py + cos (0.01f)*pz;

        particle_pos->data[i].y = py_new;
        particle_pos->data[i].z = pz_new;
      }
    }

    if(gl->button_LEFT_BUMPER)
    {
      for(i = 0; i < particle_num->data[0]; i++)
      {
        px                      = particle_pos->data[i].x;
        py                      = particle_pos->data[i].y;
        pz                      = particle_pos->data[i].z;

        px_new                  = px*0.99f;
        py_new                  = py*0.99f;
        pz_new                  = pz*0.99f;

        particle_pos->data[i].x = px_new;
        particle_pos->data[i].y = py_new;
        particle_pos->data[i].z = pz_new;
      }
    }

    if(gl->button_RIGHT_BUMPER)
    {
      for(i = 0; i < particle_num->data[0]; i++)
      {
        px                      = particle_pos->data[i].x;
        py                      = particle_pos->data[i].y;
        pz                      = particle_pos->data[i].z;

        px_new                  = px/0.99f;
        py_new                  = py/0.99f;
        pz_new                  = pz/0.99f;

        particle_pos->data[i].x = px_new;
        particle_pos->data[i].y = py_new;
        particle_pos->data[i].z = pz_new;
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
  delete color;                                                                                     // Deleting color data...
  delete position;                                                                                  // Deleting position data...
  delete position_int;                                                                              // Deleting intermediate position data...
  delete velocity;                                                                                  // Deleting velocity data...
  delete velocity_int;                                                                              // Deleting intermediate velocity data...
  delete acceleration;                                                                              // Deleting acceleration data...
  delete stiffness;                                                                                 // Deleting stiffness data...
  delete resting;                                                                                   // Deleting resting data...
  delete friction;                                                                                  // Deleting friction data...
  delete mass;                                                                                      // Deleting mass...
  delete central;                                                                                   // Deleting central...
  delete neighbour;                                                                                 // Deleting neighbours...
  delete offset;                                                                                    // Deleting offset...
  delete freedom;                                                                                   // Deleting freedom flag data...
  delete dt;                                                                                        // Deleting time step data...
  delete particle;                                                                                  // Deleting particle...
  delete particle_pos;                                                                              // Deleting particle_pos...
  delete K1;                                                                                        // Deleting OpenCL kernel...
  delete K2;                                                                                        // Deleting OpenCL kernel...
  delete S;                                                                                         // Deleting OpenGL shader...
  delete overlay;                                                                                   // Deleting OpenGL shader...
  delete spinor;                                                                                    // Deleting spinor mesh...
  delete cl;                                                                                        // Deleting OpenCL...
  delete gl;                                                                                        // Deleting OpenGL...

  return 0;
}