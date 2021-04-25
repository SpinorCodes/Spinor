/// @file     spinor_kernel_1.cl
/// @author   Erik ZORZIN
/// @date     16JAN2021
/// @brief    Does nothing.
/// @details  The nothing:
__kernel void thekernel(__global float4*    color,                                    // Color.
                        __global float4*    position,                                 // Position.
                        __global float4*    position_int,                             // Position (intermediate).
                        __global float4*    velocity,                                 // Velocity.
                        __global float4*    velocity_int,                             // Velocity (intermediate).
                        __global float4*    acceleration,                             // Acceleration.
                        __global float*     stiffness,                                // Stiffness.
                        __global float*     resting,                                  // Resting distance.
                        __global float*     friction,                                 // Friction.
                        __global float*     mass,                                     // Mass.
                        __global int*       central,                                  // Central.
                        __global int*       neighbour,                                // Neighbour.
                        __global int*       offset,                                   // Offset.
                        __global float*     dt_simulation,                            // Simulation time step.
                        __global int*       particle,                                 // Particle.
                        __global int*       particle_num,                             // Particle number.
                        __global float4*    particle_pos,                             // Particle's position.
                        __global float*     momentum_ratio,                           // Dissipative to direct momentum flow ratio.
                        __global int*       wall,                                     // Particle.
                        __global int*       wall_num,                                 // Particle number.
                        __global float4*    wall_pos)                                 // Particle's position.
{
  //////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////// INDICES /////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////
  unsigned long i = get_global_id(0);                                                 // Global index [#].
  unsigned long j;  

  //////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////// CELL VARIABLES //////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////
  float3        p                 = position[i].xyz;                                  // Central node position.
  float3        v                 = velocity[i].xyz;                                      // Central node velocity.
  float3        a                 = acceleration[i].xyz;                                  // Central node acceleration.
  float3        p_new             = (float3)(0.0f, 0.0f, 0.0f);                 // Central node position. 
  float         fr                = position[i].w;                                       // Central node freedom flag.
  float         dt                = dt_simulation[0];                                 // Simulation time step [s].
  int           p_num             = particle_num[0];                                  // Particle number.
  int           w_num             = wall_num[0];                                      // Particle number.

  // APPLYING FREEDOM CONSTRAINTS:
  if (fr == 0)
  {
    v = (float3)(0.0f, 0.0f, 0.0f);                                             // Constraining velocity...
    a = (float3)(0.0f, 0.0f, 0.0f);                                             // Constraining acceleration...
  }
/*
  for(j = 0; j < p_num; j++)
  {
    if(i == particle[j])
    {
      p = particle_pos[j];
    }
  }
*/
  for(j = 0; j < w_num; j++)
  {
    if(i == wall[j])
    {
      p = wall_pos[j].xyz;
    }
  }

  // COMPUTING NEW POSITION:
  p_new = p + v*dt + 0.5f*a*dt*dt;                                                    // Computing Taylor's approximation...
        
  // UPDATING INTERMEDIATE POSITION:
  position_int[i].xyz = p_new;                                                            // Updating intermediate position...
  velocity_int[i].xyz = v + a*dt;                                                         // Updating intermediate velocity...           
}