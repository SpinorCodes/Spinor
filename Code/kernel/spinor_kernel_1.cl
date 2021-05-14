/// @file     spinor_kernel_1.cl
/// @author   Erik ZORZIN
/// @date     16JAN2021
/// @brief    1st kernel.
/// @details  Applies freedom contraints, computes new position, uptades intermediate position.
__kernel void thekernel(__global float4*    position,                                 // vec4(position.xyz [m], freedom []).
                        __global float4*    position_int,                             // vec4(position (intermediate) [m], radiative energy [J]).
                        __global float4*    velocity,                                 // vec4(velocity.xyz [m/s], friction [N*s/m]).
                        __global float4*    velocity_int,                             // vec4(velocity (intermediate) [m/s], number of 1st + 2nd nearest neighbours []).
                        __global float4*    acceleration,                             // vec4(acceleration.xyz [m/s^2], mass [kg]).
                        __global float4*    color,                                    // vec4(color.xyz [], alpha []).
                        __global float*     stiffness,                                // Stiffness.
                        __global float*     resting,                                  // Resting distance.
                        __global int*       central,                                  // Central.
                        __global int*       neighbour,                                // Neighbour.
                        __global int*       offset,                                   // Offset.
                        __global int*       spinor,                                   // Spinor.
                        __global int*       spinor_num,                               // Spinor cells number.
                        __global float4*    spinor_pos,                               // Spinor cells position.
                        __global int*       wall,                                     // Wall.
                        __global int*       wall_num,                                 // Wall cells number.
                        __global float4*    wall_pos,                                 // Wall cells posistion.
                        __global float*     dispersion,                               // Dispersion fraction.
                        __global float*     dt_simulation                             // Simulation time step.
                        )                                 
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
  float3        v                 = velocity[i].xyz;                                  // Central node velocity.
  float3        a                 = acceleration[i].xyz;                              // Central node acceleration.
  float3        p_new             = (float3)(0.0f, 0.0f, 0.0f);                       // Central node position. 
  float         fr                = position[i].w;                                    // Central node freedom flag.
  float         dt                = dt_simulation[0];                                 // Simulation time step.
  int           s_num             = spinor_num[0];                                    // Spinor cells number.
  int           w_num             = wall_num[0];                                      // Wall cells number.

  // APPLYING FREEDOM CONSTRAINTS:
  if (fr == 0)
  {
    v = (float3)(0.0f, 0.0f, 0.0f);                                                   // Constraining velocity...
    a = (float3)(0.0f, 0.0f, 0.0f);                                                   // Constraining acceleration...
  }

  for(j = 0; j < s_num; j++)
  {
    if(i == spinor[j])
    {
      p = spinor_pos[j].xyz;
    }
  }

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
  position_int[i].xyz = p_new;                                                        // Updating intermediate position...
  velocity_int[i].xyz = v + a*dt;                                                     // Updating intermediate velocity...          
}