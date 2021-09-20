/// @file     spinor_kernel_1.cl
/// @author   Erik ZORZIN
/// @date     16JAN2021
/// @brief    1st kernel.
/// @details  Applies freedom contraints, computes new position, uptades intermediate position.
__kernel void thekernel(__global float4*    position,                                 // vec4(position.xyz [m], freedom []).
                        __global float4*    velocity,                                 // vec4(velocity.xyz [m/s], friction [N*s/m]).
                        __global float4*    velocity_int,                             // vec4(velocity (intermediate) [m/s], number of 1st + 2nd nearest neighbours []).
                        __global float4*    velocity_est,                             // vec4(velocity.xyz (estimation) [m/s], radiative energy [J]).
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
                        __global int*       frontier,                                 // Spacetime frontier.
                        __global int*       frontier_num,                             // Spacetime frontier cells number.
                        __global float4*    frontier_pos,                             // Spacetime frontier cells posistion.
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
  float3        p                 = adjzero3(position[i].xyz);                        // Central node position.
  float3        v                 = adjzero3(velocity[i].xyz);                        // Central node velocity.
  float3        a                 = adjzero3(acceleration[i].xyz);                    // Central node acceleration.
  float3        p_new             = (float3)(0.0f, 0.0f, 0.0f);                       // Central node position (new). 
  float3        v_int             = (float3)(0.0f, 0.0f, 0.0f);                       // Central node velocity (intermediate). 
  float         fr                = adjzero(position[i].w);                           // Central node freedom flag.
  float         dt                = adjzero(dt_simulation[0]);                        // Simulation time step.
  int           s_num             = spinor_num[0];                                    // Spinor cells number.
  int           f_num             = frontier_num[0];                                  // Spacetime frontier cells number.

  // APPLYING FREEDOM CONSTRAINTS:
  if (fr < FLT_EPSILON)
  {
    v = (float3)(0.0f, 0.0f, 0.0f);                                                   // Constraining velocity...
    a = (float3)(0.0f, 0.0f, 0.0f);                                                   // Constraining acceleration...
  }

  // FINDING SPINOR:
  for(j = 0; j < s_num; j++)
  {
    if(i == spinor[j])
    {
      p = spinor_pos[j].xyz;
    }
  }

  // FINDING FRONTIER:
  for(j = 0; j < f_num; j++)
  {
    if(i == frontier[j])
    {
      p = frontier_pos[j].xyz;
    }
  }

  // COMPUTING NEW POSITION:
  p_new = p + mulzero3(dt, v) + mulzero3(0.5f, mulzero3(pownzero(dt, 2), a));         // Computing new position...
  v_int = v + mulzero3(dt, a);                                                        // Computing intermediate velocity...

  // UPDATING KINEMATICS:
  position[i].xyz = p_new;                                                            // Updating new position...
  velocity_int[i].xyz = v_int;                                                        // Updating intermediate velocity...          
}