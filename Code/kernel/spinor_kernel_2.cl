/// @file     spinor_kernel_2.cl
/// @author   Erik ZORZIN
/// @date     16JAN2021
/// @brief    2nd kernel.
/// @details  Computes Verlet's integration.
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
                        __global int*       frontier,                                 // Spacetime frontier.
                        __global int*       frontier_num,                             // Spacetime frontier cells number.
                        __global float4*    frontier_pos,                             // Spacetime frontier cells posistion.
                        __global float*     dispersion,                               // Dispersion fraction.
                        __global float*     dt_simulation                             // Simulation time step.
                        )   
{
  //////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////// INDEXES /////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////
  unsigned int i = get_global_id(0);                                                  // Global index [#].
  unsigned int j = 0;                                                                 // Neighbour stride index.
  unsigned int j_min = 0;                                                             // Neighbour stride minimun index.
  unsigned int j_max = offset[i];                                                     // Neighbour stride maximum index.
  unsigned int k = 0;                                                                 // Neighbour tuple index.
  unsigned int n = central[j_max - 1];                                                // Central node index.

  //////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////// CELL VARIABLES /////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////
  float3        p_int             = position_int[n].xyz;                              // Central node position (intermediate). 
  float3        mate              = (float3)(0.0f, 0.0f, 0.0f);                       // Neighbour node position.
  float3        link              = (float3)(0.0f, 0.0f, 0.0f);                       // Neighbour link.
  float         L                 = 0.0f;                                             // Neighbour link length.
  float         R                 = 0.0f;                                             // Neighbour link resting length.
  float         S                 = 0.0f;                                             // Neighbour link strain.
  float         K                 = 0.0f;                                             // Neighbour link stiffness.
  float         D                 = dispersion[0];                                    // Dispersion.
  float         Fspring           = 0.0f;                                             // Spring force (scalar).  
  float         Jacc              = 0.0f;                                             // Central node radiated energy.
  float         b                 = 0.0f;                                             // Number of 1st + 2nd nearest neighbours.

  // COMPUTING STRIDE MINIMUM INDEX:
  if (i == 0)
  {
    j_min = 0;                                                                        // Setting stride minimum (first stride)...
  }
  else
  {
    j_min = offset[i - 1];                                                            // Setting stride minimum (all others)...
  }

  // COMPUTING ELASTIC FORCE:
  for (j = j_min; j < j_max; j++)
  {
    k = neighbour[j];                                                                 // Computing neighbour index...
    mate = position_int[k].xyz;                                                       // Getting neighbour position...
    link = p_int - mate;                                                              // Getting neighbour link vector...
    L = length(link);                                                                 // Computing neighbour link length...
    R = resting[j];                                                                   // Getting neighbour link resting length...
    S = L - R;                                                                        // Computing neighbour link strain...
    K = stiffness[j];                                                                 // Getting neighbour link stiffness...
    Fspring = -K*S;                                                                   // Computing elastic force on central node (as scalar)...

    if(K > FLT_EPSILON)
    {
      Jacc += 0.5f*D*Fspring*R;                                                       // Computing radiated energy from central node...
      b += 1.0f;                                                                      // Counting 1st and 2nd nearest neighbours around the central node...
    }
  }

  position_int[n].w = Jacc;                                                           // Accumulating central node radiative energy...
  velocity_int[n].w = b;                                                              // Setting number of 1st + 2nd nearest neighbours...
}