/// @file     spinor_kernel_3.cl
/// @author   Erik ZORZIN
/// @date     16JAN2021
/// @brief    3rd kernel.
/// @details  Computes
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
  float         freedom           = position[n].w;                                    // Central node freedom flag.
  float3        p_int             = position_int[n].xyz;                              // Central node position (intermediate).
  float3        v                 = velocity[n].xyz;                                  // Central node velocity.
  float         beta              = velocity[n].w;                                    // Central node friction.
  float3        v_int             = velocity_int[n].xyz;                              // Central node velocity (intermediate).
  float3        a                 = acceleration[n].xyz;                              // Central node acceleration.
  float         m                 = acceleration[n].w;                                // Central node mass.
  
  float3        p_new             = (float3)(0.0f, 0.0f, 0.0f);                       // Central node position (new).
  float3        v_new             = (float3)(0.0f, 0.0f, 0.0f);                       // Central node velocity (new).
  float3        a_new             = (float3)(0.0f, 0.0f, 0.0f);                       // Central node acceleration (new).
  float3        v_est             = (float3)(0.0f, 0.0f, 0.0f);                       // Central node velocity (estimation).
  float3        a_est             = (float3)(0.0f, 0.0f, 0.0f);                       // Central node acceleration (estimation).
  
  float3        Fe                = (float3)(0.0f, 0.0f, 0.0f);                       // Central node elastic force.  
  float3        Fv                = (float3)(0.0f, 0.0f, 0.0f);                       // Central node viscous force.
  float3        Fv_est            = (float3)(0.0f, 0.0f, 0.0f);                       // Central node viscous force (estimation).
  float3        F                 = (float3)(0.0f, 0.0f, 0.0f);                       // Central node total force.
  float3        F_new             = (float3)(0.0f, 0.0f, 0.0f);                       // Central node total force (new).
  int           b_central         = velocity_int[n].w;                                // Number of 1st + 2nd nearest neighbours at central node.
  int           b_mate            = 0.0f;                                             // Number of 1st + 2nd nearest neighbours at neighbour node.
  
  float3        mate              = (float3)(0.0f, 0.0f, 0.0f);                       // Neighbour node position.
  float3        link              = (float3)(0.0f, 0.0f, 0.0f);                       // Neighbour link.
  float3        direction         = (float3)(0.0f, 0.0f, 0.0f);                       // Neighbour link direction.
  float         Fspring           = 0.0f;                                             // Spring force (scalar).  
  float3        Fdirect           = (float3)(0.0f, 0.0f, 0.0f);                       // Central ndoe direct force.
  float3        Fdissipative      = (float3)(0.0f, 0.0f, 0.0f);                       // Central node dissipative force.
  float         Jacc_central      = position_int[n].w;                                // Central node radiated energy.
  float         Jacc_mate         = 0.0f;                                             // Neighbour node radiated energy.
  float         R                 = 0.0f;                                             // Neighbour link resting length.
  float         K                 = 0.0f;                                             // Neighbour link stiffness.
  float         S                 = 0.0f;                                             // Neighbour link strain.
  float         L                 = 0.0f;                                             // Neighbour link length.
  float         D                 = dispersion[0];                                    // Dispersion.
  float         dt                = dt_simulation[0];                                 // Simulation time step [s].

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
    Jacc_mate = position_int[k].w;                                                    // Radiant energy of neighbour node...
    link = p_int - mate;                                                              // Getting neighbour link vector...
    L = length(link);                                                                 // Computing neighbour link length...
    R = resting[j];                                                                   // Getting neighbour link resting length...
    S = L - R;                                                                        // Computing neighbour link strain...
    K = stiffness[j];                                                                 // Getting neighbour link stiffness...
    direction = normzero3(link);                                                      // Computing neighbour link displacement vector...
    Fspring = -K*S;                                                                   // Computing elastic force on central node (as scalar)...
    Fe = Fspring*direction;                                                           // Computing elasting force on central node (as vector)...
    Fdirect += (1.0f - fabs(D))*Fe;                                                   // Building up total elastinc force upon central node...
    b_mate = velocity_int[k].w;                                                       // Getting number of 1st + 2nd nearest neighbours...

    if(K > FLT_EPSILON)
    {
      Fdissipative += ((Jacc_central/b_central + Jacc_mate/b_mate)/R)*direction;      // Building up force from central node radiated energy...
    }
    
    if (color[j].w != 0.0f)
    {
      color[j].xyz = colormap(0.5f*(1.0f + S/R) - 0.1f);                              // Setting color...
    }
  }

  Fv = -beta*v_int;                                                                   // Computing node viscous force...
  F  = Fdirect + Fdissipative + Fv;                                                   // Computing node total force...
  a_est  = F/m;                                                                       // Computing new acceleration estimation...
  v_est = v + 0.5f*(a + a_est)*dt;                                                    // Computing new velocity estimation...
  Fv_est = -beta*v_est;                                                               // Computing new node viscous force estimation...
  F_new = Fdirect + Fdissipative + Fv_est;                                            // Computing new total node force...
  a_new = F_new/m;                                                                    // Computing acceleration...
  
  // APPLYING FREEDOM CONSTRAINTS:
  if (freedom == 0.0f)
  {
    a_new = (float3)(0.0f, 0.0f, 0.0f);                                               // Constraining acceleration...
  }

  // COMPUTING NEW VELOCITY:
  v_new = v + 0.5f*(a + a_new)*dt;                                                    // Computing velocity...

  // APPLYING FREEDOM CONSTRAINTS:
  if (freedom == 0.0f)
  {
    v_new = (float3)(0.0f, 0.0f, 0.0f);                                               // Constraining velocity...
  }

  // UPDATING KINEMATICS:
  position[n].xyz = p_int;                                                            // Updating position [m]...
  velocity[n].xyz = v_new;                                                            // Updating velocity [m/s]...
  acceleration[n].xyz = a_new;                                                        // Updating acceleration [m/s^2]...
}