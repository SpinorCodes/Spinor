/// @file     spinor_kernel_3.cl
/// @author   Erik ZORZIN
/// @date     16JAN2021
/// @brief    3rd kernel.
/// @details  Computes
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
  float         freedom           = adjzero(position[n].w);                           // Central node freedom flag.
  float3        p_new             = adjzero3(position[n].xyz);                        // Central node position (new).
  float3        v                 = adjzero3(velocity[n].xyz);                        // Central node velocity.
  float3        v_est             = adjzero3(velocity_est[n].xyz);;                   // Central node velocity (estimation).
  float3        v_new             = (float3)(0.0f, 0.0f, 0.0f);                       // Central node velocity (new).
  float3        a                 = adjzero3(acceleration[n].xyz);                    // Central node acceleration.
  float3        a_new             = (float3)(0.0f, 0.0f, 0.0f);                       // Central node acceleration (new).
  float         m                 = adjzero(acceleration[n].w);                       // Central node mass.
  float3        Fe                = (float3)(0.0f, 0.0f, 0.0f);                       // Central node elastic force.  
  float3        F_new             = (float3)(0.0f, 0.0f, 0.0f);                       // Central node total force (new).
  int           b_central         = adjzero(velocity_int[n].w);                       // Number of 1st + 2nd nearest neighbours at central node.
  int           b_mate            = 0.0f;                                             // Number of 1st + 2nd nearest neighbours at neighbour node.
  float         beta              = adjzero(velocity[n].w);                           // Central node friction.
  float3        mate              = (float3)(0.0f, 0.0f, 0.0f);                       // Neighbour node position.
  float3        pace              = (float3)(0.0f, 0.0f, 0.0f);                       // Neighbour node velocity.
  float3        link              = (float3)(0.0f, 0.0f, 0.0f);                       // Neighbour link.
  float3        rate_est          = (float3)(0.0f, 0.0f, 0.0f);                       // Neighbour rate (estimation).
  float3        direction         = (float3)(0.0f, 0.0f, 0.0f);                       // Neighbour link direction.
  float         Fspring           = 0.0f;                                             // Spring force (scalar).
  float         Fdashpot_est      = 0.0f;                                             // Dashpot force (scalar, estimation).
  float3        Fviscous_est      = (float3)(0.0f, 0.0f, 0.0f);                       // Central node viscous force (estimation).
  float3        Fdirect           = (float3)(0.0f, 0.0f, 0.0f);                       // Central node direct force.
  float3        Fdissipative      = (float3)(0.0f, 0.0f, 0.0f);                       // Central node dissipative force.
  float         Jacc_central      = adjzero(velocity_est[n].w);                       // Central node radiated energy.
  float         Jacc_mate         = 0.0f;                                             // Neighbour node radiated energy.
  float         JC                = 0.0f;                                             // Radiated energy density (central).
  float         JN                = 0.0f;                                             // Radiated energy density (neighbour).
  float         R                 = 0.0f;                                             // Neighbour link resting length.
  float         K                 = 0.0f;                                             // Neighbour link stiffness.
  float         S                 = 0.0f;                                             // Neighbour link strain.
  float         L                 = 0.0f;                                             // Neighbour link length.
  float         V_est             = 0.0f;                                             // Neighbour rate strain (estimation).
  float         D                 = adjzero(dispersion[0]);                           // Dispersion.
  float         dt                = adjzero(dt_simulation[0]);                        // Simulation time step [s].

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
    mate = adjzero3(position[k].xyz);                                                 // Getting neighbour position...
    link = adjzero3(p_new - mate);                                                    // Computing neighbour link vector...
    L = adjzero(length(link));                                                        // Computing neighbour link length...
    direction = normzero3(link);                                                      // Computing neighbour link displacement vector...
    rate_est = adjzero3(velocity_est[k].xyz);                                         // Getting neighbour velocity (estimation)...
    V_est = adjzero(dot(adjzero3(v_est - rate_est), direction));                      // Computing neighbour rate (estimation)...
    Jacc_mate = adjzero(velocity_est[k].w);                                           // Radiant energy of neighbour node...
    R = adjzero(resting[j]);                                                          // Getting neighbour link resting length...
    S = adjzero(L - R);                                                               // Computing neighbour link strain...
    K = adjzero(stiffness[j]);                                                        // Getting neighbour link stiffness...
    Fspring = mulzero(K, -S);                                                         // Computing elastic force on central node (as scalar)...
    Fe = mulzero3(Fspring, direction);                                                // Computing elasting force on central node (as vector)...
    Fdirect += mulzero3(adjzero(1.0f - fabs(D)), Fe);                                 // Building up total elastic force upon central node...
    Fdashpot_est = mulzero(beta, -V_est);                                             // Computing dashpot force on central node (as scalar)...
    Fviscous_est += mulzero3(Fdashpot_est, direction);                                // Building up total viscous force upon central node...
    b_mate = adjzero(velocity_int[k].w);                                              // Getting number of 1st + 2nd nearest neighbours...

    if(K > FLT_EPSILON)
    {
      JC = mulzero(Jacc_central, recipzero(b_central));                               // Computing radiated energy density (central)...
      JN = mulzero(Jacc_mate, recipzero(b_mate));                                     // Computing radiated energy density (neighbour)...
      Fdissipative += mulzero3(mulzero(JC + JN, recipzero(R)), direction);            // Building up force from central node radiated energy...
    }
    
    if (color[j].w != 0.0f)
    {
      color[j].xyz = colormap(0.5f*(1.0f + S/R) - 0.1f);                              // Setting color...
    }
  }

  F_new = Fdirect + Fdissipative + Fviscous_est;                                      // Computing new total node force...
  a_new = mulzero3(recipzero(m), F_new);                                              // Computing new acceleration...
  
  // APPLYING FREEDOM CONSTRAINTS:
  if (freedom < FLT_EPSILON)
  {
    a_new = (float3)(0.0f, 0.0f, 0.0f);                                               // Constraining new acceleration...
  }

  // COMPUTING NEW VELOCITY:
  v_new = v + mulzero3(0.5f, mulzero3(dt, a + a_new));                                // Computing new velocity...

  // APPLYING FREEDOM CONSTRAINTS:
  if (freedom < FLT_EPSILON)
  {
    v_new = (float3)(0.0f, 0.0f, 0.0f);                                               // Constraining new velocity...
  }

  // UPDATING KINEMATICS:
  velocity[n].xyz = v_new;                                                            // Updating velocity [m/s]...
  acceleration[n].xyz = a_new;                                                        // Updating acceleration [m/s^2]...
}