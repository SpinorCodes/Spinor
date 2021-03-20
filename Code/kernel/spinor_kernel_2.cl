/// @file     spinor_kernel_2.cl
/// @author   Erik ZORZIN
/// @date     16JAN2021
/// @brief    Does nothing.
/// @details  The nothing:
__kernel void thekernel(__global float4*    color,                              // Color.
                        __global float4*    position,                           // Position.
                        __global float4*    position_int,                       // Position (intermediate).
                        __global float4*    velocity,                           // Velocity.
                        __global float4*    velocity_int,                       // Velocity (intermediate).
                        __global float4*    acceleration,                       // Acceleration.
                        __global float*     stiffness,                          // Stiffness.
                        __global float*     resting,                            // Resting distance.
                        __global float*     friction,                           // Friction.
                        __global float*     mass,                               // Mass.
                        __global int*       central,                            // Central.
                        __global int*       neighbour,                          // Neighbour.
                        __global int*       offset,                             // Offset.
                        __global int*       freedom,                            // Freedom flag.
                        __global float*     dt_simulation,                      // Simulation time step.
                        __global int*       particle,                           // Particle.
                        __global float4*    particle_pos)                       // Paritcle's position.
{
  ////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////// INDEXES ///////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  unsigned int i = get_global_id(0);                                            // Global index [#].
  unsigned int j = 0;                                                           // Neighbour stride index.
  unsigned int j_min = 0;                                                       // Neighbour stride minimun index.
  unsigned int j_max = offset[i];                                               // Neighbour stride maximum index.
  unsigned int k = 0;                                                           // Neighbour tuple index.
  unsigned int n = central[j_max - 1];                                          // Node index.

  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////// CELL VARIABLES //////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  float4        c;                                                              // Central node color.
  float4        v                 = velocity[n];                                // Central node velocity.
  float4        a                 = acceleration[n];                            // Central node acceleration.
  float4        p_int             = position_int[n];                            // Central node position (intermediate).
  float4        v_int             = velocity_int[n];                            // Central node velocity (intermediate).
  float4        p_new             = (float4)(0.0f, 0.0f, 0.0f, 1.0f);           // Central node position (new).
  float4        v_new             = (float4)(0.0f, 0.0f, 0.0f, 1.0f);           // Central node velocity (new).
  float4        a_new             = (float4)(0.0f, 0.0f, 0.0f, 1.0f);           // Central node acceleration (new).
  float4        v_est             = (float4)(0.0f, 0.0f, 0.0f, 1.0f);           // Central node velocity (estimation).
  float4        a_est             = (float4)(0.0f, 0.0f, 0.0f, 1.0f);           // Central node acceleration (estimation).
  float         m                 = mass[n];                                    // Central node mass.
  float         B                 = friction[0];                                // Central node friction.
  int           fr                = freedom[n];                                 // Central node freedom flag.
  float4        Fe                = (float4)(0.0f, 0.0f, 0.0f, 1.0f);           // Central node elastic force.  
  float4        Fv                = (float4)(0.0f, 0.0f, 0.0f, 1.0f);           // Central node viscous force.
  float4        Fv_est            = (float4)(0.0f, 0.0f, 0.0f, 1.0f);           // Central node viscous force (estimation).
  float4        F                 = (float4)(0.0f, 0.0f, 0.0f, 1.0f);           // Central node total force.
  float4        F_new             = (float4)(0.0f, 0.0f, 0.0f, 1.0f);           // Central node total force (new).
  float4        nearest           = (float4)(0.0f, 0.0f, 0.0f, 1.0f);           // Neighbour node position.
  float4        link              = (float4)(0.0f, 0.0f, 0.0f, 1.0f);           // Neighbour link.
  float4        D                 = (float4)(0.0f, 0.0f, 0.0f, 1.0f);           // Neighbour displacement.
  float         R                 = 0.0f;                                       // Neighbour link resting length.
  float         K                 = 0.0f;                                       // Neighbour link stiffness.
  float         S                 = 0.0f;                                       // Neighbour link strain.
  float         L                 = 0.0f;                                       // Neighbour link length.
  float         dt                = dt_simulation[0];                           // Simulation time step [s].
  int           P0                = particle[0];
  int           P1                = particle[1];
  int           P2                = particle[2];
  int           P3                = particle[3];
  int           P4                = particle[4];
  int           P5                = particle[5];
  int           P6                = particle[6];
  int           P7                = particle[7];
  int           fr_spinor         = (i == P0) || (i == P1) || (i == P2) || (i == P3) || (i == P4) || (i == P5) || (i == P6) || (i == P7);

  // COMPUTING STRIDE MINIMUM INDEX:
  if (i == 0)
  {
    j_min = 0;                                                                  // Setting stride minimum (first stride)...
  }
  else
  {
    j_min = offset[i - 1];                                                      // Setting stride minimum (all others)...
  }

  // COMPUTING ELASTIC FORCE:
  for (j = j_min; j < j_max; j++)
  {
    k = neighbour[j];                                                           // Computing neighbour index...
    nearest = position_int[k];                                                  // Getting neighbour position...
    link = nearest - p_int;                                                     // Getting neighbour link vector...
    R = resting[j];                                                             // Getting neighbour link resting length...
    K = stiffness[j];                                                           // Getting neighbour link stiffness...
    L = length(link);                                                           // Computing neighbour link length...
    S = L - R;                                                                  // Computing neighbour link strain...
    
    if(L > 0.0f)
    {
      D = S*normalize(link);                                                    // Computing neighbour link displacement...
    }
    else
    {
      D = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
    }

    Fe += K*D;                                                                  // Building up elastic force on central node...
  }
  
  Fv = -B*v_int;                                                                // Computing node viscous force...

  // COMPUTING TOTAL FORCE:
  F  = Fe + Fv;                                                                 // Total force applied to the particle [N]...

  // COMPUTING NEW ACCELERATION ESTIMATION:
  a_est  = F/m;                                                                 // Computing acceleration [m/s^2]...

  // COMPUTING NEW VELOCITY ESTIMATION:
  v_est = v + 0.5f*(a + a_est)*dt;                                              // Computing velocity...

  // COMPUTING NEW VISCOUS FORCE ESTIMATION:
  Fv_est = -B*v_est;                                                            // Computing node viscous force...

  // COMPUTING NEW TOTAL FORCE:
  F_new = Fe + Fv_est;                                                          // Computing total node force...

  // COMPUTING NEW ACCELERATION:
  a_new = F_new/m;                                                              // Computing acceleration...
  
  // APPLYING FREEDOM CONSTRAINTS:
  if (fr == 0)
  {
    a_new = (float4)(0.0f, 0.0f, 0.0f, 1.0f);                                   // Constraining acceleration...
  }

  // COMPUTING NEW VELOCITY:
  v_new = v + 0.5f*(a + a_new)*dt;                                              // Computing velocity...

  // APPLYING FREEDOM CONSTRAINTS:
  if (fr == 0)
  {
    v_new = (float4)(0.0f, 0.0f, 0.0f, 1.0f);                                   // Constraining velocity...
  }

  // FIXING PROJECTIVE SPACE:
  p_int.w = 1.0f;                                                               // Adjusting projective space...
  v_new.w = 1.0f;                                                               // Adjusting projective space...
  a_new.w = 1.0f;                                                               // Adjusting projective space...

  // UPDATING KINEMATICS:
  position[n] = p_int;                                                          // Updating position [m]...
  velocity[n] = v_new;                                                          // Updating velocity [m/s]...
  acceleration[n] = a_new;                                                      // Updating acceleration [m/s^2]...
}