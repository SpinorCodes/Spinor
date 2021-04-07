/// @file

#version 460 core

uniform mat4 V_mat;                                                             // View matrix.
uniform mat4 P_mat;                                                             // Projection matrix.

in  vec4  color;                                                                // Voxel color.
in  vec2  quad;                                                                 // Billboard quad UV coordinates.
in  float depth;                                                                // z-depth.

out vec4 fragment_color;                                                        // Fragment color.

void main(void)
{
  
  float k1;                                                                     // Blooming coefficient.
  float k2;                                                                     // Smoothness coefficient.
  float k3;                                                                     // Smoothness coefficient.
  float R;                                                                      // Blooming radius.
  float z;
  float z_min;
  float z_max;
  float alpha;

  z_min = 3.0f;
  z_max = 0.1f;
  z = (clamp(depth, z_max, z_min) - z_max)/(z_min - z_max);
  R = length(quad);                                                             // Computing blooming radius.
  k1 = 1.0 - smoothstep(0.0, 0.5, R);                                           // Computing blooming coefficient...
  k2 = 1.0 - smoothstep(0.0, 0.1, R);                                           // Computing smoothing coefficient...
  k3 = 1.0 - smoothstep(0.2, 0.3, R);                                           // Computing smoothing coefficient...

  if ((abs(quad.x) < 0.46f) && (abs(quad.y) < 0.46f) ||
      ((-0.1f < quad.x) && (quad.x < 0.1f)) ||
      ((-0.1f < quad.y) && (quad.y < 0.1f)))
  {
    discard;                                                                    // Discarding fragment point...
  }

  alpha = clamp(1.0f - z, 0.6f, 0.9f);

  fragment_color = vec4(alpha, alpha, alpha, alpha);                               // Setting fragment color...  
}