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

  R = length(quad);                                                             // Computing blooming radius.
  k1 = 1.0 - smoothstep(0.0, 0.5, R);                                           // Computing blooming coefficient...
  k2 = 1.0 - smoothstep(0.0, 0.1, R);                                           // Computing smoothing coefficient...
  k3 = 1.0 - smoothstep(0.2, 0.3, R);                                           // Computing smoothing coefficient...

  if ((abs(quad.x) < 0.4f) && (abs(quad.y) < 0.4f))
  {
    discard;                                                                    // Discarding fragment point...
  }

  fragment_color = vec4(1.0f, 0.0f, 0.0f, 1.0f/(depth*depth));        // Setting fragment color...  
}