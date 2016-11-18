#version 150

in vec3 frag_vertex_normal3;
in vec2 frag_texture_coord2;

out vec3 frag_data;

uniform sampler2D dataTexture;

/**
 * computation of rainbow color map from
 * http://http.developer.nvidia.com/GPUGems/gpugems_ch08.html
 */
vec3 blend2ndOrder(vec3 f)
{
	vec3 y = 1.0-f*f;
	return max(y, vec3(0));
}


void main(void)
{
//	frag_data.xyz = vec3(0);
//	frag_data.x = frag_texture_coord2.y;
//	return;

	float red_value = texture(dataTexture, frag_texture_coord2).r;
//	red_value = frag_texture_coord2.x;

	float limited_height = max(0.0, min(1.0, red_value));
#if 1
	frag_data.xyz = blend2ndOrder(
					vec3(3.3)*(limited_height-vec3(0.75, 0.5, 0.25))	// this generates a B G R order of colours
				);
#else
	frag_data.xyz = blend2ndOrder(
					vec3(3.0)*(limited_height-vec3(1.0, 0.5, 0.0))
				);
#endif
//	frag_data.z = 0;
//	frag_data.xy = frag_texture_coord2;
	return;
}
