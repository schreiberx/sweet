#version 150

in vec2 rast_texture_coord;

out vec3 frag_data;

uniform sampler2D sampler;

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
//	vec2 texd = 1.0/textureSize(sampler, 0);

	float red_value = texture(sampler, rast_texture_coord).r;

	float limited_height = max(0.0, min(1.0, red_value));

	frag_data.xyz = blend2ndOrder(
					vec3(4.0)*(limited_height-vec3(0.75, 0.5, 0.25))
				);

	return;
}
