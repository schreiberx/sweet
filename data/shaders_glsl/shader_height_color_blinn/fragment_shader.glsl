#version 150

// TEXTURE0
uniform bool texture0_enabled = false;
uniform sampler2D texture0;

in vec2 frag_texture_coord2;

out vec4 frag_data;

in float height;

vec4 blinnShading(vec4 material_color);

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
	vec4 frag_color;
	frag_color = vec4(1, 0, 0, 0);
	
	float limited_height = max(0.05, min(0.8, height));

	frag_color.xyz = blend2ndOrder(
					vec3(4.0)*(limited_height-vec3(0.75, 0.5, 0.25))
				);

	frag_data = blinnShading(frag_color);
}
