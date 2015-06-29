#version 150

// TEXTURE0
uniform bool texture0_enabled = false;
uniform sampler2D texture0;
uniform sampler2DShadow texture_shadow_map;
uniform sampler2D texture_caustic_map;
uniform sampler2DShadow texture_caustic_depth_map;
uniform mat4 shadow_map_matrix;

in vec2 frag_texture_coord2;

out vec4 frag_data;

//////////// BLINN SHADER SKELETON //////////////
/*
 * input from vertex shader
 */
in vec3 frag_normal3;
in vec3 frag_position3;

// LIGHT0
uniform bool light0_enabled = true;
uniform vec3 light0_view_pos3 = vec3(10,10,10);	// position in view space

uniform vec3 light0_ambient_color3 = vec3(1.f,1.f,1.f);
uniform vec3 light0_diffuse_color3 = vec3(1.f,1.f,1.f);
uniform vec3 light0_specular_color3 = vec3(1.f,1.f,1.f);

uniform vec3 material_ambient_color3 = vec3(1.f,1.f,1.f);
uniform vec3 material_diffuse_color3 = vec3(1.f,1.f,1.f);
uniform vec3 material_specular_color3 = vec3(1.f,1.f,1.f);
uniform float material_specular_exponent = 20.0f;

vec4 blinnShading(	vec4 material_color,
					vec3 caustic_diffuse_color3,
					float diffuse_specular_factor
				)
{
	if (!light0_enabled)
		return vec4(0,0,0,0);

	vec3 frag_normal3_normalized = normalize(frag_normal3);

	// compute vector from surface point to light
	vec3 light_vec = normalize(light0_view_pos3 - frag_position3);

	vec3 viewer_vec = -normalize(frag_position3);
	vec3 halfway_vec = normalize(light_vec+viewer_vec);

	return 		material_color*vec4(
					(
					// ambient
					light0_ambient_color3*material_ambient_color3	+
					// direct diffuse
					max(0.0f, dot(light_vec, frag_normal3_normalized))
					*light0_diffuse_color3*material_diffuse_color3*diffuse_specular_factor +
					// caustic diffuse
					caustic_diffuse_color3*light0_diffuse_color3*material_diffuse_color3
					), 1.0f)
				+
				vec4(	(pow(max(0.0f, dot(halfway_vec, frag_normal3_normalized)), material_specular_exponent)
						*light0_specular_color3*material_specular_color3*diffuse_specular_factor
						), 1.0f);
}

////////////////// SKELETON END ////////////////////

void main(void)
{
	vec4 p = shadow_map_matrix*vec4(frag_position3, 1);
	p.xyz /= p.w;

	float shadow_depth = texture(texture_shadow_map, p.xyz);

	// set front_diffuse_surface to 1.0, if the current fragment is behind a diffuse surface
	float behind_front_diffuse_surface = float(texture(texture_caustic_depth_map, p.xyz) > p.z-0.001);

	// allow specular and diffuse lighting if the voxel is out of shadow
	float diffuse_specular_factor = (float(shadow_depth < p.z)*0.5+0.5)*behind_front_diffuse_surface;

	vec4 frag_color;
	if (texture0_enabled)	frag_color = texture(texture0, frag_texture_coord2);
	else					frag_color = vec4(1);

	// load caustic light for current voxel
	// set to 0, if the voxel is behind the 'diffuse front face' texture (determined by texture_caustic_depth_map)
	vec3 caustic_diffuse_color3 = texture(texture_caustic_map, p.xy).rgb*behind_front_diffuse_surface;

	frag_data = blinnShading(frag_color, caustic_diffuse_color3, diffuse_specular_factor);
}
