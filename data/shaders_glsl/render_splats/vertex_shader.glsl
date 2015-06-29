#version 150

/**
 * x, y: 2D position
 * z: surface height
 * w: size
 */
in vec4 splat_data;
uniform mat4 pvm_matrix;
uniform float splat_size_scalar;

out float surface_height;

// minimum surface height
uniform float surface_height_translate = -50.0;

// scalar with which the surface height after subtraction of surface_height_min is multiplied with
uniform float surface_height_scale = 2.0;

void main(void)
{
	gl_Position = vec4(vec2(pvm_matrix*vec4(splat_data.xy, 0, 1)), 0, 1);
	gl_PointSize = splat_data.w*splat_size_scalar;
	surface_height = splat_data.z;

	// subtract start height
	surface_height += surface_height_translate;
	
	// multiply 
	surface_height *= surface_height_scale;
}
