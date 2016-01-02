// use opengl 3.2 core version!
#version 150

/**
 * vertex_attrib.xy : left top origin of character to draw
 * vertex_attrib.z: ascii character number
 * vertex_attrib.w: vertex corner (0,1,2 or 4)
 */
in ivec4 vertex_attrib;
out vec2 texel_pos;

uniform ivec2 glyph_size;
uniform vec2 inv_viewport_size;

uniform sampler2D sampler;

void main(void)
{
	// compute vertex position
	ivec2 vp = vertex_attrib.xy;
	vp.x += (vertex_attrib.a & 1) * glyph_size.x;
	vp.y += (vertex_attrib.a >> 1) * glyph_size.y;

	gl_Position = vec4(vec2(vp)*inv_viewport_size*2.0-vec2(1,1), 0, 1);

	// compute texture position
	texel_pos.x =	((vertex_attrib.z & 0xf)) * glyph_size.x;
	texel_pos.y =	((vertex_attrib.z >> 4)) * glyph_size.y;

	texel_pos.x += (vertex_attrib.w & 1) * glyph_size.x;
	texel_pos.y += (vertex_attrib.w >> 1) * glyph_size.y;
}
