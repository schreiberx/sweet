#version 150

layout(lines) in;
//layout(line_strip, max_vertices = 2) out;
layout(points, max_vertices = 8) out;

uniform mat4 pvm_matrix;


// a passthrough geometry shader for color and position
void main()
{
	if (gl_PrimitiveIDIn == 0)
	{
		gl_Position = pvm_matrix*vec4(-1.0, -1.0, 0.0, 1.0);
		EmitVertex();

		gl_Position = pvm_matrix*vec4(1.0, -1.0, 0.0, 1.0);
		EmitVertex();
	}

	if (gl_PrimitiveIDIn == 0)
	{
		gl_Position = pvm_matrix*vec4(-1.0, -1.0, 0.0, 1.0);
		EmitVertex();

		gl_Position = pvm_matrix*vec4(1.0, -1.0, 0.0, 1.0);
		EmitVertex();
	}
	if (gl_PrimitiveIDIn == 1)
	{
		gl_Position = pvm_matrix*vec4(1.0, -1.0, 0.0, 1.0);
		EmitVertex();

		gl_Position = pvm_matrix*vec4(1.0, 1.0, 0.0, 1.0);
		EmitVertex();
	}
	else
	{
		gl_Position = pvm_matrix*vec4(1.0, 1.0, 0.0, 1.0);
		EmitVertex();

		gl_Position = pvm_matrix*vec4(-1.0, 1.0, 0.0, 1.0);
		EmitVertex();
	}

	EndPrimitive();

/*
    for(int i = 0; i < geom_vertex_position.length(); i++)
  {
	gl_Position = geom_vertex_position[i];
    EmitVertex();
  }
*/
}

#else
#version 120
#extension GL_EXT_geometry_shader4 : enable

void main(void)
{

	//increment variable
	int i;

	/////////////////////////////////////////////////////////////
	//This example has two parts
	//	step a) draw the primitive pushed down the pipeline
	//		 there are gl_Vertices # of vertices
	//		 put the vertex value into gl_Position
	//		 use EmitVertex => 'create' a new vertex
	// 		use EndPrimitive to signal that you are done creating a primitive!
	//	step b) create a new piece of geometry (I.E. WHY WE ARE USING A GEOMETRY SHADER!)
	//		I just do the same loop, but swizzle the x and y values
	//	result => the line we want to draw, and the same line, but along the other axis

	//Pass-thru!
	for(i=0; i< gl_VerticesIn; i++){
		gl_Position = gl_PositionIn[i];
		EmitVertex();
	}
	EndPrimitive();
	//New piece of geometry!  We just swizzle the x and y terms
	for(i=0; i< gl_VerticesIn; i++){
		gl_Position = gl_PositionIn[i];
		gl_Position.xy = gl_Position.yx;
		EmitVertex();
	}
	EndPrimitive();
																		      /////////////////////////////////////////////////////////////

}
