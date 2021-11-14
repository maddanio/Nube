////////////////////////////////////////////////////////////////////
// ONTOGENETIC MODEL FOR REAL-TIME VOLUMETRIC CLOUDS SIMULATION THESIS			
// Software Engineering and Computer Systems Deparment	
// National University for Distance Education (UNED)			    		
// Carlos Jiménez de Parga, PhD student.
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// Last revision 19/04/2019
//////////////////////////////////////////////////////////////////

#version 450 core
#pragma optimize(on)

layout (location = 0) in vec3 position;
layout (location = 1) in vec2 textCoords;

out vec2 fragCoord;

uniform mat4 MVP;

void main()
{
	fragCoord = textCoords;
	
	gl_Position = MVP * vec4(position.x,position.y, position.z, 1);
}
