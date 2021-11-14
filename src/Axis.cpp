// File: Axis.cpp
// Purpose: Implements 3D-Axis rotation


#include "Axis.h"

using namespace nimbus;


/** 
	Constructor
*/

Axis::Axis()
{
}

/**
Create 3D axis
*/
void Axis::create()
{
	glGenBuffers(1, &vertexbufferShaderAxis);
}

/**
Locate uniforms
@param shader The axis shader object
*/
void Axis::getUniforms(Shader& shader)
{
	iUniformMVPAxis = shader.getUniformLocation("MVP");
	iLineColor = shader.getUniformLocation("lineColor");
}

/**
Render 3D axis. Red = X; Green = Y; Blue = Z; White = negative axis
@param mvpAxis Model-View-Projection matrix
@param camera Camera view matrix
@param cameraPos The camera (x,y,z) position
*/
void Axis::render(Camera& cameraAxis)
{
		
	glm::mat4 cameraView = glm::lookAt(cameraAxis.getSrcLookAt(), cameraAxis.getDstLookAt(), glm::vec3(0.0, 1.0, 0.0)); // Orientate axis to user camera

	glm::mat4 mvpAxis = cameraAxis.getProjectionMatrix() * cameraView;
	
	glUniformMatrix4fv(iUniformMVPAxis, 1, GL_FALSE, glm::value_ptr(mvpAxis));

	glm::vec4 vdAxis;

	for (int i = 0; i < 6; i++)
	{

		switch (i)
		{
		case 0: // X - axis
			vdAxis =  glm::vec4(1.0, 0.0, 0.0, 1.0); 
			glUniform4f(iLineColor, 1.0, 0.0, 0.0, 1.0);
			break;
		case 1: // Y - axis
			vdAxis = glm::vec4(0.0, 1.0, 0.0,1.0);
			glUniform4f(iLineColor, 0.0, 1.0, 0.0, 1.0);
			break;
		case 2: // Z - axis
			vdAxis =  glm::vec4(0.0, 0.0, 1.0, 1.0);
			glUniform4f(iLineColor, 0.0, 0.0, 1.0, 1.0);
			break;
		case 3: // X - axis
			vdAxis = glm::vec4(-1.0, 0.0, 0.0, 1.0);
			glUniform4f(iLineColor, 1.0, 1.0, 1.0, 1.0);
			break;
		case 4: // Y - axis
			vdAxis = glm::vec4(0.0, -1.0, 0.0, 1.0);
			glUniform4f(iLineColor, 1.0, 1.0, 1.0, 1.0);
			break;
		case 5: // Z - axis
			vdAxis = glm::vec4(0.0, 0.0, -1.0, 1.0);
			glUniform4f(iLineColor, 1.0, 1.0, 1.0, 1.0);
			break;
		}


		// Axis line end position form center
		vertexBufferDataAxis[0] = cameraAxis.getPosition().x;
		vertexBufferDataAxis[1] = cameraAxis.getPosition().y;
		vertexBufferDataAxis[2] = cameraAxis.getPosition().z;
		vertexBufferDataAxis[3] = vdAxis.x + cameraAxis.getPosition().x;
		vertexBufferDataAxis[4] = vdAxis.y + cameraAxis.getPosition().y;
		vertexBufferDataAxis[5] = vdAxis.z + cameraAxis.getPosition().z;

		
		// Store it to buffer
		glBindBuffer(GL_ARRAY_BUFFER, vertexbufferShaderAxis);

		glBufferData(GL_ARRAY_BUFFER, sizeof(vertexBufferDataAxis), vertexBufferDataAxis, GL_STATIC_DRAW);

		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), nullptr);

		glDrawArrays(GL_LINES, 0, 2);
	}
}

/**
Destructor
*/

Axis::~Axis()
{
	glDeleteBuffers(1, &vertexbufferShaderAxis);
}
