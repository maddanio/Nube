// File: Canvas.cpp
// Purpose: Implement vertex shader cloud container visualization


#include "Canvas.h"

using namespace nimbus;

/**
Constructor
*/

Canvas::Canvas()
{
}

/**
Create frame buffer 2D canvas 
@param scrW screen width
@param scrH screen height
*/
void Canvas::create(GLint scrW, GLint scrH)
{
	
	GLfloat scrWp = static_cast<GLfloat>(scrW);
	GLfloat scrHp = static_cast<GLfloat>(scrH);
	// Canvas dimensions
	const GLfloat vertexBufferDataShaderFrame[] = {
	0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
	0.0f, scrHp, 0.0f, 0.0f, scrHp,
	scrWp, 0.0f, 0.0f, scrWp, 0.0f,
	scrWp, scrHp, 0.0f, scrWp, scrHp
	};

	glGenVertexArrays(1, &vertexArrayShaderFrame);
	glBindVertexArray(vertexArrayShaderFrame);

	// Generate 1 buffer, put the resulting identifier in vertexbuffer
	glGenBuffers(1, &vertexBufferShaderFrame);

	// Link the buffer
	glBindBuffer(GL_ARRAY_BUFFER, vertexBufferShaderFrame);

	// Give shader frame to buffer
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertexBufferDataShaderFrame), vertexBufferDataShaderFrame, GL_STATIC_DRAW);

}


/**
Retrieve uniforms
@param shader The canvas shader object
*/
void Canvas::getUniforms(Shader& shader)
{
	iMVPFrameCloud = shader.getUniformLocation("MVP");
	iMVPCloud = shader.getUniformLocation("iView");
}

/**
Visualize canvas
@param mvpFrame The Model-View-Projection matrix
@param cameraSky The sky camera
*/
void Canvas::render(Camera& cameraFrame, Camera& cameraSky)
{
	glm::mat4 mvpFrame = cameraFrame.getProjectionMatrix() * cameraFrame.getViewMatrix();
	
	glUniformMatrix4fv(iMVPFrameCloud, 1, GL_FALSE, glm::value_ptr(mvpFrame));
	glUniformMatrix4fv(iMVPCloud, 1, GL_FALSE, glm::value_ptr(cameraSky.getViewMatrix()));

	// Generate main shader frame form virtual texture and real coordinates

	glBindBuffer(GL_ARRAY_BUFFER, vertexBufferShaderFrame);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, (3 + 2) * sizeof(GLfloat), nullptr);

	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, (3 + 2) * sizeof(GLfloat), (void*)(3 * sizeof(float)));

	glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
}

/**
Destructor
*/
Canvas::~Canvas()
{
	glDeleteBuffers(1, &vertexBufferShaderFrame);
}