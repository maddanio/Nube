// File: Canvas.h
// Purpose: Header file for framebuffer visualization

#pragma once





#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/quaternion.hpp>

#include "Shader.h"
#include "Camera.h"

namespace nimbus
{
	/**
	Frame buffer visualization class
	*/
	class Canvas
	{
	private:
		GLint iMVPFrameCloud; // Uniforms
		GLint iMVPCloud;
		GLuint vertexBufferShaderFrame; // Vertex buffer for canvas 
		GLuint vertexArrayShaderFrame; // Vertex array for canvas
	public:
		Canvas();
		void create(GLint scrW, GLint scrH);
		void getUniforms(Shader& shader);
		void render(Camera& cameraFrame, Camera& cameraSky);
		~Canvas();
	};

}