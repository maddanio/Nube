// File: Mountain.h
// Purpose: Header file for moutains render

#pragma once

#include <iostream>
#include <vector>
#include <time.h>
#include <random>



#  include <OpenGL/gl3.h>
#  include <OpenGL/gl3ext.h>

#include "Shader.h"

namespace nimbus
{

	/**
	Class for mountainous landscape rendering
	*/
	class Mountain
	{
	private:
		float noise[256][256]; // Base noise
		GLuint textureID;	// OpenGL texture ID
		GLint iSnow;		// If snow uniforms
		GLint iNoise;		// Noise uniforms
		bool snow;			// If snow
	public:
		Mountain();
		void create(const float height, const bool snow);
		void getUniforms(Shader& shader);
		void  render();
		~Mountain();
	};
}

