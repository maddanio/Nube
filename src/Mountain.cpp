// File: Mountain.cpp
// Purpose: Implementation file for mountains render


#include "Mountain.h"

#  include <OpenGL/gl3.h>
#  include <OpenGL/gl3ext.h>


using namespace nimbus;

/**
Constructor
*/
Mountain::Mountain()
{
	iSnow = -1;
	snow = false;
}

/**
Create base 2D texture
@param height Mountain height
@param snow Snow
*/

void Mountain::create(const float height, const bool snow)
{

	std::random_device rd;
	std::default_random_engine generator;
	generator.seed(rd());
	std::uniform_int_distribution<int> distribution(0, 255);

	for (int i = 0; i < 256; i++)
		for (int j = 0; j < 256; j++)
			noise[i][j] = distribution(generator) / height;

	glActiveTexture(GL_TEXTURE0 + 2);

	glGenTextures(1, &textureID);

	glBindTexture(GL_TEXTURE_2D, textureID);

	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_REPEAT);

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, 256, 256, 0, GL_RED, GL_FLOAT, noise);

	this->snow = snow;

}

/**
Retrieve uniforms
@param The cloud shader object
*/

void Mountain::getUniforms(Shader& shader)
{
	iSnow = shader.getUniformLocation("iSnow");
	iNoise = shader.getUniformLocation("iChannel0");
}

/**
Pass to shader and render
*/

void  Mountain::render()
{

	glActiveTexture(GL_TEXTURE0 + 1);
	// Pass 2D texture
	glBindTexture(GL_TEXTURE_2D, textureID);
	glUniform1i(iNoise, 1);
	glUniform1i(iSnow, snow);


}

/**
Destructor
*/
Mountain::~Mountain()
{
}
