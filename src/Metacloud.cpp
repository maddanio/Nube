// File: Metacloud.cpp
// Purpose: Implementation file for metaballs based clouds class

#include "Metacloud.h"

using namespace nimbus;

GLint Metacloud::iMean;
GLfloat Metacloud::mean[MAXCLOUDS];

/**
Constructor. Assume metecloud as a cumulus cloud.
*/

Metacloud::Metacloud() :Cumulus()
{
}

/**
Retrieve uniforms
@param shader The metacloud shader object
*/

void Metacloud::getUniforms(Shader& shader)
{
	iMean = shader.getUniformLocation("iMean");
}

/**
Create metacloud
@param spheres Number of spheres
@param siz Size of cloud
@param center Cloud 3D position
@param nuX Gaussian X-mean
@param sigX Gaussian X-standard deviation
@param nuY Gaussian Y-mean
@param sigY Gaussian Y-standard deviation
@param nuZ Gaussian Z-mean
@param sigZ Gaussian Z-standard deviation
@param isFlat Cloud with level of condensation
@param optimze Optimze cloud spheres
*/

void Metacloud::create(int spheres, GLfloat siz, glm::vec3& center, GLfloat nuX, GLfloat sigX, GLfloat nuY, GLfloat sigY, GLfloat nuZ, GLfloat sigZ, bool isFlat, bool optimize)
{
	Cumulus::create(spheres, siz, center, nuX, sigX, nuY, sigY, nuZ, sigZ, isFlat, optimize);

	mean[id] = 0.0f;

	for (int i = 0; i < spheres; i++)
		mean[id] += sphRads[i];
	mean[id] = (mean[id] * mean[id]) / spheres; // Update mean for this cloud

	std::cout << "MEAN = " << mean[id] << std::endl;
}

/**
Retrieve sphere mean
@param i Sphere index
@return mean
*/

float Metacloud::getMean(int i)
{
	return mean[i];
}

/**
Render metacloud
*/

void Metacloud::render()
{
	glUniform1fv(iMean, cloudNum, (GLfloat*)mean);
}