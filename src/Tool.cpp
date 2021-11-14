// File: Tool.cpp
// Purpose: Implementation file for math assistant class

#include "Tool.h"

using namespace nimbus;

/**
Constructor
*/

Tool::Tool()
{}

/** 
Olinde Rodrigues' equations custom
@param a Source direction vector
@param b Destination direction vector
@param angle Rotation angle
@return The Rodrigues' rotation matrix
*/
glm::mat3 Tool::rotAtoB(glm::vec3& a, glm::vec3& b, float& angle)
{
	glm::vec3 x = glm::normalize(glm::cross(a, b));
	float theta = acos(glm::dot(a, b) / (glm::length(a)*length(b)));
	angle = theta;
	glm::mat3 A = glm::mat3(glm::vec3(0, -x.z, x.y), glm::vec3(x.z, 0, -x.x), glm::vec3(-x.y, x.x, 0));
	glm::mat3 R = glm::mat3(1) + sin(theta)*A + (1 - cos(theta))*(A*A);
	return R;
}

/**
Olinde Rodrigues' equations for morphing
@param a Source direction vector
@param b Destination direction vector
@param angle Rotation angle
@return The Rodrigues' rotation matrix
*/
glm::mat3 Tool::rotAtoBMorph(glm::vec3& a, glm::vec3& b, float theta)
{
	glm::vec3 x = glm::normalize(glm::cross(a, b));
	glm::mat3 A = glm::mat3(glm::vec3(0, -x.z, x.y), glm::vec3(x.z, 0, -x.x), glm::vec3(-x.y, x.x, 0));
	glm::mat3 R = glm::mat3(1) + sin(theta)*A + (1 - cos(theta))*(A*A);
	return R;
}

/**
Destructor
*/

Tool::~Tool()
{
}