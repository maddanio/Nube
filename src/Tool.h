// File: Tool.h
// Purpose: Header file for math assistance class

#pragma once

#define GLM_FORCE_RADIANS

#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/quaternion.hpp>

namespace nimbus
{
	/**
	Utility class for math assistance
	*/
	class Tool
	{
	private:
		int w;
		int h;
	public:
		Tool();
		static glm::mat3 rotAtoB(glm::vec3& a, glm::vec3& b, float& angle);
		static glm::mat3 rotAtoBMorph(glm::vec3& a, glm::vec3& b, float theta);
		~Tool();
	};
}