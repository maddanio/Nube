// File: OBJLoader.h
// Purpose: Header file for .obj file reader

#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/quaternion.hpp>

#include "Exception.h"


namespace nimbus
{
	// Info of ellipsoid

	/**
	Struct for sphere storage
	*/
	struct SPH
	{
		glm::vec3 pos;  // 3D position
		float siz;		// Size
	};

	/**
	.obj file reader class
	*/
	class OBJLoader
	{
	private:
		std::ifstream file; // File object
	public:
		bool openObjFile(const std::string& filename);
		void readVertices(std::vector<glm::vec3>& vertices);
		void readFaces(std::vector<glm::ivec3>& faces);
		void generateEllipsoids(std::vector<glm::vec3>& vertices, std::vector<glm::ivec3>& faces, std::vector<SPH>& out);
		void generateBarycenter(std::vector<glm::vec3>& vertices, std::vector<glm::ivec3>& faces, std::vector<glm::vec3>& out);
		~OBJLoader();
	};
}

