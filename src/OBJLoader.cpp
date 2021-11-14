// File: OBJLoader.cpp
// Purpose: Implementation file for .obj file reader


#include "OBJLoader.h"

using namespace nimbus;

/**
Open file
@param filename OBJ path file
@return True if not exception is raised
*/
bool OBJLoader::openObjFile(const std::string& filename)
{
	file.open(filename, std::fstream::in);
	
	if (!file)
		throw nimbus::NimbusException("Error opening OBJ mesh file", __FILE__, __FUNCTION__, __LINE__);

	return true;
}

/**
Read vertices from OBJ file
@param vertices A vector of 3D vertices
*/

void OBJLoader::readVertices(std::vector<glm::vec3>& vertices)
{
	if (file.is_open()) {
		std::string line, line_aux;

		while (line.find("s off") == -1)
		{
			std::getline(file, line);
			int pos;

			if (static_cast<int>(line.find("#")) == -1 && (pos = static_cast<int>(line.find("v"))) == 0)
			{

				line_aux = line.substr(pos + 1);
				float a;
				float b;
				float c;
				std::istringstream iss(line_aux);
				iss >> a >> b >> c;
				vertices.push_back(glm::vec3(a, b, c));

			}
		}
	}

}

/**
Read faces from OBJ file
@param faces A vector of 3D faces
*/

void OBJLoader::readFaces(std::vector<glm::ivec3>& faces)
{
	if (file.is_open()) {
		std::string line, line_aux;

		while (!file.eof())
		{
			std::getline(file, line);
			int pos;

			if (static_cast<int>(line.find("#")) == -1 && (pos = static_cast<int>(line.find("f"))) == 0)
			{

				line_aux = line.substr(pos + 1);
				int a;
				int b;
				int c;
				std::istringstream iss(line_aux);
				iss >> a >> b >> c;
				faces.push_back(glm::ivec3(a, b, c));

			}
		}
	}
}

/**
Generate ellipsoids
@param vertices A vector of 3D vertices
@param faces A vector of 3D faces
@param out A vector of spheres
*/

void OBJLoader::generateEllipsoids(std::vector<glm::vec3>& vertices, std::vector<glm::ivec3>& faces, std::vector<SPH>& out)
{
	std::vector<glm::ivec3>::iterator it;

	for (it = faces.begin(); it != faces.end(); it++)
	{
		glm::vec3 vert1 = vertices.at((*it).x - 1);
		glm::vec3 vert2 = vertices.at((*it).y - 1);
		glm::vec3 vert3 = vertices.at((*it).z - 1);

		glm::vec3 bar = glm::vec3(((vert1.x + vert2.x + vert3.x) / 3.0), ((vert1.y + vert2.y + vert3.y) / 3.0), ((vert1.z + vert2.z + vert3.z) / 3.0));

		float d1 = glm::distance(bar, vert1);
		float d2 = glm::distance(bar, vert2);
		float d3 = glm::distance(bar, vert3);

		float siz = glm::min(d3, glm::min(d1, d2));


		out.push_back(SPH{ bar, siz });
	}
}

/**
Calculate barycenter
@param vertices A vector of 3D vertices
@param faces A vector of 3D faces
@param out A vector of spheres
*/

void OBJLoader::generateBarycenter(std::vector<glm::vec3>& vertices, std::vector<glm::ivec3>& faces, std::vector<glm::vec3>& out)
{
	std::vector<glm::ivec3>::iterator it;

	for (it = faces.begin(); it != faces.end(); it++)
	{
		glm::vec3 vert1 = vertices.at((*it).x - 1);
		glm::vec3 vert2 = vertices.at((*it).y - 1);
		glm::vec3 vert3 = vertices.at((*it).z - 1);

		out.push_back(glm::vec3(((vert1.x + vert2.x + vert3.x) / 3.0), ((vert1.y + vert2.y + vert3.y) / 3.0), ((vert1.z + vert2.z + vert3.z) / 3.0)));
	}
}

/**
Destructor
*/

OBJLoader::~OBJLoader()
{
	file.close();
}
