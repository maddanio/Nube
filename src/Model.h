// File: Model.h
// Purpose: Header file for 3D meshes based clouds

#pragma once

#include "Cloud.h"
#include "Shader.h"
#include "OBJLoader.h"
#include "Tool.h"
#include "Precompute.h"

namespace nimbus
{
	/**
	Class for 3D mesh model rendering
	*/
	class Model : public Cloud
	{
	private:
		// Uniforms
		static GLint iAlpha; // Linear interpolation factor
		static GLint iEvolute; // If evolute direction
		static GLint iCloudPosBlockR;	// OpenGL buffers for rotation matrix
		static GLuint cloudPosUBOR;
		OBJLoader loader;		// OBJ mesh file loader
		static glm::mat4 sphPos[MAXSPHERES]; // Ellipsoids position
		static glm::mat4 cloudPosR[MAXSPHERES]; // Rodrigues rotational matrix array
		glm::vec3 dirTriaArray[MAXSPHERES / 2]; // Direction of triangles array
		int numSph;
	public:
		Model();
		static void getUniforms(Shader& shader);
		static void precomputeLight(PrecomputeLight& precomp, glm::vec3& sunDir, float sunDistance, float darkLevel, int model1NumSph, int model2NumSph);
		void create(glm::vec3& center, const std::string& file, float scale);
		glm::vec3 getDirTriangle(int i);
		static glm::mat4 getEllipsoid(int i);
		static glm::mat4 getEllipsoidR(int i);
		static glm::mat4* getCloudPosR();
		static void setEllipsoidR(int i, glm::mat4 r);
		static void renderFirstTime(int offset, bool evolute);
		static void render(Shader& shader, glm::mat4* cloudR, glm::mat4* cloudPos, float alpha);
		int getNumEllipsoids();
		~Model();
	};
}