// File: Precompute.h
// Purpose: Header file for precomputed light interface

#pragma once
//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"

#define GLM_FORCE_INLINE 

#include <glm/glm.hpp>

#  include <OpenGL/gl3.h>
#  include <OpenGL/gl3ext.h>



#include "Defines.h"
#include "Exception.h"

#define NV 40 // Size of voxel grid

typedef float precArray[NV][NV];

namespace nimbus
{

	/**
	Precompute light interface
	*/
	class PrecomputeLight
	{
	public:
		virtual void setTotalSPH(int totalSPH) = 0;
		virtual int getVoxelsSize() = 0;
		virtual void precomputeCloud(glm::vec4* pos, int numSph, int numClouds, glm::vec3 sundir, float sunDistance, float darkLevel, glm::vec3* vmin, glm::vec3* vmax, GLuint* voxelTextureID) = 0;
		virtual void precomputeModel(glm::mat4* pos, int numSph, int numClouds, glm::vec3 sundir, float sunDistance, float darkLevel, glm::vec3* vmin, glm::vec3* vmax, GLuint* voxelTextureID) = 0;
	};

}