// File: PrecomputeCUDA.h
// Purpose: Header file for precomputed light CUDA class

#pragma once

#include "precompute.h"

// Number of threads per block

#define THREADS_X 16
#define THREADS_Y 8
#define THREADS_Z 8
namespace nimbus
{
	/**
	GPGPU-CUDA implementation of precomputed light
	*/
	class PrecomputeCUDA : public PrecomputeLight
	{
	private:
		precArray* precomp;
		precArray* devPrecomp;
		int numBlocksX; // Number of blocks
		int numBlocksY;
		int numBlocksZ;
#ifdef CUMULUS
		glm::vec4* devPos;
#else
		glm::mat4* devPos;
#endif
		int totalSph; // Total spheres
		bool allocated; // If shperes/ellipsoids allocated
	public:
		PrecomputeCUDA();
		void setTotalSPH(int totalSPH);
		int getVoxelsSize();
		void precomputeCloud(glm::vec4* pos, int numSph, int numClouds, glm::vec3& sundir, float sunDistance, float darkLevel, glm::vec3* vmin, glm::vec3* vmax, GLuint* voxelTextureID);
		void precomputeModel(glm::mat4* pos, int numSph, int numClouds, glm::vec3& sundir, float sunDistance, float darkLevel, glm::vec3* vmin, glm::vec3* vmax, GLuint* voxelTextureID);
		~PrecomputeCUDA();
	};
}