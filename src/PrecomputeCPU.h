#pragma once

#include "Precompute.h"
namespace nimbus
{
	/**
	CPU implementation of precomputed light class
	*/
	class PrecomputeCPU : public PrecomputeLight
	{
	private:
		precArray* precomp;
		int totalSph;
		float darkLevel;
	private:
		inline int orderCumulusCPU(glm::vec3 org, glm::vec4* pos, int numSph, glm::vec3 dir, glm::vec3* candidates);
		inline int orderModelCPU(glm::vec3 org, glm::mat4* pos, int src, int dst, glm::vec3 dir, glm::vec3* candidates);
		inline float traceRayCumulusCPU(glm::vec3 org, glm::vec4* cen, int numSph, glm::vec3 dir);
		inline float traceRayModelCPU(glm::vec3 org, glm::mat4* cen, int src, int dst, glm::vec3 dir);
		inline glm::vec3 indexToCoordCPU(glm::vec3 index, glm::vec3 vmin, glm::vec3 cellSiz);
		void precomputeCumulusCPU(glm::vec4* pos, int numSph, glm::vec3 sunpos, float cellSizX, float cellSizY, float cellSizZ, glm::vec3 min);
		void precomputeModelCPU(glm::mat4* pos, int numSph, int totalSph, int bound, glm::vec3 sunpos, float cellSizX, float cellSizY, float cellSizZ, glm::vec3 min);
	public:
		PrecomputeCPU();
		void setTotalSPH(int totalSPH);
		int getVoxelsSize();
		void precomputeCloud(glm::vec4* pos, int numSph, int numClouds, glm::vec3 sundir, float sunDistance, float darkLevel, glm::vec3* vmin, glm::vec3* vmax, GLuint* voxelTextureID);
		void precomputeModel(glm::mat4* pos, int numSph, int numClouds, glm::vec3 sundir, float sunDistance, float darkLevel, glm::vec3* vmin, glm::vec3* vmax, GLuint* voxelTextureID);
		~PrecomputeCPU();
	};
}