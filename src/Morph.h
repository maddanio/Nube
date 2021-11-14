#pragma once
// File: Morph.h
// Purpose: Header file for morphing

#include "Model.h"

namespace nimbus
{
	/**
	Class for morphing performance
	*/
	class Morph
	{
	private:
		Model* modelA; // Pointer to initial 3D mesh
		Model* modelB;	// Pointer to final 3D mesh
		glm::mat4 cloudPosSrc[MAXSPHERES / 2]; // Source ellipsoids array
		glm::mat4 cloudPosDst[MAXSPHERES / 2]; // Final ellipsois array
		glm::mat4 cloudPosRDst[MAXSPHERES / 2]; // Final rotation matrix array
		glm::vec3 dirTriaSrc[MAXSPHERES / 2]; // Source triangle direction array
		glm::vec3 dirTriaDst[MAXSPHERES / 2]; // Final triangle direction array
		float thetaSrc[MAXSPHERES / 2];		// Source rotation angle
		float thetaDst[MAXSPHERES / 2];		// Final rotation angle
		int offsetSrc;						// Source cloud ellipsoid offset
		int offsetDst;						// Final cloud ellipsoid offset
	public:
		Morph();
		void setModels(Model* modelA, Model* modelB, bool evolute);
		void prepareMorphEvolute();
		void prepareMorphInvolute();
		void morphEvolute(float a);
		void morphInvolute(float a);
		void morph(float step);
		glm::mat4* getCloudPosRDst();
		glm::mat4* getCloudPosDst();
		~Morph();
	};
}