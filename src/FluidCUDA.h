// File: FluidCUDA.h
// Purpose: Header file for fluid dynamics in CUDA
// Based on the article of Jos Stam

#pragma once

// CUDA programming required headers

#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define GLM_FORCE_RADIANS
#define GLM_FORCE_INLINE 
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/quaternion.hpp>

#include "Fluid.h"

// Block size

#define THREADS_X 16
#define THREADS_Y 8
#define THREADS_Z 8

// Define array 3D

typedef float nRarray[N + 2][O + 2];

namespace nimbus
{

	/**
	Class for GPGPU-CUDA fluid simulation
	*/
	class FluidCUDA : public Fluid
	{
	private:
		//fluid field information
		float dt; // time delta
		float diff; // diffuse
		float visc; // viscosity
		int size;
		nRarray *u, *v, *w;
		nRarray *devU, *devV, *devW, *devUPrev, *devVPrev, *devWPrev;  // storage for result computed on device
		nRarray *devUy, *devVy, *devWy, *devUPrevy;

		int numBlocksX;
		int numBlocksY;
		int numBlocksZ;

	private:
		inline int  IX(int i, int j, int k);
		inline void diffuse(float x[][N + 2][O + 2], float x0[][N + 2][O + 2], float diff, float dt, float dev_y[][N + 2][O + 2]);
		inline void advect(float d[][N + 2][O + 2], float d0[][N + 2][O + 2], float u[][N + 2][O + 2], float v[][N + 2][O + 2], float w[][N + 2][O + 2], float dt);
		inline void project(float u[][N + 2][O + 2], float v[][N + 2][O + 2], float w[][N + 2][O + 2], float p[][N + 2][O + 2], float div[][N + 2][O + 2], float dev_y[][N + 2][O + 2]);
		inline void velStep();
	public:
		FluidCUDA(float dt, float diff, float visc);
		void  setUForce(float force, int i, int j, int k);
		void  setVForce(float force, int i, int j, int k);
		void  setWForce(float force, int i, int j, int k);
		float getUForce(int i, int j, int k);
		float getVForce(int i, int j, int k);
		float getWForce(int i, int j, int k);
		glm::ivec3 getDimensions();
		void clearUVW();
		void freeData();
		void clearData();
		int allocateData();
		void sendData();
		void receiveData();
		void sim();
		~FluidCUDA();
	};

}