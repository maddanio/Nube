// File: Cloud.h
// Purpose: Header file for all cloud types base class

#pragma once





#include <iostream>
#include <vector>
#include <time.h>
#include <random>
#include <stdlib.h>   

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/quaternion.hpp>

#include "Shader.h" 
#include "Precompute.h"
#include "Exception.h"

const int MAXTEXSIZE = 512; // Cloud noise size
const int MAXCLOUDS = 10;   // Total clouds in scene
const int MAXSPHERES = 1000; // Total spheres in scene

namespace nimbus
{
	/**
	Base class for all cloud types
	*/
	class Cloud
	{
	private:
		static float* noise; // For texture memory allocation
		static GLuint textureID;
		static int dim; // Dimension of texture
		// Uniforms
		static GLint iNoiseTextCloud;   // Cloud texture
		static GLint iResolutionCloud;  // Resolution of screen
		static GLint iMouseCloud;		// Mouse position
		static GLint iCloudDepthCloud;  // Camera distance to cloud
		static GLint iSkyTurnCloud;		// Scene day time
		static GLint iPosCloud;			// Cloud sphere positions
		static GLint iNumSph;			// Number of spheres
		static GLint iDebug;			// For debug purposes
		static GLint iVoxelText[MAXCLOUDS];	// Precompute light grids
		static GLint iTime;				// Time for rendering
		static GLint iSunDir;			// Sun direction
	protected:
		int id;							// Cloud ID
		static int cloudNum;
		static int cloudOffset;			// Spheres offset for all clouds
		static glm::vec3 vmin[MAXCLOUDS];	// Bounding box 
		static glm::vec3 vmax[MAXCLOUDS];
		static GLuint lowLimits[MAXCLOUDS]; // Downlimit of each cloud spheres in the array
		static GLuint upLimits[MAXCLOUDS];	// Uplimit of each cloud spheres in the array
		static GLuint voxelTextureID[MAXCLOUDS]; // OpenGL texture ID
		static GLuint cloudPosUBO;				 // OpenGL UBO
		// Uniforms
		static GLint iNumClouds;				// Number of clouds
		static GLint iLowLimits;				// Lowlimit of each cloud spheres in the array
		static GLint iUpLimits;					// Uplimit of each cloud spheres in the array
		static GLint iVmin;						// Bounding box
		static GLint iVmax;
		static float timeDir;					// Evolution of time
		static GLint iCloudPosBlock;			// OpenGL-GLSL great buffer
		static int voxelSize;					// Grid size
		static glm::vec3 sunDir;				// Sun direction

	public:
		Cloud();
		static int IX(int i, int j, int k);
		static void createTexture(int dim);
		static void getUniforms(Shader& shader);
		static void renderTexture();
		static void renderFirstTime(GLint scrW, GLint scrH);
		static void render(glm::vec2& mousePos, float time, float cloudDepth, int skyTurn, glm::vec3 cameraPos, bool debug);
		static int getNumClouds();
		static void resetNumClouds();
		static void renderVoxelTexture(int cloudIndex);
		static void freeTexture();
		static int getNumSpheres();
		virtual ~Cloud();
	};
}