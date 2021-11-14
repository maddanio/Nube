// File: Cumulus.h
// Purpose: Header file for cumulus rendering

#pragma once

#include "Cloud.h"
#include "Fluid.h"

namespace nimbus
{
	enum Winds { NORTH, SOUTH, EAST, WEST }; // Wind direction

	/**
	Class for cumulus rendering
	*/
	class Cumulus : public Cloud
	{
	private:
		static int totalSph;
	protected:
		static float maxRad; // Maximum radius for all spheres
		int sphPivot;		 // Cloud guide point	
		int numSph;			 // Number of spheres
		float sphRads[MAXSPHERES]; // Final radius of spheres for animation
		GLuint lowLimit; // Spheres lowlimit in the array
		GLuint upLimit;	 // Spheres uplimit in the array
		static GLint iWindDirection; // Windirection uniform
		static GLint iFlat;			 // Castellanus flat cumulus flag
		static int dimM;			 // Dimension of fluid simulator
		static int dimN;
		static int dimO;
		static Winds windDirection;	// Wind direction object
		bool extinctedCloud;		// Checks if cloud is extincted
		bool createdCloud;			// Checks if cloud is created
		float radIndex;				// Animation radius
		static glm::vec3 extinctionLimits; // Scene cloud extinction limits
		static glm::vec4 sphPos[MAXSPHERES];	// Position of all clouds spheres
		static GLint flat[MAXCLOUDS];			// Castellanus cumulus activation
	public:
		Cumulus();
		void create(int spheres, GLfloat siz, glm::vec3 center, GLfloat nuX, GLfloat sigX, GLfloat nuY, GLfloat sigY, GLfloat nuZ, GLfloat sigZ, bool isFlat, bool optimize);
		static void getUniforms(Shader& shader);
		void setGuidePoint(Winds windDirection);
		int getGuidePoint();
		static void setLimits(glm::vec3 limits);
		bool checkLimits();
		void applyWind(float force, Fluid* windGrid);
		void bornCloud(float step);
		void extinctCloud(float step);
		void computeWind(Fluid* windGrid);
		static void render(Shader& shader);
		static void setWind(Winds windDirection);
		static void precomputeLight(PrecomputeLight& precomp, glm::vec3 sunDir, float sunDistance, float darkLevel);
		bool getExtincted();
		bool getCreated();
		void setCreated(bool created);
		static float getMaxRadius();
		void incRadIndex(float step);
		float getRadIndex();
		static void setTimeDir(float timeDir);
		static float getTimeDir();
		virtual ~Cumulus();
	};
}
