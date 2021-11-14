// File: PrecomputeCPU.cpp
// Purpose: Implementation file for CPU precompute light

#define NOMINMAX

#include <iostream>
#include "PrecomputeCPU.h"

using namespace nimbus;

/**
Constructor
*/
PrecomputeCPU::PrecomputeCPU()
{

	precomp = nullptr;

	precomp = (precArray*)malloc(NV*NV*NV * sizeof(float));

	if (!precomp)
		throw nimbus::NimbusException("Can't allocate data in host for CPU light precomputation", __FILE__, __FUNCTION__, __LINE__);

}

/**
Order the spheres from front to back for cumulus
@param org Grid 3D position
@param cen Array of sphere positions
@param numSph Number of spheres
@param dir Sun direction
@param candidates Array of selected candidates
@return Number of candidates
*/
int PrecomputeCPU::orderCumulusCPU(glm::vec3 org, glm::vec4* pos, int numSph, glm::vec3 dir, glm::vec3* candidates)
{
	int n = 0;
	glm::vec3 aux;

	for (int j = 0; j < numSph; j++)
	{
		glm::vec3 cloudPos = glm::vec3(pos[j]);
		float radius = pos[j].w;

		glm::vec3 temp = org - cloudPos;
		float b = 2.0f*glm::dot(dir, temp);
		float c = glm::dot(temp, temp) - radius * radius;

		float disc = b * b - 4.0f*c;

		if (disc > 0.0)
		{
			disc = glm::sqrt(disc);
			float t1 = ((-b - disc) / 2.0f);
			float t2 = ((-b + disc) / 2.0f);

			if (t1 > 0.0 && t2 > 0.0)
			{
				candidates[n] = glm::vec3(t1, t2, j);
				n++;
			}
			else if (t1 <= 0.0  && t2 > 0.0)
			{
				candidates[n] = glm::vec3(0, t2, j);
				n++;
			}
		}
	}

	// Insertion-sort algorithm

	int h;
	for (int i = 1; i < n; i++)
	{
		aux = candidates[i];
		h = i - 1;
		while ((h >= 0) && (aux.x < candidates[h].x))
		{
			candidates[h + 1] = candidates[h];
			h--;
		}
		candidates[h + 1] = aux;

	}


	return n;
}

/**
Order ellipsoids from front to back for 3D mesh based clouds
@param org Grid 3D position
@param cen Array of ellipsoid positions
@param src Source model index
@param dst Final model index
@param dir Sun direction
@param candidates Array of selected candidates
@return Number of candidates
*/
int PrecomputeCPU::orderModelCPU(glm::vec3 org, glm::mat4* pos, int src, int dst, glm::vec3 dir, glm::vec3* candidates)
{
	int n = 0;
	glm::vec3 aux;

	for (int j = src; j < dst; j++)
	{


		glm::vec3 cloudPos = glm::vec3(pos[j][0]);
		float radius = glm::max(pos[j][1][0], glm::max(pos[j][1][1], pos[j][1][2]));

		glm::vec3 temp = org - cloudPos;
		float b = 2.0f*glm::dot(dir, temp);
		float c = glm::dot(temp, temp) - radius * radius;

		float disc = b * b - 4.0f*c;

		if (disc > 0.0)
		{
			disc = glm::sqrt(disc);
			float t1 = ((-b - disc) / 2.0f);
			float t2 = ((-b + disc) / 2.0f);

			if (t1 > 0.0 && t2 > 0.0)
			{
				candidates[n] = glm::vec3(t1, t2, j);
				n++;
			}
			else if (t1 <= 0.0  && t2 > 0.0)
			{
				candidates[n] = glm::vec3(0, t2, j);
				n++;
			}
		}
	}

	// Insertion-sort algorithm

	int h;
	for (int i = 1; i < n; i++)
	{
		aux = candidates[i];
		h = i - 1;
		while ((h >= 0) && (aux.x < candidates[h].x))
		{
			candidates[h + 1] = candidates[h];
			h--;
		}
		candidates[h + 1] = aux;

	}

	return n;
}

/**
The Henyey-Greenstein phase function
@param g Asymmetry factor
@return The phase
*/
float phaseCPU(float g)
{

	return 0.0795f * ((1.0f - g * g) / pow(1.0f + g * g - 2.0f *g, 1.5f));

}

/**
Trace ray for each cumulus sphere to other spheres in the same cloud and others towards the sun direction
@param org Grid 3D position
@param cen Array of sphere positions
@param numSph Number of spheres
@param dir Sun direction
@return Shadow factor
*/
float PrecomputeCPU::traceRayCumulusCPU(glm::vec3 org, glm::vec4* cen, int numSph, glm::vec3 dir)
{
	glm::vec3 candidates[100]; // Allocate space for candidates

	int n = orderCumulusCPU(org, cen, numSph, dir, candidates);

	if (n == 0) return 1.0f; // No candidates

	float  scatterLight = 0.0;

	float T = 1.0f;

	float tOut;

	float t = candidates[0].x;

	float ph = phaseCPU(0.3f);

	float totalLight = 0.0f;

	for (int i = 0; i < n; i++)
	{

		if (t > candidates[i].y) // No-duplicate-tracing algorithm
			continue;
		else if (t < candidates[i].x)
			t = candidates[i].x;

		tOut = candidates[i].y;

		glm::vec3 centers = glm::vec3(cen[(int)candidates[i].z]);
		float radius = cen[(int)candidates[i].z].w;

		while (t < tOut) // Iterate spheroid
		{
			glm::vec3  pos = org + t * dir;
			float den = 1.0f - glm::distance(pos, centers) / radius;
			float deltaT = exp(-0.02f * den);
			// Scattering                                    
			scatterLight += ph * T* ((darkLevel == 0.2f) ? 0.001f : 0.0001f);
			// Absorted light
			float absortLight = T;
			totalLight = absortLight + scatterLight;
			T *= deltaT;
			if (T < darkLevel)
				return totalLight;
			t += 0.1f;
		}


	}
	return  totalLight;
}

/**
Trace ray for each mesh ellipsoid to other ellipsoids in the same cloud towards the sun direction
@param org Grid 3D position
@param cen Array of ellipsoid positions
@param src Source model index 
@param dst Final model index
@param dir Sun direction
@return Shadow factor
*/

float PrecomputeCPU::traceRayModelCPU(glm::vec3 org, glm::mat4* cen, int src, int dst, glm::vec3 dir)
{
	glm::vec3 candidates[100]; // Allocate space for candidates

	int n = orderModelCPU(org, cen, src, dst, dir, candidates);

	if (n == 0) return 1.0f; // No candidates found

	float  T = 1.0f, scatterLight = 0.0, totalLight;


	float tIn, tOut;

	tIn = candidates[0].x;
	tOut = candidates[n - 1].y;

	float t = tIn;

	float ph = phaseCPU(0.9f);

	while (t <= tOut) // Iterate thorugh the external ellipsoids
	{
		glm::vec3  pos = org + t * dir;
		float den = 1 - (t - tIn) / (tOut - tIn);
		float deltaT = exp(-0.03f*den);
		// Scattering                                    
		scatterLight += ph * T*0.0001f;
		// Absorted light
		float absortLight = T;
		totalLight = absortLight + scatterLight;
		T *= deltaT;
		if (T < darkLevel)
			return totalLight;
		t += 0.1f;

	}
	return totalLight;
}

/**
Convert grid index to 3D world coordinates
@param index Voxel index
@param vmin Bounding-box min
@param cellSiz Voxel size
@return Grid 3D position
*/

glm::vec3 PrecomputeCPU::indexToCoordCPU(glm::vec3 index, glm::vec3 vmin, glm::vec3 cellSiz)
{
	return (index + glm::vec3(0.5))*cellSiz + vmin;
}

/**
Iterate all the grid to perform light precomputation in cumulus
@param pos Array of sphere positions
@param numSph Number of ellipsoids
@param sunpos Sun 3D position
@param cellSizX X voxel size
@param cellSizY Y voxel size
@param cellSizZ Z voxel size
@param min Bounding box min
*/

void PrecomputeCPU::precomputeCumulusCPU(glm::vec4* pos, int numSph, glm::vec3 sunpos, float cellSizX, float cellSizY, float cellSizZ, glm::vec3 min)
{
	for (int i = 0; i < NV; i++)
		for (int j = 0; j < NV; j++)
			for (int k = 0; k < NV; k++)
			{
				glm::vec3 coord = indexToCoordCPU(glm::vec3(i, j, k), min, glm::vec3(cellSizX, cellSizY, cellSizZ));
				glm::vec3 dir = glm::normalize(sunpos - coord); // Direction to the sun
				precomp[k][j][i] = traceRayCumulusCPU(coord, pos, numSph, dir);
			}
}



/**
Iterate all the grid to perform light precomputation in 3D meshe
@param pos Array of ellipsoid positions
@param numSph Number of ellipsoids
@param totalSph Total number of ellipsoids
@param bound Index of model for morphing
@param sunpos Sun 3D position
@param cellSizX X voxel size
@param cellSizY Y voxel size
@param cellSizZ Z voxel size
@param min Bounding box min
*/
void PrecomputeCPU::precomputeModelCPU(glm::mat4* pos, int numSph, int totalSph, int bound, glm::vec3 sunpos, float cellSizX, float cellSizY, float cellSizZ, glm::vec3 min)
{
	for (int i = 0; i < NV; i++)
		for (int j = 0; j < NV; j++)
			for (int k = 0; k < NV; k++)
			{
				glm::vec3 coord = indexToCoordCPU(glm::vec3(i, j, k), min, glm::vec3(cellSizX, cellSizY, cellSizZ));
				glm::vec3 dir = glm::normalize(sunpos - coord); // Sun direction
				if (bound == 0) // Condition for morphing
					precomp[k][j][i] = traceRayModelCPU(coord, pos, 0, numSph, dir);
				else
					precomp[k][j][i] = traceRayModelCPU(coord, pos, numSph, totalSph, dir);
			}


}

/**
Set number of spheres/ellipsoids
@param totalSPH Total number of spheres/ellipsoids
*/
void PrecomputeCPU::setTotalSPH(int totalSPH)
{
	this->totalSph = totalSPH;
}

/**
Retrieve grid size
@return 3D grid size
*/
int PrecomputeCPU::getVoxelsSize()
{
	return NV;
}

/**
Precompute light for all cumulus clouds
@param pos Array of spheres
@param numSph Number of spheres
@param sundDir Sun direction
@param sunDistance Distance to sun
@param darkLevel Shadow level
@param vmin 3D vector of min bounding box coordinates
@param vmax 3D vector of max bounding box coordinates
@param voxelTextureID OpenGL texture ID
*/

void PrecomputeCPU::precomputeCloud(glm::vec4* pos, int numSph, int numClouds, glm::vec3 sundir, float sunDistance, float darkLevel, glm::vec3* vmin, glm::vec3* vmax, GLuint* voxelTextureID)
{

	glm::vec3 sunpos = sundir * glm::vec3(-sunDistance);

	this->darkLevel = darkLevel;

	for (int bound = 0; bound < numClouds; bound++)
	{

		glm::vec3  min = vmin[bound]; // Bounding box limits
		glm::vec3  max = vmax[bound];

		float cellSizX = ceil((max.x - min.x)) / NV; // Calculate cell size in 3D units
		float cellSizY = ceil((max.y - min.y)) / NV;
		float cellSizZ = ceil((max.z - min.z)) / NV;

		precomputeCumulusCPU(pos, numSph, sunpos, cellSizX, cellSizY, cellSizZ, min);


		glBindTexture(GL_TEXTURE_3D, voxelTextureID[bound]);

		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);

		glTexImage3D(GL_TEXTURE_3D, 0, GL_RED, NV, NV, NV, 0, GL_RED, GL_FLOAT, precomp);
	}
}

/**
Precompute light for all mesh clouds
@param pos Array of ellipsoids
@param numSph Number of ellipsoids
@param numClouds Number of morphing model
@param sundDir Sun direction
@param sunDistance Distance to sun
@param darkLevel Shadow level
@param vmin 3D vector of bounding box coordinates
@param vmax 3D vector of bounding box coordinates
@param voxelTextureID OpenGL texture ID
*/
void PrecomputeCPU::precomputeModel(glm::mat4* pos, int numSph, int numClouds, glm::vec3 sundir, float sunDistance, float darkLevel, glm::vec3* vmin, glm::vec3* vmax, GLuint* voxelTextureID)
{

	glm::vec3 sunpos = sundir * glm::vec3(-sunDistance);

	this->darkLevel = darkLevel;

	for (int bound = 0; bound < numClouds; bound++)
	{

		glm::vec3  min = vmin[bound]; // 3D mesh bounding box
		glm::vec3  max = vmax[bound];

		float cellSizX = ceil((max.x - min.x)) / NV;
		float cellSizY = ceil((max.y - min.y)) / NV;
		float cellSizZ = ceil((max.z - min.z)) / NV;


		precomputeModelCPU(pos, numSph, totalSph, bound, sunpos, cellSizX, cellSizY, cellSizZ, min);

		// Create OpenGL texture

		glGenTextures(1, &voxelTextureID[bound]);

		glBindTexture(GL_TEXTURE_3D, voxelTextureID[bound]);

		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);

		glTexImage3D(GL_TEXTURE_3D, 0, GL_RED, NV, NV, NV, 0, GL_RED, GL_FLOAT, precomp);
	}
}

/**
Destructor
*/
PrecomputeCPU::~PrecomputeCPU()
{
	if (precomp)
		free(precomp);
}

