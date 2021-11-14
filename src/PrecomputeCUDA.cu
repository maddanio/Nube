// File: PrecomputeCUDA.cpp
// Purpose: Implementation file for precomputed light for CUDA class

#include <iostream>

#include "precomputeCUDA.h"

using namespace nimbus;

/**
Order the spheres from front to back for cumulus
@param org Grid 3D position
@param cen Array of sphere positions
@param numSph Number of spheres
@param dir Sun direction
@param candidates Array of selected candidates
@return Number of candidates
*/

__device__  int orderCumulusCUDA(glm::vec3& org, glm::vec4* pos, int numSph, glm::vec3& dir, glm::vec3* candidates)
{
	int n = 0;
	glm::vec3 aux;

	for (int j = 0; j < numSph; j++)
	{

		glm::vec3 cloudPos = glm::vec3(pos[j]);
		float radius = pos[j].w;

		glm::vec3 temp = org - cloudPos;
		float b = 2.0*glm::dot(dir, temp);
		float c = glm::dot(temp, temp) - radius * radius;

		float disc = b * b - 4.0*c;

		if (disc > 0.0)
		{
			disc = glm::sqrt(disc);
			float t1 = ((-b - disc) / 2.0);
			float t2 = ((-b + disc) / 2.0);

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

__device__ int orderModelCUDA(glm::vec3& org, glm::mat4* pos, int src, int dst, glm::vec3& dir, glm::vec3* candidates)
{
	int n = 0;
	glm::vec3 aux;

	for (int j = src; j < dst; j++)
	{

		glm::vec3 cloudPos = glm::vec3(pos[j][0]);
		float radius = glm::max(pos[j][1][0], glm::max(pos[j][1][1], pos[j][1][2]));

		glm::vec3 temp = org - cloudPos;
		float b = 2.0*glm::dot(dir, temp);
		float c = glm::dot(temp, temp) - radius * radius;

		float disc = b * b - 4.0*c;

		if (disc > 0.0)
		{
			disc = glm::sqrt(disc);
			float t1 = ((-b - disc) / 2.0);
			float t2 = ((-b + disc) / 2.0);

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

__device__ float phaseCUDA(float g)
{

	return 0.0795 * ((1.0 - g * g) / pow(1.0 + g * g - 2.0 *g, 1.5));

}

/**
Trace ray for each cumulus sphere to other spheres in the same cloud and others towards the sun direction
@param org Grid 3D position
@param cen Array of sphere positions
@param numSph Number of spheres
@param dir Sun direction
@return Shadow factor
*/

__device__  float traceRayCumulusCUDA(glm::vec3& org, glm::vec4* cen, int numSph, glm::vec3& dir, float darkLevel)
{
	glm::vec3 candidates[100]; // Allocate space for candidates

	int n = orderCumulusCUDA(org, cen, numSph, dir, candidates);

	if (n == 0) return 1.0f; // No candidates

	float  scatterLight = 0.0;

	float T = 1.0f;

	float tOut;

	float t = candidates[0].x;

	float ph = phaseCUDA(0.3f);

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
			float deltaT = exp(-0.01*den);
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

__device__ float traceRayModelCUDA(glm::vec3& org, glm::mat4* cen, int src, int dst, glm::vec3& dir, float darkLevel)
{
	glm::vec3 candidates[100]; // Allocate space for candidates

	int n = orderModelCUDA(org, cen, src, dst, dir, candidates);

	if (n == 0) return 1.0f; // No candidates

	float  T = 1.0f, scatterLight = 0.0, totalLight;

	float tIn, tOut;

	tIn = candidates[0].x;
	tOut = candidates[n - 1].y;

	float t = tIn;

	float ph = phaseCUDA(0.9);

	while (t <= tOut) // Iterate thorugh the external ellipsoids
	{
		glm::vec3  pos = org + t * dir;
		float den = 1 - (t - tIn) / (tOut - tIn);
		float deltaT = exp(-0.03*den);

		// Scattering                                    
		scatterLight += ph * T*0.0001;
		// Absroted light
		float absortLight = T;
		totalLight = absortLight + scatterLight;
		T *= deltaT;
		if (T < darkLevel)
			return totalLight;
		t += 0.1;

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

__device__ glm::vec3 indexToCoordCUDA(const glm::vec3& index, const glm::vec3& vmin, const glm::vec3& cellSiz)
{
	return (index + glm::vec3(0.5))*cellSiz + vmin;
}

/**
Iterate all the grid to perform light precomputation in cumulus
@param dev_precomp CUDA 3D grid
@param pos Array of sphere positions
@param numSph Number of ellipsoids
@param sunpos Sun 3D position
@param cellSizX X voxel size
@param cellSizY Y voxel size
@param cellSizZ Z voxel size
@param min Bounding box min coordinates
*/

__global__ void precomputeCumulus(float dev_precomp[][NV][NV], glm::vec4* pos, int numSph, glm::vec3 sunpos, float darkLevel, float cellSizX, float cellSizY, float cellSizZ, glm::vec3 min)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	int j = blockIdx.y*blockDim.y + threadIdx.y;
	int k = blockIdx.z*blockDim.z + threadIdx.z;

	if ((i < NV) && (j < NV) && (k < NV)) // Avoid limits ovwerflow
	{
		glm::vec3 coord = indexToCoordCUDA(glm::vec3(i, j, k), min, glm::vec3(cellSizX, cellSizY, cellSizZ));
		glm::vec3 dir = glm::normalize(sunpos - coord);
		dev_precomp[k][j][i] = traceRayCumulusCUDA(coord, pos, numSph, dir, darkLevel);
	}

}

/**
Iterate all the grid to perform light precomputation in 3D meshe
@param dev_precomp CUDA 3D grid
@param pos Array of ellipsoid positions
@param numSph Number of ellipsoids
@param totalSph Total number of ellipsoids
@param bound Index of model for morphing
@param sunpos Sun 3D position
@param cellSizX X voxel size
@param cellSizY Y voxel size
@param cellSizZ Z voxel size
@param min Bounding box min coordinates
*/

__global__ void precomputeModelCUDA(float dev_precomp[][NV][NV], glm::mat4* pos, int numSph, int totalSph, int bound, glm::vec3 sunpos, float darkLevel, float cellSizX, float cellSizY, float cellSizZ, glm::vec3 min)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	int j = blockIdx.y*blockDim.y + threadIdx.y;
	int k = blockIdx.z*blockDim.z + threadIdx.z;

	if ((i < NV) && (j < NV) && (k < NV)) // Avoid limits ovwerflow
	{
		glm::vec3 coord = indexToCoordCUDA(glm::vec3(i, j, k), min, glm::vec3(cellSizX, cellSizY, cellSizZ));
		glm::vec3 dir = glm::normalize(sunpos - coord);
		if (bound == 0) // Condition for morphing
			dev_precomp[k][j][i] = traceRayModelCUDA(coord, pos, 0, numSph, dir, darkLevel);
		else dev_precomp[k][j][i] = traceRayModelCUDA(coord, pos, numSph, totalSph, dir, darkLevel);
	}

}

/**
Constructor
*/

PrecomputeCUDA::PrecomputeCUDA()
{

	precomp = nullptr;
	devPrecomp = nullptr;
	devPos = nullptr;

	precomp = (precArray*)malloc(NV*NV*NV * sizeof(float));

	if (!precomp)
		throw nimbus::NimbusException("Can't allocate data in host for CUDA light precomputation", __FILE__, __FUNCTION__, __LINE__);

	cudaMalloc((void**)&devPrecomp, NV*NV*NV * sizeof(float));

	if (!devPrecomp)
		throw nimbus::NimbusException("Can't allocate data in device for CUDA light precomputation", __FILE__, __FUNCTION__, __LINE__);


	// Calculate block dimensions in grid

	numBlocksX = (int)ceil((double)(NV / (double)THREADS_X));
	numBlocksY = (int)ceil((double)(NV / (double)THREADS_Y));
	numBlocksZ = (int)ceil((double)(NV / (double)THREADS_Z));

	allocated = false;
}

/**
Set number of spheres/ellipsoids
@param totalSPH Total number of spheres/ellipsoids
*/
void PrecomputeCUDA::setTotalSPH(int totalSPH)
{
	if (allocated)
	{
		cudaFree(devPos);
		allocated = false;
	}
		
	if (!allocated)
	{
		std::cout << " LLAMANDO" << std::endl;
#ifdef CUMULUS
		cudaMalloc((void**)&devPos, totalSPH * sizeof(glm::vec4));
#else
		cudaMalloc((void**)&devPos, totalSPH * sizeof(glm::mat4));
#endif
		this->totalSph = totalSPH;
		allocated = true;
	}
}

/** Retrieve grid size
@return grid size
*/

int PrecomputeCUDA::getVoxelsSize()
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

void PrecomputeCUDA::precomputeCloud(glm::vec4* pos, int numSph, int numClouds, glm::vec3& sundir, float sunDistance, float darkLevel, glm::vec3* vmin, glm::vec3* vmax, GLuint* voxelTextureID)
{

	dim3 block(numBlocksX, numBlocksY, numBlocksZ);
	dim3 thread(THREADS_X, THREADS_Y, THREADS_Z);

	glm::vec3 sunpos = sundir * glm::vec3(-sunDistance);

	cudaMemcpy(devPos, pos, numSph * sizeof(glm::vec4), cudaMemcpyHostToDevice);

	for (int bound = 0; bound < numClouds; bound++)
	{

		glm::vec3  min = vmin[bound]; // Bounding box limits
		glm::vec3  max = vmax[bound];

		float cellSizX = ceil((max.x - min.x)) / NV; // Cell size calculation
		float cellSizY = ceil((max.y - min.y)) / NV;
		float cellSizZ = ceil((max.z - min.z)) / NV;

#ifdef CUMULUS

		precomputeCumulus << <block, thread >> > (devPrecomp, devPos, numSph, sunpos, darkLevel, cellSizX, cellSizY, cellSizZ, min);
#endif
		cudaMemcpy(precomp, devPrecomp, NV*NV*NV * sizeof(float), cudaMemcpyDeviceToHost);

		// Create OpenGL texture


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
@param vmin 3D vector of min bounding box coordinates
@param vmax 3D vector of max bounding box coordinates
@param voxelTextureID OpenGL texture ID
*/
void PrecomputeCUDA::precomputeModel(glm::mat4* pos, int numSph, int numClouds, glm::vec3& sundir, float sunDistance, float darkLevel, glm::vec3* vmin, glm::vec3* vmax, GLuint* voxelTextureID)
{


	dim3 block(numBlocksX, numBlocksY, numBlocksZ);
	dim3 thread(THREADS_X, THREADS_Y, THREADS_Z);

	glm::vec3 sunpos = sundir * glm::vec3(-sunDistance);

	cudaMemcpy(devPos, pos, totalSph * sizeof(glm::mat4), cudaMemcpyHostToDevice);

	for (int bound = 0; bound < numClouds; bound++)
	{

		glm::vec3  min = vmin[bound]; // Bounding box limits
		glm::vec3  max = vmax[bound];

		float cellSizX = ceil((max.x - min.x)) / NV; // Cell size calculation
		float cellSizY = ceil((max.y - min.y)) / NV;
		float cellSizZ = ceil((max.z - min.z)) / NV;

#ifdef MODEL
		precomputeModelCUDA << <block, thread >> > (devPrecomp, devPos, numSph, totalSph, bound, sunpos, darkLevel, cellSizX, cellSizY, cellSizZ, min);
#endif
		cudaMemcpy(precomp, devPrecomp, NV*NV*NV * sizeof(float), cudaMemcpyDeviceToHost);

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
PrecomputeCUDA::~PrecomputeCUDA()
{
	if (allocated)
		cudaFree(devPos);
	if (precomp)
		free(precomp);
	if (devPrecomp)
		cudaFree(devPrecomp);

}
