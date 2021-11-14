// File: FluidCUDA.cpp
// Purpose: Implementation file for fluid dynamics in CUDA
// Based on the article of Jos Stam

#include "FluidCUDA.h"
#include <iostream>

#define MAX(a,b)            (((a) > (b)) ? (a) : (b))
#define LINEARSOLVERTIMES 10
#define SWAP(x0,x) {nRarray *tmp=x0;x0=x;x=tmp;} // Exchange pointers

using namespace nimbus;

/**
Constructor
@param dt Time spacing between the snapshots
@param visc Fluid viscosity
@param diff Diffusion rate
*/
FluidCUDA::FluidCUDA(float dt, float diff, float visc)
{
	this->dt = dt;
	this->diff = diff;
	this->visc = visc;

	allocateData();

	clearData();

	numBlocksX = (int)ceil((float)(M + 2) / (float)THREADS_X);
	numBlocksY = (int)ceil((float)(N + 2) / (float)THREADS_Y);
	numBlocksZ = (int)ceil((float)(O + 2) / (float)THREADS_Z);

	clearUVW();
}


/**
Set wind U component force (X,Y,Z)
@param force Wind force
@param i-index
@param j-index
@param k-index
*/
void FluidCUDA::setUForce(float force, int i, int j, int k)
{
	u[i][j][k] = force;
}

/**
Set wind V component force (X,Y,Z)
@param force Wind force
@param i-index
@param j-index
@param k-index
*/
void FluidCUDA::setVForce(float force, int i, int j, int k)
{
	v[i][j][k] = force;
}

/**
Set wind W component force (X,Y,Z)
@param force Wind force
@param i-index
@param j-index
@param k-index
*/
void FluidCUDA::setWForce(float force, int i, int j, int k)
{
	w[i][j][k] = force;
}

/**
Retrieve wind U component force (X,Y,Z)
@param i-index
@param j-index
@param k-index
@return U force
*/

float FluidCUDA::getUForce(int i, int j, int k)
{
	return u[i][j][k];
}

/**
Retrieve wind V component force (X,Y,Z)
@param i-index
@param j-index
@param k-index
@return V force
*/
float FluidCUDA::getVForce(int i, int j, int k)
{
	return v[i][j][k];
}

/**
Retrieve wind W component force (X,Y,Z)
@param i-index
@param j-index
@param k-index
@return W force
*/

float FluidCUDA::getWForce(int i, int j, int k)
{
	return w[i][j][k];
}

/**
Retrieve grid dimensions
@param 3D grid dimensions
*/
glm::ivec3 FluidCUDA::getDimensions()
{
	return glm::ivec3(M, N, O);
}

/**
Clear FluidCUDA data
*/
void FluidCUDA::clearData()
{

	for (int i = 0; i < M + 2; i++)
		for (int j = 0; j < N + 2; j++)
			for (int k = 0; k < O + 2; k++)
				u[i][j][k] = v[i][j][k] = w[i][j][k] = 0.0;
}

/**
Clear pre-data in device
*/
__global__ void clearUVWDev(float x[][N + 2][O + 2], float y[][N + 2][O + 2], float z[][N + 2][O + 2])
{


	int i = blockIdx.x*blockDim.x + threadIdx.x;
	int j = blockIdx.y*blockDim.y + threadIdx.y;
	int k = blockIdx.z*blockDim.z + threadIdx.z;

	if ((i >= 1) && (j >= 1) && (k >= 1) && (i <= M) && (j <= N) && (k <= O))
		x[i][j][k] = y[i][j][k] = z[i][j][k] = 0.0;
}

/**
Clear pre-data
*/
void FluidCUDA::clearUVW()
{
	dim3 block(numBlocksX, numBlocksY, numBlocksZ);
	dim3 thread(THREADS_X, THREADS_Y, THREADS_Z);
	clearUVWDev << <block, thread >> > (devUPrev, devVPrev, devWPrev);
}

/**
Allocate Fluid CUDA data
*/
int FluidCUDA::allocateData()
{
	size = (M + 2)*(N + 2)*(O + 2);

	u = (nRarray*)malloc(size * sizeof(float));
	v = (nRarray*)malloc(size * sizeof(float));
	w = (nRarray*)malloc(size * sizeof(float));

	if (!u || !v || !w)
		throw nimbus::NimbusException("Can't allocate data in host for CUDA fluids", __FILE__, __FUNCTION__, __LINE__);


	cudaMalloc((void**)&devU, size * sizeof(float));
	cudaMalloc((void**)&devV, size * sizeof(float));
	cudaMalloc((void**)&devW, size * sizeof(float));
	cudaMalloc((void**)&devUPrev, size * sizeof(float));
	cudaMalloc((void**)&devVPrev, size * sizeof(float));
	cudaMalloc((void**)&devWPrev, size * sizeof(float));

	cudaMalloc((void**)&devUy, size * sizeof(float));
	cudaMalloc((void**)&devVy, size * sizeof(float));
	cudaMalloc((void**)&devWy, size * sizeof(float));
	cudaMalloc((void**)&devUPrevy, size * sizeof(float));


	if (!devU || !devV || !devW || !devUPrev || !devVPrev || !devWPrev || !devUy || !devVy || !devWy || !devUPrevy)
		throw nimbus::NimbusException("Can't allocate data in device for CUDA fluids", __FILE__, __FUNCTION__, __LINE__);

	return 1;
}

/**
Free FluidCUDA data
*/
void FluidCUDA::freeData()
{
	if (u) free(u);
	if (v) free(v);
	if (w) free(w);

	if (devU) cudaFree(devU);
	if (devV) cudaFree(devV);
	if (devW) cudaFree(devW);
	if (devUPrev) cudaFree(devUPrev);
	if (devVPrev) cudaFree(devVPrev);
	if (devWPrev) cudaFree(devWPrev);

	if (devUy) cudaFree(devUy);
	if (devVy) cudaFree(devVy);
	if (devWy) cudaFree(devWy);
	if (devUPrevy) cudaFree(devUPrevy);

}

/**
Send data to device
*/
void FluidCUDA::sendData()
{
	cudaMemcpy(devU, u, size * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(devV, v, size * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(devW, w, size * sizeof(float), cudaMemcpyHostToDevice);
}

/**
Received calculated data
*/
void FluidCUDA::receiveData()
{
	cudaMemcpy(u, devU, size * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(v, devV, size * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(w, devW, size * sizeof(float), cudaMemcpyDeviceToHost);
}

/**
Simulate FluidCUDA
*/
void FluidCUDA::sim()
{
	velStep();
}


/**
Add source
*/

__global__ void addSource(float  x[][N + 2][O + 2], float s[][N + 2][O + 2], float dt)
{

	int i = blockIdx.x*blockDim.x + threadIdx.x;
	int j = blockIdx.y*blockDim.y + threadIdx.y;
	int k = blockIdx.z*blockDim.z + threadIdx.z;

	if ((i >= 1) && (j >= 1) && (k >= 1) && (i <= M) && (j <= N) && (k <= O))
		x[i][j][k] += dt * s[i][j][k];
}


/**
Linear solve
*/
__global__ void kernelLinSolve(float  x[][N + 2][O + 2], float x0[][N + 2][O + 2], float y[][N + 2][O + 2], float a, float c)
{


	int i = blockIdx.x*blockDim.x + threadIdx.x;
	int j = blockIdx.y*blockDim.y + threadIdx.y;
	int k = blockIdx.z*blockDim.z + threadIdx.z;


	if ((i >= 1) && (j >= 1) && (k >= 1) && (i <= M) && (j <= N) && (k <= O))
	{
		// update for each cell
		y[i][j][k] = (x0[i][j][k] + a * (x[i - 1][j][k] + x[i + 1][j][k] + x[i][j - 1][k] + x[i][j + 1][k] + x[i][j][k - 1] + x[i][j][k + 1])) / c;
	}

}

/**
Diffuse
*/
void FluidCUDA::diffuse(float x[][N + 2][O + 2], float x0[][N + 2][O + 2], float diff, float dt, float dev_y[][N + 2][O + 2])
{
	int max = MAX(MAX(M, N), MAX(N, O));
	float a = dt * diff*max*max*max;

	dim3 block(numBlocksX, numBlocksY, numBlocksZ);
	dim3 thread(THREADS_X, THREADS_Y, THREADS_Z);

	for (int k = 0; k < LINEARSOLVERTIMES; k++)
	{
		kernelLinSolve << <block, thread >> > (x, x0, dev_y, a, 1.0f + 6.0f * a);
		SWAP(x, dev_y);
	}

}

/**
Advect
*/

__global__ void kernelAdvect(float d[][N + 2][O + 2], float d0[][N + 2][O + 2], float u[][N + 2][O + 2], float v[][N + 2][O + 2], float w[][N + 2][O + 2], float dt, float dtx, float dty, float dtz)
{


	int i = blockIdx.x*blockDim.x + threadIdx.x;
	int j = blockIdx.y*blockDim.y + threadIdx.y;
	int k = blockIdx.z*blockDim.z + threadIdx.z;

	int i0, j0, k0, i1, j1, k1;
	float x, y, z, s0, t0, s1, t1, u1, u0;



	if ((i >= 1) && (j >= 1) && (k >= 1) && (i <= M) && (j <= N) && (k <= O))
	{

		x = i - dtx * u[i][j][k]; y = j - dty * v[i][j][k]; z = k - dtz * w[i][j][k];
		if (x < 0.5) x = 0.5; if (x > M + 0.5) x = M + 0.5; i0 = (int)x; i1 = i0 + 1;
		if (y < 0.5) y = 0.5; if (y > N + 0.5) y = N + 0.5; j0 = (int)y; j1 = j0 + 1;
		if (z < 0.5) z = 0.5; if (z > O + 0.5) z = O + 0.5; k0 = (int)z; k1 = k0 + 1;

		s1 = x - i0; s0 = 1 - s1; t1 = y - j0; t0 = 1 - t1; u1 = z - k0; u0 = 1 - u1;
		d[i][j][k] = s0 * (t0*u0*d0[i0][j0][k0] + t1 * u0*d0[i0][j1][k0] + t0 * u1*d0[i0][j0][k1] + t1 * u1*d0[i0][j1][k1]) +
			s1 * (t0*u0*d0[i1][j0][k0] + t1 * u0*d0[i1][j1][k0] + t0 * u1*d0[i1][j0][k1] + t1 * u1*d0[i1][j1][k1]);
	}



}

/**
Basic idea behind advection step: look for particles which end up exactly at the cell centers by tracking backwards in time from the cell centers (with a linear backtrace)
*/

void FluidCUDA::advect(float d[][N + 2][O + 2], float d0[][N + 2][O + 2], float u[][N + 2][O + 2], float v[][N + 2][O + 2], float w[][N + 2][O + 2], float dt)
{


	dim3 block(numBlocksX, numBlocksY, numBlocksZ);
	dim3 thread(THREADS_X, THREADS_Y, THREADS_Z);

	float dtx, dty, dtz;

	dtx = dty = dtz = dt * MAX(MAX(M, N), MAX(N, O));

	kernelAdvect << < block, thread >> > (d, d0, u, v, w, dt, dtx, dty, dtz);


}

/**
Project part 1
*/

__global__ void kernelProject1(float u[][N + 2][O + 2], float v[][N + 2][O + 2], float w[][N + 2][O + 2], float p[][N + 2][O + 2], float div[][N + 2][O + 2])
{
	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;
	int j = blockIdx.y*blockDim.y + threadIdx.y + 1;
	int k = blockIdx.z*blockDim.z + threadIdx.z + 1;


	if ((i >= 1) && (j >= 1) && (k >= 1) && (i <= M) && (j <= N) && (k <= O))
	{
		div[i][j][k] = -1.0 / 3.0*((u[i + 1][j][k] - u[i - 1][j][k]) / M + (v[i][j + 1][k] - v[i][j - 1][k]) / M + (w[i][j][k + 1] - w[i][j][k - 1]) / M);
		p[i][j][k] = 0.0;
	}
}

/**
Project part 2
*/
__global__ void kernelProject2(float u[][N + 2][O + 2], float v[][N + 2][O + 2], float w[][N + 2][O + 2], float p[][N + 2][O + 2], float div[][N + 2][O + 2])
{
	int i = blockIdx.x*blockDim.x + threadIdx.x + 1;
	int j = blockIdx.y*blockDim.y + threadIdx.y + 1;
	int k = blockIdx.z*blockDim.z + threadIdx.z + 1;



	if ((i >= 1) && (j >= 1) && (k >= 1) && (i <= M) && (j <= N) && (k <= O))
	{

		u[i][j][k] -= 0.5*M*(p[i + 1][j][k] - p[i - 1][j][k]);
		v[i][j][k] -= 0.5*M*(p[i][j + 1][k] - p[i][j - 1][k]);
		w[i][j][k] -= 0.5*M*(p[i][j][k + 1] - p[i][j][k - 1]);
	}

}


/**
Computing the height field involves the solution of some linear system called a Poisson equation - reuse Gauss-Seidel relaxation code from the density diffusion function.
*/
void FluidCUDA::project(float u[][N + 2][O + 2], float v[][N + 2][O + 2], float w[][N + 2][O + 2], float p[][N + 2][O + 2], float div[][N + 2][O + 2], float dev_y[][N + 2][O + 2])
{

	dim3 block(numBlocksX, numBlocksY, numBlocksZ);
	dim3 thread(THREADS_X, THREADS_Y, THREADS_Z);

	kernelProject1 << < block, thread >> > (u, v, w, p, div);

	nRarray* aux = p;

	for (int k = 0; k < LINEARSOLVERTIMES; k++)
	{
		kernelLinSolve << < block, thread >> > (p, div, dev_y, 1, 6); SWAP(p, dev_y);
	}

	p = aux;

	kernelProject2 << < block, thread >> > (u, v, w, p, div);

}

/**
Simulation step
*/
void FluidCUDA::velStep()
{
	dim3 block(numBlocksX, numBlocksY, numBlocksZ);
	dim3 thread(THREADS_X, THREADS_Y, THREADS_Z);

	addSource << <block, thread >> > (devU, devUPrev, dt); addSource << <block, thread >> > (devV, devVPrev, dt); addSource << <block, thread >> > (devW, devWPrev, dt);
	SWAP(devUPrev, devU);
	diffuse(devU, devUPrev, visc, dt, devUy);
	SWAP(devVPrev, devV);
	diffuse(devV, devVPrev, visc, dt, devVy);
	SWAP(devWPrev, devW);
	diffuse(devW, devWPrev, visc, dt, devWy);
	project(devU, devV, devW, devUPrev, devVPrev, devUPrevy);
	SWAP(devUPrev, devU); SWAP(devVPrev, devV); SWAP(devWPrev, devW);
	advect(devU, devUPrev, devUPrev, devVPrev, devWPrev, dt); advect(devV, devVPrev, devUPrev, devVPrev, devWPrev, dt); advect(devW, devWPrev, devUPrev, devVPrev, devWPrev, dt);
	project(devU, devV, devW, devUPrev, devVPrev, devUPrevy);
}

/**
Destructor
*/
FluidCUDA::~FluidCUDA()
{
	freeData();
}

