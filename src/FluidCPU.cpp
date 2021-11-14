// File: FluidCPU.cpp
// Purpose: Implementation file for fluid dynamics for CPU class

#include "FluidCPU.h"
#include "time.h"

using namespace nimbus;

#define _USE_MATH_DEFINES

#define MAX(a,b)            (((a) > (b)) ? (a) : (b))
#define SWAP(x0, x) { float * tmp = x0; x0 = x; x = tmp; } // Exchange pointers
#define LINEARSOLVERTIMES 10

/** 
Constructor
@param dt Time spacing between the snapshots
@param visc Fluid viscosity
@param diff Diffusion rate
*/

FluidCPU::FluidCPU(float dt, float diff, float visc)
{
	this->dimM = M;
	this->dimN = N;
	this->dimO = O;
	this->dt = dt;
	this->diff = diff;
	this->visc = visc;

	if (!allocateData())
		throw nimbus::NimbusException("Can't allocate data for CPU fluid", __FILE__, __FUNCTION__, __LINE__);
	clearData();
}

/**
Fast index retrieval
@param i i-index
@param j j-index
@param k k-index
*/
int FluidCPU::IX(int i, int j, int k)
{
	return ((i)+(dimM + 2)*(j)+(dimM + 2)*(dimN + 2)*(k));
}

/**
Set U wind forces (X,Y,Z)
@param force Wind force
@param i i-index
@param j j-index
@param k k-index
*/

void FluidCPU::setUForce(float force, int i, int j, int k)
{
	u[IX(i, j, k)] = force;
}

/** Set V wind forces(X, Y, Z)
@param force Wind force
@param i i-index
@param j j-index
@param k k-index
*/

void FluidCPU::setVForce(float force, int i, int j, int k)
{
	v[IX(i, j, k)] = force;
}

/** Set W wind forces (X, Y, Z)
@param force Wind force
@param i i-index
@param j j-index
@param k k-index
*/

void FluidCPU::setWForce(float force, int i, int j, int k)
{
	w[IX(i, j, k)] = force;
}

/**
Get wind forces (X,Y,Z)
@param i i-index
@param j j-index
@param k k-index
@return The U force
*/

float FluidCPU::getUForce(int i, int j, int k)
{
	return u[IX(i, j, k)];
}

/**
Get wind forces (X,Y,Z)
@param i i-index
@param j j-index
@param k k-index
@return The V force
*/

float FluidCPU::getVForce(int i, int j, int k)
{
	return v[IX(i, j, k)];
}

/**
Get wind forces (X,Y,Z)
@param i i-index
@param j j-index
@param k k-index
@return The W force
*/

float FluidCPU::getWForce(int i, int j, int k)
{
	return w[IX(i, j, k)];
}

/**
Reinitialize grid
*/

void FluidCPU::clearData()
{
	int i;

	for (i = 0; i < size; i++) {
		u[i] = v[i] = w[i] = uPrev[i] = vPrev[i] = wPrev[i] = 0.0f;
	}

}

/**
Reinitialize prevs
*/
void FluidCPU::clearUVW()
{
	for (int i = 0; i < size; i++) {
		uPrev[i] = vPrev[i] = wPrev[i] = 0.0f;
	}
}

/**
Allocate 3D grid
*/

int FluidCPU::allocateData()
{
	size = (dimM + 2)*(dimN + 2)*(dimO + 2);

	u = (float *)malloc(size * sizeof(float));
	v = (float *)malloc(size * sizeof(float));
	w = (float *)malloc(size * sizeof(float));
	uPrev = (float *)malloc(size * sizeof(float));
	vPrev = (float *)malloc(size * sizeof(float));
	wPrev = (float *)malloc(size * sizeof(float));

	if (!u || !v || !w || !uPrev || !vPrev || !wPrev) {
		return 0;
	}

	return 1;
}

/**
Free 3D grid
*/

void FluidCPU::freeData()
{
	if (u) free(u);
	if (v) free(v);
	if (w) free(w);
	if (uPrev) free(uPrev);
	if (vPrev) free(vPrev);
	if (wPrev) free(wPrev);
}

/**
Simulate fluid
*/

void FluidCPU::sim()
{
	velStep();
}

/**
Add source
*/

void FluidCPU::addSource(float* x, float* s, float dt)
{
	int i;
	for (i = 0; i < size; i++) x[i] += dt * s[i];
}

/**
Linear solver
*/

void FluidCPU::linSolve(int b, float* x, float* x0, float a, float c)
{
	int i, j, k, l;

	// iterate the solver
	for (l = 0; l < LINEARSOLVERTIMES; l++) {
		// update for each cell
		for (i = 1; i <= dimM; i++) {
			for (j = 1; j <= dimN; j++) {
				for (k = 1; k <= dimO; k++) {
					x[IX(i, j, k)] = (x0[IX(i, j, k)] + a * (x[IX(i - 1, j, k)] + x[IX(i + 1, j, k)] + x[IX(i, j - 1, k)] + x[IX(i, j + 1, k)] + x[IX(i, j, k - 1)] + x[IX(i, j, k + 1)])) / c;
				}
			}
		}
	}

}

/**
Diffuse
*/
void FluidCPU::diffuse(int b, float* x, float* x0, float diff, float dt)
{
	int max = MAX(MAX(dimM, dimN), MAX(dimN, dimO));
	float a = dt * diff*max*max*max;
	linSolve(b, x, x0, a, 1 + 6 * a);
}



/**
Basic idea behind advection step: look for particles which end up exactly at the cell centers by tracking backwards in time from the cell centers (with a linear backtrace)
*/

void FluidCPU::advect(int b, float* d, float* d0, float* u, float* v, float* w, float dt)
{
	int i, j, k, i0, j0, k0, i1, j1, k1;
	float x, y, z, s0, t0, s1, t1, u1, u0, dtx, dty, dtz;

	dtx = dty = dtz = dt * MAX(MAX(dimM, dimN), MAX(dimN, dimO));

	for (i = 1; i <= dimM; i++) {
		for (j = 1; j <= dimN; j++) {
			for (k = 1; k <= dimO; k++) {
				x = i - dtx * u[IX(i, j, k)]; y = j - dty * v[IX(i, j, k)]; z = k - dtz * w[IX(i, j, k)];
				if (x < 0.5f) x = 0.5f; if (x > dimM + 0.5f) x = dimM + 0.5f; i0 = (int)x; i1 = i0 + 1;
				if (y < 0.5f) y = 0.5f; if (y > dimN + 0.5f) y = dimN + 0.5f; j0 = (int)y; j1 = j0 + 1;
				if (z < 0.5f) z = 0.5f; if (z > dimO + 0.5f) z = dimO + 0.5f; k0 = (int)z; k1 = k0 + 1;

				s1 = x - i0; s0 = 1 - s1; t1 = y - j0; t0 = 1 - t1; u1 = z - k0; u0 = 1 - u1;
				d[IX(i, j, k)] = s0 * (t0*u0*d0[IX(i0, j0, k0)] + t1 * u0*d0[IX(i0, j1, k0)] + t0 * u1*d0[IX(i0, j0, k1)] + t1 * u1*d0[IX(i0, j1, k1)]) +
					s1 * (t0*u0*d0[IX(i1, j0, k0)] + t1 * u0*d0[IX(i1, j1, k0)] + t0 * u1*d0[IX(i1, j0, k1)] + t1 * u1*d0[IX(i1, j1, k1)]);
			}
		}
	}
}

/**
Computing the height field involves the solution of some linear system called a Poisson equation - reuse Gauss-Seidel relaxation code from the density diffusion function.
*/
void FluidCPU::project(float * u, float * v, float * w, float * p, float * div)
{
	int i, j, k;

	for (i = 1; i <= dimM; i++) {
		for (j = 1; j <= dimN; j++) {
			for (k = 1; k <= dimO; k++) {
				div[IX(i, j, k)] = -1.0f / 3.0f*((u[IX(i + 1, j, k)] - u[IX(i - 1, j, k)]) / dimM + (v[IX(i, j + 1, k)] - v[IX(i, j - 1, k)]) / dimM + (w[IX(i, j, k + 1)] - w[IX(i, j, k - 1)]) / dimM);
				p[IX(i, j, k)] = 0;
			}
		}
	}

	linSolve(0, p, div, 1, 6);

	for (i = 1; i <= dimM; i++) {
		for (j = 1; j <= dimN; j++) {
			for (k = 1; k <= dimO; k++) {
				u[IX(i, j, k)] -= 0.5f*dimM*(p[IX(i + 1, j, k)] - p[IX(i - 1, j, k)]);
				v[IX(i, j, k)] -= 0.5f*dimM*(p[IX(i, j + 1, k)] - p[IX(i, j - 1, k)]);
				w[IX(i, j, k)] -= 0.5f*dimM*(p[IX(i, j, k + 1)] - p[IX(i, j, k - 1)]);
			}
		}
	}

}

/**
Velstep
*/

void FluidCPU::velStep()
{
	addSource(u, uPrev, dt); addSource(v, vPrev, dt); addSource(w, wPrev, dt);
	SWAP(uPrev, u); diffuse(1, u, uPrev, visc, dt);
	SWAP(vPrev, v); diffuse(2, v, vPrev, visc, dt);
	SWAP(wPrev, w); diffuse(3, w, wPrev, visc, dt);
	project(u, v, w, uPrev, vPrev);
	SWAP(uPrev, u); SWAP(vPrev, v); SWAP(wPrev, w);
	advect(1, u, uPrev, uPrev, vPrev, wPrev, dt); advect(2, v, vPrev, uPrev, vPrev, wPrev, dt); advect(3, w, wPrev, uPrev, vPrev, wPrev, dt);
	project(u, v, w, uPrev, vPrev);
}

/**
Destructor
*/
FluidCPU::~FluidCPU()
{
	freeData();
}

