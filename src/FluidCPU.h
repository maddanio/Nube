// File: FluidCPU.h
// Purpose: Header file for fluid dynamic class for CPU
// Based on the article of Jos Stam

#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "Fluid.h"

namespace nimbus
{
	/**
	 Fluid dynamic class for CPU
	*/
	class FluidCPU : public Fluid
	{
	private:
		//fluid field information
		int dimM; // grid x
		int dimN; // grid y
		int dimO; // grid z
		float dt; // time delta
		float diff; // diffuse
		float visc; // viscosity
		int size;
		float *u, *v, *w, *uPrev, *vPrev, *wPrev;
	private:
		int  IX(int i, int j, int k);
		void addSource(float* x, float* s, float dt);
		void linSolve(int b, float* x, float* x0, float a, float c);
		void diffuse(int b, float* x, float* x0, float diff, float dt);
		void advect(int b, float* d, float* d0, float* u, float* v, float* w, float dt);
		void project(float * u, float * v, float * w, float * p, float * div);
		void velStep();
	public:
		FluidCPU(float dt, float diff, float visc);
		void  setUForce(float force, int i, int j, int k);
		void  setVForce(float force, int i, int j, int k);
		void  setWForce(float force, int i, int j, int k);
		float getUForce(int i, int j, int k);
		float getVForce(int i, int j, int k);
		float getWForce(int i, int j, int k);
		void clearUVW();
		void freeData();
		void clearData();
		int  allocateData();
		void sim();
		~FluidCPU();
	};
}
