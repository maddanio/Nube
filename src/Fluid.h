// File: Fluid.h
// Purpose: Header file for fluid interface

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "Exception.h"

// Fluid grid dimensions

const int  M = 300;// grid x
const int  N = 50; // grid y
const int  O = 100; // grid z

namespace nimbus
{

	/**
	Fluid interface handler
	*/
	class Fluid
	{
	public:
		virtual void  setUForce(float force, int i, int j, int k) = 0;
		virtual void  setVForce(float force, int i, int j, int k) = 0;
		virtual void  setWForce(float force, int i, int j, int k) = 0;
		virtual float getUForce(int i, int j, int k) = 0;
		virtual float getVForce(int i, int j, int k) = 0;
		virtual float getWForce(int i, int j, int k) = 0;
		virtual void clearUVW() = 0;
		virtual void freeData() = 0;
		virtual void clearData() = 0;
		virtual int allocateData() = 0;
		virtual void sim() = 0;
	};
}