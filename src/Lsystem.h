// File: LSystem.h
// Purpose: Header file for L-System based clouds
// Based on the MATLAB code of the Department of Computing for Neurobiology at Cornell University

#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <stack>

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/quaternion.hpp>

#include "Cumulus.h"

namespace nimbus
{

	/**
	BNF grammar rule class
	*/
	class Rule
	{
	private:
		std::string axiom; // Grammar axiom
		std::string before; // Rule precondition
		std::string after;  // Rule postcondition
	public:
		void setAxiom(const std::string& axiom);
		void setBefore(const std::string& before);
		void setAfter(const std::string& after);

		std::string& getAxiom();
		std::string& getBefore();
		std::string& getAfter();

	};

	/**
	The LSystem fractal cloud generator class
	*/
	class LSystem : public Cumulus
	{
	private:
		std::string Lsystem;
		std::stack<glm::vec3> stack;

		float xT;
		float yT;
		float zT;
		float alpha;

		glm::vec3 hT;
		glm::vec3 lT;
		glm::vec3 uT;

		glm::mat3 Rup;
		glm::mat3 Rum;
		glm::mat3 Rlp;
		glm::mat3 Rlm;
		glm::mat3 Rhp;
		glm::mat3 Rhm;

	private:
		void generate(const std::string& axiom, std::vector<Rule>& rules, int nReps, float a);
		void interpret(std::vector<glm::vec4>& result);

	public:
		LSystem();
		void create(glm::vec3& center, std::vector<Rule>& rules, int iter, float alpha, float radiusFactor);
	};
}
