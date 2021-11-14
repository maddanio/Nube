// File: LSystem.cpp
// Purpose: Implementation file for L-System based clouds
// Based on the MATLAB code of the Department of Computing for Neurobiology at Cornell University

#include "Lsystem.h"

using namespace nimbus;


/**
Set grammar axiom
@param axiom Axiom
*/
void Rule::setAxiom(const std::string& axiom)
{
	this->axiom = axiom;
}

/**
Rule antecedent
@param before The rule antecedent 
*/
void Rule::setBefore(const std::string& before)
{
	this->before = before;
}

/**
Rule consequent
@param after The rule consequent
*/
void Rule::setAfter(const std::string& after)
{
	this->after = after;
}

/**
Retrieve axiom
@return The grammar axiom
*/

std::string& Rule::getAxiom()
{
	return axiom;
}

/**
Retrieve grammar antecedent
@return The grammar antecedent
*/

std::string& Rule::getBefore()
{
	return before;
}

/**
Retrieve grammar consequent
@return The grammar consequent
*/

std::string& Rule::getAfter()
{
	return after;
}

/**
L-System constructor. Assume as cumulus cloud.
*/

LSystem::LSystem() : Cumulus::Cumulus()
{
}

/**
Generate derivation
@param axiom The grammar axiom
@param rules A set of rules
@param nReps iterations
@param a The angle
*/

void LSystem::generate(const std::string& axiom, std::vector<Rule>& rules, int nReps, float a)
{

	hT = glm::vec3(1, 0, 0);
	lT = glm::vec3(0, 1, 0);
	uT = glm::vec3(0, 0, 1);

	xT = 0.0f;
	yT = 0.0f;
	zT = 0.0f;


	alpha = glm::radians(a);

	Rup = glm::mat3(glm::vec3(cos(alpha), sin(alpha), 0),
		glm::vec3(-sin(alpha), cos(alpha), 0),
		glm::vec3(0, 0, 1));

	Rup = glm::transpose(Rup);


	Rum = glm::mat3(glm::vec3(cos(-alpha), sin(-alpha), 0),
		glm::vec3(-sin(-alpha), cos(-alpha), 0),
		glm::vec3(0, 0, 1));

	Rum = glm::transpose(Rum);

	Rlp = glm::mat3(glm::vec3(cos(alpha), 0, -sin(alpha)),
		glm::vec3(0, 1, 0),
		glm::vec3(sin(alpha), 0, cos(alpha)));

	Rlp = glm::transpose(Rlp);

	Rlm = glm::mat3(glm::vec3(cos(-alpha), 0, -sin(-alpha)),
		glm::vec3(0, 1, 0),
		glm::vec3(sin(-alpha), 0, cos(-alpha)));

	Rlm = glm::transpose(Rlm);

	Rhp = glm::mat3(glm::vec3(1, 0, 0),
		glm::vec3(0, cos(alpha), -sin(alpha)),
		glm::vec3(0, sin(alpha), cos(alpha)));

	Rhp = glm::transpose(Rhp);

	Rhm = glm::mat3(glm::vec3(1, 0, 0),
		glm::vec3(0, cos(-alpha), -sin(-alpha)),
		glm::vec3(0, sin(-alpha), cos(-alpha)));

	Rhm = glm::transpose(Rhm);

	std::string aux;

	Lsystem = axiom;


	for (int i = 0; i < nReps; i++)
	{
		int j = 0;
		while (j < Lsystem.length())
		{
			bool isRule = false;

			for (int k = 0; k < rules.size(); k++)
			{
				if (Lsystem[j] == rules[k].getBefore()[0])
				{
					aux += rules[k].getAfter();
					isRule = true;
				}
			}
			if (!isRule)
				aux += Lsystem[j];
			j++;
		}

		Lsystem = aux;
		aux.clear();

	}

	std::string s(Lsystem.begin(), Lsystem.end());
	std::cout << s << std::endl;

}

/**
Interpret language string
@param result A vector of spheres
*/

void LSystem::interpret(std::vector<glm::vec4>& result)
{
	srand(static_cast<int>(time(nullptr)));

	for (int i = 0; i < Lsystem.length(); i++)
	{
		char cmdT = Lsystem[i];
		float radius = rand() / static_cast<float>(RAND_MAX) + 0.1f;
		float lenF = radius;

		switch (cmdT)
		{
		case 'X':
		{
			xT += lenF * hT.x;
			yT += lenF * hT.y;
			zT += lenF * hT.z;
			result.push_back(glm::vec4(xT, yT, zT, radius));
			break;
		}
		case 'Y':
		{
			xT += lenF * hT.x;
			yT += lenF * hT.y;
			zT += lenF * hT.z;
			result.push_back(glm::vec4(xT, yT, zT, radius));
			break;
		}
		case 'Z':
		{
			xT += lenF * hT.x;
			yT += lenF * hT.y;
			zT += lenF * hT.z;
			result.push_back(glm::vec4(xT, yT, zT, radius));
			break;
		}
		case '+':
		{
			hT = hT * Rup;
			lT = lT * Rup;
			uT = uT * Rup;
			break;
		}
		case '-':
		{
			hT = hT * Rum;
			lT = lT * Rum;
			uT = uT * Rum;
			break;
		}
		case '&':
		{
			hT = hT * Rlp;
			lT = lT * Rlp;
			uT = uT * Rlp;
			break;
		}
		case '^':
		{
			hT = hT * Rlm;
			lT = lT * Rlm;
			uT = uT * Rlm;
			break;
		}
		case '\\':
		{
			hT = hT * Rhp;
			lT = lT * Rhp;
			uT = uT * Rhp;
			break;
		}
		case '/':
		{
			hT = hT * Rhm;
			lT = lT * Rhm;
			uT = uT * Rhm;
			break;
		}
		case '[':
		{
			stack.push(glm::vec3(xT, yT, zT));
			break;
		}

		case ']':
		{
			glm::vec3 s = stack.top();
			xT = s.x;
			yT = s.y;
			zT = s.z;
			stack.pop();
			break;
		}
		default:
			std::cout << "ERROR in L-System grammar interpretation" << std::endl;
		}
	}
}


/**
Create L-System cloud
@param center The 3D position
@param rules The L-system rules
@param iter Number of iterations
@param alpha The angle
@param radiusFactor The sphere radius factor
*/

void LSystem::create(glm::vec3& center, std::vector<Rule>& rules, int iter, float alpha, float radiusFactor)
{
	std::vector<glm::vec4> result;
	generate(rules[0].getAxiom(), rules, iter, alpha);
	interpret(result);

	std::vector<glm::vec4>::iterator it;



	glm::vec3 max = glm::vec3(-999999999.0f);
	glm::vec3 min = glm::vec3(99999999.0f);

	int i = 0;
	for (it = result.begin(); it != result.end(); it++)
	{
		sphPos[i + cloudOffset].x = (*it).x + center.x;
		sphPos[i + cloudOffset].y = (*it).y + center.y;
		sphPos[i + cloudOffset].z = (*it).z + center.z;
		sphPos[i + cloudOffset].w = 0;

		float rad = sphRads[i] = (*it).w*radiusFactor;

		if (sphPos[i + cloudOffset].x - rad < min.x)
			min.x = sphPos[i + cloudOffset].x - rad;
		if (sphPos[i + cloudOffset].y - rad < min.y)
			min.y = sphPos[i + cloudOffset].y - rad;
		if (sphPos[i + cloudOffset].z - rad < min.z)
			min.z = sphPos[i + cloudOffset].z - rad;


		if (sphPos[i + cloudOffset].x + rad > max.x)
			max.x = sphPos[i + cloudOffset].x + rad;

		if (sphPos[i + cloudOffset].y + rad > max.y)
			max.y = sphPos[i + cloudOffset].y + rad;

		if (sphPos[i + cloudOffset].z + rad > max.z)
			max.z = sphPos[i + cloudOffset].z + rad;

		std::cout << sphPos[i + cloudOffset].x << "  " << sphPos[i + cloudOffset].y << "  " << sphPos[i + cloudOffset].z << " " << sphRads[i] << std::endl;

		if (rad > maxRad)
			maxRad = rad;


		i++;
	}

	numSph = i;

	std::cout << "LSYSTEM-BALLS = " << i << std::endl;

	std::cout << "MAX RADIUS = " << maxRad << std::endl;

	this->lowLimit = lowLimits[cloudNum] = cloudOffset;
	this->upLimit = upLimits[cloudNum] = cloudOffset + i;

	cloudOffset += i;

	vmin[cloudNum] = min;
	vmax[cloudNum] = max;

	std::cout << "BOUNDING BOX(MIN) = " << vmin[cloudNum].x << "  " << vmin[cloudNum].y << "  " << vmin[cloudNum].z << std::endl;
	std::cout << "BOUNDING BOX(MAX) = " << vmax[cloudNum].x << "  " << vmax[cloudNum].y << "  " << vmax[cloudNum].z << std::endl;

	id = cloudNum;

	cloudNum++;

}

