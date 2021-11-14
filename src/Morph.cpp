// File: Morph.cpp
// Purpose: Implementation file for morphing

#include "Morph.h"

using namespace nimbus;

/**
Constructor
*/

Morph::Morph()
{}

/**
Establish models
@param modelA Source model
@param modelB Destination model
@param evolute Set for evolution or involution
*/

void Morph::setModels(Model* modelA, Model* modelB, bool evolute)
{
	this->modelA = modelA;
	this->modelB = modelB;

	for (int i = 0; i < modelA->getNumEllipsoids(); i++)
	{

		if (evolute) // In case of evolution
			dirTriaSrc[i] = modelA->getDirTriangle(i);
		else // In case of involution
			dirTriaDst[i] = modelA->getDirTriangle(i);
	}

	for (int i = 0; i < modelB->getNumEllipsoids(); i++)
	{

		if (evolute) // In case of evolution
			dirTriaDst[i] = modelB->getDirTriangle(i);
		else // In case of involution
			dirTriaSrc[i] = modelB->getDirTriangle(i);
	}

	offsetSrc = modelA->getNumEllipsoids();
	offsetDst = modelB->getNumEllipsoids();

}

/**
Prepare for evolution
*/
void Morph::prepareMorphEvolute()
{

	for (int i = 0; i < offsetDst; i++)
	{
		if (i >= offsetSrc)
		{
			cloudPosDst[i] = cloudPosSrc[i] = Model::getEllipsoid(i % offsetSrc);
			dirTriaSrc[i] = dirTriaSrc[i % offsetSrc];
			Model::setEllipsoidR(i, Model::getEllipsoidR(i % offsetSrc));
		}
		else
			cloudPosDst[i] = cloudPosSrc[i] = Model::getEllipsoid(i);

		thetaSrc[i] = 0.0f;
		thetaDst[i] = acos(glm::dot(dirTriaSrc[i], dirTriaDst[i]) / (glm::length(dirTriaSrc[i])*length(dirTriaDst[i])));

	}
}

/**
Prepare for involution
*/

void Morph::prepareMorphInvolute()
{


	for (int i = 0; i < offsetDst; i++)
	{

		cloudPosDst[i] = cloudPosSrc[i] = Model::getEllipsoid(i + offsetSrc);

		thetaSrc[i] = 0.0f;
		thetaDst[i] = acos(glm::dot(dirTriaSrc[i], dirTriaDst[i % offsetSrc]) / (glm::length(dirTriaSrc[i])*length(dirTriaDst[i % offsetSrc])));

	}
}

/**
Render morphing for evolute
@param alpha Linear interpolation alpha
*/
void Morph::morphEvolute(float alpha)
{
	for (int i = 0; i < offsetDst; i++)
		cloudPosDst[i] = glm::mix(cloudPosSrc[i], Model::getEllipsoid(i + offsetSrc), alpha);
}

/**
Render morphing for involute
@param alpha Linear interpolation alpha
*/

void Morph::morphInvolute(float alpha)
{
	for (int i = 0; i < offsetDst; i++)
		cloudPosDst[i] = glm::mix(cloudPosSrc[i], Model::getEllipsoid(i % offsetSrc), alpha);
}

/**
Main render
@param Morphing step
*/

void Morph::morph(float step)
{
	for (int i = 0; i < offsetDst; i++)
	{
		if (thetaSrc[i] < thetaDst[i])
		{
			thetaSrc[i] += step;

			glm::mat3 Mat = Tool::rotAtoBMorph(dirTriaSrc[i], dirTriaDst[i % offsetSrc], thetaSrc[i]);
			cloudPosRDst[i] = glm::mat4(glm::vec4(Mat[0][0], Mat[0][1], Mat[0][2], 0),
				glm::vec4(Mat[1][0], Mat[1][1], Mat[1][2], 0),
				glm::vec4(Mat[2][0], Mat[2][1], Mat[2][2], 0),
				glm::vec4(0)) * Model::getEllipsoidR(i);
		}
	}
}

/**
Retrieve rotation matrix for final mesh
@return The destination Rodrigues' rotation matrix
*/

glm::mat4* Morph::getCloudPosRDst()
{
	return cloudPosRDst;
}

/** Retrieve ellipsoid positions for final mesh
@return The ellipsoid destination position
*/
glm::mat4* Morph::getCloudPosDst()
{
	return cloudPosDst;
}

/**
Destructor
*/
Morph::~Morph()
{}