// File: Cumulus.cpp
// Purpose: Implementation file for cumulus class

#include "Cumulus.h"

using namespace nimbus;

GLint Cumulus::iWindDirection;
GLint Cumulus::iFlat;
int Cumulus::dimM;
int Cumulus::dimN;
int Cumulus::dimO;
float Cumulus::maxRad;
glm::vec3 Cumulus::extinctionLimits;
Winds Cumulus::windDirection;
glm::vec4 Cumulus::sphPos[MAXSPHERES];
GLint Cumulus::flat[MAXCLOUDS];
int Cumulus::totalSph = 0;

/**
Constructor
*/

Cumulus::Cumulus()
{
	maxRad = -99999.9f;
	createdCloud = false;
	extinctedCloud = false;
	radIndex = 0.0f;
	timeDir = -0.06f;
	dimM = M;
	dimN = N;
	dimO = O;
}

/**
Create cumulus
@param spheres Number of spheres to generate cloud
@param siz Size of Gaussian cloud
@param center Cloud position
@param nuX Gaussian X-mean
@param sigX Gaussian X-standard deviation
@param nuY Gaussian Y-mean
@param sigY Gaussian Y-standard deviation
@param nuZ Gaussian Z-mean
@param sigZ Gaussian Z-standard deviation
@param isFlat Cloud with level of condensation
@param optimze Optimze cloud spheres
*/

void Cumulus::create(int spheres, GLfloat siz, glm::vec3 center, GLfloat nuX, GLfloat sigX, GLfloat nuY, GLfloat sigY, GLfloat nuZ, GLfloat sigZ, bool isFlat, bool optimize)
{

	std::random_device rd;

	std::default_random_engine generator;

	generator.seed(rd());

	std::normal_distribution<GLfloat> distribution_x(nuX, sigX);
	std::normal_distribution<GLfloat> distribution_y(nuY, sigY);
	std::normal_distribution<GLfloat> distribution_z(nuZ, sigZ);

	glm::vec3 max = glm::vec3(-999999999.0f);
	glm::vec3 min = glm::vec3(99999999.0f);

	std::cout << "X     Y       Z        R" << std::endl;

	int valid = 0;

	numSph = spheres;


	for (int i = 0; i < spheres; i++)
	{

		float a = glm::clamp(distribution_x(generator), -2.0f * sigX, 2.0f * sigX);
		float b = glm::clamp(distribution_y(generator), 0.0f, 2.0f * sigY);
		float c = glm::clamp(distribution_z(generator), -2.0f * sigZ, 2.0f * sigZ);


		if ((optimize) && !((abs(a) >= 1.5*sigX) || (b <= 0.5) || (b >= 1.5*sigY) || (abs(c) >= 1.5*sigZ))) continue;
		else
		{

			sphPos[valid + cloudOffset].x = center.x + a;
			sphPos[valid + cloudOffset].y = center.y + b;
			sphPos[valid + cloudOffset].z = center.z + c;



			sphPos[valid + cloudOffset].w = sphRads[i] = siz * (1.0f - 0.1f*glm::sqrt(glm::pow((sphPos[valid + cloudOffset].x - center.x) / (2.0f*sigX), 2.0f) + glm::pow((sphPos[valid + cloudOffset].y - center.y) / (2.0f*sigY), 2.0f) + glm::pow((sphPos[valid + cloudOffset].z - center.z) / (2.0f*sigZ), 2.0f)));

			std::cout << " sphPos = " << sphPos[valid + cloudOffset].x << "   " << sphPos[valid + cloudOffset].y << "   " << sphPos[valid + cloudOffset].z << "   " << sphRads[i] << std::endl;


			if (sphRads[i] > maxRad)
				maxRad = sphRads[i];


			if (sphPos[valid + cloudOffset].x - sphRads[i] < min.x)
				min.x = sphPos[valid + cloudOffset].x - sphRads[i];
			if (sphPos[valid + cloudOffset].y - sphRads[i] < min.y)
				min.y = sphPos[valid + cloudOffset].y - sphRads[i];
			if (sphPos[valid + cloudOffset].z - sphRads[i] < min.z)
				min.z = sphPos[valid + cloudOffset].z - sphRads[i];


			if (sphPos[valid + cloudOffset].x + sphRads[i] > max.x)
				max.x = sphPos[valid + cloudOffset].x + sphRads[i];

			if (sphPos[valid + cloudOffset].y + sphRads[i] > max.y)
				max.y = sphPos[valid + cloudOffset].y + sphRads[i];

			if (sphPos[valid + cloudOffset].z + sphRads[i] > max.z)
				max.z = sphPos[valid + cloudOffset].z + sphRads[i];

			valid++;
		}

	}

	std::cout << "MAX RADIUS = " << maxRad << std::endl;


	this->lowLimit = lowLimits[cloudNum] = cloudOffset;
	this->upLimit = upLimits[cloudNum] = cloudOffset + valid;


	cloudOffset += valid;

	vmin[cloudNum] = min;
	vmax[cloudNum] = max;

	if (isFlat)
	{
		vmin[cloudNum].y+=0.5f;
		flat[cloudNum] = true;
	}
	else flat[cloudNum] = false;

	std::cout << "TOTAL SPH = " << cloudOffset << std::endl;
	std::cout << "BOUNDING BOX(MIN) = " << vmin[cloudNum].x << "  " << vmin[cloudNum].y << "  " << vmin[cloudNum].z << std::endl;
	std::cout << "BOUNDING BOX(MAX) = " << vmax[cloudNum].x << "  " << vmax[cloudNum].y << "  " << vmax[cloudNum].z << std::endl;

	id = cloudNum;

	cloudNum++;


}

/**
Retrieve uniforms
@param shader The cloud shader object
*/

void Cumulus::getUniforms(Shader& shader)
{
	iWindDirection = shader.getUniformLocation("iWindDirection");
	iFlat = shader.getUniformLocation("isFlat");
}


/**
Calculate guide point
@param windDirection The wind direction
*/
void Cumulus::setGuidePoint(Winds windDirection)
{
	float maxX = -99999.9f;
	float minX = 99999.9f;
	float maxZ = -999999.9f;
	float minZ = 99999.9f;


	for (GLuint i = lowLimit; i < upLimit; i++)
	{

		if ((sphPos[i].x > maxX) && (windDirection == WEST))
		{
			maxX = sphPos[i].x;
			sphPivot = i;
		}
		else if ((sphPos[i].x < minX) && (windDirection == EAST))
		{
			minX = sphPos[i].x;
			sphPivot = i;
		}
		else if ((sphPos[i].z > maxZ) && (windDirection == SOUTH))
		{
			maxZ = sphPos[i].z;
			sphPivot = i;
		}
		else if ((sphPos[i].z < minZ) && (windDirection == NORTH))
		{
			minZ = sphPos[i].z;
			sphPivot = i;
		}
	}

	Cumulus::windDirection = windDirection;

	std::cout << "FINAL PIVOT = " << sphPos[sphPivot].x << "   " << sphPos[sphPivot].y << "   " << sphPos[sphPivot].z << "   " << sphPos[sphPivot].w << std::endl;

}

/**
Check for cloud guide point location and apply wind forces to grid tunnel
@param force Wind force
@param windGrid Pointer to wind 3D grid
*/

void Cumulus::applyWind(float force, Fluid* windGrid)
{

	glm::vec4 p = sphPos[sphPivot] + glm::vec4(dimM / 2, dimN / 2, dimO / 2, 0);

/*
	if ((p.x < 1 || p.x > dimM || p.y < 1 || p.y > dimN || p.z < 1 || p.z > dimO) && !extinctedCloud)
	{
		std::cout << "****CLOUD " << id << "OUT OF BOUND****" << std::endl;
		extinctedCloud = true;
		std::cout << "PIVOT APPLY WIND = " << p.x << "   " << p.y << "   " << p.z << "   " << std::endl;
		return;
	}
*/
	for (int i = 1; i < dimM; i++)
		for (int j = 1; j < dimN; j++)
			for (int k = 1; k < dimO; k++)
				switch (windDirection)
				{
				case NORTH:
					windGrid->setWForce(-force, i, j, k);
					break;
				case SOUTH:
					windGrid->setWForce(force, i, j, k);
					break;
				case EAST:
					windGrid->setUForce(-force, i, j, k);
					break;
				case WEST:
					windGrid->setUForce(force, i, j, k);
					break;
				}
}

/**
Borning cloud animation
@param step Borning step
*/
void Cumulus::bornCloud(float step)
{
	for (int i = 0; i < numSph; i++)
		if (sphPos[i + lowLimit].w < sphRads[i])
			sphPos[i + lowLimit].w += step;

}

/**
Extinct cloud animation
@param step Extinction step
*/
void Cumulus::extinctCloud(float step)
{
	for (int i = 0; i < numSph; i++)
		if (sphPos[i + lowLimit].w > 0.0)
		{
			sphPos[i + lowLimit].w -= step;
			if (sphPos[i + lowLimit].w <= 0.1)
				sphPos[i + lowLimit].w = 0.0;
		}

}

/**
Calculate cloud guide point position according to fluid grid state
@param windGrid Pointer to 3D wind grid
*/
void Cumulus::computeWind(Fluid* windGrid)
{
	glm::ivec4 p = sphPos[sphPivot] + glm::vec4(dimM / 2, dimN / 2, dimO / 2, 0);


	float auxU = windGrid->getUForce(p.x, p.y, p.z);
	float auxV = windGrid->getVForce(p.x, p.y, p.z);
	float auxW = windGrid->getWForce(p.x, p.y, p.z);

	for (int k = 0; k < numSph; k++)
	{

		sphPos[k + lowLimit].x += auxU;

		sphPos[k + lowLimit].y += auxV;

		sphPos[k + lowLimit].z += auxW;
	}

	vmin[id] += glm::vec3(auxU, auxV, auxW);
	vmax[id] += glm::vec3(auxU, auxV, auxW);

}

/**
Pass winds to shader
@param winDirection wind direction
*/
void Cumulus::setWind(Winds winDirection)
{

	switch (windDirection)
	{

	case SOUTH:
		timeDir = -0.1f;
		glUniform3f(iWindDirection, 0.0f, 0.0f, 1.0f);
		break;
	case NORTH:
		timeDir = 0.1f;
		glUniform3f(iWindDirection, 0.0f, 0.0f, 1.0f);
		break;
	case EAST:
		timeDir = 0.06f;
		glUniform3f(iWindDirection, 1.0f, 0.0f, 0.0f);
		break;
	case WEST:
		timeDir = -0.06f;
		glUniform3f(iWindDirection, 1.0f, 0.0f, 0.0f);
		break;
	}
}

/**
Render cumulus
@param shader The cloud shader object
*/
void Cumulus::render(Shader& shader)
{
	glUniform1iv(iFlat, cloudNum, (GLint*)flat);

	glBindBuffer(GL_UNIFORM_BUFFER, cloudPosUBO);
	glBufferData(GL_UNIFORM_BUFFER, sizeof(glm::vec4) * cloudOffset, sphPos, GL_DYNAMIC_DRAW);
	glUniformBlockBinding(shader.getProgram(), iCloudPosBlock, 0);

	glBindBufferBase(GL_UNIFORM_BUFFER, 0, cloudPosUBO);

	glUniform3fv(iVmin, cloudNum, (GLfloat*)vmin);
	glUniform3fv(iVmax, cloudNum, (GLfloat*)vmax);

}


/**
Precompute light for cumulus
@param precomp The light precomputer concrete object
@param sundDir The sun direction vector
@param sunDistance The distance to sun
@param darkLevel Shadow level

*/

void Cumulus::precomputeLight(PrecomputeLight& precomp, glm::vec3 sunDir, float sunDistance, float darkLevel)
{
	voxelSize = precomp.getVoxelsSize();
	
	if (cloudOffset != totalSph)
	{
		precomp.setTotalSPH(cloudOffset);
		totalSph = cloudOffset;
	}

	precomp.precomputeCloud(sphPos, cloudOffset, cloudNum, sunDir, sunDistance, darkLevel, vmin, vmax, voxelTextureID);

	Cloud::sunDir = sunDir;
}

/**
Set scene extinction limits (horizons)
@param limits Set wind grid limits
*/

void Cumulus::setLimits(glm::vec3 limits)
{
	extinctionLimits = limits;
}

/**
Check extinction limits
@return True if out of limits
*/

bool Cumulus::checkLimits()
{
	return ((sphPos[sphPivot].x < extinctionLimits.x) &&
		(sphPos[sphPivot].x > -extinctionLimits.x) &&
		(sphPos[sphPivot].z < extinctionLimits.z) &&
		(sphPos[sphPivot].z > -extinctionLimits.z));
}

/**
Get maximum cloud radius for animation
@return The maximum radius
*/
float Cumulus::getMaxRadius()
{
	return maxRad;
}

/**
Get animation cloud born index status
@return Borning radius
*/
float Cumulus::getRadIndex()
{
	return radIndex;
}

/**
Increment animation cloud born index
*/
void Cumulus::incRadIndex(float step)
{
	radIndex += step;
}

/**
Checks if cloud extincted
@return True if cloud extincted
*/

bool Cumulus::getExtincted()
{
	return extinctedCloud;
}

/**
Check if cloud created
@return True if cloud created
*/

bool Cumulus::getCreated()
{
	return createdCloud;
}

/**
Sets if cloud is created
*/

void Cumulus::setCreated(bool created)
{
	createdCloud = created;
}

/**
Set time evolution
@param Time direction
*/
void Cumulus::setTimeDir(float timeDir)
{
	Cumulus::timeDir = timeDir;
}

/**
Get time evolution
@return The time direction
*/
float Cumulus::getTimeDir()
{
	return Cumulus::timeDir;
}

/**
Retrieve guide point index in the cloud
@return Cloud sphere pivot
*/
int Cumulus::getGuidePoint()
{
	return sphPivot;
}

/**
Destructor
*/
Cumulus::~Cumulus()
{
}