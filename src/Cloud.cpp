// File: Cloud.cpp
// Purpose: Implementation file for all cloud types base class

#include "Cloud.h"

using namespace nimbus;

float* Cloud::noise = nullptr;
GLuint Cloud::textureID;
int Cloud::dim;
GLint Cloud::iNoiseTextCloud;
GLint Cloud::iResolutionCloud;
GLint Cloud::iMouseCloud;
GLint Cloud::iCloudDepthCloud;
GLint Cloud::iSkyTurnCloud;
GLint Cloud::iPosCloud;
GLint Cloud::iCloudPosBlock;
GLint Cloud::iNumSph;
GLint Cloud::iVmin;
GLint Cloud::iVmax;
GLint Cloud::iDebug;
GLint Cloud::iNumClouds;
GLint Cloud::iLowLimits;
GLint Cloud::iUpLimits;
GLint Cloud::iTime;
GLuint Cloud::cloudPosUBO;
GLint Cloud::iVoxelText[MAXCLOUDS];
GLint Cloud::iSunDir;
glm::vec3 Cloud::sunDir;


int Cloud::cloudNum = 0;
int Cloud::cloudOffset = 0;
int Cloud::voxelSize;

glm::vec3 Cloud::vmin[MAXCLOUDS];
glm::vec3 Cloud::vmax[MAXCLOUDS];
GLuint Cloud::lowLimits[MAXCLOUDS];
GLuint Cloud::upLimits[MAXCLOUDS];
GLuint Cloud::voxelTextureID[MAXCLOUDS];

float Cloud::timeDir;

/**
Constructor
*/

Cloud::Cloud()
{
}

/**
Location in 1D array
@param i index i
@param j index j
@param k index k
*/
int Cloud::IX(int i, int j, int k)
{
	return i + dim * j + dim * dim*k;
}

/**
Create cloud base 3D texture
@param dim Texture size
*/
void Cloud::createTexture(int dim)
{

	Cloud::dim = dim;

	noise = new float[dim * dim * dim];

	if (noise == nullptr || dim > MAXTEXSIZE)
		throw NimbusException("Can't allocate data for cloud texture", __FILE__, __FUNCTION__, __LINE__);


	std::random_device rd;
	std::default_random_engine generator;
	generator.seed(rd());
	std::uniform_int_distribution<int> distribution(0, 255);

	// Normalize
	for (int i = 0; i < dim; i++)
		for (int j = 0; j < dim; j++)
			for (int k = 0; k < dim; k++)
				noise[IX(i, j, k)] = distribution(generator) / 255.0f;


	// Pass to OpenGL
	glActiveTexture(GL_TEXTURE0);
	glGenTextures(1, &textureID);
	glBindTexture(GL_TEXTURE_3D, textureID);

	glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_REPEAT);

	glTexImage3D(GL_TEXTURE_3D, 0, GL_RED, dim, dim, dim, 0, GL_RED, GL_FLOAT, noise);
}


/**
Retrieve uniforms
@param shader The cloud shader
*/
void Cloud::getUniforms(Shader& shader)
{
	for (int i = 0; i < MAXCLOUDS; i++)
	{
		std::string id_str = "iVoxel[" + std::to_string(i) + "]";
		iVoxelText[i] = shader.getUniformLocation(id_str.c_str());
	}
	iNoiseTextCloud = shader.getUniformLocation("iNoise");
	iResolutionCloud = shader.getUniformLocation("iResolution");
	iCloudDepthCloud = shader.getUniformLocation("iDepth");
	iSkyTurnCloud = shader.getUniformLocation("iTurn");
	iPosCloud = shader.getUniformLocation("iPos");
	iCloudPosBlock = shader.getUniformBlockIndex("iCloudPosBlock");
	iNumSph = shader.getUniformLocation("iNumSph");
	iVmin = shader.getUniformLocation("iVmin");
	iVmax = shader.getUniformLocation("iVmax");
	iDebug = shader.getUniformLocation("iDebug");
	iNumClouds = shader.getUniformLocation("iNumClouds");
	iLowLimits = shader.getUniformLocation("iLowLimits");
	iUpLimits = shader.getUniformLocation("iUpLimits");
	iTime = shader.getUniformLocation("iTime");
	iSunDir = shader.getUniformLocation("iSunDir");
}

/**
Pass cloud base 3D texture to shader
*/
void Cloud::renderTexture()
{
	// Pass 3D texture
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_3D, textureID);
	glUniform1i(iNoiseTextCloud, 0);
}

/**
Pass this uniforms only once to shader
@param scrW screen width
@param scrH screen height
*/
void Cloud::renderFirstTime(GLint scrW, GLint scrH)
{
	glUniform2f(iResolutionCloud,static_cast<GLfloat>(scrW),static_cast<GLfloat>(scrH));
	glUniform3fv(iVmin, cloudNum, (GLfloat*)vmin);
	glUniform3fv(iVmax, cloudNum, (GLfloat*)vmax);
	glUniform1i(iNumClouds, cloudNum);
	glUniform1iv(iLowLimits, cloudNum, (GLint*)lowLimits);
	glUniform1iv(iUpLimits, cloudNum, (GLint*)upLimits);
	glGenBuffers(1, &cloudPosUBO);
	// Create OpenGL texture
	glGenTextures(cloudNum, voxelTextureID);
		
}

/**
Pass this uniforms each frame to shader
@param mousePose The mouse position
@param time Time increment
@param cloudDepth Cloud distance to camera
@param skyTurn Time of day
@param cameraPos The camera position
@param debud For debug pourposes
*/
void Cloud::render(glm::vec2& mousePos, float time, float cloudDepth, int skyTurn, glm::vec3 cameraPos, bool debug)
{
	glUniform1i(iNumSph, cloudOffset);
	glUniform1f(iTime, time);
	glUniform2f(iMouseCloud, mousePos.x, mousePos.y);
	glUniform1f(iCloudDepthCloud, cloudDepth);
	glUniform1i(iSkyTurnCloud, skyTurn);
	glUniform3f(iPosCloud, cameraPos.x, cameraPos.y, cameraPos.z);
	glUniform1i(iDebug, debug);
	glUniform3f(iSunDir, sunDir.x, sunDir.y, sunDir.z);

}

/**
Render voxel grid texture for precompute light (shading)
@param cloudIndex Cloud index
*/

void Cloud::renderVoxelTexture(int cloudIndex)
{
	glActiveTexture(GL_TEXTURE0 + 2 + cloudIndex);
	// Pass 3D texture
	glBindTexture(GL_TEXTURE_3D, voxelTextureID[cloudIndex]);
	glUniform1i(iVoxelText[cloudIndex], 2 + cloudIndex);

}

/**
Retrieve number of clouds
@return Total number of clouds
*/

int Cloud::getNumClouds()
{
	return cloudNum;
}

/**
Reset cloud counter
*/
void Cloud::resetNumClouds()
{
	glDeleteTextures(cloudNum, voxelTextureID);
	cloudNum = 0;
	cloudOffset = 0;

}

/**
Retrieve number of spheres
@return The spheres offset
*/
int Cloud::getNumSpheres()
{
	return cloudOffset;
}

/** 
Destroy 3D texture
*/
void Cloud::freeTexture()
{
	if (noise != nullptr)
		delete[] noise;
}

/** 
Destructor 
*/
Cloud::~Cloud()
{}