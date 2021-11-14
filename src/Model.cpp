// File: Model.cpp
// Purpose: Implementation file for 3D meshes based clouds

#include "Model.h"

using namespace nimbus;

GLint Model::iAlpha;
GLint Model::iEvolute;
GLint Model::iCloudPosBlockR;
GLuint Model::cloudPosUBOR;
glm::mat4 Model::sphPos[MAXSPHERES];
glm::mat4 Model::cloudPosR[MAXSPHERES];

/**
Constructor
*/

Model::Model()
{
	cloudOffset = 0;
	numSph = 0;
}

/**
Retrieve uniforms
@param shader The model shader object
*/

void Model::getUniforms(Shader& shader)
{
	iAlpha = shader.getUniformLocation("iAlpha");
	iEvolute = shader.getUniformLocation("iEvolute");
	iCloudPosBlockR = shader.getUniformBlockIndex("iCloudPosBlockR");
}


/**
Create model spheres
@param center Model 3D position
@param file OBJ decimated mesh file
@param scl Scale
*/
void Model::create(glm::vec3& center, const std::string& file, float scl)
{

	OBJLoader objLoader;
	std::vector<glm::vec3> vertices;
	std::vector<glm::ivec3> faces;
	std::vector<SPH> spheres;
	std::vector<glm::vec3> barsList;

	glm::vec3 max = glm::vec3(-999999999.0f);
	glm::vec3 min = glm::vec3(99999999.0f);

	objLoader.openObjFile(file);
	objLoader.readVertices(vertices);
	objLoader.readFaces(faces);

	objLoader.generateEllipsoids(vertices, faces, spheres);

	objLoader.generateBarycenter(vertices, faces, barsList);

	std::vector<glm::ivec3>::iterator itf;
	std::vector<SPH>::iterator its;

	glm::vec3 scale = glm::vec3(scl);

	int i = 0;

	for (itf = faces.begin(); itf != faces.end(); itf++)
	{
		glm::vec3 bari = barsList[i];

		// Scale triangle

		glm::vec3 r1 = (vertices.at((*itf).x - 1) - bari)*scale + bari;
		glm::vec3 r2 = (vertices.at((*itf).y - 1) - bari)*scale + bari;
		glm::vec3 r3 = (vertices.at((*itf).z - 1) - bari)*scale + bari;

		float dist1, dist2, dist3;

		// Distance to barycenter of each triangle vertex

		dist1 = glm::distance(bari, r1);
		dist2 = glm::distance(bari, r2);
		dist3 = glm::distance(bari, r3);

		sphPos[i + cloudOffset][0] = glm::vec4(bari + center, 0);

		// Radius

		float rad = glm::max(dist1, glm::max(dist2, dist3));

		sphPos[i + cloudOffset][1] = glm::vec4(dist1, dist2, dist3, 0);
		sphPos[i + cloudOffset][2].x = rad;

		glm::vec3 dirTria, dirElip;

		// Maximum distance from barycenter to each ellipsoid radius (a,b,c)

		if (rad == dist1)
		{
			dirTria = r1 - bari;
			dirElip = bari - glm::vec3(dist1 + bari.x, bari.y, bari.z);
		}
		else if (rad == dist2)
		{
			dirTria = r2 - bari;
			dirElip = bari - glm::vec3(bari.x, dist2 + bari.y, bari.z);
		}
		else
		{
			dirTria = r3 - bari;
			dirElip = bari - glm::vec3(bari.x, bari.y, dist3 + bari.z);
		}

		float theta;

		dirTriaArray[i] = dirTria;

		// Rodrigues' transformation

		glm::mat3 Mat = Tool::rotAtoB(dirElip, dirTria, theta);

		glm::mat4 R = glm::mat4(glm::vec4(Mat[0][0], Mat[0][1], Mat[0][2], 0),
			glm::vec4(Mat[1][0], Mat[1][1], Mat[1][2], 0),
			glm::vec4(Mat[2][0], Mat[2][1], Mat[2][2], 0),
			glm::vec4(0));


		cloudPosR[i + cloudOffset] = R;

		if (sphPos[i + cloudOffset][0].x - rad < min.x)
			min.x = sphPos[i + cloudOffset][0].x - rad;
		if (sphPos[i + cloudOffset][0].y - rad < min.y)
			min.y = sphPos[i + cloudOffset][0].y - rad;
		if (sphPos[i + cloudOffset][0].z - rad < min.z)
			min.z = sphPos[i + cloudOffset][0].z - rad;


		if (sphPos[i + cloudOffset][0].x + rad > max.x)
			max.x = sphPos[i + cloudOffset][0].x + rad;

		if (sphPos[i + cloudOffset][0].y + rad > max.y)
			max.y = sphPos[i + cloudOffset][0].y + rad;

		if (sphPos[i + cloudOffset][0].z + rad > max.z)
			max.z = sphPos[i + cloudOffset][0].z + rad;

		i++;

	}

	numSph = i;

	lowLimits[cloudNum] = cloudOffset;
	upLimits[cloudNum] = cloudOffset + i;

	cloudOffset += i;

	vmin[cloudNum] = min;
	vmax[cloudNum] = max;

	std::cout << "TOTAL POLYGONS=" << i << std::endl;
	std::cout << "BOUNDING BOX(MIN) = " << vmin[cloudNum].x << "  " << vmin[cloudNum].y << "  " << vmin[cloudNum].z << std::endl;
	std::cout << "BOUNDING BOX(MAX) = " << vmax[cloudNum].x << "  " << vmax[cloudNum].y << "  " << vmax[cloudNum].z << std::endl;

	cloudNum++;

}

/**
Precompute light for 3D mesh (shading)
@param precomp The precomputer concrete object
@param sunDir The sun direction vector
@param sunDistance The sun distance to cloud
@param darkLevel Shadow level
@param model1NumSph Number of ellipsoids of model1
@param model2NumSph Number of ellipsoids of model2
*/

void Model::precomputeLight(PrecomputeLight& precomp, glm::vec3& sunDir, float sunDistance, float darkLevel, int model1NumSph, int model2NumSph)
{
	voxelSize = precomp.getVoxelsSize();

	precomp.setTotalSPH(model1NumSph + model2NumSph);
	precomp.precomputeModel(sphPos, model1NumSph, cloudNum, sunDir, sunDistance, darkLevel, vmin, vmax, voxelTextureID);

	Cloud::sunDir = sunDir;
}


/**
Pass uniforms to shader only once
@param offset Elliposids offset
@param evolute Set evolute or involute
*/

void Model::renderFirstTime(int offset, bool evolute)
{
	glUniform1i(iEvolute, evolute);
	vmin[2] = glm::min(vmin[0], vmin[1]);
	vmax[2] = glm::max(vmax[0], vmax[1]);
	glUniform3fv(iVmin, 3, (GLfloat*)vmin);
	glUniform3fv(iVmax, 3, (GLfloat*)vmax);
	glUniform1i(iNumClouds, cloudNum);
	lowLimits[0] = 0;
	upLimits[0] = offset;
	glUniform1iv(iLowLimits, cloudNum, (GLint*)lowLimits);
	glUniform1iv(iUpLimits, cloudNum, (GLint*)upLimits);
	glGenBuffers(1, &cloudPosUBOR);
}

/**
Render 3D mesh cloud
@param shader The shader object
@param cloudR The Rodrigues' rotation matrix
@param cloudPos The ellipsoids vector
@param alpha Linear interpolation alpha
*/
void Model::render(Shader& shader, glm::mat4* cloudR, glm::mat4* cloudPos, float alpha)
{
	glUniform1f(iAlpha, alpha);

	glBindBuffer(GL_UNIFORM_BUFFER, cloudPosUBO);
	glBufferData(GL_UNIFORM_BUFFER, sizeof(glm::mat4) * cloudOffset, cloudPos, GL_DYNAMIC_DRAW);
	glUniformBlockBinding(shader.getProgram(), iCloudPosBlock, 0);
	glBindBufferBase(GL_UNIFORM_BUFFER, 0, cloudPosUBO);

	glBindBuffer(GL_UNIFORM_BUFFER, cloudPosUBOR);
	glBufferData(GL_UNIFORM_BUFFER, sizeof(glm::mat4) * cloudOffset, cloudR, GL_DYNAMIC_DRAW);
	glUniformBlockBinding(shader.getProgram(), iCloudPosBlockR, 1);
	glBindBufferBase(GL_UNIFORM_BUFFER, 1, cloudPosUBOR);
}

/**
Retrieve rotation matrix
@return The Rodrigues' roation matrix
*/

glm::mat4* Model::getCloudPosR()
{
	return cloudPosR;
}

/**
Get number of ellipsoids
@return Number of ellipsoids
*/
int Model::getNumEllipsoids()
{
	return numSph;
}


/**
Get triangle direction
@param i Ellipsoid index
@return The ellipsoid triangle direction
*/
glm::vec3 Model::getDirTriangle(int i)
{
	return dirTriaArray[i];
}

/**
Get ellipsoid data
@param i Ellipsoid index
@return Ellipsoid data
*/

glm::mat4 Model::getEllipsoid(int i)
{
	return sphPos[i];
}

/**
Set ellipsoid rotation matrix
@param i Ellipsoid index
@param r Rodrigues' rotation matrix
*/

void Model::setEllipsoidR(int i, glm::mat4 r)
{
	cloudPosR[i] = r;
}

/**
Get ellipsod rotation matrix
@param i Ellipsoid index
@return The ellipsoid Rodrigues' rotation matrix
*/

glm::mat4 Model::getEllipsoidR(int i)
{
	return cloudPosR[i];
}

/**
Destructor
*/

Model::~Model()
{}