// File: Shader.cpp
// Purpose: Implementation file for OpenGL shader handling

#include "Shader.h"

using namespace nimbus;

/**
Constructor
*/

Shader::Shader()
{
	shaderProgram = 0;
}

/** Loads a shader and returns the compiled shader object.
 If the shader source file could not be opened or compiling the 
 shader fails, then this function raises a NimbusException.
@param shaderType The shader type
@param shaderFile The shader file
 */

void Shader::loadShader(GLenum shaderType, const std::string& shaderFile)
{
	std::cerr << " loading " << shaderFile << std::endl;

	std::ifstream ifs;

	// Load the shader.
	ifs.open(shaderFile);

	if (!ifs)
	{
		std::string errMsg = "Can not open shader file: \"" + shaderFile + "\"";
		throw NimbusException(errMsg.c_str(), __FILE__, __FUNCTION__, __LINE__);
	}

	std::string source(std::istreambuf_iterator<char>(ifs), (std::istreambuf_iterator<char>()));
	ifs.close();

	// Create a shader object.
	GLuint shader = glCreateShader(shaderType);

	// Load the shader source for each shader object.
	const GLchar* sources[] = { source.c_str() };
	glShaderSource(shader, 1, sources, nullptr);

	// Compile the shader.
	glCompileShader(shader);

	// Check for errors
	GLint compileStatus;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &compileStatus);
	if (compileStatus != GL_TRUE)
	{
		GLint logLength;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &logLength);
		GLchar* infoLog = new GLchar[logLength];
		glGetShaderInfoLog(shader, logLength, nullptr, infoLog);
		std::string errMsg = infoLog;
		delete infoLog;
		throw NimbusException(errMsg.c_str(), __FILE__, __FUNCTION__, __LINE__);
	}

	shadersList.push_back(shader);
}

/**
Create a shader program from a set of compiled shader objects
*/
void Shader::createShaderProgram()
{
	// Create a shader program.
	GLuint program = glCreateProgram();

	// Attach the appropriate shader objects.

	std::vector<GLuint>::iterator it;

	for (it = shadersList.begin(); it != shadersList.end(); it++)
		glAttachShader(program, *it);

	// Link the program
	glLinkProgram(program);

	// Check the link status.
	GLint linkStatus;
	glGetProgramiv(program, GL_LINK_STATUS, &linkStatus);
	if (linkStatus != GL_TRUE)
	{
		GLint logLength;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &logLength);
		GLchar* infoLog = new GLchar[logLength];
		glGetProgramInfoLog(program, logLength, nullptr, infoLog);
		std::string errMsg = infoLog;
		delete infoLog;
		throw NimbusException(errMsg.c_str(), __FILE__, __FUNCTION__, __LINE__);
	}

	shaderProgram = program;
	assert(shaderProgram);
}

/**
Get uniform from shader
@param name Uniform variable name
@return Error code
*/
GLint Shader::getUniformLocation(const GLchar* name)
{
	return glGetUniformLocation(shaderProgram, name);
}

/**
Get uniform block index from shader
@param uniformBlockName Uniform block variable name
@return Error code
*/
GLuint Shader::getUniformBlockIndex(const GLchar* uniformBlockName)
{
	return glGetUniformBlockIndex(shaderProgram, uniformBlockName);
}

/**
Retrieve program
@return Shader ID
*/
GLuint Shader::getProgram()
{
	return shaderProgram;
}

/**
Use program
*/
void Shader::useProgram()
{
	glUseProgram(shaderProgram);
}

/** 
Destructor
*/

Shader::~Shader()
{
	glDeleteProgram(shaderProgram);
}
