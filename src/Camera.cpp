// File: Camera.cpp
// Purpose: Implements a basic camera class

#include "Camera.h"

using namespace nimbus;

/**
Default constructor
*/

Camera::Camera()
	: viewport(0)
	, position(0)
	, rotation()
	, projectionMatrix(1)
	, viewMatrix(1)
	, viewDirty(false)
{}

/**
Viewport Camera constructor
@param screenWidth The screen width
@param screenHeight The screen height
*/
Camera::Camera(int screenWidth, int screenHeight)
	: viewport(0, 0, screenWidth, screenHeight)
	, position(0)
	, rotation()
	, projectionMatrix(1)
	, viewMatrix(1)
	, viewDirty(false)
{

}

/**
Set viewport method
@param x X position of viewport
@param y Y position of viewport
@param width Width of viewport
@param height Height of viewport
*/
void Camera::setViewport(GLint x, GLint y, GLsizei width, GLsizei height)
{
	viewport = glm::vec4(x, y, width, height);
	glViewport(x, y, 2 * width, 2 * height);
}

/**
Retrieve viewport size
@return The viewport size
*/
glm::vec4 Camera::getViewport() const
{
	return viewport;
}

/**
Set projection matrix
@param fov FOV
@param aspectoRation Frustum perspective aspect ration
@param zNear Front plane
@param zFar Far plane
*/
void Camera::setProjectionRH(float fov, float aspectRatio, float zNear, float zFar)
{
	projectionMatrix = glm::perspective(glm::radians(fov), aspectRatio, zNear, zFar);
}

/**
Apply view matrix
*/

void Camera::applyViewMatrix()
{
	updateViewMatrix();
}

/**
Set camera look at
@param srcPos Camera source position
@param dstPos Camera target look at
*/

void Camera::setLookAt(const glm::vec3& srcPos, const glm::vec3& dstPos)
{
	this->srcPos = srcPos;
	this->lookAt = dstPos;
}

/**
Return camera look-at source position
@return The source look at position
*/
glm::vec3 Camera::getSrcLookAt()
{
	return srcPos;
}


/**
Return camera look-at position
@return The look at position
*/
glm::vec3 Camera::getDstLookAt()
{
	return lookAt;
}



/**
Retrieve camera direction
@return Camera direction vector
*/

glm::vec3 Camera::getDirection()
{
	return direction;
}

/**
Set camera position
@param pos The 3D position
*/

void Camera::setPosition(const glm::vec3& pos)
{
	position = pos;
}


/** Retrieve camera position
@return Retrieve camera position
*/

glm::vec3 Camera::getPosition() const
{
	return position;
}

/** Translate camera
@param translate Translation vector
*/

void Camera::translate(const glm::vec3& translate)
{
	translation = translate;
	position = translate;
	viewDirty = true;
}

/**
Set camera rotation
@param Rotational quaternion
*/

void Camera::setRotation(const glm::quat& rot)
{
	rotation = rot;
	viewDirty = true;
}

/**
Retrieve camera rotation
@return Camera rotational quaternion
*/

glm::quat Camera::getRotation() const
{
	return rotation;
}

/**
Rotate camera
@param rot Rotation quaternion
*/
void Camera::rotate(const glm::quat& rot)
{
	rotation = rotation * rot;
	viewDirty = true;
}

/**
Retrieve projection matrix
@return The projection matrix
*/
glm::mat4 Camera::getProjectionMatrix()
{
	return projectionMatrix;
}

/**
Retrieve view matrix
@return The view matrix
*/
glm::mat4 Camera::getViewMatrix()
{
	updateViewMatrix();
	return viewMatrix;
}

/** 
Set view Matrix
@param The view matrix
*/
void Camera::setViewMatrix(glm::mat4& viewMatrix)
{
	viewMatrix = viewMatrix;
	viewDirty = true;
	updateViewMatrix();
}

/**
Update view matrix
*/
void Camera::updateViewMatrix()
{
	if (viewDirty)
	{
		direction = glm::normalize(lookAt - srcPos);

		glm::vec3 cameraRight = glm::normalize(glm::cross(glm::vec3(0.0, 1.0, 0.0), direction));

		glm::vec3 cameraUp = glm::cross(direction, cameraRight);

		glm::mat4 mat = glm::mat4(glm::vec4(cameraRight, 0.0), glm::vec4(cameraUp, 0.0), glm::vec4(direction, 0.0), glm::vec4(0.0, 0.0, 0.0, 1.0));

		glm::mat4 translate = glm::translate(translation);

		glm::mat4 rotate = glm::transpose(glm::toMat4(rotation));

		viewMatrix = mat * rotate * translate;

		viewDirty = false;
	}
}
