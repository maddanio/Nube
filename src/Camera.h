// File: Camera.h
// Purpose: Basic camera class

#pragma once




#  include <OpenGL/gl3.h>
#  include <OpenGL/gl3ext.h>


#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/quaternion.hpp>
namespace nimbus
{
	/**
	Scene camera class
	*/
	class Camera
	{
	protected:
		void updateViewMatrix();
		glm::vec4 viewport;
		glm::vec3 srcPos;
		glm::vec3 lookAt;
		glm::vec3 position;
		glm::vec3 direction;
		glm::quat rotation;
		glm::mat4 rot;
		glm::vec3 translation;
		glm::mat4 viewMatrix;
		glm::mat4 projectionMatrix;
	private:
		bool viewDirty;
	public:
		Camera();
		Camera(int screenWidth, int screenHeight);
		void setViewport(GLint x, GLint y, GLsizei width, GLsizei height);
		glm::vec4 getViewport() const;
		void setProjectionRH(float fov, float aspectRatio, float zNear, float zFar);
		void applyViewMatrix();
		void setLookAt(const glm::vec3& srcPos, const glm::vec3& dstPos);
		glm::vec3 getSrcLookAt();
		void setPosition(const glm::vec3& pos);
		glm::vec3 getPosition() const;
				glm::vec3 getDstLookAt();
		glm::vec3 getDirection();
		void translate(const glm::vec3& translate);
		void setRotation(const glm::quat& rot);
		glm::quat getRotation() const;
		void rotate(const glm::quat& rot);
		glm::mat4 getProjectionMatrix();
		glm::mat4 getViewMatrix();
		void setViewMatrix(glm::mat4& viewMatrix);
	};
}