////////////////////////////////////////////////////////////////////
// ONTOGENETIC MODEL FOR REAL-TIME VOLUMETRIC CLOUDS SIMULATION THESIS			
// Software Engineering and Computer Systems Deparment	
// National University for Distance Education (UNED)			    		
// (c) Carlos Jiménez de Parga, PhD student.
// Last revision 18/01/2019
// Version 1.0
//////////////////////////////////////////////////////////////////

#include <GLFW/glfw3.h>
#include <stdexcept>
#define _CRTDBMAP_ALLOC // For debug purposes
#include <stdlib.h>

#include <time.h>
#include <random>

#include "Defines.h"
#include "Main.h"
#include "Camera.h"
#include "FluidCPU.h"
#ifdef CUDA
#include "FluidCUDA.h"
#include "PrecomputeCUDA.h"
#endif
#ifdef CPU
#include "PrecomputeCPU.h"
#endif

struct AtExit // For debug purposes
{
	~AtExit() {}
} doAtExit;

// Sun/moon light direction
glm::vec3 sunDir = glm::normalize(glm::vec3(0.0, -0.1, -1.0));
const float windForce = 0.009f;
const int FLUIDLIMIT = 100;
bool parallel = true;
bool mean = false;

nimbus::Shader shaderCloud;
nimbus::Shader shaderAxis;
nimbus::Shader shaderSky;
nimbus::Canvas canvas;
#ifdef MOUNT
nimbus::Mountain mountain;
#endif
nimbus::Axis axis;
nimbus::Model model1;
nimbus::Model model2;
nimbus::Morph morphing;


nimbus::Cumulus myCloud[3]; // Create cumulus clouds

glm::vec3 EXTINCTION_LIMIT(100.0, 15.0, 30.0); // Limit of scenary
const GLint SCR_W = 1200;  // Screen dimensions
const GLint SCR_H = 600;
const GLfloat SCR_Z = 480.0;  // Depth of camera //AXIS_Z-SCR_Z 890.0
const int TEXTSIZ = 128; // Size of cloud base texture for a later fBm
bool EVOLUTE = true; // If morphing evolution

float alpha; // Linear interpolation increment for morphing
float alphaDir; // Linear interpolation direction for morphing

// Wind direction
nimbus::Winds windDirection = nimbus::Winds::EAST;

// 30 FPS minimum real-time
const int32_t FPS = 1000 / 35;


int windowWidth =  SCR_W;
int windowHeight = SCR_H;
int windowHandle = 0;

// Required cameras

nimbus::Camera cameraFrame;
nimbus::Camera cameraSky;
nimbus::Camera cameraAxis;

// Main camera positions
glm::vec3 userCameraPos;

// Camera initial position
glm::vec3 initialCameraFramePosition;
glm::vec2 mousePos;


//////////////////////////////////////// FLUIDS ////////////////////////////////////////////

float dt = 0.4f; // time delta
float diff = 0.0f; // diffuse
float visc = 0.00001f; // viscosity

#ifdef CUDA
nimbus::FluidCUDA windGridCUDA(dt, diff, visc);
nimbus::PrecomputeCUDA precompCUDA;
#endif
#ifdef CPU
nimbus::FluidCPU windGridCPU(dt, diff, visc);
nimbus::PrecomputeCPU precompCPU;
#endif

float cloudDepth = 20.0f; // Distance from viewer
float frameCount = 0.0; //FPS measurement
int previousTime;

int  skyTurn = 0; // Day hour
float timeDelta = 99999.0f; // Timers 
float timeDir = -0.06f;
int simcont = 1000;

const int TOTALTIME = 5;
const int PRECOMPTIMEOUT = 150;

int totalTime = 0; // Time to check born/extinction
int precomputeTimeOut = PRECOMPTIMEOUT; // Time for regular precompute light

bool debug = false; // If debugging
bool onPlay = true; // The loop is idle
bool firstPass = true; // First pass in loop

// Function prototypes

void reshapeGL(GLFWwindow*, int w, int h);
void displayGL();
void idleGL();
void keyboardGL(GLFWwindow*, int key, int scancode, int, int);
void keyboardUpGL(unsigned char c);
void specialGL(int key, int x, int y);
void specialUpGL(int key, int x, int y);
void mouseGL(int button, int state, int x, int y);
void mouseWheel(GLFWwindow* window, double x, double y);
void motionGL(GLFWwindow* window, double x, double y);
void applyWind();
void syncFPS();


// Initialize the OpenGL context and create a render window.

void GLFWerrorCallback(int code, const char *message)
{
	throw std::runtime_error{"glfw error " + std::string(message)};
}

GLFWwindow* initGL(int argc, char* argv[])
{
	std::cerr << "Initialize glut..." << std::endl;

    if (!glfwInit())
    {
        throw std::runtime_error("glfw won't start");
    }

    glfwSetErrorCallback(GLFWerrorCallback);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    glfwWindowHint(GLFW_RED_BITS, 8);
    glfwWindowHint(GLFW_GREEN_BITS, 8);
    glfwWindowHint(GLFW_BLUE_BITS, 8);
    glfwWindowHint(GLFW_ALPHA_BITS, 8);
    glfwWindowHint(GLFW_DEPTH_BITS, 8);

    auto window = glfwCreateWindow(windowWidth, windowHeight, "nube", NULL, NULL);

	// Register GLUT callbacks.
    glfwSetKeyCallback(window, keyboardGL);
    glfwSetScrollCallback(window, mouseWheel);
    glfwSetCursorPosCallback(window, motionGL);
    glfwSetWindowSizeCallback(window, reshapeGL);

	//glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);

	glfwMakeContextCurrent(window);
    GLint framebuffer_id;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &framebuffer_id);
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffer_id);
	glClearColor(0.46f, 0.78f, 0.97f, 0.0f);
	glClearDepth(1.0f);
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_TEXTURE_2D);
	glEnable(GL_TEXTURE_3D);
	glFrontFace(GL_CW);
	std::cerr << "Initialize OpenGL Success!" << std::endl;
	return window;
}


// Initialize Glew

// Entry point for cumulus rendering
#ifdef CUMULUS

int main(int argc, char* argv[])
{
	// Initialize FreeGLUT and Glew
	auto window = initGL(argc, argv);

	std::cerr << "run..." << std::endl;

	try {
		// Cameras setup
		cameraFrame.setProjectionRH(30.0f, SCR_W / SCR_H, 0.1f, 2000.0f);
		cameraAxis.setProjectionRH(30.0f, SCR_W / SCR_H, 0.1f, 2000.0f);
		cameraFrame.setViewport(0, 0, SCR_W, SCR_H);
		cameraFrame.setLookAt(glm::vec3(0, 0, -SCR_Z), glm::vec3(0, 0, SCR_Z));
		cameraFrame.translate(glm::vec3(-SCR_W / 2.0,-SCR_H / 2.0, -SCR_Z));
	
		userCameraPos = glm::vec3(0.0, 0.4, 0.0);

		// Create framgent shader canvas
		canvas.create(SCR_W, SCR_H);
		// Create cloud base texture
		nimbus::Cloud::createTexture(TEXTSIZ);

#ifdef MOUNT
		mountain.create(800.0, false); // Create mountain
#endif
		axis.create(); // Create 3D axis

		//Create cumulus clouds
		myCloud[0].create(35, 2.8f, glm::vec3(0.0, 5.0, 0.0), 0.0f, 3.0f, 0.0f, 1.9f, 0.0f, 3.0f, true, false);

		// Calculate guide points for cumulus
		myCloud[0].setGuidePoint(nimbus::Winds::EAST);

		// Load shaders
		// Main shader
		shaderCloud.loadShader(GL_VERTEX_SHADER, "../data/shaders/canvasCloud.vert");
#ifdef MOUNT
		// Mountains shader for cumulus
		shaderCloud.loadShader(GL_FRAGMENT_SHADER, "../data/shaders/clouds_CUMULUS_MOUNT.frag");
#endif
#ifdef SEA
		// Sea shader for cumulus
		shaderCloud.loadShader(GL_FRAGMENT_SHADER, "../data/shaders/clouds_CUMULUS_SEA.frag");
#endif
		// Axis shaders
		shaderAxis.loadShader(GL_VERTEX_SHADER, "../data/shaders/axis.vert");
		shaderAxis.loadShader(GL_FRAGMENT_SHADER, "../data/shaders/axis.frag");

		// Create shader programs
		shaderCloud.createShaderProgram();
		shaderAxis.createShaderProgram();

#ifdef MOUNT
		mountain.getUniforms(shaderCloud);
#endif
		canvas.getUniforms(shaderCloud);
		nimbus::Cloud::getUniforms(shaderCloud);
		nimbus::Cumulus::getUniforms(shaderCloud);
		axis.getUniforms(shaderAxis);

		// Start main loop
		for (;;)
		{
			idleGL();
			displayGL();
			std::cerr << " swap " << std::endl;
            glfwWaitEventsTimeout(0);
			glfwSwapBuffers(window);
		}

	}
	catch (nimbus::NimbusException& exception)
	{
		exception.printError();
		system("pause");
	}

	// Free texture
	nimbus::Cloud::freeTexture();
	return 0;
}

#endif


// Entry point for mesh morphing
#ifdef MODEL

void main(int argc, char* argv[])
{
	// Initialize FreeGLUT and Glew

	initGL(argc, argv);

	try {
		// Cameras setup
		cameraFrame.setProjectionRH(30.0f, SCR_W / SCR_H, 0.1f, 2000.0f);
		cameraAxis.setProjectionRH(30.0f, SCR_W / SCR_H, 0.1f, 2000.0f);
		cameraFrame.setViewport(0, 0, SCR_W, SCR_H);
		cameraFrame.setLookAt(glm::vec3(0, 0, -SCR_Z), glm::vec3(0, 0, SCR_Z));
		cameraFrame.translate(glm::vec3(-SCR_W / 2.0, -SCR_H / 2.0, -SCR_Z));
		
		userCameraPos = glm::vec3(0.0, 0.4, 0.0);

		// Create framgent shader canvas
		canvas.create(SCR_W, SCR_H);
		// Create cloud base texture
		nimbus::Cloud::createTexture(TEXTSIZ);

#ifdef MOUNT
		mountain.create(300.0, false); // Create mountain
#endif
		axis.create(); // Create 3D axis
		model1.create(glm::vec3(-1.0, 7.0, 0.0), MESH1, 1.1f); // Create mesh 1
		model2.create(glm::vec3(1.0, 7.0, -3.0), MESH2, 1.1f);  // Create mesh 2
		morphing.setModels(&model1, &model2, EVOLUTE); // Setup modes for morphing

		// Load shaders
		// Main shader
		shaderCloud.loadShader(GL_VERTEX_SHADER, "../Nube/data/shaders/canvasCloud.vert");
#ifdef MOUNT
		// Mountains shader for 3D meshes based clouds
		shaderCloud.loadShader(GL_FRAGMENT_SHADER, "../Nube/data/shaders/clouds_MORPH_MOUNT.frag");
#endif
#ifdef SEA
		// Sea shader for 3D meshes based clouds
		shaderCloud.loadShader(GL_FRAGMENT_SHADER, "../Nube/data/shaders/clouds_MORPH_SEA.frag");
#endif

		// Axis shaders
		shaderAxis.loadShader(GL_VERTEX_SHADER, "../Nube/data/shaders/axis.vert");
		shaderAxis.loadShader(GL_FRAGMENT_SHADER, "../Nube/data/shaders/axis.frag");

		// Create shader programs
		shaderCloud.createShaderProgram();
		shaderAxis.createShaderProgram();

		// Locate uniforms
		nimbus::Cloud::getUniforms(shaderCloud);

#ifdef MOUNT

		mountain.getUniforms(shaderCloud);

#endif

		canvas.getUniforms(shaderCloud);

		nimbus::Model::getUniforms(shaderCloud);

		axis.getUniforms(shaderAxis);

		// Start main loop
		for (;;)
		{
			idleGL();
			displayGL();
			std::cerr << " swap " << std::endl;
            glfwWaitEventsTimeout(0);
			glfwSwapBuffers(window);
		}

	}
	catch (nimbus::NimbusException& exception)
	{
		exception.printError();
		system("pause");
	}

	// Free texture
	nimbus::Cloud::freeTexture();
}

#endif


// Reshape window

void reshapeGL(GLFWwindow*, int w, int h)
{
	if (h == 0)
	{
		h = 1;
	}

	windowWidth = w;
	windowHeight = h;
}

// Application speed synchronization

void syncFPS()
{
	/*
	static int32_t dwLastTime = 0;

	int32_t dwCurrentTime = GetTickCount(); // Get milliseconds from the system start up

	int32_t dwElapsed = dwCurrentTime - dwLastTime; // Calculates elapsed time

	if (dwElapsed < FPS) // The frame loop lasted less than the defined time 				
	{
		Sleep(FPS - dwElapsed);	// Sleeps the application
		dwLastTime = dwCurrentTime + FPS - dwElapsed; // Adds the sleeped time

	}
	else dwLastTime = dwCurrentTime;	// The frame loop exceeded the time
*/
}

// Calculate Frame-per-second

void calculateFPS()
{
	static bool firstTime = true;
	float fps;
	int currentTime;

	//  Increase frame count
	frameCount++;

	if (firstTime)
	{
		previousTime = currentTime;
		firstTime = false;
	}

	//  Calculate time passed
	int timeInterval = currentTime - previousTime;

	if (timeInterval > 1000)
	{
		//  calculate the number of frames per second
		fps = frameCount / (timeInterval / 1000.0f);
		std::cerr << "FPS = " << fps << std::endl;

		static float samples = 0.0f;
		static float sumFPS = 0.0f;

		if (mean)
		{
			sumFPS += fps;
			samples++;
			if (samples > 15)
			{
				std::cerr << "====================MEAN FPS for" << ((parallel) ? " CUDA = " : " CPU = ") << sumFPS / samples << std::endl;
				sumFPS = samples = 0.0f;
				mean = false;
			}
		}
		else
		{
			samples = 0.0f;
			sumFPS = 0.0f;
		}

		//  Set time
		previousTime = currentTime;

		//  Reset frame count
		frameCount = 0;

	}
}

#ifdef CUMULUS
void applyWind() // Apply wind force
{

#ifdef CUDA
	if (parallel)
	{
		windGridCUDA.clearUVW(); // Clear fluid internal data
		for (int i = 0; i < nimbus::Cumulus::getNumClouds(); i++)
			myCloud[i].applyWind(windForce, &windGridCUDA);
	} else
#endif
	{
		windGridCPU.clearUVW(); // Clear fluid internal data
		for (int i = 0; i < nimbus::Cumulus::getNumClouds(); i++)
			myCloud[i].applyWind(windForce, &windGridCPU);
	}
	
}

#endif

#ifdef CUMULUS

// Render function
void displayGL()
{
	if (onPlay)
	{

		// Clear back-buffer
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		//////////////// CLOUDS ///////////////////

		// Change coordinate system

		glm::vec2 mouseScale = mousePos / glm::vec2(SCR_W, SCR_H);

		// Rotate camera
		glm::vec3 userCameraRotatePos = glm::vec3(sin(mouseScale.x*3.0), mouseScale.y, cos(mouseScale.x*3.0));
		
		glDisable(GL_DEPTH_TEST);

		shaderCloud.useProgram();

		if (firstPass) // If first iteration
		{

			// Render cloud base texture
			nimbus::Cloud::renderTexture();

#ifdef MOUNT
			// Render mountain
			mountain.render();
#else
			glActiveTexture(GL_TEXTURE0 + 1);
#endif
			nimbus::Cloud::renderFirstTime(SCR_W, SCR_H);
		}
		
			// Render cloud base class uniforms
			nimbus::Cloud::render(mousePos, timeDelta, cloudDepth, skyTurn, cloudDepth  * userCameraRotatePos, debug);
		
			// Calculate and apply wind
			for (int i = 0; i < nimbus::Cumulus::getNumClouds(); i++)
#ifdef CUDA
				(parallel) ?
				myCloud[i].computeWind(&windGridCUDA) :
#endif
				myCloud[i].computeWind(&windGridCPU);
			
			// Wind setup
			nimbus::Cumulus::setWind(windDirection);
			
			// Render cumulus class
			nimbus::Cumulus::render(shaderCloud);

			if (precomputeTimeOut >= PRECOMPTIMEOUT) // Check for regular precompute light (shading)
			{
				clock_t start = clock();

#ifdef CUDA
				if (parallel)
				{
					if (skyTurn == 0) // If morning sun is near else sun is far (sunset)
						nimbus::Cumulus::precomputeLight(precompCUDA, sunDir, 100.0f, 0.2f);
					else nimbus::Cumulus::precomputeLight(precompCUDA, sunDir, 10000.0f, 1e-6f);
					cudaDeviceSynchronize();
				}
				else
#endif
				{
					if (skyTurn == 0)
						nimbus::Cumulus::precomputeLight(precompCPU, sunDir, 100.0f, 0.2f);
					else nimbus::Cumulus::precomputeLight(precompCPU, sunDir, 10000.0f, 1e-6f);
				}
				clock_t end = clock();
				
				for (int i = 0; i < nimbus::Cumulus::getNumClouds(); i++)
					nimbus::Cloud::renderVoxelTexture(i);
			
				precomputeTimeOut = 0;
				float msElapsed = static_cast<float>(end - start);
				std::cerr << "PRECOMPUTING LIGHT TIME =" << msElapsed / CLOCKS_PER_SEC << std::endl;
			}


			if (totalTime > TOTALTIME)
			{
				timeDelta += nimbus::Cumulus::getTimeDir();
				totalTime = 0;
			}			

			totalTime++;
			precomputeTimeOut++;

			// User camera setup
			cameraSky.setLookAt(cloudDepth * userCameraRotatePos, userCameraPos);
			cameraSky.translate(userCameraPos);
				
			canvas.render(cameraFrame, cameraSky);

			/////////////// AXIS ////////////////////

			shaderAxis.useProgram();

			// Render axis
			cameraAxis.setViewport(SCR_W / 3, SCR_H / 3, SCR_W, SCR_H);
			cameraAxis.setLookAt(10.0f * userCameraRotatePos, userCameraPos);
			cameraAxis.setPosition(userCameraPos);

			axis.render(cameraAxis);

			// Restore landscape viewport
			cameraFrame.setViewport(0, 0, SCR_W, SCR_H);

			glEnable(GL_DEPTH_TEST);
			//  Calculate FPS
			calculateFPS();
			firstPass = false; // First pass ended
		
	}
}

#endif

#ifdef MODEL

void displayGL()
{
	if (onPlay)
	{

		// Clear back-buffer
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Change coordinate system

		glm::vec2 mouseScale = mousePos / glm::vec2(SCR_W, SCR_H);

		// Rotate camera

		glm::vec3 userCameraRotatePos =  glm::vec3(sin(mouseScale.x*3.0), mouseScale.y, cos(mouseScale.x*3.0));
		
		//////////////// MORPHING ///////////////////

		glDisable(GL_DEPTH_TEST);

		shaderCloud.useProgram();
		if (firstPass) // If first iteration
		{

			nimbus::Cloud::renderFirstTime(SCR_W, SCR_H);

			clock_t start = clock();
#ifdef CUDA	
			// Precompute light for meshes
			nimbus::Model::precomputeLight(precompCUDA, sunDir, 100.0f, 1e-6f, model1.getNumEllipsoids(), model2.getNumEllipsoids());
#else
			nimbus::Model::precomputeLight(precompCPU, sunDir, 100.0f, 1e-6f, model1.getNumEllipsoids(), model2.getNumEllipsoids());
#endif
			clock_t end = clock();
			std::cerr << "PRECOMPUTE TIME LIGHT = " << end - start << std::endl;

			// Prepare for morphing
			(EVOLUTE) ? morphing.prepareMorphEvolute() : morphing.prepareMorphInvolute();
			alpha = alphaDir = 0.01f;
			// First morphing render
			nimbus::Model::renderFirstTime(model2.getNumEllipsoids(), EVOLUTE);
			
			// Render cloud base texture
			nimbus::Cloud::renderTexture();
#ifdef MOUNT
			// Render mountain
			mountain.render();
#endif
		
			// Render clouds precomputed light textures
			for (int i = 0; i < nimbus::Cloud::getNumClouds(); i++)
				nimbus::Cloud::renderVoxelTexture(i);
		}

		// Render cloud base class uniforms
		nimbus::Cloud::render(mousePos, timeDelta, cloudDepth, skyTurn, cloudDepth * userCameraRotatePos, debug);
			
		static bool totalTimePass = false;

		if (totalTime > TOTALTIME) // Check time for morphing animation
		{
			totalTimePass = true;
			timeDelta += timeDir;
			if (alpha < 1.0 && alpha > 0.0)
			{
				alpha += alphaDir; // Animate morphing
				(EVOLUTE) ? morphing.morphEvolute(alpha) : morphing.morphInvolute(alpha);
				morphing.morph(0.1f); // Animation speed
			}

			totalTime = 0;
		}

		totalTime++;

		// Mesh renderer
		nimbus::Model::render(shaderCloud, (totalTimePass) ? morphing.getCloudPosRDst() : nimbus::Model::getCloudPosR(), morphing.getCloudPosDst(), alpha);


		// User camera setup
		cameraSky.setLookAt(cloudDepth * userCameraRotatePos, userCameraPos);
		cameraSky.translate(userCameraPos);

		canvas.render(cameraFrame, cameraSky);

		/////////////// AXIS ////////////////////

		shaderAxis.useProgram();

		// Render axis
		cameraAxis.setViewport(SCR_W / 3, SCR_H / 3, SCR_W, SCR_H);
		cameraAxis.setLookAt(10.0f * userCameraRotatePos, userCameraPos);
		cameraAxis.setPosition(userCameraPos);

		axis.render(cameraAxis);

		// Restore landscape viewport
		cameraFrame.setViewport(0, 0, SCR_W, SCR_H);

		glEnable(GL_DEPTH_TEST);
		//  Calculate FPS
		calculateFPS();
		firstPass = false; // First pass ended
	}
}

#endif


// Idle function

void idleGL()
{
	if (!onPlay) return;

	syncFPS();

#ifdef CUMULUS

	simcont++;
	if (simcont > FLUIDLIMIT) // Simulate fluid
	{
		applyWind();
		clock_t start = clock();		
#ifdef CUDA
		if (parallel)
		{
			windGridCUDA.sendData();
			windGridCUDA.sim();
			windGridCUDA.receiveData();
			cudaDeviceSynchronize(); // For clock_t usage
		} else 
#endif
			windGridCPU.sim();

		clock_t end = clock();
		float msElapsed = static_cast<float>(end - start);

		std::cerr << "FLUID SIMULATION TIME = " << msElapsed / CLOCKS_PER_SEC << std::endl;
		simcont = 0;
	}
#endif
}

// Key press function



void keyboardGL(GLFWwindow*, int key, int scancode, int action, int mods)
{
	char c = key;
	if (action == GLFW_RELEASE)
		keyboardUpGL(c);
	switch (c)
	{
	case 'w':
	case 'W':
		break;
	case 'a':
	case 'A':
		break;
	case 's':
		break;
	case 'S':
		break;
	case 't':
		parallel = !parallel;
		(parallel) ? std::cerr << "......CUDA MODE......." << std::endl : std::cerr << "......CPU MODE......." << std::endl;
		break;
	case 'm':
		mean = !mean;
		(mean) ? std::cerr << "MEAN ENABLED for " : std::cerr << "MEAN DISABLED for ";
		std::cerr << ((parallel) ? "CUDA" : "CPU") << std::endl;

		break;
	case 'd':
	case 'D':
		break;
	case 'q':
	case 'Q':
		break;
	case 'e':
	case 'E':
		break;
#ifdef CUMULUS
	case '8': // The wind blows from the south
		windDirection = nimbus::Winds::SOUTH;
		for (int i = 0; i < nimbus::Cloud::getNumClouds(); i++)
			myCloud[i].setGuidePoint(windDirection);
		break;
	case '2': // The wind blows from the north
		windDirection = nimbus::Winds::NORTH;
		for (int i = 0; i < nimbus::Cloud::getNumClouds(); i++)
			myCloud[i].setGuidePoint(windDirection);

		break;
	case '4': // The wind blows from the east
		windDirection = nimbus::Winds::EAST;
		for (int i = 0; i < nimbus::Cloud::getNumClouds(); i++)
			myCloud[i].setGuidePoint(windDirection);

		break;
	case '6': // the wind blows from the west
		windDirection = nimbus::Winds::WEST;
		for (int i = 0; i < nimbus::Cloud::getNumClouds(); i++)
			myCloud[i].setGuidePoint(windDirection);
		break;
	case 'z':
		break;
	case 'n':
		break;
	case 'r': // Recreate new cloud
		nimbus::Cloud::resetNumClouds();
		precomputeTimeOut = PRECOMPTIMEOUT;
		myCloud[0].create(35, 2.1f, glm::vec3(10.0, 11.0, 5.0), 0.0f, 3.0f, 0.0f, 1.2f, 0.0f, 5.0f, true, false);
		myCloud[0].setGuidePoint(nimbus::Winds::EAST);
		myCloud[1].create(35, 2.1f, glm::vec3(-10.0, 12, 10.0), 0.0f, 3.0f, 0.0f, 1.2f, 0.0f, 5.0f, true, false);
		myCloud[1].setGuidePoint(nimbus::Winds::EAST);
		myCloud[2].create(35, 2.3f, glm::vec3(0.0, 11.0, 20.0), 0.0f, 3.0f, 0.0f, 1.2f, 0.0f, 5.0f, true, false);
		myCloud[2].setGuidePoint(nimbus::Winds::EAST);
		firstPass = true;
		break;
#endif
	case 27:
		//glutLeaveMainLoop();
		break;
	case 32:
		break;
	case 13: // Change daytime
		skyTurn = (skyTurn + 1) % 3;
	}
}

// Change linear interpolation direction back/front<->front/back

void changeADir()
{
	if (alpha <= 0.0)
	{
		alpha = 0.01f;
		alphaDir = 0.01f;
	}
	else if (alpha >= 1.0)
	{
		alpha = 0.99f;
		alphaDir = -0.01f;
	}
}

// Keyboard press

void keyboardUpGL(unsigned char c)
{
	switch (c)
	{
	case 'w':
		break;
	case 'W':
		break;
	case 'a':
	case 'A':
		break;
	case 's':
		break;
	case 'S':
		break;
	case 'd':
		debug = !debug; // Activating debug
		std::cerr << "DEBUG = " << debug << std::endl;
		break;
	case 'D':
		break;
	case 'q':
		changeADir(); // Setup evolution/involution direction
		break;
	case 'Q':
		break;
	case 'e':
		alphaDir = -alphaDir; // Change evolution/involution direction
		break;
	case 'E':
		break;
	case 'p': // Pause 
		onPlay = !onPlay;
		break;
	case GLFW_KEY_UP:
		userCameraPos.y += 0.1f;
		break;
	case GLFW_KEY_DOWN:
		userCameraPos.y -= 0.1f;
		break;
	case GLFW_KEY_RIGHT:
		userCameraPos.x += 0.1f;
		break;
	case GLFW_KEY_LEFT:
		userCameraPos.x -= 0.1f;

	default:
		break;
	}
}



// GLUT mouse wheel handler


void mouseWheel(GLFWwindow* window, double x, double y)
{

	if (x > 0)
	{
		// Camera zoom in
		cloudDepth -= 0.4f;
	}
	else
	{
		// Camera zoom out
		cloudDepth += 0.4f;
	}
	std::cerr << "Depth-1 = " << cloudDepth << std::endl;
}

// Mouse motion function

void motionGL(GLFWwindow*, double x, double y)
{
	mousePos = glm::vec2(x * 4, y * 4);
}
