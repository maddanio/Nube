////////////////////////////////////////////////////////////////////
// ONTOGENETIC MODEL FOR REAL-TIME VOLUMETRIC CLOUDS SIMULATION THESIS			
// Software Engineering and Computer Systems Deparment	
// National University for Distance Education (UNED)			    		
// Carlos Jiménez de Parga, PhD student.
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// Last revision 19/04/2019
//////////////////////////////////////////////////////////////////

#version 450 core
#pragma optimize(on)

//#define FULLMOON

#define SC (250.0)

// Declare uniforms

uniform sampler3D iNoise;
uniform sampler3D iVoxel[10];
uniform sampler2D iChannel0;
uniform float iDepth;
uniform int iNumSph;
uniform int iNumClouds;
uniform bool isFlat[10];
uniform int iTurn;
uniform vec2 iResolution;
uniform mat4 iView;
uniform vec3 iPos;
uniform vec3 iVmin[10];
uniform vec3 iVmax[10];
uniform int iLowLimits[10];
uniform int iUpLimits[10];
uniform int iDebug;
uniform float iTime;
uniform vec3 iWindDirection;
uniform float iMean;
uniform vec3 iSunDir;
uniform bool iSnow;

// Global variables

vec3 cloudColor;         // Color of the cloud
vec3 sunColor;           // Light color
vec4 candidates[100];    // List of candidates
float kPhase;           // Phase function constant
float T;                // Light threshold
int n = 0;              // Number of candidates
vec3 rdNorm;            // Normalized ray-direction
vec3 skyColor;          // Color of the sky
float sunSize = 0.5;   // Sun radius
float tmin, tmax, tymin, tymax, tzmin, tzmax; // For Smits' algorithm

in vec2 fragCoord; // Fragment shader 2D coordinates
out vec4 FBColor;  // Returned frame buffer color

// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// Created by S. Guillitte 2015
// Galaxy morphology based on http://iopscience.iop.org/0004-637X/783/2/138/pdf/0004-637X_783_2_138.pdf

const mat2 starM = mat2(0.8, 0.6, -0.6, 0.8);

float noise(in vec3 p)
{

	p *= 2.;
	float res = 0.;
	float f = 1.;
	for (int i = 0; i < 3; i++)
	{
		p.xy = starM * p.xy;
		p = p.zxy*f + .6;
		f *= 1.15;
		res += sin(p.y + 1.3*sin(1.2*p.x) + 1.7*sin(1.7*p.z));
	}
	return res / 3.;

}

float fbmgal(vec3 p) {

	p = p * 10.0;

	float f = 1.;
	float r = 0.0;
	for (int i = 1; i < 5; i++) {
		r += noise(p*(20. + 3.*f)) / f;
		p.xz *= starM;
		f += 1.;

	}
	return pow(abs(r), 4.);
}


/*

	Non physical based atmospheric scattering made by robobo1221
	Site: http://www.robobo1221.net/shaders
	Shadertoy: http://www.shadertoy.com/user/robobo1221

*/

const float pi = 3.14159265359;
const float invPi = 1.0 / pi;

const float zenithOffset = 0.16;
const float multiScatterPhase = 0.4;
const float density = 0.7;

const float anisotropicIntensity = 0.1; //Higher numbers result in more anisotropic scattering

const vec3 sunSetColor = vec3(0.39, 0.57, 1.0) * (1.0 + anisotropicIntensity); //Make sure one of the conponents is never 0.0

#define smooth(x) x*x*(3.0-2.0*x)
#define zenithDensity(x) density / pow(max(x - zenithOffset, 0.35e-2), 0.75)


vec3 getSkyAbsorption(vec3 x, float y) {

	vec3 absorption = x * -y;
	absorption = exp2(absorption) * 2.0;

	return absorption;
}



float getSunPoint(vec3 rd, vec3 sd) {
	return smoothstep(0.03, 0.026, 1.0 + dot(rd, sd) + 0.028)*50.0;
}

float getRayleigMultiplier(vec3 rd, vec3 sd) {
	return 1.0 + pow(1.0 - clamp(1.0 + dot(rd, sd) + 0.028, 0.0, 1.0), 2.0) * pi * 0.5;
}

float getMie(vec3 rd, vec3 sd) {
	float disk = clamp(1.0 - pow(1.0 + dot(rd, sd) + 0.028, 0.1), 0.0, 1.0);

	return disk * disk*(3.0 - 2.0 * disk) * 2.0 * pi;
}

vec3 getAtmosphericScattering(vec2 p, vec3 lp, vec3 rd) {
	vec3 correctedLp = lp;

	float zenith = zenithDensity(p.y);
	float sunPointDistMult = clamp(length(max(correctedLp.y + multiScatterPhase - zenithOffset, 0.0)), 0.0, 1.0);

	float rayleighMult = getRayleigMultiplier(rd, correctedLp);

	vec3 absorption = getSkyAbsorption(sunSetColor, zenith);
	vec3 sunAbsorption = getSkyAbsorption(sunSetColor, zenithDensity(correctedLp.y + multiScatterPhase));
	vec3 sky = sunSetColor * zenith * rayleighMult;
	vec3 sun = getSunPoint(rd, correctedLp) * absorption;
	vec3 mie = getMie(rd, correctedLp) * sunAbsorption;

	float grad = 1.0 - exp(-p.y * 3.0);

	vec3 totalSky = mix(sky * absorption, vec3(0.0), grad);

	totalSky += sun + mie;
	totalSky *= sunAbsorption * 0.5 + 0.5 * length(sunAbsorption);

	return totalSky;
}

vec3 jodieReinhardTonemap(vec3 c) {
	float l = dot(c, vec3(0.2126, 0.7152, 0.0722));
	vec3 tc = c / (c + 1.0);

	return mix(c / (l + 1.0), tc, tc);
}



vec3 getSunset(in vec3 rd, in vec2 px)
{
  vec3 color = getAtmosphericScattering(px, iSunDir, rd) * pi;
	color = jodieReinhardTonemap(color);
	color = pow(color, vec3(2.2)); //Back to linear

	return color;


}


// sky
vec3 getSkyColor(vec3 e) {

	e.y = max(e.y, 0.0) + 0.2;
	return vec3(pow(1.0 - e.y, 2.0), 1.0 - e.y, 0.6 + (1.0 - e.y)*0.4);

}


float diffuseSphere(vec3 rd, vec3 c, float r, vec3 l)
{
	float px = rd.x - c.x;
	float py = rd.y - c.y;
	float pz = rd.z - c.z;
	float sq = r * r - px * px - py * py - pz * pz;
	if (sq < 0)
	{
		return 0;
	}

	float z = sqrt(sq);
	vec3 normal = normalize(vec3(px, -py, z));
	float diffuse = max(0.0, dot(normal, l));
	return diffuse;
}

// Spheres array

layout(std140) uniform iCloudPosBlock
{
	vec4 iCloudPos[150];
};


// Calculate fBm noise

float fbm(in vec3 q)
{

	q = q + iWindDirection * iTime;

	float f = 0.0;
	float a = 0.5;

	for (int i = 0; i < 5; i++) {
		f += a * textureLod(iNoise, q / 64.0, -100.0).r;
		q = q * 2.0 + 0.03;
		a *= 0.5;
	}

	return f;


}

// Collision detection routine

void addCandidates(in vec3 rayOrg, in vec3 rayDir, int boundIdx)
{
	// Limits of the cloud spheres
	
	int sphStart = iLowLimits[boundIdx];
	int sphEnd = iUpLimits[boundIdx];

	// Iterate over spheres
	
	for (int j = sphStart; j < sphEnd; j++)
	{

		vec3 cloudPos = iCloudPos[j].xyz; // Sphere (x,y,z) position
		float radius = iCloudPos[j].w;   // Sphere radius 

		vec3 temp = rayOrg - cloudPos;
		float a = dot(rayDir, rayDir);
		float b = 2.0*dot(rayDir, temp);
		float c = dot(temp, temp) - radius * radius;

		float disc = b * b - 4.0*a*c;

		if (disc > 0.0) // There is a collision
		{

			disc = sqrt(disc);
			float t = ((-b - disc) / (2.0*a));
			float limit = ((-b + disc) / (2.0*a));
			candidates[n] = vec4(t, limit, j, boundIdx); // Add to candidates array
			n++;

		}
	}
}

// Insertion-sort algorithm

void order()
{
	int h;
	vec4 aux;

	for (int i = 1; i < n; i++)
	{
		aux = candidates[i];
		h = i - 1;
		while ((h >= 0) && (aux.x < candidates[h].x))
		{
			candidates[h + 1] = candidates[h];
			h--;
		}
		candidates[h + 1] = aux;
	}
}


// Simplified Henyey-Greenstein phase function

float phase(in float g)
{

	return 0.0795 * ((1.0 - g * g) / pow(1.0 + g * g - 2.0*g*dot(rdNorm, iSunDir), 1.5));

}

// Ray-trace pseudo-sphere

vec4 trace(in vec3 rayOrg, in vec3 rayDir, in vec3 color)
{

	float localDen, den;

	float tIn, tOut;
  
	for (int i = 0; i < n; i++)
	{

		int sphIdx = int(candidates[i].z);    // Sphere index
		int boundIdx = int(candidates[i].w);  // Cloud index
		vec3 cloudPos = iCloudPos[sphIdx].xyz; // Sphere 3D position
		float radius = iCloudPos[sphIdx].w;   // Sphere radius
		vec3 iVoxMin = iVmin[boundIdx];        // Bounding box 3D min values
		vec3 iVoxMax = iVmax[boundIdx];        // Bounding box 3D max values
		vec3 difVox = iVoxMax - iVoxMin;       // Voxel distance 

		tIn = candidates[i].x;                 // Lambda-in
		tOut = candidates[i].y;                // Lambda-out
		  

		bool bFlat = isFlat[boundIdx];        // Is the cloud flat?

		// LOD (Level-of-Detail)
		float lambda = clamp(10.0*exp(-distance(rayOrg, cloudPos)*0.23), 0.1, 10.0);

		while (tIn <= tOut) // Iterate sphere
		{
			vec3  pos = rayOrg + tIn * rayDir; // Point in the straight line

			vec3 index = (pos - iVoxMin) / difVox; // Voxel index

			// Pre-compute light value
			float precLight = texture(iVoxel[boundIdx], index, -100.0).r;

			if (!bFlat || (pos.y > iVoxMin.y + precLight * 1.4)) // Condition for flat cloud
			{
				den = fbm(pos); // Density

				// Create pseudo-sphere
				localDen = exp(-distance(pos, cloudPos) / (radius*(0.7 + 0.4*den)));

				// Render pseudo-sphere with lighting
				if (den < localDen)
				{

					float deltaT; // Calculate light

					deltaT = exp(-0.1 * den);

					vec3 absorpLight = sunColor * precLight; // Absorption                                                           

					vec3 scatterLight = sunColor * phase(kPhase)*precLight; // Scattering

					vec3 totalLight = absorpLight * 0.9 + scatterLight; // Total light

					color += (1.0 - deltaT) * totalLight * T; // Coud color

					T *= deltaT; // Increment threshold

					if (T < 1e-6) // Exit condition
						return vec4(color, 1 - T);

				}
			}

			tIn+=0.1; // Step ray-marching

		}

	}
 return vec4(color, 1 - T); // Return color with transparency
}

// Render cloud

vec4 render(in vec3 rayOrg, in vec3 rayDir, in vec2 px, in vec4 resM)
{

	rdNorm = normalize(rayDir);

	// background sky     
	vec4 skycol;
	vec4 res = vec4(0);

	float sun = 1.0 + dot(iSunDir, rdNorm);

	if (iTurn == 0) // Morning scene
	{
		float fexp = exp(-sun * 1000.0 / sunSize);
		float fexp2 = exp(-sun * 1000.0 / 0.1);

		skycol = vec4(0.2*sunColor*fexp + 0.6*sunColor*fexp2, 1.0);
		// sun glare    
		skycol += vec4(getSkyColor(rdNorm.xyz) + 0.2*sunColor*exp(-sun), 0);

	}

	else if (iTurn == 1) // Sunset scene
	{

		vec3 stars;

		if ((rdNorm.y > 0.1) && (dot(iSunDir, rdNorm) > -0.995))
		{
			float g = 0.2*fbmgal(rdNorm); // Generate stars
			stars = 0.3*vec3(g*g*g, g*g*1.3, 8.5*g);
			if (length(stars) < 2.4)
				stars = vec3(0);
		}

		// sun glare    
		skycol = vec4(getSunset(rdNorm, px) + 0.2*sunColor*exp(-sun) + stars, 1.0);

		resM.xyz = resM.xyz / 5.0; // Mountains darkness

	}
	else // Night scene
	{
		vec3  stars = vec3(0);

		if ((dot(iSunDir, rdNorm) > -0.9988))
		{
			float g = 0.2*fbmgal(rdNorm); // Generate stars
			stars = 0.3*vec3(g*g*g, g*g*1.3, 1.5*g);
		}


#ifdef FULLMOON
    // Black moon by DeMaCia
    // License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

		// Render fullmoon
		float fexp = exp(-sun * 1000.0 / sunSize);
		float fexp2 = exp(-sun * 1000.0 / 0.2);


		skycol = vec4(0.2*sunColor*fexp + sunColor * fexp2, 1.0);
		// moon glare    
		skycol += vec4(0.2*sunColor*exp(-sun) + stars, 0);

#else

		// Render quarter moon
		float time = 812.81;

		//Diffuse
		float r = 0.05;
		vec3 vp = vec3(sin(time*0.2), cos(time*0.2), sin(time*0.2));
		vec3 vl = normalize(vp);
		float diffuse = diffuseSphere(rdNorm, -iSunDir, r, vl);

		skycol = vec4(vec3(diffuse) + stars, 1);
#endif

		resM.xyz = resM.xyz / 8.0; // Mountains darkness
	}

	// Smits' algorithm   

	for (int i = 0; i < iNumClouds; i++)
	{

		float vminX = iVmin[i].x;
		float vmaxX = iVmax[i].x;
		float vminY = iVmin[i].y;
		float vmaxY = iVmax[i].y;
		float vminZ = iVmin[i].z;
		float vmaxZ = iVmax[i].z;

		bool flag = true;

		if (rayDir.x >= 0)
		{
			tmin = (vminX - rayOrg.x) / rayDir.x;
			tmax = (vmaxX - rayOrg.x) / rayDir.x;
		}
		else
		{
			tmin = (vmaxX - rayOrg.x) / rayDir.x;
			tmax = (vminX - rayOrg.x) / rayDir.x;
		}
		if (rayDir.y >= 0)
		{
			tymin = (vminY - rayOrg.y) / rayDir.y;
			tymax = (vmaxY - rayOrg.y) / rayDir.y;
		}
		else
		{
			tymin = (vmaxY - rayOrg.y) / rayDir.y;
			tymax = (vminY - rayOrg.y) / rayDir.y;
		}

		if ((tmin > tymax) || (tymin > tmax))
			flag = false;

		if (tymin > tmin)
			tmin = tymin;

		if (tymax < tmax)
			tmax = tymax;

		if (rayDir.z >= 0)
		{
			tzmin = (vminZ - rayOrg.z) / rayDir.z;
			tzmax = (vmaxZ - rayOrg.z) / rayDir.z;
		}
		else
		{
			tzmin = (vmaxZ - rayOrg.z) / rayDir.z;
			tzmax = (vminZ - rayOrg.z) / rayDir.z;
		}

		if ((tmin > tzmax) || (tzmin > tmax))
			flag = false;

		// In bounding-box    
		if (flag)
			addCandidates(rayOrg, rayDir, i);
			
		

	}

	T = 1.0; // Initialize T

	if (n > 0 && resM.w < 0.0) // No mountain collision
	{
		order();
		vec3 color = vec3(0.0);
		res = trace(rayOrg, rayDir, color);

		return res + skycol * T;
	}
	else if (n > 0 && resM.w > 0.0) // Mountain collision
	{
		order();
		vec3 color = vec3(0.0);
		res = trace(rayOrg, rayDir, color);

		return resM * T + res;
	}
	else if (n == 0 && resM.w < 0.0)
		return skycol;
	else if (n == 0 && resM.w > 0.0)
		return resM;

}



void main()
{

	// Change coordinate system
	vec2 p = (-iResolution.xy + 2.0*fragCoord.xy) / iResolution.y;

	vec2 position = fragCoord.xy / iResolution.x;

	switch (iTurn)
	{
	case 0:
		// MORNING	
		sunColor = vec3(1.0, 1.0, 1.0);
		kPhase = -0.4;
		break;
	case 1:
		// SUNSET	
    sunColor = vec3(1.0, 0.8, 0.2);
	  kPhase = -0.5;
    break;
	case 2:
		// NIGHT
		sunColor = vec3(0.9, 0.9, 0.9);
		kPhase = -0.5;
	}

	// Camera ray    
	vec4 rayDir = iView * normalize(vec4(p, 1.5, 1));

	vec4 res = renderMountains(iPos*iDepth * 300 + terrainL(iPos.xz) + 19.0*SC, rayDir.xyz);
   
  // Frame-buffer return color
	FBColor = render(iPos, rayDir.xyz, position, res);
 
  if (iCloudPos[2].x == 2)
    FBColor = vec4(0,1,0,1);
 
}

