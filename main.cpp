#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
 
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
// Needed on MsWindows
#define NOMINMAX
#include <windows.h>
#endif // Win32 platform
 
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
// Download glut from: http://www.opengl.org/resources/libraries/glut/
#include <GLUT/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <float.h>

#include "float2.h"
#include "float3.h"
#include "float4.h"
#include "float4x4.h"
#include <vector>
#include <algorithm>
// #include "perlin.h"

// For procedural texturing
// Perlin perlin;
// for quadrics, so that we do not need a float4.cpp
const float4x4 float4x4::identity(
	1, 0, 0, 0,
	0, 1, 0, 0,
	0, 0, 1, 0,
	0, 0, 0, 1);

// Abstract base class for light sources
class LightSource
{
public:
	virtual float3 getPowerDensityAt(float3 x)=0;
	virtual float3 getLightDirAt(float3 x)=0;
	virtual float  getDistanceFrom(float3 x)=0;

};

// TO BE CREATED AT PRACTICAL
class DirectionalLightSource : public LightSource
{
	float3 powerDensity;
	float3 direction;
public:
	DirectionalLightSource():powerDensity(float3(.5,.5,.5)), direction(float3(1,1,1)) {}
	DirectionalLightSource(float3 pd, float3 dir):powerDensity(pd), direction(dir) {}

	float3 getPowerDensityAt(float3 x) {
		return powerDensity;
	}

	float3 getLightDirAt(float3 x) {
		return direction;
	}

	float getDistanceFrom(float3 x) {
		return FLT_MAX;
	}
};

// TO BE CREATED AT PRACTICAL
class PointLightSource : public LightSource
{
	float3 powerDensity;
	float3 origin;

	float square(float x) {
		return x * x;
	}

	float distSq(float3 x) {
		return square(x.x - origin.x) + square(x.y - origin.y) + square(x.z - origin.z);
	}
public:
	// PointLightSource() {}
	PointLightSource(float3 pd, float3 origin):powerDensity(pd), origin(origin) {}

	float3 getPowerDensityAt(float3 x) {
		return powerDensity/(4 * 3.14 * distSq(x));
	}

	float3 getLightDirAt(float3 x) {
		return (x - origin).normalize();
	}

	float getDistanceFrom(float3 x) {
		return sqrt(distSq(x));
	}
};

// Skeletal Material class. Feel free to add methods e.g. for illumination computation (shading).
class Material
{
public:
	bool reflective;
	bool refractive;
	bool textured;
	float3 minReflectance;		// Fresnel coefficient
	float refractiveIndex;			// index of refraction
	float3 kd;			// diffuse reflection coefficient - color
	float3 ks;			// specular reflection coefficient
	float shininess;	// specular exponent
	Material()
	{
		reflective = false;
		refractive = false;
		textured = false;
		minReflectance = float3(0.93, 0.85, 0.4);
		refractiveIndex = 1;
		kd = float3(0.5, 0.5, 0.5) + kd * 0.5;
		ks = float3(1, 1, 1);
		shininess = 15;
	}

	float3 reflect(float3 inDir, float3 normal) {
		return inDir - normal * normal.dot(inDir) * 2;
	}
	float3 refract(float3 inDir, float3 normal) {
		float ri = refractiveIndex;
		float cosa = -normal.dot(inDir);
		if(cosa < 0) { cosa = -cosa; normal = -normal; ri = 1 / ri; }
		float disc = 1 - (1 - cosa * cosa) / ri / ri;
		if(disc < 0) return reflect(inDir, normal);
		return inDir * (1.0 / ri) + normal * (cosa / ri - sqrt(disc));
	}
	float3 getReflectance(float3 inDir, float3 normal) {
		float cosa = fabs(normal.dot(inDir));
		return minReflectance + (float3(1, 1, 1) - minReflectance) * pow(1 - cosa, 5);
	}

	virtual float3 shade(
		float3 position,
		float3 normal,
		float3 viewDir,
		float3 lightDir,
		float3 lightPowerDensity)
	{
		// return normal * lightPowerDensity;
	    float cosTheta = normal.dot(lightDir);
	    if(cosTheta < 0) return float3(0,0,0);
	    float3 halfway = (viewDir + lightDir).normalize(); // this is H/halfway vector
	    float cosDelta = normal.dot(halfway);
	    float3 diffuse = kd * lightPowerDensity * cosTheta;
	    if(cosDelta < 0) return diffuse; // return float3(0,0,0);
	    return diffuse + lightPowerDensity * ks 
	                    * pow(cosDelta, shininess);
	    // float cosTheta = normal.dot(lightDir);
	    // return lightDir * lightDir;          
	}
};

// Marble
// float cosTheta = normal.dot(lightDir);
// if(cosTheta < 0) return float3(0,0,0);
// float3 halfway = (viewDir + lightDir).normalize(); // this is H/halfway vector
// float cosDelta = normal.dot(halfway);
// float nc = perlin.marble(position);

// Skeletal Camera class. Feel free to add custom initialization, set aspect ratio to fit viewport dimensions, or animation.
class Camera
{
	float3 eye;

	float3 lookAt;
	float3 right;
	float3 up;

public:
	float3 getEye()
	{
		return eye;
	}
	Camera()
	{
		eye = float3(0, 0, 3);
		lookAt = float3(0, 0, 2);
		right = float3(1, 0, 0);
		up = float3(0, 1, 0);
	}

	float3 rayDirFromNdc(const float2 ndc) {
		return (lookAt - eye
			+ right * ndc.x
			+ up    * ndc.y
			).normalize();
	}
};

// Ray structure.
class Ray
{
public:
    float3 origin;
    float3 dir;
    Ray(float3 o, float3 d)
    {
        origin = o;
        dir = d;
    }
};

// Hit record structure. Contains all data that describes a ray-object intersection point.
class Hit
{
public:
	Hit()
	{
		t = -1;
	}
	float t;
	float3 position;
	float3 normal;
	Material* material;
};

// Object abstract base class.
class Intersectable
{
protected:
	Material* material;
public:
	Intersectable(Material* material):material(material) {}
    virtual Hit intersect(const Ray& ray)=0;
};

// Object realization.
class Sphere : public Intersectable
{
	float3 center;
	float radius;
public:
    Sphere(const float3& center, float radius, Material* material):
		Intersectable(material),
		center(center),
		radius(radius)
    {
    }
    Hit intersect(const Ray& ray)
    {
        float3 diff = ray.origin - center;
        double a = ray.dir.dot(ray.dir);
        double b = diff.dot(ray.dir) * 2.0;
        double c = diff.dot(diff) - radius * radius;
 
        double discr = b * b - 4.0 * a * c;
        if ( discr < 0 ) 
            return Hit();
        double sqrt_discr = sqrt( discr );
        double t1 = (-b + sqrt_discr)/2.0/a;
        double t2 = (-b - sqrt_discr)/2.0/a;
 
		float t = (t1<t2)?t1:t2;
		if(t < 0)
			t = (t1<t2)?t2:t1;
		if (t < 0)
            return Hit();

		Hit h;
		h.t = t;
		h.material = material;
		h.position = ray.origin + ray.dir * t;
		h.normal = h.position - center;
		h.normal.normalize();

		return h;
    }
}; 

// TO BE CREATED AT PRACTICAL
class Plane : public Intersectable 
{
	float3 normal;
	float3 x0;
public:
	Plane(const float3& normal, const float3& x0, Material* material):
		Intersectable(material),
		normal(normal),
		x0(x0){}

    Hit intersect(const Ray& ray)
    {
        
 		float denom = ray.dir.dot(normal);

        if ( denom == 0 ) 
            return Hit();
  //       double sqrt_discr = sqrt( discr );
  //       double t1 = (-b + sqrt_discr)/2.0/a;
  //       double t2 = (-b - sqrt_discr)/2.0/a;
 
		// float t = (t1<t2)?t1:t2;
		// if(t < 0)
		// 	t = (t1<t2)?t2:t1;
		// if (t < 0)
  //           return Hit();

        float3 diff = x0 - ray.origin;
        float d = diff.dot(normal)/denom;
		Hit h;
		h.t = d;
		h.material = material;
		h.position = ray.origin + ray.dir * d;
		h.normal = normal;
		h.normal.normalize();

		return h;
    }


};

// TO BE CREATED AT PRACTICAL
class Quadric : public Intersectable
{
	float4x4 A;
public:
    Quadric(float3 center, float4 size, Material* material):
      Intersectable(material)
    { // ellipsoid hardwired here
    // you should add methods
    // to set up different quadrics
    	A = float4x4::identity;
    	A._00 = size.x;
    	A._11 = size.y;
    	A._22 = size.z;
    	A._33 = size.w;

    	A._03 = center.x * size.x;
    	A._13 = center.y * size.y;
    	A._23 = center.z * size.z;
    }

    Hit intersect(const Ray& ray)
    {
	        // ray in homo coords
	    float4 e = float4(ray.origin.x,
	      ray.origin.y, ray.origin.z, 1);
	    float4 d = float4(ray.dir.x,
	      ray.dir.y, ray.dir.z, 0);
	    // quadratic coeffs.
	    double a = d.dot( A * d );
	    double b = e.dot( A * d ) 
	            + d.dot( A * e );
	    double c = e.dot( A * e );
	    // from here on identical to Sphere

        double discr = b * b - 4.0 * a * c;
        if ( discr < 0 ) 
            return Hit();
        double sqrt_discr = sqrt( discr );
        double t1 = (-b + sqrt_discr)/2.0/a;
        double t2 = (-b - sqrt_discr)/2.0/a;
 
		float t = (t1<t2)?t1:t2;
		if(t < 0)
			t = (t1<t2)?t2:t1;
		if (t < 0)
            return Hit();

		Hit h;
		h.t = t;
		h.material = material;
		h.position = ray.origin + ray.dir * t;
		float4 hPos = float4(h.position.x,
		  h.position.y, h.position.z, 1);
		// homo normal per quadric normal formula
		float4 hNormal = A * hPos +  hPos * A;
		// Cartesian normal
		h.normal = float3(hNormal.x, hNormal.y, hNormal.z).normalize();
		// h.normal.normalize();

		return h;
    }

};


class Scene
{
	Camera camera;
	std::vector<LightSource*> lightSources;
	std::vector<Intersectable*> objects;
	std::vector<Material*> materials;
	int maxDepth = 5;
public:
	Scene()
	{
		// ADD LIGHT SOURCES HERE
		lightSources.push_back(new DirectionalLightSource(float3(.8,.8,.8), float3(.7,0,0)));
		lightSources.push_back(new DirectionalLightSource(float3(.8,.8,.8), float3(0,.7,0)));

		// lightSources.push_back(new DirectionalLightSource(float3(.15,.15,.15), float3(-1,-1,1)));
		// lightSources.push_back(new PointLightSource(float3(.8,0,.8), float3(0,-.7,.7)));
		// lightSources.push_back(new PointLightSource(float3(.6, .6, .6), float3(.5,.5,.5)));
		// ADD MATERIALS HERE
		materials.push_back(new Material());
		materials.at(0)->reflective = true;
		materials.push_back(new Material());
		materials.at(1)->kd = float3(1,1,1);
		materials.at(1)->ks = float3(1,1,1);

		// ADD OBJECTS HERE
		// objects.push_back(new Sphere(float3(0, 0, 0), 1, materials.at(0)));
		objects.push_back(new Quadric(float3(0,-1,0), float4(6,8,9,-1), materials.at(1)));
		objects.push_back(new Sphere(float3(.75,.75,-.75), 1, materials.at(0)));
		objects.push_back(new Quadric(float3(0,2.5,1), float4(4,6,7,-1), materials.at(1)));
		objects.push_back(new Plane(float3(0,1,0), float3(0,-2.2, 0), materials.at(1)));

	}
	~Scene()
	{
		for (std::vector<LightSource*>::iterator iLightSource = lightSources.begin(); iLightSource != lightSources.end(); ++iLightSource)
			delete *iLightSource;
		for (std::vector<Material*>::iterator iMaterial = materials.begin(); iMaterial != materials.end(); ++iMaterial)
			delete *iMaterial;
		for (std::vector<Intersectable*>::iterator iObject = objects.begin(); iObject != objects.end(); ++iObject)
			delete *iObject;		
	}

	float max(float a, float b) {
		if (a < b) 
			return b;
		else
			return a;
	}

	float F(float3 position, float3 dir) {

	}

public:
	Camera& getCamera()
	{
		return camera;
	}

	// IMPLEMENTED FOR YOUR CONVENIENCE, CALL THIS WHEN APPROPRIATE
	Hit firstIntersect(const Ray& ray)
	{
		Hit bestHit;
		bestHit.t = FLT_MAX;
		for(int oi=0; oi < objects.size(); oi++)
		{
			Hit hit = objects[oi]->intersect(ray);
			if(hit.t > 0 && hit.t < bestHit.t)
				bestHit = hit;
		}
		if(bestHit.t == FLT_MAX)
			return Hit();
		return bestHit;
	}

	float3 trace(const Ray& ray)
	{
		if (objects.size() < 1)
			return ray.dir * ray.dir;

		// Hit hit = firstIntersect(ray);

		// if(hit.t < 0)
		// 	return ray.dir * ray.dir;

		// float3 lightSum = float3(0,0,0);
		// float3 shade = float3(0,0,0);
		// for (int i = 0; i < lightSources.size(); i++) {
		// 	lightSum += lightSources.at(i)->getPowerDensityAt(hit.position) 
		// 		* max(0.0, hit.normal.dot(lightSources.at(i) ->getLightDirAt(hit.position)));
		// 	if (hit.material != NULL)
		// 		shade += hit.material->shade(hit.position, hit.normal, camera.getEye(),
		// 			lightSources.at(i)->getLightDirAt(hit.position),
		// 			lightSources.at(i)->getPowerDensityAt(hit.position));
		// }

		// return lightSum;//shade;

		Hit hit = firstIntersect(ray);

		if(hit.t < 0) return ray.dir * ray.dir; // nothing

		float3 outRadiance (0, 0, 0);
		for(int i = 0; i < lightSources.size(); i++) {
			Ray shadowRay(hit.position + hit.normal * .0001, lightSources.at(i)->getLightDirAt(hit.position));
			Hit shadowHit = firstIntersect(shadowRay);
			if(shadowHit.t < 0 || shadowHit.t > lightSources.at(i)->getDistanceFrom(hit.position) ) {
			    outRadiance += hit.material->shade(hit.position, hit.normal, -ray.dir, lightSources.at(i)->getLightDirAt(hit.position),
			    	lightSources.at(i)->getPowerDensityAt(hit.position));
			}
		}
		if(hit.material->reflective){
		    float3 reflectionDir = hit.material->reflect(hit.position, hit.normal);
		    Ray reflectedRay(hit.position + hit.normal*.0001, reflectionDir);
		    outRadiance += trace(reflectedRay); //* .4 * //F(hit.position,hit.normal);
		}
		if(hit.material->refractive) {
		    float3 refractionDir = hit.material->refract(hit.position, hit.normal);
		    Ray refractedRay(hit.position + hit.normal * .0001, refractionDir );
		    outRadiance += trace(refractedRay)*(float3(1,1,1));//-F(V,N));
		}
		return outRadiance;


		// return lightSum;// + hit.normal; 
		//return hit.normal;
	}

	float3 trace(Ray ray, int depth) {
	    if(depth > maxDepth) return ray.dir * ray.dir;
	    Hit hit = firstIntersect(ray);
	    if(hit.t < 0) return ray.dir * ray.dir; // nothing

	    float3 outRadiance(0, 0, 0);
	    for(int i = 0; i < lightSources.size(); i++) {
		    Ray shadowRay(hit.position + hit.normal * .0001, lightSources.at(i)->getLightDirAt(hit.position));
		    Hit shadowHit = firstIntersect(shadowRay);
		    if(shadowHit.t < 0 || shadowHit.t > lightSources.at(i)->getDistanceFrom(hit.position) ) {
			    outRadiance += hit.material->shade(hit.position, hit.normal, -ray.dir, lightSources.at(i)->getLightDirAt(hit.position),
			    	lightSources.at(i)->getPowerDensityAt(hit.position));
			}
	    }
	    if(hit.material->reflective){
	    	if (depth != 0){
			    float3 reflectionDir = hit.material->reflect(hit.position, hit.normal);
			    Ray reflectedRay(hit.position + hit.normal*.0001, reflectionDir);
			    outRadiance += trace(reflectedRay, depth+1) * hit.normal * .4;//F(hit.position,hit.normal);
			}
	    }
	    if(hit.material->refractive) {
		    float3 refractionDir = hit.material->refract(hit.position, hit.normal);
		    Ray refractedRay(hit.position + hit.normal * .0001, refractionDir );
		    outRadiance += trace(refractedRay, depth+1)*(float3(1,1,1)-hit.normal*.4);
	    }
	    return outRadiance;
	}

};

////////////////////////////////////////////////////////////////////////////////////////////////////////
// global application data

// screen resolution
const int screenWidth = 600;
const int screenHeight = 600;
// image to be computed by ray tracing
float3 image[screenWidth*screenHeight];

Scene scene;

bool computeImage()
{
	static unsigned int iPart = 0;

	if(iPart >= 64)
		return false;
    for(int j = iPart; j < screenHeight; j+=64)
	{
        for(int i = 0; i < screenWidth; i++)
		{
			float3 pixelColor = float3(0, 0, 0);
			float2 ndcPixelCentre( (2.0 * i - screenWidth) / screenWidth, (2.0 * j - screenHeight) / screenHeight );

			Camera& camera = scene.getCamera();
			Ray ray = Ray(camera.getEye(), camera.rayDirFromNdc(ndcPixelCentre));
			
			image[j*screenWidth + i] = scene.trace(ray,0);
		}
	}
	iPart++;
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OpenGL starts here. In the ray tracing example, OpenGL just outputs the image computed to the array.

// display callback invoked when window needs to be redrawn
void onDisplay( ) {
    glClearColor(0.1f, 0.2f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // clear screen

	if(computeImage())
		glutPostRedisplay();
    glDrawPixels(screenWidth, screenHeight, GL_RGB, GL_FLOAT, image);
 
    glutSwapBuffers(); // drawing finished
}

int main(int argc, char **argv) {
    glutInit(&argc, argv);						// initialize GLUT
    glutInitWindowSize(600, 600);				// startup window size 
    glutInitWindowPosition(100, 100);           // where to put window on screen
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);    // 8 bit R,G,B,A + double buffer + depth buffer
 
    glutCreateWindow("Ray caster");				// application window is created and displayed
 
    glViewport(0, 0, screenWidth, screenHeight);

    glutDisplayFunc(onDisplay);					// register callback
 
    glutMainLoop();								// launch event handling loop
    
    return 0;
}

