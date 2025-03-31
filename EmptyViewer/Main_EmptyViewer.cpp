#include <Windows.h>
#include <iostream>
#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/freeglut.h>

#define GLFW_INCLUDE_GLU
#define GLFW_DLL
#include <GLFW/glfw3.h>
#include <vector>
#include <utility>

#define GLM_SWIZZLE
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>

using namespace glm;
using namespace std;

// -------------------------------------------------
// Global Variables
// -------------------------------------------------
int Width = 512;
int Height = 512;
std::vector<float> OutputImage;
// -------------------------------------------------

class Scene;


class Ray {
private:
	vec3 origin;
	vec3 direction;
	vec3 intensity;
public:
	Ray(vec3 origin, vec3 dir) : origin(origin), direction(dir)
	{
		intensity = vec3(0.0f, 0.0f, 0.0f);
	}

	Ray(vec3 origin, vec3 dir, vec3 intensity) : origin(origin), direction(dir), intensity(intensity) {}

	vec3 getOrigin() const
	{
		return origin;
	}

	vec3 getDirection() const
	{
		return direction;
	}

	vec3 getIntensity() const
	{
		return intensity;
	}
};

class Surface {
private:
	vec3 position;
public:
	Surface(vec3 pos) : position(pos) {}

	vec3 getPosition() const
	{
		return position;
	}

	virtual std::pair<Surface*, float> intersect(const Ray& ray, float tMin, float tMax) = 0;

	virtual vec3 shade(vec3& targetPoint, vec3& n, vec3& v, Ray& light, Scene* scene) = 0;

	vec3 phongShading(const vec3& targetPoint, const vec3& n, const vec3& v,
		const vec3& ka, const vec3& kd, const vec3& ks, float& specularPower, const Ray& light)
	{
		vec3 l = normalize(light.getOrigin() - targetPoint);
		vec3 diffuse = kd * light.getIntensity() * glm::max(dot(n, l), 0.0f);

		vec3 h = normalize(l + v);
		vec3 specular = ks * light.getIntensity() * glm::pow(glm::max(dot(n, h), 0.0f), specularPower);

		return ka * light.getIntensity() + diffuse + specular;
	}
};

class Camera {
private:
	vec3 eye;
	vec3 u, v, w;
	float l, r, b, t, d;
	int nx, ny;
public:
	Camera() {
		eye = vec3(0.0f, 0.0f, 0.0f);
		u = vec3(1.0f, 0.0f, 0.0f); //x
		v = vec3(0.0f, 1.0f, 0.0f); //y
		w = vec3(0.0f, 0.0f, 1.0f); //z, -direction
		l = -0.1f, r = 0.1f, b = -0.1f, t = 0.1f;
		d = 0.1f; //distance
		nx = 512, ny = 512;
	}

	Ray getRay(int x, int y)
	{
		float u_coord = l + (r - l) * (x + 0.5f) / nx;
		float v_coord = b + (t - b) * (y + 0.5f) / ny;

		vec3 direction = normalize(u_coord * u + v_coord * v - d * w);

		return Ray(eye, direction);
	}
};

class Scene {
private:
	vector<Surface*> objects; // sphere 3 and plane 1
	vector<Ray> lights;
	Camera camera;
	vec3 background_color; // 0.0f 0.0f 0.0f
public:
	Scene() {
		background_color = vec3(0.0f, 0.0f, 0.0f);
	}

	void addObject(Surface* obj)
	{
		objects.push_back(obj);
	}

	void addLight(Ray light)
	{
		lights.push_back(light);
	}

	void setCamera(Camera c)
	{
		camera = c;
	}

	vector<Surface*> getObjects()
	{
		return objects;
	}

	vec3 trace(const Ray& ray, float tMin, float tMax)
	{
		Surface* closestSurface = nullptr;
		float closestT = FLT_MAX;
		vec3 intersectionPoint, n;

		for (Surface* object : objects)
		{
			pair<Surface*, float> result = object->intersect(ray, tMin, tMax);
			Surface* surface = result.first;
			float t = result.second;

			if (surface && t < closestT)
			{
				closestSurface = surface;
				closestT = t;
				intersectionPoint = ray.getOrigin() + t * ray.getDirection();
				n = normalize(intersectionPoint - surface->getPosition());
			}
		}

		if (!closestSurface)
			return background_color;

		vec3 v = normalize(-ray.getDirection());

		return closestSurface->shade(intersectionPoint, n, v, lights[0], this);
	}

	void render()
	{
		//Create our image. We don't want to do this in 
		//the main loop since this may be too slow and we 
		//want a responsive display of our beautiful image.
		//Instead we draw to another buffer and copy this to the 
		//framebuffer using glDrawPixels(...) every refresh

		OutputImage.clear();
		OutputImage.reserve(Width * Height * 3);

		for (int j = 0; j < Height; ++j)
		{
			for (int i = 0; i < Width; ++i)
			{
				// ---------------------------------------------------
				// --- Implement your code here to generate the image
				// ---------------------------------------------------

				Ray ray = camera.getRay(i, j);
				vec3 color = trace(ray, 0.001f, FLT_MAX);

				float gamma = 2.2f;
				vec3 corrected_color = pow(color, vec3(1.0f / gamma));

				OutputImage.push_back(corrected_color.r);
				OutputImage.push_back(corrected_color.g);
				OutputImage.push_back(corrected_color.b);

			}
		}
	}
};

class Plane : public Surface {
private:
	vec3 normal;
	float d;
public:
	Plane(vec3 pos, vec3 normal, float d) : Surface(pos), normal(normal), d(d) {}

	std::pair<Surface*, float> intersect(const Ray& ray, float tMin, float tMax) override
	{
		float denom = glm::dot(ray.getDirection(), normal);

		if (abs(denom) < 1e-6) {
			return { nullptr, FLT_MAX };
		}

		float t = glm::dot(getPosition() - ray.getOrigin(), normal) / denom;

		if (t < tMin || t > tMax) {
			return { nullptr, FLT_MAX };
		}

		return { this, t };
	}

	vec3 shade(vec3& targetPoint, vec3& n, vec3& v, Ray& light, Scene* scene) override
	{
		vec3 ka, kd, ks;
		float specularPower;

		ka = vec3(0.2f, 0.2f, 0.2f);
		kd = vec3(1.0f, 1.0f, 1.0f);
		ks = vec3(0.0f, 0.0f, 0.0f);
		specularPower = 0.0f;


		vec3 lightDir = normalize(light.getOrigin() - targetPoint);
		Ray shadowRay(targetPoint + n * 0.001f, lightDir);

		bool isInShadow = false;

		for (Surface* object : scene->getObjects())
		{
			if (object == this)
				continue;

			pair<Surface*, float> intersection = object->intersect(shadowRay, 0.001f, FLT_MAX);

			if (intersection.first != nullptr)
			{
				isInShadow = true;
				break;
			}
		}
		
		if (isInShadow)
			return vec3(0.0f, 0.0f, 0.0f);

		return phongShading(targetPoint, n, v, ka, kd, ks, specularPower, light);
	}
};

class Sphere : public Surface {
private:
	float r;
public:
	Sphere(vec3 pos, float r) : Surface(pos), r(r) {}

	std::pair<Surface*, float> intersect(const Ray& ray, float tMin, float tMax) override
	{
		vec3 oc = ray.getOrigin() - getPosition();
		vec3 d = ray.getDirection();

		float p = glm::length(oc);

		float pd = glm::dot(oc, d);

		float tm = -pd;

		float lm_square = p * p - (pd * pd);

		if (lm_square < 0.0f) {
			return { nullptr, FLT_MAX };
		}

		float delta_t = glm::sqrt(r * r - lm_square);

		float t0 = tm - delta_t;
		float t1 = tm + delta_t;

		if (t0 >= tMin && t0 <= tMax) {
			return { this, t0 };
		}
		if (t1 >= tMin && t1 <= tMax) {
			return { this, t1 };
		}

		return { nullptr, FLT_MAX };
	}

	vec3 shade(vec3& targetPoint, vec3& n, vec3& v, Ray& light, Scene* scene) override
	{
		vec3 ka, kd, ks;
		float specularPower;

		if (getPosition() == vec3(-4.0f, 0.0f, -7.0f)) {
			ka = vec3(0.2f, 0.0f, 0.0f);
			kd = vec3(1.0f, 0.0f, 0.0f);
			ks = vec3(0.0f, 0.0f, 0.0f);
			specularPower = 0.0f;
		}
		else if (getPosition() == vec3(0.0f, 0.0f, -7.0f)) {
			ka = vec3(0.0f, 0.2f, 0.0f);
			kd = vec3(0.0f, 0.5f, 0.0f);
			ks = vec3(0.5f, 0.5f, 0.5f);
			specularPower = 32.0f;
		}
		else {
			ka = vec3(0.0f, 0.0f, 0.2f);
			kd = vec3(0.0f, 0.0f, 1.0f);
			ks = vec3(0.0f, 0.0f, 0.0f);
			specularPower = 0.0f;
		}

		vec3 lightDir = normalize(light.getOrigin() - targetPoint);
		Ray shadowRay(targetPoint + n * 0.001f, lightDir);

		bool isInShadow = false;

		for (Surface* object : scene->getObjects())
		{
			if (object == this)
				continue;

			pair<Surface*, float> intersection = object->intersect(shadowRay, 0.001f, FLT_MAX);

			if (intersection.first != nullptr)
			{
				isInShadow = true;
				break;
			}
		}

		if (isInShadow)
			return vec3(0.0f, 0.0f, 0.0f);

		return phongShading(targetPoint, n, v, ka, kd, ks, specularPower, light);
	}
};




void resize_callback(GLFWwindow* window, int nw, int nh) 
{
	//This is called in response to the window resizing.
	//The new width and height are passed in so we make 
	//any necessary changes:
	Width = nw;
	Height = nh;
	//Tell the viewport to use all of our screen estate
	glViewport(0, 0, nw, nh);

	//This is not necessary, we're just working in 2d so
	//why not let our spaces reflect it?
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glOrtho(0.0, static_cast<double>(Width)
		, 0.0, static_cast<double>(Height)
		, 1.0, -1.0);

	//Reserve memory for our render so that we don't do 
	//excessive allocations and render the image
	OutputImage.reserve(Width * Height * 3);

	Scene* scene = static_cast<Scene*>(glfwGetWindowUserPointer(window));

	if (scene)
		scene->render();
}


int main(int argc, char* argv[])
{
	// -------------------------------------------------
	// Initialize Window
	// -------------------------------------------------

	GLFWwindow* window;

	Scene* scene = new Scene();
	Surface* p1 = new Plane(vec3(0.0f, -2.0f, 0.0f), vec3(0.0f, 1.0f, 0.0f), 2);
	Surface* s1 = new Sphere(vec3(-4.0f, 0.0f, -7.0f), 1.0f);
	Surface* s2 = new Sphere(vec3(0.0f, 0.0f, -7.0f), 2.0f);
	Surface* s3 = new Sphere(vec3(4.0f, 0.0f, -7.0f), 1.0f);
	Camera camera = Camera();
	Ray light = Ray(vec3(-4.0f, 4.0f, -3.0f), vec3(0.0f, 0.0f, 0.0f), vec3(1.0f, 1.0f, 1.0f));

	scene->addObject(p1);
	scene->addObject(s1);
	scene->addObject(s2);
	scene->addObject(s3);
	scene->setCamera(camera);
	scene->addLight(light);

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(Width, Height, "OpenGL Viewer", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);

	//We have an opengl context now. Everything from here on out 
	//is just managing our window or opengl directly.

	//Tell the opengl state machine we don't want it to make 
	//any assumptions about how pixels are aligned in memory 
	//during transfers between host and device (like glDrawPixels(...) )
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);

	//We call our resize function once to set everything up initially
	//after registering it as a callback with glfw
	glfwSetWindowUserPointer(window, scene);
	glfwSetFramebufferSizeCallback(window, resize_callback);
	resize_callback(window, Width, Height);

	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		//Clear the screen
		glClear(GL_COLOR_BUFFER_BIT);

		// -------------------------------------------------------------
		//Rendering begins!
		glDrawPixels(Width, Height, GL_RGB, GL_FLOAT, &OutputImage[0]);
		//and ends.
		// -------------------------------------------------------------

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();

		//Close when the user hits 'q' or escape
		if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS
			|| glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
		{
			glfwSetWindowShouldClose(window, GL_TRUE);
		}
	}


	delete p1;
	delete s1;
	delete s2;
	delete s3;
	delete scene;


	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}
