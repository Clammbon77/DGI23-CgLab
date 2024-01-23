#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::mat3;
using glm::vec3;

/*struct Triangle
{
	vec3 v0 ;
	vec3 v1 ;
	vec3 v2 ;
	vec3 normal ;
	vec3 color ;
} ;
*/
struct Intersection
{
	vec3 position;
	float distance;
	int triangleIndex;
};

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface *screen;
int t;

vector<Triangle> triangles;
float focalLength = SCREEN_HEIGHT;
vec3 cameraPos(0, 0, -3);
vec3 lightPos(0, -0.5, -0.7);
vec3 lightColor = 14.f * vec3(1, 1, 1);
vec3 indirectLight = 0.5f * vec3(1, 1, 1);

mat3 R = mat3(
	1, 0, 0,
	0, 1, 0,
	0, 0, 1); // Rotation Matrix

const float pi = 3.14159265358979323846f;
float yaw = 0.0; // Angle that the camera rotate around y-axis

// ----------------------------------------------------------------------------
// FUNCTIONS

void Update();
void Draw();
bool ClosestIntersection(
	vec3 start,
	vec3 dir,
	const vector<Triangle> &triangles,
	Intersection &closestIntersection);
vec3 DirectLight(const Intersection &i);

int main(int argc, char *argv[])
{
	screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT);
	t = SDL_GetTicks(); // Set start value for timer.

	LoadTestModel(triangles);

	while (NoQuitMessageSDL())
	{
		Update();
		Draw();
	}

	SDL_SaveBMP(screen, "screenshot.bmp");
	return 0;
}

void Rotate()
{
	R = mat3(
		cos(yaw), 0, sin(yaw),
		0, 1, 0,
		-sin(yaw), 0, cos(yaw));
}
void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2 - t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

	Uint8 *keystate = SDL_GetKeyState(0);
	float weight = (2.0f * pi / 1000.0f) * dt;
	// vectors representing the right (x-axis), down (y-axis) and forward (z-axis) directions
	vec3 right(R[0][0], R[0][1], R[0][2]);	 // x axis
	vec3 down(R[1][0], R[1][1], R[1][2]);	 // y axis
	vec3 forward(R[2][0], R[2][1], R[2][2]); // z axis

	// Control Camera Position
	if (keystate[SDLK_UP])
	{
		// Move camera forward
		cameraPos += weight * forward;
	}
	if (keystate[SDLK_DOWN])
	{
		// Move camera backward
		cameraPos -= weight * forward;
	}
	if (keystate[SDLK_LEFT])
	{
		// Move camera to the left
		cameraPos -= weight * right;
	}
	if (keystate[SDLK_RIGHT])
	{
		// Move camera to the right
		cameraPos += weight * right;
	}

	// Rotate Camera
	if (keystate[SDLK_q])
	{
		// Move camera to the right
		yaw -= weight * 1;
		Rotate();
	}

	if (keystate[SDLK_e])
	{
		// Move camera to the right
		yaw += weight * 1;
		Rotate();
	}

	// Control Light Position
	if (keystate[SDLK_w])
	{
		lightPos += weight * forward;
	}

	if (keystate[SDLK_a])
	{
		lightPos -= weight * right;
	}

	if (keystate[SDLK_s])
	{
		lightPos -= weight * forward;
	}

	if (keystate[SDLK_d])
	{
		lightPos += weight * right;
	}
}

void Draw()
{
	if (SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);

	for (int y = 0; y < SCREEN_HEIGHT; ++y)
	{
		for (int x = 0; x < SCREEN_WIDTH; ++x)
		{
			vec3 direction((x - SCREEN_WIDTH / 2), (y - SCREEN_HEIGHT / 2), focalLength);
			Intersection intersection;
			vec3 black(0, 0, 0);
			direction = R * glm::normalize(direction);
			if (ClosestIntersection(cameraPos, direction, triangles, intersection))
			{
				vec3 colour = triangles[intersection.triangleIndex].color * (DirectLight(intersection) + indirectLight);
				PutPixelSDL(screen, x, y, colour);
			}
			else
			{
				PutPixelSDL(screen, x, y, black);
			}
		}
	}

	if (SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect(screen, 0, 0, 0, 0);
}

bool ClosestIntersection(
	vec3 start,
	vec3 dir,
	const vector<Triangle> &triangles,
	Intersection &closestIntersection)
{
	bool exist = false;
	float m = std::numeric_limits<float>::max();

	for (int i = 0; i < triangles.size(); i++)
	{
		Triangle triangle = triangles[i];
		vec3 v0 = triangle.v0;
		vec3 v1 = triangle.v1;
		vec3 v2 = triangle.v2;
		vec3 e1 = v1 - v0;
		vec3 e2 = v2 - v0;
		vec3 b = start - v0;
		mat3 A(-dir, e1, e2);
		vec3 X = glm::inverse(A) * b;

		float t = X.x; // scalar coordinate for the intersection point position on the ray
		float u = X.y; // scalar coordinate for the intersection point position on e1 direction
		float v = X.z; // scalar coordinate for the intersection point position on e2 direction

		if (u >= 0 && v >= 0 && u + v <= 1)
		{
			if (t >= 0 && t <= m)
			{
				m = t; // Set the maximum to the current ray position to see if is there much closer intersection exists in the scenario
				closestIntersection.distance = t;
				closestIntersection.position = start + (dir * closestIntersection.distance);
				closestIntersection.triangleIndex = i;
			}
			exist = true;
		}
	}
	return exist;
}

vec3 DirectLight(const Intersection &i)
{
	// Let n̂ be a unit vector describing the normal pointing out from the surface
	// Let r̂ be a unit vector describing the direction from the surface point to the light source.

	vec3 n = triangles[i.triangleIndex].normal;
	vec3 r = glm::normalize(i.position - lightPos);
	vec3 D;

	float r_length = glm::distance(lightPos, i.position);

	Intersection intersection;

	// Cast another ray from surface to the light source using ClosetIntersection
	// Check if distance to the closest intersecting surface for surface point is closer than to the light source
	if (ClosestIntersection(lightPos, r, triangles, intersection))
	{
		if (intersection.distance < r_length)
		{
			return vec3(0, 0, 0);
		}
	}

	// lightColor describes the power P
	// i.e. energy per time unit E/t of the emitted light for each color component

	// D = B max(r̂ . n̂ , 0) = (P max (r̂ . n̂ , 0))/4πr^2
	D = vec3((lightColor * glm::max(glm::dot(-r, n), 0.0f)) / (4.0f * glm::pow(r_length, 2.0f) * pi));

	return D;
};
