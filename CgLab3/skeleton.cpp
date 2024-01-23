#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::mat3;
using glm::vec2;
using glm::vec3;

//----------------------------------------------------------------------------
// Structures
struct Pixel
{
	int x;
	int y;
	float zinv;
	vec3 pos3d;
};

struct Vertex
{
	vec3 position;
};

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface *screen;
int t;
vector<Triangle> triangles;
vec3 cameraPos(0, 0, -3.001);

float yaw = 0;	   // Angle controlling camera rotation around y-axis
float cameraAngle; // Angle controlling camera rotation around x-axis
float tilt;		   // Angle controlling camera rotation around z-axis

float focalLength = SCREEN_HEIGHT;
vec3 currentColor;
float depthBuffer[SCREEN_WIDTH][SCREEN_HEIGHT];

vec3 lightPos(0, -0.5, -0.7);
vec3 lightPower = 1.1f * vec3(1, 1, 1);
vec3 indirectLightPowerPerArea = 0.5f * vec3(1, 1, 1);

vec3 currentNormal;
vec3 currentReflectance;

// ----------------------------------------------------------------------------
// FUNCTIONS

void Update();
void Draw();
void VertexShader(const vec3 &v, Pixel &p);
void PixelShader(const Pixel &p);
void Interpolate(Pixel a, Pixel b, vector<Pixel> &result);
void DrawLineSDL(SDL_Surface *surface, Pixel a, Pixel b, vec3 color);
void DrawPolygonEdges(const vector<Vertex> &vertices);
void Rotate(Vertex &in);

void ComputePolygonRows(const vector<Pixel> &vertexPixels,
						vector<Pixel> &leftPixels,
						vector<Pixel> &rightPixels);
void DrawRows(const vector<Pixel> &leftPixels, const vector<Pixel> &rightPixels);
void DrawPolygon(const vector<Vertex> &vertices);
void DrawPolygonRows(vector<Pixel> leftPixels, vector<Pixel> rightPixels);

int main(int argc, char *argv[])
{
	LoadTestModel(triangles);
	screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT);
	t = SDL_GetTicks(); // Set start value for timer.

	while (NoQuitMessageSDL())
	{
		Update();
		Draw();
	}
	SDL_SaveBMP(screen, "screenshot_extra.bmp");
	return 0;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2 - t);
	t = t2;
	// cout << "Render time: " << dt << " ms." << endl;

	Uint8 *keystate = SDL_GetKeyState(0);

	if (keystate[SDLK_UP])
		yaw += 0.01;

	if (keystate[SDLK_DOWN])
		yaw -= 0.01;

	if (keystate[SDLK_RIGHT])
		cameraAngle += 0.01;

	if (keystate[SDLK_LEFT])
		cameraAngle -= 0.01;

	if (keystate[SDLK_RSHIFT])
		tilt += 0.01;

	if (keystate[SDLK_RCTRL])
		tilt -= 0.01;

	if (keystate[SDLK_w])
		cameraPos.y -= 0.01;

	if (keystate[SDLK_s])
		cameraPos.y += 0.01;

	if (keystate[SDLK_d])
		cameraPos.x += 0.01;

	if (keystate[SDLK_a])
		cameraPos.x -= 0.01;

	if (keystate[SDLK_e])
		cameraPos.z += 0.01;

	if (keystate[SDLK_q])
		cameraPos.z -= 0.01;

	if (keystate[SDLK_z])
		lightPos.y += 0.01;

	if (keystate[SDLK_x])
		lightPos.y -= 0.01;

	if (keystate[SDLK_c])
		lightPos.x += 0.01;

	if (keystate[SDLK_v])
		lightPos.x -= 0.01;

	if (keystate[SDLK_b])
		lightPos.z += 0.01;

	if (keystate[SDLK_n])
		lightPos.z -= 0.01;
}

void Draw()
{
	for (int x = 0; x < SCREEN_WIDTH; ++x)
	{
		for (int y = 0; y < SCREEN_HEIGHT; ++y)
		{
			depthBuffer[x][y] = 0;
		}
	}

	SDL_FillRect(screen, 0, 0);
	if (SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);

	for (int i = 0; i < triangles.size(); ++i)
	{
		vector<Vertex> vertices(3);
		vertices[0].position = triangles[i].v0;
		vertices[1].position = triangles[i].v1;
		vertices[2].position = triangles[i].v2;
        
		Rotate(vertices[0]);
		Rotate(vertices[1]);
		Rotate(vertices[2]);

		currentNormal = triangles[i].normal;
		currentReflectance = vec3(1, 1, 1) * 10.f;
		currentColor = triangles[i].color;

		DrawPolygon(vertices);
	}
	if (SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);
	SDL_UpdateRect(screen, 0, 0, 0, 0);
}

void VertexShader(const Vertex &v, Pixel &p)
{
	vec3 pos = (v.position - cameraPos);
	p.zinv = 1 / pos.z;
	p.x = int(focalLength * pos.x * p.zinv) + SCREEN_WIDTH / 2;
	p.y = int(focalLength * pos.y * p.zinv) + SCREEN_HEIGHT / 2;
	p.pos3d = v.position;
}

void PixelShader(const Pixel &p)
{
	int x = p.x;
	int y = p.y;
	if (p.zinv >= depthBuffer[x][y])
	{
		depthBuffer[x][y] = p.zinv;

		float r = glm::length(p.pos3d - lightPos);
		float d = glm::dot(currentNormal, glm::normalize(lightPos - p.pos3d)) / (4 * 3 * r * r);
		vec3 light = (lightPower * SDL_max(d, 0)) * currentReflectance + indirectLightPowerPerArea;

		PutPixelSDL(screen, x, y, light * currentColor);
	}
}

void Interpolate(Pixel a, Pixel b, vector<Pixel> &result)
{
	int N = result.size();
	float stepX = float(b.x - a.x) / float(glm::max(N - 1, 1));
	float stepY = float(b.y - a.y) / float(glm::max(N - 1, 1));
	float stepZinv = float(b.zinv - a.zinv) / float(glm::max(N - 1, 1));
	vec3 step3d = (b.pos3d - a.pos3d) / float(glm::max(N - 1, 1));

	Pixel current(a);
	for (int i = 0; i < N; ++i)
	{
		result[i] = current;
		current.x = a.x + i * stepX;
		current.y = a.y + i * stepY;
		current.zinv = a.zinv + i * stepZinv;
		current.pos3d += step3d;
	}
}

void DrawLineSDL(SDL_Surface *surface, Pixel a, Pixel b, vec3 color)
{
	Pixel delta;
	delta.x = glm::abs(a.x - b.x);
	delta.y = glm::abs(a.y - b.y);

	int pixels = glm::max(delta.x, delta.y) + 1;

	vector<Pixel> line(pixels);
	Interpolate(a, b, line);

	for (int i = 0; i < line.size(); i++)
	{
		if (line[i].x >= 0 && line[i].x < SCREEN_WIDTH && line[i].y >= 0 && line[i].y < SCREEN_HEIGHT)
		{
			PixelShader(line[i]);
		}
	}
}

void DrawPolygonEdges(const vector<Vertex> &vertices)
{
	int V = vertices.size();
	// Transform each vertex from 3D world position to 2D image position:
	vector<Pixel> projectedVertices(V);
	for (int i = 0; i < V; ++i)
	{
		VertexShader(vertices[i], projectedVertices[i]);
	}
	// Loop over all vertices and draw the edge from it to the next vertex:
	for (int i = 0; i < V; ++i)
	{
		int j = (i + 1) % V; // The next vertex
		vec3 color(1, 1, 1);
		DrawLineSDL(screen, projectedVertices[i], projectedVertices[j], color);
	}
}

void Rotate(Vertex &in)
{
	// Rotate around y axis
	vec3 row1(cos(cameraAngle), 0, sin(cameraAngle));
	vec3 row2(0, 1, 0);
	vec3 row3(-sin(cameraAngle), 0, cos(cameraAngle));

	in.position.x = glm::dot(row1, in.position);
	in.position.y = glm::dot(row2, in.position);
	in.position.z = glm::dot(row3, in.position);

	// Rotate around x axis
	row1 = vec3(1, 0, 0);
	row2 = vec3(0, cos(yaw), -sin(yaw));
	row3 = vec3(0, sin(yaw), cos(yaw));

	in.position.x = glm::dot(row1, in.position);
	in.position.y = glm::dot(row2, in.position);
	in.position.z = glm::dot(row3, in.position);

	// Rotate around z axis
	row1 = vec3(cos(tilt), -sin(tilt), 0);
	row2 = vec3(sin(tilt), cos(tilt), 0);
	row3 = vec3(0, 0, 1);

	in.position.x = glm::dot(row1, in.position);
	in.position.y = glm::dot(row2, in.position);
	in.position.z = glm::dot(row3, in.position);
}

void DrawPolygon(const vector<Vertex> &vertices)
{
	// cout << "DrawPolygon" << endl;
	int V = vertices.size();
	vector<Pixel> vertexPixels(V);
	for (int i = 0; i < V; ++i)
		VertexShader(vertices[i], vertexPixels[i]);
	vector<Pixel> leftPixels;
	vector<Pixel> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	DrawPolygonRows(leftPixels, rightPixels);
}

void ComputePolygonRows(const vector<Pixel> &vertexPixels,
						vector<Pixel> &leftPixels,
						vector<Pixel> &rightPixels)
{
	// 1. Find max and min y-value of the polygon
	// and compute the number of rows it occupies.
	int max = numeric_limits<int>::min();
	int min = numeric_limits<int>::max();

	for (Pixel v : vertexPixels)
	{
		if (v.y < min)
		{
			min = v.y;
		}
		if (v.y > max)
		{
			max = v.y;
		}
	}

	// 2. Resize leftPixels and rightPixels
	// so that they have an element for each row.
	leftPixels.resize(max - min + 1);
	rightPixels.resize(max - min + 1);

	// 3. Initialize the x-	coordinates in leftPixels
	// to some really large value and the x-coordinates
	// in rightPixels to some really small value.
	for (int i = 0; i < leftPixels.size(); i++)
	{
		leftPixels[i].x = SCREEN_WIDTH;
		leftPixels[i].y = min + i;

		rightPixels[i].x = 0;
		rightPixels[i].y = min + i;
	}

	// 4. Loop through all edges of	the polygon and use
	// linear interpolation to find the x-coordinate for
	// each row it occupies. Update the corresponding
	// values in rightPixels and leftPixels.

	for (int j = 0; j < vertexPixels.size(); j++)
	{
		Pixel p1 = vertexPixels[j];
		Pixel p2 = vertexPixels[(j + 1) % vertexPixels.size()];

		int steps = abs(p1.y - p2.y) + 1;
		vector<Pixel> line = vector<Pixel>(steps);
		Interpolate(p1, p2, line);

		for (Pixel l : line)
		{
			int row = l.y - min;

			if (l.x < leftPixels[row].x)
			{
				leftPixels[row].x = l.x;
				leftPixels[row].zinv = l.zinv;
				leftPixels[row].pos3d = l.pos3d;
			}
			if (l.x >= rightPixels[row].x)
			{
				rightPixels[row].x = l.x;
				rightPixels[row].zinv = l.zinv;
				rightPixels[row].pos3d = l.pos3d;
			}
		}
	}
}

void DrawPolygonRows(vector<Pixel> leftPixels, vector<Pixel> rightPixels)
{
	for (int i = 0; i < leftPixels.size(); i++)
	{
		DrawLineSDL(screen, leftPixels[i], rightPixels[i], currentColor);
	}
}