#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::ivec2;
using glm::mat3;
using glm::vec3;

struct Pixel
{
	int x;
	int y;
	float zinv;
};

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface *screen;
int t;
vector<Triangle> triangles;

vec3 cameraPos(0, 0, -3.001);
mat3 R;
float yaw = 0;	   // Yaw angle controlling camera rotation around y-axis
float cameraAngle; // angle controlling camera rotation around x-axis
float tilt;		   // angle controlling camera rotation around z-axis
const int focalLength = SCREEN_HEIGHT;

vec3 currentColor;

// ----------------------------------------------------------------------------
// DEPTH BUFFER
float depthBuffer[SCREEN_WIDTH][SCREEN_HEIGHT];

// ----------------------------------------------------------------------------
// FUNCTIONS

void Update();
void Draw();
void VertexShader(const vec3 &v, Pixel &p);
void Interpolate(Pixel a, Pixel b, vector<Pixel> &result);
void DrawLineSDL(SDL_Surface *surface, Pixel a, Pixel b, vec3 color);
void DrawPolygonEdges(const vector<vec3> &vertices);
void ComputePolygonRows(const vector<Pixel> &vertexPixels, vector<Pixel> &leftPixels, vector<Pixel> &rightPixels);
void DrawPolygonRows(const vector<Pixel> &leftPixels, const vector<Pixel> &rightPixels);
void DrawPolygon(const vector<vec3> &vertices);

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

	SDL_SaveBMP(screen, "screenshot.bmp");

	return 0;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2 - t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

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
}

void Draw()
{
	for (int y = 0; y < SCREEN_HEIGHT; ++y)
	{
		for (int x = 0; x < SCREEN_WIDTH; ++x)
		{
			depthBuffer[x][y] = numeric_limits<float>::max();
		}
	}

	SDL_FillRect(screen, 0, 0);

	if (SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);

	for (int i = 0; i < triangles.size(); ++i)
	{
		currentColor = triangles[i].color;
		vector<vec3> vertices(3);
		vertices[0] = triangles[i].v0;
		vertices[1] = triangles[i].v1;
		vertices[2] = triangles[i].v2;
		DrawPolygon(vertices);
	}

	if (SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);
	SDL_UpdateRect(screen, 0, 0, 0, 0);
}

void VertexShader(const vec3 &v, Pixel &p)
{
	vec3 c = (cameraPos - v);
	p.x = int(focalLength * c.x / c.z) + SCREEN_WIDTH / 2;
	p.y = int(focalLength * c.y / c.z) + SCREEN_HEIGHT / 2;
	p.zinv = glm::length(c);
}

void Interpolate(Pixel a, Pixel b, vector<Pixel> &result)
{
	int N = result.size();
	float stepX = float(b.x - a.x) / float(glm::max(N - 1, 1));
	float stepY = float(b.y - a.y) / float(glm::max(N - 1, 1));
	float stepZinv = float(b.zinv - a.zinv) / float(glm::max(N - 1, 1));
	Pixel current(a);
	for (int i = 0; i < N; ++i)
	{
		result[i] = current;
		current.x = a.x + i * stepX;
		current.y = a.y + i * stepY;
		current.zinv = a.zinv + i * stepZinv;
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
			if (line[i].zinv <= depthBuffer[line[i].x][line[i].y])
			{
				PutPixelSDL(screen, line[i].x, line[i].y, color);
				depthBuffer[line[i].x][line[i].y] = line[i].zinv;
			}
		}
	}
}

void DrawPolygonEdges(const vector<vec3> &vertices)
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

void ComputePolygonRows(const vector<Pixel> &vertexPixels, vector<Pixel> &leftPixels,
						vector<Pixel> &rightPixels)
{
	// 1. Find max and min y-value of the polygon
	// and compute the number of rows it occupies.
	int max = numeric_limits<int>::min();
	int min = numeric_limits<int>::max();

	for (Pixel vertex : vertexPixels)
	{
		if (vertex.y < min)
		{
			min = vertex.y;
		}
		if (vertex.y > max)
		{
			max = vertex.y;
		}
	}
	// 2. Resize leftPixels and rightPixels
	// so that they have an element for each row.
	leftPixels.resize(max - min + 1);
	rightPixels.resize(max - min + 1);

	// 3. Initialize the x-coordinates in leftPixels
	// to some really large value and the x-coordinates
	// in rightPixels to some really small value.
	for (int i = 0; i < leftPixels.size(); i++)
	{
		leftPixels[i].x = SCREEN_WIDTH;
		leftPixels[i].y = min + i;

		rightPixels[i].x = 0;
		rightPixels[i].y = min + i;
	}

	// 4. Loop through all edges of the polygon and use
	// linear interpolation to find the x-coordinate for
	// each row it occupies. Update the corresponding
	// values in rightPixels and leftPixels.
	for (int i = 0; i < vertexPixels.size(); i++)
	{
		Pixel currentVertex = vertexPixels[i];
		Pixel nextVertex = vertexPixels[(i + 1) % vertexPixels.size()];

		int steps = glm::abs(nextVertex.y - currentVertex.y) + 1;
		vector<Pixel> line(steps);
		Interpolate(currentVertex, nextVertex, line);

		for (Pixel l : line)
		{
			int row = l.y - min;

			if (l.x < leftPixels[row].x)
			{
				leftPixels[row].x = l.x;
				leftPixels[row].zinv = l.zinv;
			}
			if (l.x > rightPixels[row].x)
			{
				rightPixels[row].x = l.x;
				rightPixels[row].zinv = l.zinv;
			}
		}
	}
}

void DrawPolygonRows(const vector<Pixel> &leftPixels, const vector<Pixel> &rightPixels)
{
	for (int i = 0; i < leftPixels.size(); i++)
	{
		DrawLineSDL(screen, leftPixels[i], rightPixels[i], currentColor);
	}
}

void DrawPolygon(const vector<vec3> &vertices)
{
	int V = vertices.size();
	vector<Pixel> vertexPixels(V);
	for (int i = 0; i < V; ++i)
	{
		VertexShader(vertices[i], vertexPixels[i]);
	}
	vector<Pixel> leftPixels;
	vector<Pixel> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	DrawPolygonRows(leftPixels, rightPixels);
}