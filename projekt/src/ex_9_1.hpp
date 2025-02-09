#include "glew.h"
#include <GLFW/glfw3.h>
#include "glm.hpp"
#include "ext.hpp"
#include <iostream>
#include <cmath>
#include <vector>

#include "Shader_Loader.h"
#include "Render_Utils.h"
#include "Texture.h"

#include "Box.cpp"
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <string>
#include "SOIL/SOIL.h"
#include <utility>
#include <limits>

#define AMPLITUDE 150.0f
#define OCTAVES 4
#define ROUGHNESS 0.01f

#define GRID_SIZE 250
#define HEIGHT_SCALE 0.1f

float generateHeight(int x, int z, int seed, int xOffset, int zOffset);

const unsigned int SHADOW_WIDTH = 1024, SHADOW_HEIGHT = 1024;

int WIDTH = 500, HEIGHT = 500;

namespace models {
	Core::RenderContext bedContext;
	Core::RenderContext chairContext;
	Core::RenderContext deskContext;
	Core::RenderContext doorContext;
	Core::RenderContext drawerContext;
	Core::RenderContext marbleBustContext;
	Core::RenderContext materaceContext;
	Core::RenderContext pencilsContext;
	Core::RenderContext planeContext;
	Core::RenderContext roomContext;
	Core::RenderContext spaceshipContext;
	Core::RenderContext sphereContext;
	Core::RenderContext windowContext;
	Core::RenderContext testContext;
}

namespace texture {
	GLuint earth;
	GLuint earth_map;
	GLuint moon;
	GLuint moon_map;
	GLuint chair;
	GLuint chair_map;
	GLuint table;
	GLuint table_map;
	GLuint bird1_map;
	GLuint bird1;
	GLuint bird2_map;
	GLuint bird2;
	GLuint bird3_map;
	GLuint bird3;
}
GLuint programTex;

GLuint depthMapFBO;
GLuint depthMap;

GLuint program;
GLuint programSun;
GLuint programTest;

Core::Shader_Loader shaderLoader;

Core::RenderContext shipContext;
Core::RenderContext sphereContext;
Core::RenderContext birdContext;

glm::vec3 sunPos = glm::vec3(-4.740971f, 2.149999f, 0.369280f);
glm::vec3 sunDir = glm::vec3(-0.93633f, 0.351106, 0.003226f);
glm::vec3 sunColor = glm::vec3(0.9f, 0.9f, 0.7f)*5;

glm::vec3 cameraPos = glm::vec3(0.479490f, 1.250000f, -2.124680f);
glm::vec3 cameraDir = glm::vec3(-0.354510f, 0.000000f, 0.935054f);


glm::vec3 spaceshipPos = glm::vec3(0.065808f, 1.250000f, -2.189549f);
glm::vec3 spaceshipDir = glm::vec3(-0.490263f, 0.000000f, 0.871578f);
GLuint VAO,VBO;

float aspectRatio = 1.f;

float exposition = 1.f;

glm::vec3 pointlightPos = glm::vec3(0, 2, 0);
glm::vec3 pointlightColor = glm::vec3(0.9, 0.6, 0.6);

glm::vec3 spotlightPos = glm::vec3(0, 0, 0);
glm::vec3 spotlightConeDir = glm::vec3(0, 0, 0);
glm::vec3 spotlightColor = glm::vec3(0.4, 0.4, 0.9)*3;
float spotlightPhi = 3.14 / 4;


GLuint programDepth;
bool shadowsEnabled = true;
bool boidsPause = false;


float lastTime = -1.f;
float deltaTime = 0.f;

struct AABB {
	float xMin, yMin, zMin;
	float xMax, yMax, zMax;

	AABB(float xMin, float yMin, float zMin, float xMax, float yMax, float zMax)
		: xMin(xMin), yMin(yMin), zMin(zMin), xMax(xMax), yMax(yMax), zMax(zMax) {}

};

struct Boid {
	glm::vec3 position;
	glm::vec3 velocity;
	glm::vec3 direction;
	int groupId;
	AABB boundingBox;
	Boid()
		: position(0.0f, 0.0f, 0.0f), velocity(0.0f, 0.0f, 0.0f), direction(0.0f, 0.0f, 0.0f), groupId(0), boundingBox(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f) {}

};

struct Obstacle {
	glm::vec3 position;
	AABB boundingBox;
	Obstacle()
		: position(0.0f, 0.0f, 0.0f), boundingBox(0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f) {}
};

void updateDeltaTime(float time) {
	if (lastTime < 0) {
		lastTime = time;
		return;
	}

	deltaTime = time - lastTime;
	if (deltaTime > 0.1) deltaTime = 0.1;
	lastTime = time;
}
glm::mat4 createCameraMatrix()
{
	glm::vec3 cameraSide = glm::normalize(glm::cross(cameraDir,glm::vec3(0.f,1.f,0.f)));
	glm::vec3 cameraUp = glm::normalize(glm::cross(cameraSide,cameraDir));
	glm::mat4 cameraRotrationMatrix = glm::mat4({
		cameraSide.x,cameraSide.y,cameraSide.z,0,
		cameraUp.x,cameraUp.y,cameraUp.z ,0,
		-cameraDir.x,-cameraDir.y,-cameraDir.z,0,
		0.,0.,0.,1.,
		});
	cameraRotrationMatrix = glm::transpose(cameraRotrationMatrix);
	glm::mat4 cameraMatrix = cameraRotrationMatrix * glm::translate(-cameraPos);

	return cameraMatrix;
}

glm::mat4 createPerspectiveMatrix()
{
	
	glm::mat4 perspectiveMatrix;
	float n = 0.05;
	float f = 20.;
	float a1 = glm::min(aspectRatio, 1.f);
	float a2 = glm::min(1 / aspectRatio, 1.f);
	perspectiveMatrix = glm::mat4({
		1,0.,0.,0.,
		0.,aspectRatio,0.,0.,
		0.,0.,(f+n) / (n - f),2*f * n / (n - f),
		0.,0.,-1.,0.,
		});

	
	perspectiveMatrix=glm::transpose(perspectiveMatrix);

	return perspectiveMatrix;
}

void drawObjectPBR(Core::RenderContext& context, glm::mat4 modelMatrix, glm::vec3 color, float roughness, float metallic) {

	glm::mat4 viewProjectionMatrix = createPerspectiveMatrix() * createCameraMatrix();
	glm::mat4 transformation = viewProjectionMatrix * modelMatrix;
	glUniformMatrix4fv(glGetUniformLocation(program, "transformation"), 1, GL_FALSE, (float*)&transformation);
	glUniformMatrix4fv(glGetUniformLocation(program, "modelMatrix"), 1, GL_FALSE, (float*)&modelMatrix);

	glUniform1f(glGetUniformLocation(program, "exposition"), exposition);

	glUniform1f(glGetUniformLocation(program, "roughness"), roughness);
	glUniform1f(glGetUniformLocation(program, "metallic"), metallic);

	glUniform3f(glGetUniformLocation(program, "color"), color.x, color.y, color.z);

	glUniform3f(glGetUniformLocation(program, "cameraPos"), cameraPos.x, cameraPos.y, cameraPos.z);

	glUniform3f(glGetUniformLocation(program, "sunDir"), sunDir.x, sunDir.y, sunDir.z);
	glUniform3f(glGetUniformLocation(program, "sunColor"), sunColor.x, sunColor.y, sunColor.z);

	glUniform3f(glGetUniformLocation(program, "lightPos"), pointlightPos.x, pointlightPos.y, pointlightPos.z);
	glUniform3f(glGetUniformLocation(program, "lightColor"), pointlightColor.x, pointlightColor.y, pointlightColor.z);

	glUniform3f(glGetUniformLocation(program, "spotlightConeDir"), spotlightConeDir.x, spotlightConeDir.y, spotlightConeDir.z);
	glUniform3f(glGetUniformLocation(program, "spotlightPos"), spotlightPos.x, spotlightPos.y, spotlightPos.z);
	glUniform3f(glGetUniformLocation(program, "spotlightColor"), spotlightColor.x, spotlightColor.y, spotlightColor.z);
	glUniform1f(glGetUniformLocation(program, "spotlightPhi"), spotlightPhi);
	glUniform1i(glGetUniformLocation(program, "shadowsEnabled"), shadowsEnabled);
	

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, depthMap);  // powiązanie tekstury depthMap z jednostką tekstury 0
	//glUniform1i(glGetUniformLocation(program, "shadowMap"), 0);  // Przypisanie tekstury do uniformu w shaderze

	glm::mat4 LightVP = glm::ortho(-1.f, 1.f, -1.f, 1.f, 1.0f, 30.0f) * glm::lookAt(sunPos, sunPos - sunDir, glm::vec3(0, 1, 0));
	glUniformMatrix4fv(glGetUniformLocation(program, "LightVP"), 1, GL_FALSE, (float*)&LightVP);

	
	
	Core::DrawContext(context);

}

void drawObjectTexture(Core::RenderContext& context, glm::mat4 modelMatrix, GLuint textureID, GLuint shaderProgram, GLuint normalmapId) {

	glUseProgram(shaderProgram);
	glm::mat4 viewProjectionMatrix = createPerspectiveMatrix() * createCameraMatrix();
	glm::mat4 transformation = viewProjectionMatrix * modelMatrix;

	glm::mat4 LightVP = glm::ortho(-1.f, 1.f, -1.f, 1.f, 1.0f, 30.0f) * glm::lookAt(sunPos, sunPos - sunDir, glm::vec3(0, 1, 0));


	glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "transformation"), 1, GL_FALSE, (float*)&transformation);
	glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "modelMatrix"), 1, GL_FALSE, (float*)&modelMatrix);
	glUniform3f(glGetUniformLocation(shaderProgram, "lightPos"), spotlightConeDir.x, spotlightConeDir.y, spotlightConeDir.z);
	Core::SetActiveTexture(textureID, "colorTexture", shaderProgram, 0);
	Core::SetActiveTexture(normalmapId, "normalSampler", shaderProgram, 1);
	glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "LightVP"), 1, GL_FALSE, (float*)&LightVP);
	glUniform1i(glGetUniformLocation(shaderProgram, "shadowsEnabled"), shadowsEnabled);
	glUniform1f(glGetUniformLocation(shaderProgram, "exposition"), exposition);
	glUniform3f(glGetUniformLocation(shaderProgram, "cameraPos"), cameraPos.x, cameraPos.y, cameraPos.z);
	glUniform3f(glGetUniformLocation(shaderProgram, "lightColor"), pointlightColor.x, pointlightColor.y, pointlightColor.z);



	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, depthMap);  // powiązanie tekstury depthMap z jednostką tekstury 0
	glUniform1i(glGetUniformLocation(shaderProgram, "shadowMap"), 2);  // Przypisanie tekstury do uniformu w shaderze

	Core::DrawContext(context);
}

void initDepthMap()
{
	glGenFramebuffers(1, &depthMapFBO);

	glGenTextures(1, &depthMap);
	glBindTexture(GL_TEXTURE_2D, depthMap);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT,
		SHADOW_WIDTH, SHADOW_HEIGHT, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthMap, 0);
	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void drawObjectDepth(Core::RenderContext& context, const glm::mat4& viewProjection, const glm::mat4& model)
{
	glUseProgram(programDepth);

	GLuint modelLoc = glGetUniformLocation(programDepth, "modelMatrix");
	GLuint viewProjectionLoc = glGetUniformLocation(programDepth, "viewProjectionMatrix");

	//wysłanie macierzy modelu i widoku-projekcji do GPU
	glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
	glUniformMatrix4fv(viewProjectionLoc, 1, GL_FALSE, glm::value_ptr(viewProjection));

	Core::DrawContext(context);
}

std::vector<Boid> boids;
std::vector<Obstacle> obstacles;

void renderShadowmapSun() {
	float time = glfwGetTime();

	// ustawienie przestrzeni rysowania
	glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT);

	glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);

	glClear(GL_DEPTH_BUFFER_BIT);

	glUseProgram(programDepth);

	glm::mat4 lightVP = glm::ortho(-1.f, 1.f, -1.f, 1.f, 1.0f, 30.0f) * glm::lookAt(sunPos, sunPos - sunDir, glm::vec3(0, 1, 0));

	drawObjectDepth(sphereContext, lightVP, glm::translate(pointlightPos) * glm::scale(glm::vec3(0.1)) * glm::eulerAngleY(time / 3) * glm::translate(glm::vec3(4.f, 0, 0)) * glm::scale(glm::vec3(0.3f)));

	drawObjectDepth(sphereContext, lightVP,
		glm::translate(pointlightPos) * glm::scale(glm::vec3(0.1)) * glm::eulerAngleY(time / 3) * glm::translate(glm::vec3(4.f, 0, 0)) * glm::eulerAngleY(time) * glm::translate(glm::vec3(1.f, 0, 0)) * glm::scale(glm::vec3(0.1f)));

	drawObjectDepth(models::bedContext, lightVP, glm::mat4());
	drawObjectDepth(models::chairContext, lightVP, glm::mat4());
	drawObjectDepth(models::deskContext, lightVP, glm::mat4());
	drawObjectDepth(models::doorContext, lightVP, glm::mat4());
	drawObjectDepth(models::drawerContext, lightVP, glm::mat4());
	drawObjectDepth(models::marbleBustContext, lightVP, glm::mat4());
	drawObjectDepth(models::materaceContext, lightVP, glm::mat4());
	drawObjectDepth(models::pencilsContext, lightVP, glm::mat4());
	drawObjectDepth(models::planeContext, lightVP, glm::mat4());
	drawObjectDepth(models::roomContext, lightVP, glm::mat4());
	drawObjectDepth(models::windowContext, lightVP, glm::mat4());

	glm::vec3 spaceshipSide = glm::normalize(glm::cross(spaceshipDir, glm::vec3(0.f, 1.f, 0.f)));
	glm::vec3 spaceshipUp = glm::normalize(glm::cross(spaceshipSide, spaceshipDir));
	glm::mat4 spaceshipCameraRotationMatrix = glm::mat4({
		spaceshipSide.x, spaceshipSide.y, spaceshipSide.z, 0,
		spaceshipUp.x, spaceshipUp.y, spaceshipUp.z, 0,
		-spaceshipDir.x, -spaceshipDir.y, -spaceshipDir.z, 0,
		0.f, 0.f, 0.f, 1.f
		});

	drawObjectDepth(shipContext, lightVP,
		glm::translate(spaceshipPos) * spaceshipCameraRotationMatrix * glm::eulerAngleY(glm::pi<float>()) * glm::scale(glm::vec3(0.03f))
	);

	for (auto& boid : boids) {
		glm::mat4 model = glm::translate(glm::mat4(1.0f), boid.position);
		model = glm::scale(model, glm::vec3(0.3f));

		// Rysowanie boidów 
		if (boid.groupId == 0) {
			drawObjectDepth(birdContext, lightVP, model);
		}
		else if (boid.groupId == 1) {
			drawObjectDepth(birdContext, lightVP, model);
		}
		else if (boid.groupId == 2) {
			drawObjectDepth(birdContext, lightVP, model);
		}
	}

	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glViewport(0, 0, WIDTH, HEIGHT);
}

std::pair<glm::vec3, glm::vec3> getMinMaxPosition(std::string path, Core::RenderContext& context) {
	return { context.minVertexPos, context.maxVertexPos };
}

void updateBoundingBox(Boid& boid) {
	std::pair<glm::vec3, glm::vec3> minMax = getMinMaxPosition("./models/bird.obj", birdContext);

	float modelSizeX = (minMax.second.x - minMax.first.x) / 2.0f;
	float modelSizeY = (minMax.second.y - minMax.first.y) / 2.0f;
	float modelSizeZ = (minMax.second.z - minMax.first.z) / 2.0f;

	boid.boundingBox = AABB(boid.position.x - modelSizeX,
		boid.position.y - modelSizeY,
		boid.position.z - modelSizeZ,
		boid.position.x + modelSizeX,
		boid.position.y + modelSizeY,
		boid.position.z + modelSizeZ);
}

bool checkCollision(const AABB& box1, const AABB& box2) {
	return !(box1.xMax < box2.xMin || box1.xMin > box2.xMax ||
		box1.yMax < box2.yMin || box1.yMin > box2.yMax ||
		box1.zMax < box2.zMin || box1.zMin > box2.zMax);
}

glm::vec3 chcekColisionBoidsWithObstacles(const Boid& self, const std::vector<Obstacle>& obstacles) {
	glm::vec3 avoidanceVector(0.0f);
	bool collisionDetected = false;

	for (const auto& obstacle : obstacles) {
		if (checkCollision(self.boundingBox, obstacle.boundingBox)) {
			collisionDetected = true;
			glm::vec3 difference = self.position - obstacle.position;
			float distance = glm::length(difference);

			if (distance > 0.0f) {
				avoidanceVector += glm::normalize(difference) / distance;
			}
		}
	}

	return avoidanceVector;
}


void initBoids(int numBoids) {
	std::pair<glm::vec3, glm::vec3> minMax = getMinMaxPosition("./models/bird.obj", birdContext);

	float modelSizeX = (minMax.second.x - minMax.first.x) / 2.0f;
	float modelSizeY = (minMax.second.y - minMax.first.y) / 2.0f;
	float modelSizeZ = (minMax.second.z - minMax.first.z) / 2.0f;


	for (int i = 0; i < numBoids; ++i) {
		Boid newBoid;

		newBoid.position = glm::vec3(rand() % 50 - 25, rand() % 50 - 25, rand() % 50 - 25); // -25 do 24.
		newBoid.velocity = glm::vec3(rand() % 2 - 1, rand() % 2 - 1, rand() % 2 - 1);// od -1 do 1
		newBoid.groupId = i % 3;

		newBoid.boundingBox = AABB(
			newBoid.position.x - modelSizeX,
			newBoid.position.y - modelSizeY,
			newBoid.position.z - modelSizeZ,
			newBoid.position.x + modelSizeX,
			newBoid.position.y + modelSizeY,
			newBoid.position.z + modelSizeZ
		);

		boids.push_back(newBoid);
	}
}

void addObstacle(glm::vec3 position, Core::RenderContext& context, std::string path) {
	std::pair<glm::vec3, glm::vec3> minMax = getMinMaxPosition(path, context);

	float modelSizeX = (minMax.second.x - minMax.first.x) / 2.0f;
	float modelSizeY = (minMax.second.y - minMax.first.y) / 2.0f;
	float modelSizeZ = (minMax.second.z - minMax.first.z) / 2.0f;

	Obstacle newObstacle;
	newObstacle.position = position;

	float offsetX = minMax.first.x + modelSizeX - position.x;
	float offsetY = minMax.first.y + modelSizeY - position.y;
	float offsetZ = minMax.first.z + modelSizeZ - position.z;

	newObstacle.boundingBox = AABB(
		position.x - modelSizeX + offsetX, 
		position.y - modelSizeY + offsetY,
		position.z - modelSizeZ + offsetZ,
		position.x + modelSizeX + offsetX, 
		position.y + modelSizeY + offsetY,
		position.z + modelSizeZ + offsetZ
	);

	obstacles.push_back(newObstacle);
}

glm::vec3 separation(const Boid& self, const std::vector<Boid>& boids, float separationRadius) {
	glm::vec3 avoidVector(0.0f, 0.0f, 0.0f);

	for (const auto& boid : boids) {
		if (&boid == &self || boid.groupId != self.groupId) continue;

		float distance = glm::distance(self.position, boid.position);

		if (distance < separationRadius && distance > 0.0f) {
			glm::vec3 difference = self.position - boid.position;
			avoidVector += glm::normalize(difference) / distance;
		}
	}

	return avoidVector;
}


glm::vec3 align(const Boid& currentBoid, const std::vector<Boid>& boids, float neighborRadius, float maxSpeed) {
	glm::vec3 avgBoidVector(0.0f, 0.0f, 0.0f);
	int count = 0;

	for (const auto& boid : boids) {
		if (&boid == &currentBoid || boid.groupId != currentBoid.groupId) {
			continue;
		}

		float distance = glm::distance(currentBoid.position, boid.position);

		if (distance < neighborRadius) {
			avgBoidVector += boid.velocity;
			count++;
		}
	}

	if (count > 0) {
		avgBoidVector /= static_cast<float>(count);

		if (glm::length(avgBoidVector) > 0.0f) {
			avgBoidVector = glm::normalize(avgBoidVector) * maxSpeed;
		}

		return avgBoidVector - currentBoid.velocity;
	}

	return glm::vec3(0.0f, 0.0f, 0.0f);
}



glm::vec3 cohesion(const Boid& currentBoid, const std::vector<Boid>& boids, float neighborRadius) {
	glm::vec3 centerOfGroup(0.0f, 0.0f, 0.0f);
	int count = 0;

	for (const auto& other : boids) {
		if (&currentBoid == &other || other.groupId != currentBoid.groupId) continue;

		float distance = glm::distance(currentBoid.position, other.position);
		if (distance < neighborRadius) {
			centerOfGroup += other.position;
			count++;
		}
	}

	if (count > 0) {
		centerOfGroup /= static_cast<float>(count);
		return glm::normalize(centerOfGroup - currentBoid.position);
	}

	return glm::vec3(0.0f, 0.0f, 0.0f);
}


void updateBoids() {
	float maxSpeed = 4.0f;
	if (boidsPause) {
		maxSpeed = 0.0f;
	}
	else {
		maxSpeed = 4.0f;
	}
	
	float limit = 4.0f;
	float separationRadius = 1.5f;
	float alignmentRadius = 3.0f;
	float cohesionRadius = 2.0f;

	float separationWeight = 1.5f;
	float alignmentWeight = 1.0f;
	float cohesionWeight = 1.0f;
	float avoidanceWeight = 40.0f;//1.5f;

	for (auto& boid : boids) {
		glm::vec3 separationForce = separation(boid, boids, separationRadius) * separationWeight;
		glm::vec3 alignmentForce = align(boid, boids, alignmentRadius, maxSpeed) * alignmentWeight;
		glm::vec3 cohesionForce = cohesion(boid, boids, cohesionRadius) * cohesionWeight;
		glm::vec3 avoidanceForce = chcekColisionBoidsWithObstacles(boid, obstacles) * avoidanceWeight;

		boid.velocity += (separationForce + alignmentForce + cohesionForce + avoidanceForce) * deltaTime;

		if (glm::length(boid.velocity) > maxSpeed) {
			boid.velocity = glm::normalize(boid.velocity) * maxSpeed;
		}

		boid.position += boid.velocity * deltaTime;

		for (int i = 0; i < 3; ++i) {
			if (boid.position[i] < -limit) {
				boid.position[i] = -limit;
				boid.velocity[i] = -boid.velocity[i];
			}
			if (boid.position[i] > limit) {
				boid.position[i] = limit;
				boid.velocity[i] = -boid.velocity[i];
			}
		}
	}
}




/* terrain generation */

GLuint terrainVAO, terrainVBO;
std::vector<float> terrainVertices;
std::vector<unsigned int> terrainIndices;
GLuint terrainTexture, terrainTexture_map;


// Function to get random noise based on x, z coordinates and seed
float getNoise(int x, int z, int seed) {
	srand(x * 49632 + z * 325176 + seed);
	return (float)(rand()) / (float)(RAND_MAX / 2.0) - 1.0f;
}


// Function to get smooth noise by averaging surrounding values
float getSmoothNoise(int x, int z, int seed) {
	float corners = (getNoise(x - 1, z - 1, seed) + getNoise(x + 1, z - 1, seed) + getNoise(x - 1, z + 1, seed) + getNoise(x + 1, z + 1, seed)) / 16.0f;
	float sides = (getNoise(x - 1, z, seed) + getNoise(x + 1, z, seed) + getNoise(x, z - 1, seed) + getNoise(x, z + 1, seed)) / 8.0f;
	float center = getNoise(x, z, seed) / 4.0f;
	return corners + sides + center;
}

// Function to interpolate between two values using cosine interpolation
float interpolate(float a, float b, float blend) {
	float theta = blend * 3.1415927f;
	float f = (1.0f - cos(theta)) * 0.5f;
	return a * (1.0f - f) + b * f;
}

// Function to get interpolated noise for smooth terrain
float getInterpolatedNoise(float x, float z, int seed, int xOffset, int zOffset) {
	int intX = (int)x;
	int intZ = (int)z;
	float fracX = x - intX;
	float fracZ = z - intZ;

	float v1 = getSmoothNoise(intX, intZ, seed);
	float v2 = getSmoothNoise(intX + 1, intZ, seed);
	float v3 = getSmoothNoise(intX, intZ + 1, seed);
	float v4 = getSmoothNoise(intX + 1, intZ + 1, seed);

	float i1 = interpolate(v1, v2, fracX);
	float i2 = interpolate(v3, v4, fracX);

	return interpolate(i1, i2, fracZ);
}

// Function to generate the height at a specific grid position (x, z)
float generateHeight(int x, int z, int seed, int xOffset, int zOffset) {
	float total = 0.0f;
	float d = std::pow(2, OCTAVES - 1);
	for (int i = 0; i < OCTAVES; ++i) {
		float freq = std::pow(2, i) / d;
		float amp = std::pow(ROUGHNESS, i) * AMPLITUDE;
		total += getInterpolatedNoise((x + xOffset) * freq, (z + zOffset) * freq, seed, xOffset, zOffset) * amp;
	}
	return total;
}

void generateTerrain(int gridSize, float heightScale) {
	terrainVertices.clear();
	terrainIndices.clear();

	float halfSize = gridSize / 2.0f;

	for (int z = 0; z < gridSize; z++) {
		for (int x = 0; x < gridSize; x++) {
			float worldX = x - halfSize; 
			float worldZ = z - halfSize; 
			//float height = sin(x * 0.1f) * cos(z * 0.1f) * heightScale;
			float height = generateHeight(x, z, 1, 0, 0) * heightScale;

			// pos
			terrainVertices.push_back(worldX);
			terrainVertices.push_back(height);
			terrainVertices.push_back(worldZ);

			// normal
			terrainVertices.push_back(0.0f);
			terrainVertices.push_back(1.0f);
			terrainVertices.push_back(0.0f);

			// tex
			float tileScale = 10.0f; 

			terrainVertices.push_back(x * tileScale / (float)gridSize); // U
			terrainVertices.push_back(z * tileScale / (float)gridSize); // V


			// tan, bitan
			terrainVertices.push_back(1.0f); terrainVertices.push_back(0.0f); terrainVertices.push_back(0.0f); // Tangent
			terrainVertices.push_back(0.0f); terrainVertices.push_back(0.0f); terrainVertices.push_back(1.0f); // Bitangent
		}
	}

	for (int z = 0; z < gridSize - 1; z++) {
		for (int x = 0; x < gridSize - 1; x++) {
			int topLeft = z * gridSize + x;
			int topRight = topLeft + 1;
			int bottomLeft = (z + 1) * gridSize + x;
			int bottomRight = bottomLeft + 1;

			terrainIndices.push_back(topLeft);
			terrainIndices.push_back(bottomLeft);
			terrainIndices.push_back(topRight);

			terrainIndices.push_back(topRight);
			terrainIndices.push_back(bottomLeft);
			terrainIndices.push_back(bottomRight);
		}
	}

	glGenVertexArrays(1, &terrainVAO);
	glGenBuffers(1, &terrainVBO);
	GLuint terrainEBO;
	glGenBuffers(1, &terrainEBO);

	glBindVertexArray(terrainVAO);

	glBindBuffer(GL_ARRAY_BUFFER, terrainVBO);
	glBufferData(GL_ARRAY_BUFFER, terrainVertices.size() * sizeof(float), terrainVertices.data(), GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, terrainEBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, terrainIndices.size() * sizeof(unsigned int), terrainIndices.data(), GL_STATIC_DRAW);

	// pos
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 14 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	// normal
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 14 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);

	// tex
	glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 14 * sizeof(float), (void*)(6 * sizeof(float)));
	glEnableVertexAttribArray(2);

	// tan
	glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 14 * sizeof(float), (void*)(8 * sizeof(float)));
	glEnableVertexAttribArray(3);

	// bitan
	glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, 14 * sizeof(float), (void*)(11 * sizeof(float)));
	glEnableVertexAttribArray(4);

	glBindVertexArray(0);
}




void drawTerrain() {
	glUseProgram(programTex);

	glm::mat4 modelMatrix = glm::scale(glm::mat4(1.0f), glm::vec3(0.1f));
	glm::mat4 mvp = createPerspectiveMatrix() * createCameraMatrix() * modelMatrix;

	glUniformMatrix4fv(glGetUniformLocation(programTex, "transformation"), 1, GL_FALSE, &mvp[0][0]);

	Core::SetActiveTexture(terrainTexture, "colorTexture", programTex, 0);
	Core::SetActiveTexture(terrainTexture_map, "normalSampler", programTex, 1); 

	glBindVertexArray(terrainVAO);
	glDrawElements(GL_TRIANGLES, terrainIndices.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
}





void renderScene(GLFWwindow* window)
{
	glClearColor(0.4f, 0.4f, 0.8f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	float time = glfwGetTime();
	updateDeltaTime(time);



	if (shadowsEnabled == true)
	{
		renderShadowmapSun();
	}


	//space lamp
	glUseProgram(programSun);
	glm::mat4 viewProjectionMatrix = createPerspectiveMatrix() * createCameraMatrix();
	glm::mat4 transformation = viewProjectionMatrix * glm::translate(pointlightPos) * glm::scale(glm::vec3(0.1));
	glUniformMatrix4fv(glGetUniformLocation(programSun, "transformation"), 1, GL_FALSE, (float*)&transformation);
	glUniform3f(glGetUniformLocation(programSun, "color"), sunColor.x / 2, sunColor.y / 2, sunColor.z / 2);
	glUniform1f(glGetUniformLocation(programSun, "exposition"), exposition);
	Core::DrawContext(sphereContext);

	glUseProgram(programTex);
	drawTerrain();
	drawObjectTexture(sphereContext, glm::translate(pointlightPos) * glm::scale(glm::vec3(0.1)) * glm::eulerAngleY(time / 3) * glm::translate(glm::vec3(4.f, 0, 0)) * glm::scale(glm::vec3(0.3f)),
		texture::earth, programTex, texture::earth_map);

	drawObjectTexture(sphereContext, glm::translate(pointlightPos) * glm::scale(glm::vec3(0.1)) * glm::eulerAngleY(time / 3) * glm::translate(glm::vec3(4.f, 0, 0)) * glm::eulerAngleY(time) * glm::translate(glm::vec3(1.f, 0, 0)) * glm::scale(glm::vec3(0.1f)),
		texture::moon, programTex, texture::moon_map);

	//drawObjectTexture(models::chairContext, glm::mat4(), texture::chair, programTex, texture::chair_map);
	drawObjectTexture(models::deskContext, glm::mat4(), texture::table, programTex, texture::table_map);

	glUseProgram(program);

	//dwa poniższe to planety wogół lampy
	//drawObjectPBR(sphereContext, glm::translate(pointlightPos) * glm::scale(glm::vec3(0.1)) * glm::eulerAngleY(time / 3) * glm::translate(glm::vec3(4.f, 0, 0)) * glm::scale(glm::vec3(0.3f)), glm::vec3(0.2, 0.7, 0.3), 0.3, 0.0);

	//drawObjectPBR(sphereContext,
	//	glm::translate(pointlightPos) * glm::scale(glm::vec3(0.1)) * glm::eulerAngleY(time / 3) * glm::translate(glm::vec3(4.f, 0, 0)) * glm::eulerAngleY(time) * glm::translate(glm::vec3(1.f, 0, 0)) * glm::scale(glm::vec3(0.1f)),
	//	glm::vec3(0.5, 0.5, 0.5), 0.7, 0.0);




	drawObjectPBR(models::bedContext, glm::mat4(), glm::vec3(0.03f, 0.03f, 0.03f), 0.2f, 0.0f);
	drawObjectPBR(models::chairContext, glm::mat4(), glm::vec3(0.195239f, 0.37728f, 0.8f), 0.4f, 0.0f);
	//drawObjectPBR(models::deskContext, glm::mat4(), glm::vec3(0.428691f, 0.08022f, 0.036889f), 0.2f, 0.0f);
	drawObjectPBR(models::doorContext, glm::mat4(), glm::vec3(0.402978f, 0.120509f, 0.057729f), 0.2f, 0.0f);
	drawObjectPBR(models::drawerContext, glm::mat4(), glm::vec3(0.428691f, 0.08022f, 0.036889f), 0.2f, 0.0f);
	drawObjectPBR(models::marbleBustContext, glm::mat4(), glm::vec3(1.f, 1.f, 1.f), 0.5f, 1.0f);
	drawObjectPBR(models::materaceContext, glm::mat4(), glm::vec3(0.9f, 0.9f, 0.9f), 0.8f, 0.0f);
	drawObjectPBR(models::pencilsContext, glm::mat4(), glm::vec3(0.10039f, 0.018356f, 0.001935f), 0.1f, 0.0f);
	drawObjectPBR(models::planeContext, glm::mat4(), glm::vec3(0.402978f, 0.120509f, 0.057729f), 0.2f, 0.0f);


	glm::vec3 objectPosition(0.0f, 0.0f, 0.0f);

	// Renderowanie obiektu w tej pozycji
	glm::mat4 modelMatrix = glm::translate(glm::mat4(1.0f), objectPosition);
	drawObjectPBR(models::roomContext, modelMatrix, glm::vec3(0.9f, 0.9f, 0.9f), 0.8f, 0.0f);

	//addObstacle(objectPosition, models::roomContext,"./models/room.obj");

	drawObjectPBR(models::windowContext, glm::mat4(), glm::vec3(0.402978f, 0.120509f, 0.057729f), 0.2f, 0.0f);

	

	glm::vec3 spaceshipSide = glm::normalize(glm::cross(spaceshipDir, glm::vec3(0.f, 1.f, 0.f)));
	glm::vec3 spaceshipUp = glm::normalize(glm::cross(spaceshipSide, spaceshipDir));
	glm::mat4 specshipCameraRotrationMatrix = glm::mat4({
		spaceshipSide.x,spaceshipSide.y,spaceshipSide.z,0,
		spaceshipUp.x,spaceshipUp.y,spaceshipUp.z ,0,
		-spaceshipDir.x,-spaceshipDir.y,-spaceshipDir.z,0,
		0.,0.,0.,1.,
		});


	//drawObjectColor(shipContext,
	//	glm::translate(cameraPos + 1.5 * cameraDir + cameraUp * -0.5f) * inveseCameraRotrationMatrix * glm::eulerAngleY(glm::pi<float>()),
	//	glm::vec3(0.3, 0.3, 0.5)
	//	);
	drawObjectPBR(shipContext,
		glm::translate(spaceshipPos) * specshipCameraRotrationMatrix * glm::eulerAngleY(glm::pi<float>()) * glm::scale(glm::vec3(0.03f)),
		glm::vec3(0.3, 0.3, 0.5),
		0.2,1.0
	);

	spotlightPos = spaceshipPos + 0.2 * spaceshipDir;
	spotlightConeDir = spaceshipDir;


	//wyświetlenie shadow map - czy inaczej depth map cieni
	//test depth buffer
	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//glUseProgram(programTest);
	//glActiveTexture(GL_TEXTURE0);
	//glBindTexture(GL_TEXTURE_2D, depthMap);
	//Core::DrawContext(models::testContext);

	for (auto& boid : boids) {
	

		glm::mat4 model = glm::translate(glm::mat4(1.0f), boid.position);
		model = glm::scale(model, glm::vec3(0.3f));
		//model = glm::rotate(model, angleXY, glm::vec3(0.0f, 1.0f, 0.0f));

		// Rysowanie boidów (obiektów 3D)
		if (boid.groupId == 0) {
			drawObjectTexture(birdContext, model, texture::bird1, programTex, texture::bird1_map);
		}
		else if (boid.groupId == 1) {
			drawObjectTexture(birdContext, model, texture::bird2, programTex, texture::bird2_map);
		}
		else if (boid.groupId == 2) {
			drawObjectTexture(birdContext, model, texture::bird3, programTex, texture::bird3_map);
		}

		// Teraz rysowanie AABB
		updateBoundingBox(boid);  
	}

	updateBoids();
	 
	

	glUseProgram(0);
	glfwSwapBuffers(window);
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	aspectRatio = width / float(height);
	glViewport(0, 0, width, height);
	WIDTH = width;
	HEIGHT = height;
}


void loadModelToContext(std::string path, Core::RenderContext& context)
{
	Assimp::Importer import;
	const aiScene* scene = import.ReadFile(path, aiProcess_Triangulate | aiProcess_CalcTangentSpace);

	if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode)
	{
		std::cout << "ERROR::ASSIMP::" << import.GetErrorString() << std::endl;
		return;
	}
	context.initFromAssimpMesh(scene->mMeshes[0]);
}



void init(GLFWwindow* window)
{	
	initDepthMap();

	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	glEnable(GL_DEPTH_TEST);
	program = shaderLoader.CreateProgram("shaders/shader_9_1.vert", "shaders/shader_9_1.frag");
	programTest = shaderLoader.CreateProgram("shaders/test.vert", "shaders/test.frag");
	programSun = shaderLoader.CreateProgram("shaders/shader_8_sun.vert", "shaders/shader_8_sun.frag");
	programDepth = shaderLoader.CreateProgram("shaders/shadowSun.vert", "shaders/shader_8_sun.frag");
	programTex = shaderLoader.CreateProgram("shaders/shader_tex.vert", "shaders/shader_tex.frag");



	loadModelToContext("./models/sphere.obj", sphereContext);
	loadModelToContext("./models/spaceship.obj", shipContext);


	loadModelToContext("./models/bed.obj", models::bedContext);
	loadModelToContext("./models/chair.obj", models::chairContext);
	loadModelToContext("./models/desk.obj", models::deskContext);
	loadModelToContext("./models/door.obj", models::doorContext);
	loadModelToContext("./models/drawer.obj", models::drawerContext);
	loadModelToContext("./models/marbleBust.obj", models::marbleBustContext);
	loadModelToContext("./models/materace.obj", models::materaceContext);
	loadModelToContext("./models/pencils.obj", models::pencilsContext);
	loadModelToContext("./models/plane.obj", models::planeContext);
	loadModelToContext("./models/room.obj", models::roomContext);
	loadModelToContext("./models/spaceship.obj", models::spaceshipContext);
	loadModelToContext("./models/sphere.obj", models::sphereContext);
	loadModelToContext("./models/window.obj", models::windowContext);
	loadModelToContext("./models/test.obj", models::testContext);
	loadModelToContext("./models/bird.obj", birdContext);

	texture::earth = Core::LoadTexture("textures/earth.png");
	texture::earth_map = Core::LoadTexture("textures/earth_normalmap.png");
	texture::moon = Core::LoadTexture("textures/moon.jpg");
	texture::moon_map = Core::LoadTexture("textures/moon_normal.jpg");
	texture::chair = Core::LoadTexture("textures/green_rough_planks_diff_4k.jpg");
	texture::chair_map = Core::LoadTexture("textures/green_rough_planks_nor_gl_4k.png");
	texture::table = Core::LoadTexture("textures/wood_table_001_diff_4k.jpg");
	texture::table_map = Core::LoadTexture("textures/wood_table_001_nor_gl_4k.png");
	texture::bird1 = Core::LoadTexture("textures/bird1.png");
	texture::bird1_map = Core::LoadTexture("textures/bird1_normal.png");
	texture::bird2 = Core::LoadTexture("textures/bird2.png");
	texture::bird2_map = Core::LoadTexture("textures/bird2_normal.png");
	texture::bird3 = Core::LoadTexture("textures/bird3.png");
	texture::bird3_map = Core::LoadTexture("textures/bird3_normal.png");
	
	terrainTexture = Core::LoadTexture("textures/grass1.png");
	terrainTexture_map = Core::LoadTexture("textures/grass1_normal.png");

	generateTerrain(GRID_SIZE, HEIGHT_SCALE); 

	initBoids(100);
	addObstacle(glm::vec3(0.0f, 0.0f, 0.0f), models::roomContext, "./models/room.obj");
	//initObstacles(1);

}

void shutdown(GLFWwindow* window)
{
	shaderLoader.DeleteProgram(program);
}

//obsluga wejscia
void processInput(GLFWwindow* window)
{
	glm::vec3 spaceshipSide = glm::normalize(glm::cross(spaceshipDir, glm::vec3(0.f,1.f,0.f)));
	glm::vec3 spaceshipUp = glm::vec3(0.f, 1.f, 0.f);
	float angleSpeed = 0.05f * deltaTime * 60;
	float moveSpeed = 0.05f * deltaTime * 60;
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
		glfwSetWindowShouldClose(window, true);
	}
	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
		spaceshipPos += spaceshipDir * moveSpeed;
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
		spaceshipPos -= spaceshipDir * moveSpeed;
	if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS)
		spaceshipPos += spaceshipSide * moveSpeed;
	if (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS)
		spaceshipPos -= spaceshipSide * moveSpeed;
	if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
		spaceshipPos += spaceshipUp * moveSpeed;
	if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
		spaceshipPos -= spaceshipUp * moveSpeed;
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
		spaceshipDir = glm::vec3(glm::eulerAngleY(angleSpeed) * glm::vec4(spaceshipDir, 0));
	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
		spaceshipDir = glm::vec3(glm::eulerAngleY(-angleSpeed) * glm::vec4(spaceshipDir, 0));

	if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS)
	{
		shadowsEnabled = true;

	}
	if (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS)
	{
		shadowsEnabled = false;
	}
	if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) {
		boidsPause = true;
	}
	if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS) {
		boidsPause = false;
	}


	cameraPos = spaceshipPos - 0.5 * spaceshipDir + glm::vec3(0, 1, 0) * 0.2f;
	cameraDir = spaceshipDir;

	if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
		exposition -= 0.05;
	if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
		exposition += 0.05;

	if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS) {
		printf("spaceshipPos = glm::vec3(%ff, %ff, %ff);\n", spaceshipPos.x, spaceshipPos.y, spaceshipPos.z);
		printf("spaceshipDir = glm::vec3(%ff, %ff, %ff);\n", spaceshipDir.x, spaceshipDir.y, spaceshipDir.z);
	}

	//cameraDir = glm::normalize(-cameraPos);

}

// funkcja jest glowna petla
void renderLoop(GLFWwindow* window) {
	while (!glfwWindowShouldClose(window))
	{
		processInput(window);

		renderScene(window);
		glfwPollEvents();
	}
}
//}