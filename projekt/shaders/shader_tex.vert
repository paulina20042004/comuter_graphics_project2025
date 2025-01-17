#version 430 core

layout(location = 0) in vec3 vertexPosition;
layout(location = 1) in vec3 vertexNormal;
layout(location = 2) in vec2 vertexTexCoord;
layout(location = 3) in vec3 vertexTangent;
layout(location = 4) in vec3 vertexBitangent;

uniform mat4 transformation;
uniform mat4 modelMatrix;
uniform mat4 LightVP;

//out vec3 vecNormal;
//out vec3 worldPos;

uniform vec3 lightPos;
uniform vec3 spotlightPos;
uniform vec3 cameraPos;
uniform vec3 sunDir;


out vec2 vecTex;


out vec3 viewDirTS;
out vec3 lightDirTS;
out vec3 spotlightDirTS;
out vec3 sunDirTS;

out vec4 sunSpacePos;


void main()
{
	vec3 worldPos = (modelMatrix* vec4(vertexPosition,1)).xyz;
	vec3 vecNormal = (modelMatrix* vec4(vertexNormal,0)).xyz;

	gl_Position = transformation * vec4(vertexPosition, 1.0);
	sunSpacePos=LightVP*modelMatrix*vec4(vertexPosition,1);
	
	vec3 normal = normalize((modelMatrix*vec4(vertexNormal, 0)).xyz);
	vec3 tangent = normalize((modelMatrix*vec4(vertexTangent, 0)).xyz);
	vec3 bitangent = normalize((modelMatrix*vec4(vertexBitangent, 0)).xyz); 

	mat3 TBN = transpose(mat3(tangent, bitangent, normal));
	
	vec3 viewDir = normalize(cameraPos - worldPos);
	vec3 lightDir = normalize(lightPos - worldPos);

	viewDirTS = TBN * viewDir;
	lightDirTS = TBN * lightDir;

	vecTex = vertexTexCoord;
	vecTex.y = 1.0 - vecTex.y;

}
