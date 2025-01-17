#version 430 core

float AMBIENT = 0.1;

uniform vec3 color;

uniform sampler2D colorTexture;
uniform sampler2D normalSampler;
uniform sampler2D shadowMap;
uniform bool shadowsEnabled;
uniform float exposition;


in vec3 vecNormal;
in vec3 worldPos;
in vec2 vecTex;

in vec3 viewDirTS;
in vec3 lightDirTS;
in vec4 sunSpacePos;


out vec4 outColor;

float calculateShadow(vec4 sunSpacePos) {
    
    // perform perspective divide
    vec3 sunSpacePos3 = sunSpacePos.xyz / sunSpacePos.w;
    // transform to [0,1] range
    vec3 sunSpacePosNormalized = sunSpacePos3 * 0.5 + 0.5;
    
    // get closest depth value from light's perspective (using [0,1] range fragPosLight as coords)
    float closestDepth = texture(shadowMap, sunSpacePosNormalized.xy).r; 
    
    // get depth of current fragment from light's perspective
    float currentDepth = sunSpacePosNormalized.z;
    
    // check whether current frag pos is in shadow
    float bias = 0.005;
    float shadow = currentDepth > closestDepth + bias  ? 0.0 : 1.0;

    return shadow;
}


void main()
{	

	vec3 lightDir = normalize(lightDirTS);
	vec3 viewDir = normalize(viewDirTS);
	vec3 normal = vec3(0.0, 0.0, 1.0);
    vec3 N = texture(normalSampler, vecTex).xyz;

	vec3 textureColor = texture2D(colorTexture, vecTex).xyz;

    float shadow = calculateShadow(sunSpacePos);

    N = N * 2.0 - 1.0;
	N = normalize(N);
    float diffuse=max(0,dot(N,lightDir));

    if (shadowsEnabled) {
        diffuse *= shadow;
    } 
    vec3 finalColor = textureColor * min(1.0, AMBIENT + diffuse);
    finalColor *= exposition;
    outColor = vec4(finalColor, 1.0);

}
