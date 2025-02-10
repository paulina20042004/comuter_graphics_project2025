#version 330 core
out vec4 FragColor;
in vec3 TexCoords;

uniform samplerCube skybox;

void main() {
    vec3 flippedCoords = TexCoords;

    if (flippedCoords.y > 0.99) {
        flippedCoords.x = -flippedCoords.x;
        flippedCoords.z = -flippedCoords.z;
    }

    FragColor = texture(skybox, flippedCoords);
}
