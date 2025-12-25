#version 430 core
layout(location = 0) in vec4 aPosTime; // xyz = 位置, w = spawn time (秒)

out float vSpawnTime;

// 直接把 world position 塞進 gl_Position，讓 GS 再去做 viewProj
void main() {
    gl_Position = vec4(aPosTime.xyz, 1.0);
    vSpawnTime = aPosTime.w;
}
