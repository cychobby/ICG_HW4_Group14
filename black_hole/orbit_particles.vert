#version 430 core
layout(location = 0) in vec3 aPos;

// 直接把 world position 塞進 gl_Position，讓 GS 再去做 viewProj
void main() {
    gl_Position = vec4(aPos, 1.0);
}
