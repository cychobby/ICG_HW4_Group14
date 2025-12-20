#version 430 core

layout(points) in;
layout(triangle_strip, max_vertices = 4) out;

uniform mat4 viewProj;
uniform float uSizeWorld;     // 粒子在 world space 的大小（會依年齡縮放）
uniform int   uPointCount;    // trail 長度（用來算 fade）
uniform vec3  uCamPos;        // 用來做 billboard

out vec2  vUV;
out float vT;                 // 0(舊) -> 1(新)

void main() {
    vec3 center = gl_in[0].gl_Position.xyz;

    // Billboard basis（面向相機）
    vec3 forward = normalize(uCamPos - center);
    vec3 worldUp = vec3(0.0, 1.0, 0.0);
    vec3 right   = normalize(cross(worldUp, forward));
    vec3 up      = normalize(cross(forward, right));

    float t = 0.0;
    if (uPointCount > 1) {
        t = float(gl_PrimitiveIDIn) / float(uPointCount - 1);
    }
    float age = 1.0 - t;                 // 新的點 age 越大
    float size = uSizeWorld * mix(0.25, 1.0, age);

    vec3 p0 = center + (-right - up) * size;
    vec3 p1 = center + ( right - up) * size;
    vec3 p2 = center + (-right + up) * size;
    vec3 p3 = center + ( right + up) * size;

    vT  = t;

    vUV = vec2(0, 0); gl_Position = viewProj * vec4(p0, 1.0); EmitVertex();
    vUV = vec2(1, 0); gl_Position = viewProj * vec4(p1, 1.0); EmitVertex();
    vUV = vec2(0, 1); gl_Position = viewProj * vec4(p2, 1.0); EmitVertex();
    vUV = vec2(1, 1); gl_Position = viewProj * vec4(p3, 1.0); EmitVertex();
    EndPrimitive();
}
