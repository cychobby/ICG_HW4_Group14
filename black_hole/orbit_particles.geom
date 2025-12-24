#version 430 core

layout(points) in;
// 每顆粒子 1 個 quad = 4 vertices；SAT_COUNT 顆就是 4*SAT_COUNT
layout(triangle_strip, max_vertices = 96) out;

uniform mat4  viewProj;
uniform float uSizeWorld;   // 基礎粒子大小（quad 的半徑，世界座標）
uniform int   uPointCount;  // trail 長度（用來算 t）
uniform vec3  uCamPos;      // billboard 面向相機
uniform float uTime;        // 秒（你 C++ 端每 frame 更新）

// 新增：0 = 單顆（彗星 / meteor），1 = 多顆（軌跡粒子群）
uniform int   uMode;

in float vSpawnTime[];

out vec2  vUV;
out float vT;
out float vAgeSec;

// 你要的「至少 5 顆」：這裡用 10 顆更明顯
const int   SAT_COUNT = 10;

// 粒子本體大小變大
const float SIZE_SCALE = 2.2;

// 粒子群「繞紅光」的半徑（距離紅光中心）變大
const float SUB_ORBIT_RADIUS_FACTOR = 6.0;

// 粒子群繞紅光自轉角速度（rad/s）
const float SUB_ORBIT_OMEGA = 2.5;

// 每顆粒子半徑做一點點差異（看起來更自然）
float radiusJitter(int i) {
    // 0.85 ~ 1.15
    return 0.85 + 0.30 * fract(sin(float(i) * 12.9898) * 43758.5453);
}

void EmitQuad(vec3 center, vec3 camRight, vec3 camUp, float size)
{
    vec3 p0 = center + (-camRight - camUp) * size;
    vec3 p1 = center + ( camRight - camUp) * size;
    vec3 p2 = center + (-camRight + camUp) * size;
    vec3 p3 = center + ( camRight + camUp) * size;

    vUV = vec2(0, 0); gl_Position = viewProj * vec4(p0, 1.0); EmitVertex();
    vUV = vec2(1, 0); gl_Position = viewProj * vec4(p1, 1.0); EmitVertex();
    vUV = vec2(0, 1); gl_Position = viewProj * vec4(p2, 1.0); EmitVertex();
    vUV = vec2(1, 1); gl_Position = viewProj * vec4(p3, 1.0); EmitVertex();
    EndPrimitive();
}

void main()
{
    // base：這個點本身就在「主軌道」上（跟紅光一起繞黑洞）
    vec3 base = gl_in[0].gl_Position.xyz;

    // billboard 的相機座標基底（決定 quad 面向相機）
    vec3 forward = normalize(uCamPos - base);
    vec3 worldUp = vec3(0.0, 1.0, 0.0);
    vec3 camRight = normalize(cross(worldUp, forward));
    vec3 camUp    = normalize(cross(forward, camRight));

    // trail 進度（0 舊 -> 1 新）
    float t = 0.0;
    if (uPointCount > 1) {
        t = float(gl_PrimitiveIDIn) / float(uPointCount - 1);
    }
    vT = t;

    // 讓新的點比較亮/大，舊的點比較淡/小
    float age = 1.0 - t;

    float size = (uSizeWorld * SIZE_SCALE) * mix(0.35, 1.0, age);

    float ageSec = max(0.0, uTime - vSpawnTime[0]);
    vAgeSec = ageSec;

    // === 關鍵：彗星模式只吐 1 顆，避免吃到 SAT_COUNT ===
    if (uMode == 0) {
        EmitQuad(base, camRight, camUp, size);
        return;
    }

    // 粒子群「繞紅光」的半徑（距離紅光中心）
    float R0 = uSizeWorld * SUB_ORBIT_RADIUS_FACTOR;

    // 「繞紅光公轉」：在世界 XZ 平面繞 base 轉
    vec3 planeX = vec3(1.0, 0.0, 0.0);
    vec3 planeZ = vec3(0.0, 0.0, 1.0);

    // 給每顆粒子一個時間相位：確保它們在 base 周圍「自轉」
    float spin = uTime * SUB_ORBIT_OMEGA;

    for (int i = 0; i < SAT_COUNT; ++i) {
        float phase = 6.28318530718 * (float(i) / float(SAT_COUNT));
        float ang   = phase + spin;
        float Ri    = R0 * radiusJitter(i);

        vec3 offset = (cos(ang) * planeX + sin(ang) * planeZ) * Ri;
        EmitQuad(base + offset, camRight, camUp, size);
    }
}
