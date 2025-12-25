#version 430 core

in vec2  vUV;
in float vT;
in float vAgeSec;

out vec4 FragColor;

// Allow external control of particle colors (meteor vs default)
uniform vec3 uHeadColor;
uniform vec3 uTailColor;

void main() {
    vec2 p  = vUV * 2.0 - 1.0;
    float r2 = dot(p, p);

    // 邊緣柔化
    float mask = exp(-2.2 * r2);

    float t = clamp(vT, 0.0, 1.0);
    float age = 1.0 - t;

    // 依時間衰減：k 越大淡得越快
    const float k = 0.5; // 每秒衰減係數，可調
    float timeFade = exp(-k * vAgeSec);

    // 讓尾巴不要直接歸零：至少保留 0.25 的亮度，然後再乘時間衰減
    float alpha = mask * (0.35 + 0.90 * age) * timeFade;

    vec3 headCol = uHeadColor;
    vec3 tailCol = uTailColor;

    // 讓顏色過渡更「慢」
    vec3 col = mix(tailCol, headCol, sqrt(age));

    FragColor = vec4(col, clamp(alpha, 0.0, 1.0));
}
