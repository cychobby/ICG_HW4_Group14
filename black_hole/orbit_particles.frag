#version 430 core

in vec2  vUV;
in float vT;

out vec4 FragColor;

void main() {
    vec2 p  = vUV * 2.0 - 1.0;
    float r2 = dot(p, p);

    // 邊緣柔化
    float mask = exp(-4.0 * r2);

    float t = clamp(vT, 0.0, 1.0);
    float age = 1.0 - t;

    // 讓尾巴不要直接歸零：至少保留 0.25 的亮度
    float alpha = mask * (0.25 + 0.75 * age);

    vec3 headCol = vec3(1.0, 0.95, 0.2);
    vec3 tailCol = vec3(1.0, 0.25, 0.05);

    // 讓顏色過渡更「慢」
    vec3 col = mix(tailCol, headCol, sqrt(age));

    FragColor = vec4(col, clamp(alpha, 0.0, 1.0));
}
