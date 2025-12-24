#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <chrono>
#include <fstream>
#include <sstream>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace glm;
using namespace std;
using Clock = std::chrono::high_resolution_clock;

// VARS
double lastPrintTime = 0.0;
int    framesCount   = 0;
double c = 299792458.0;
double G = 6.67430e-11;
struct Ray;
bool Gravity = false;
struct Camera;
void spawnMeteor(const Camera& cam, bool fast);
void updateMeteor(double dt);
void deactivateObject(int index);
void clearMeteorObjects();
int allocateDynamicSlot();
bool isDynamicIndex(int index);
double computeTidalRatio(double distance, double mass);
void disruptMeteor(const dvec3& pos, const dvec3& vel, const dvec3& rVec);
bool updateDynamicBody(int index, double dt);

constexpr float  METEOR_RADIUS = 3e9f;
constexpr double METEOR_TIDAL_RADIUS = 1e8;
constexpr double METEOR_MASS = 1e30;
constexpr double METEOR_DESPAWN_RADIUS = 2e12;
constexpr int    METEOR_SUBSTEPS = 128;
constexpr double METEOR_TIME_SCALE = 1500.0;
constexpr double METEOR_SPEED_SLOW = 0.9;
constexpr double METEOR_SPEED_FAST = 1.1;
constexpr double METEOR_INWARD_BIAS = 0.12;
constexpr double METEOR_DRAG_BASE = 0.02;
constexpr double METEOR_DRAG_NEAR = 0.6;
constexpr int    MAX_OBJECTS = 16;
constexpr int    BASE_OBJECTS = 3;
constexpr float  GRID_BASE_OFFSET = -3e10f * BASE_OBJECTS;
constexpr int    METEOR_DEBRIS_COUNT = 4;
constexpr double METEOR_DEBRIS_RADIUS_SCALE = 0.35;
constexpr double METEOR_DEBRIS_SPREAD = 4.0;
constexpr double METEOR_DEBRIS_VEL_RADIAL = 0.08;
constexpr double METEOR_DEBRIS_VEL_TANGENT = 0.05;
constexpr double METEOR_TIDAL_THRESHOLD = 1.0;
constexpr double METEOR_STRETCH_GAIN = 1.2;
constexpr double METEOR_STRETCH_MAX = 3.0;

int meteorIndex = -1;
bool meteorActive = false;
bool meteorDisrupted = false;
vector<int> meteorDebrisIndices;

struct Camera {
    // Center the camera orbit on the black hole at (0, 0, 0)
    vec3 target = vec3(0.0f, 0.0f, 0.0f); // Always look at the black hole center
    float radius = 6.34194e10f;
    float minRadius = 1e10f, maxRadius = 1e12f;

    float azimuth = 0.0f;
    float elevation = M_PI / 2.0f;

    float orbitSpeed = 0.01f;
    float panSpeed = 0.01f;
    double zoomSpeed = 25e9f;

    bool dragging = false;
    bool panning = false;
    bool moving = false; // For compute shader optimization
    double lastX = 0.0, lastY = 0.0;

    // Calculate camera position in world space
    vec3 position() const {
        float clampedElevation = glm::clamp(elevation, 0.01f, float(M_PI) - 0.01f);
        // Orbit around (0,0,0) always
        return vec3(
            radius * sin(clampedElevation) * cos(azimuth),
            radius * cos(clampedElevation),
            radius * sin(clampedElevation) * sin(azimuth)
        );
    }
    void update() {
        // Always keep target at black hole center
        target = vec3(0.0f, 0.0f, 0.0f);
        if(dragging | panning) {
            moving = true;
        } else {
            moving = false;
        }
    }

    void processMouseMove(double x, double y) {
        float dx = float(x - lastX);
        float dy = float(y - lastY);

        if (dragging && panning) {
            // Pan: Shift + Left or Middle Mouse
            // Disable panning to keep camera centered on black hole
        }
        else if (dragging && !panning) {
            // Orbit: Left mouse only
            azimuth   += dx * orbitSpeed;
            elevation -= dy * orbitSpeed;
            elevation = glm::clamp(elevation, 0.01f, float(M_PI) - 0.01f);
        }

        lastX = x;
        lastY = y;
        update();
    }
    void processMouseButton(int button, int action, int mods, GLFWwindow* win) {
        if (button == GLFW_MOUSE_BUTTON_LEFT || button == GLFW_MOUSE_BUTTON_MIDDLE) {
            if (action == GLFW_PRESS) {
                dragging = true;
                // Disable panning so camera always orbits center
                panning = false;
                glfwGetCursorPos(win, &lastX, &lastY);
            } else if (action == GLFW_RELEASE) {
                dragging = false;
                panning = false;
            }
        }
        if (button == GLFW_MOUSE_BUTTON_RIGHT) {
            if (action == GLFW_PRESS) {
                Gravity = true;
            } else if (action == GLFW_RELEASE) {
                Gravity = false;
            }
        }
    }
    void processScroll(double xoffset, double yoffset) {
        radius -= yoffset * zoomSpeed;
        radius = glm::clamp(radius, minRadius, maxRadius);
        update();
    }
    void processKey(int key, int scancode, int action, int mods) {
        if (action == GLFW_PRESS && key == GLFW_KEY_G) {
            Gravity = !Gravity;
            cout << "[INFO] Gravity turned " << (Gravity ? "ON" : "OFF") << endl;
        }
        if (action == GLFW_PRESS && key == GLFW_KEY_M) {
            bool fast = (mods & GLFW_MOD_SHIFT) != 0;
            spawnMeteor(*this, fast);
        }
    }
};
Camera camera;

struct BlackHole {
    vec3 position;
    double mass;
    double radius;
    double r_s;

    BlackHole(vec3 pos, float m) : position(pos), mass(m) {r_s = 2.0 * G * mass / (c*c);}
    bool Intercept(float px, float py, float pz) const {
        double dx = double(px) - double(position.x);
        double dy = double(py) - double(position.y);
        double dz = double(pz) - double(position.z);
        double dist2 = dx * dx + dy * dy + dz * dz;
        return dist2 < r_s * r_s;
    }
};
BlackHole SagA(vec3(0.0f, 0.0f, 0.0f), 8.54e36); // Sagittarius A black hole
struct ObjectData {
    vec4 posRadius; // xyz = position, w = radius
    vec4 color;     // rgb = color, a = unused
    float  mass;
    vec3 velocity = vec3(0.0f, 0.0f, 0.0f); // Initial velocity
};
vector<ObjectData> objects = {
    { vec4(4e11f, 0.0f, 0.0f, 7e10f)   , vec4(1,1,0,1), 1.98892e30 },
    { vec4(0.0f, 0.0f, 4e11f, 4e10f)   , vec4(1,0,0,1), 1.98892e30 },
    { vec4(0.0f, 0.0f, 0.0f, SagA.r_s) , vec4(0,0,0,1), static_cast<float>(SagA.mass)  },
    //{ vec4(6e10f, 0.0f, 0.0f, 5e10f), vec4(0,1,0,1) }
};
void deactivateObject(int index) {
    if (index < 0 || index >= static_cast<int>(objects.size())) {
        return;
    }
    ObjectData& obj = objects[index];
    obj.mass = 0.0f;
    obj.posRadius.w = 0.0f;
    obj.color.w = 0.0f;
    obj.velocity = vec3(0.0f);
}
void clearMeteorObjects() {
    if (meteorIndex >= 0) {
        deactivateObject(meteorIndex);
    }
    for (int idx : meteorDebrisIndices) {
        if (idx != meteorIndex) {
            deactivateObject(idx);
        }
    }
    meteorDebrisIndices.clear();
    meteorActive = false;
    meteorDisrupted = false;
}
int allocateDynamicSlot() {
    for (int i = BASE_OBJECTS; i < static_cast<int>(objects.size()); ++i) {
        if (objects[i].mass == 0.0f && objects[i].posRadius.w == 0.0f) {
            return i;
        }
    }
    if (static_cast<int>(objects.size()) >= MAX_OBJECTS) {
        return -1;
    }
    objects.push_back({ vec4(0.0f), vec4(0.0f), 0.0f, vec3(0.0f) });
    return static_cast<int>(objects.size() - 1);
}
bool isDynamicIndex(int index) {
    if (index == meteorIndex) {
        return true;
    }
    for (int idx : meteorDebrisIndices) {
        if (idx == index) {
            return true;
        }
    }
    return false;
}
double computeTidalRatio(double distance, double mass) {
    if (distance <= 0.0 || mass <= 0.0) {
        return 0.0;
    }
    double r = METEOR_TIDAL_RADIUS;
    double r3 = r * r * r;
    double d3 = distance * distance * distance;
    return (2.0 * SagA.mass * r3) / (mass * d3);
}
void spawnMeteor(const Camera& cam, bool fast) {
    vec3 camPos = cam.position();
    double camR = glm::length(camPos);
    if (camR <= 0.0) {
        return;
    }

    clearMeteorObjects();
    if (meteorIndex < BASE_OBJECTS || meteorIndex >= static_cast<int>(objects.size())) {
        meteorIndex = -1;
    }
    if (meteorIndex < 0) {
        meteorIndex = allocateDynamicSlot();
        if (meteorIndex < 0) {
            cout << "[WARN] No free slot for meteor\n";
            return;
        }
    }

    double spawnR = camR * 0.8;
    if (spawnR < SagA.r_s * 3.0) {
        spawnR = SagA.r_s * 3.0;
    }
    vec3 radial = normalize(camPos);
    vec3 forward = normalize(cam.target - camPos);
    vec3 up = vec3(0.0f, 1.0f, 0.0f);
    vec3 right = normalize(cross(forward, up));
    if (glm::length(right) < 1e-5f) {
        right = vec3(1.0f, 0.0f, 0.0f);
    }
    vec3 spawnPos = radial * float(spawnR) + right * float(camR * 0.05);

    dvec3 rVec = dvec3(spawnPos) - dvec3(SagA.position);
    double rLen = glm::length(rVec);
    if (rLen <= SagA.r_s) {
        rVec = dvec3(radial) * (SagA.r_s * 3.0);
        spawnPos = vec3(rVec);
        rLen = glm::length(rVec);
    }
    dvec3 rHat = rVec / rLen;
    dvec3 tangent = cross(dvec3(0.0, 1.0, 0.0), rHat);
    double tLen = glm::length(tangent);
    if (tLen < 1e-6) {
        tangent = cross(dvec3(1.0, 0.0, 0.0), rHat);
        tLen = glm::length(tangent);
    }
    tangent /= tLen;

    double vCirc = sqrt(G * SagA.mass / rLen);
    double speedFactor = fast ? METEOR_SPEED_FAST : METEOR_SPEED_SLOW;
    dvec3 v = tangent * (vCirc * speedFactor);
    if (!fast) {
        v -= rHat * (vCirc * METEOR_INWARD_BIAS);
    }

    ObjectData meteor;
    meteor.posRadius = vec4(vec3(spawnPos), METEOR_RADIUS);
    meteor.color = vec4(0.9f, 0.9f, 1.0f, 1.0f);
    meteor.mass = static_cast<float>(METEOR_MASS);
    meteor.velocity = vec3(v);

    objects[meteorIndex] = meteor;
    meteorActive = true;
    meteorDisrupted = false;
    cout << "[INFO] Meteor spawned (" << (fast ? "fast" : "slow") << ")\n";
}
void disruptMeteor(const dvec3& pos, const dvec3& vel, const dvec3& rVec) {
    if (meteorIndex < 0) {
        return;
    }
    double totalMass = objects[meteorIndex].mass;
    if (totalMass <= 0.0) {
        return;
    }

    meteorDisrupted = true;
    meteorDebrisIndices.clear();

    dvec3 rHat = normalize(rVec);
    dvec3 tHat = cross(dvec3(0.0, 1.0, 0.0), rHat);
    if (glm::length(tHat) < 1e-6) {
        tHat = cross(dvec3(1.0, 0.0, 0.0), rHat);
    }
    tHat = normalize(tHat);
    dvec3 bHat = normalize(cross(rHat, tHat));

    double rLen = glm::length(rVec);
    double vCirc = sqrt(G * SagA.mass / rLen);
    double spread = METEOR_DEBRIS_SPREAD * METEOR_RADIUS;
    double velRad = vCirc * METEOR_DEBRIS_VEL_RADIAL;
    double velTan = vCirc * METEOR_DEBRIS_VEL_TANGENT;
    int count = METEOR_DEBRIS_COUNT;
    if (count < 1) {
        count = 1;
    }
    double fragMass = totalMass / double(count);

    for (int i = 0; i < count; ++i) {
        double t = (count == 1) ? 0.0 : (double(i) / double(count - 1) - 0.5);
        dvec3 offset = rHat * (t * spread) + bHat * (t * spread * 0.2);
        dvec3 fragPos = pos + offset;
        dvec3 fragVel = vel + rHat * (t * velRad) + tHat * (t * velTan);

        ObjectData frag;
        frag.posRadius = vec4(vec3(fragPos), float(METEOR_RADIUS * METEOR_DEBRIS_RADIUS_SCALE));
        frag.color = vec4(1.0f, 0.85f, 0.6f, 1.0f);
        frag.mass = static_cast<float>(fragMass);
        frag.velocity = vec3(fragVel);

        int slot = (i == 0) ? meteorIndex : allocateDynamicSlot();
        if (slot < 0) {
            continue;
        }
        objects[slot] = frag;
        meteorDebrisIndices.push_back(slot);
    }
    meteorActive = true;
}
bool updateDynamicBody(int index, double dt) {
    if (index < 0 || index >= static_cast<int>(objects.size())) {
        return false;
    }
    ObjectData& body = objects[index];
    if (body.posRadius.w <= 0.0f || body.mass <= 0.0f) {
        return false;
    }

    double scaledDt = dt * METEOR_TIME_SCALE;
    if (scaledDt <= 0.0) {
        return true;
    }

    dvec3 center = dvec3(SagA.position);
    dvec3 pos = dvec3(vec3(body.posRadius));
    dvec3 vel = dvec3(body.velocity);

    double h = scaledDt / double(METEOR_SUBSTEPS);
    for (int step = 0; step < METEOR_SUBSTEPS; ++step) {
        dvec3 rVec = pos - center;
        double rLen = glm::length(rVec);
        if (rLen <= SagA.r_s) {
            deactivateObject(index);
            return false;
        }

        double invR = 1.0 / rLen;
        dvec3 grav = -(G * SagA.mass) * (rVec * (invR * invR * invR));
        double rScale = SagA.r_s * invR;
        double drag = METEOR_DRAG_BASE + METEOR_DRAG_NEAR * (rScale * rScale);
        double damp = 1.0 - drag * h;
        if (damp < 0.0) {
            damp = 0.0;
        }

        vel = vel * damp + grav * h;
        pos += vel * h;

        double nextR = glm::length(pos - center);
        if (nextR > METEOR_DESPAWN_RADIUS) {
            deactivateObject(index);
            return false;
        }
    }

    body.velocity = vec3(vel);
    body.posRadius = vec4(vec3(pos), body.posRadius.w);
    return true;
}
void updateMeteor(double dt) {
    if (!meteorActive || meteorIndex < 0) {
        return;
    }

    if (meteorDisrupted) {
        bool anyActive = false;
        for (int idx : meteorDebrisIndices) {
            if (updateDynamicBody(idx, dt)) {
                anyActive = true;
            }
        }
        if (!anyActive) {
            clearMeteorObjects();
        }
        return;
    }

    ObjectData& meteor = objects[meteorIndex];
    if (meteor.posRadius.w <= 0.0f || meteor.mass <= 0.0f) {
        meteorActive = false;
        return;
    }

    double scaledDt = dt * METEOR_TIME_SCALE;
    if (scaledDt <= 0.0) {
        return;
    }

    dvec3 center = dvec3(SagA.position);
    dvec3 pos = dvec3(vec3(meteor.posRadius));
    dvec3 vel = dvec3(meteor.velocity);

    double h = scaledDt / double(METEOR_SUBSTEPS);
    double stretchScale = 1.0;
    for (int step = 0; step < METEOR_SUBSTEPS; ++step) {
        dvec3 rVec = pos - center;
        double rLen = glm::length(rVec);
        if (rLen <= SagA.r_s) {
            clearMeteorObjects();
            return;
        }

        double ratio = computeTidalRatio(rLen, meteor.mass);
        stretchScale = 1.0 + METEOR_STRETCH_GAIN * ratio;
        if (stretchScale > METEOR_STRETCH_MAX) {
            stretchScale = METEOR_STRETCH_MAX;
        }
        if (ratio >= METEOR_TIDAL_THRESHOLD) {
            disruptMeteor(pos, vel, rVec);
            return;
        }

        double invR = 1.0 / rLen;
        dvec3 grav = -(G * SagA.mass) * (rVec * (invR * invR * invR));
        double rScale = SagA.r_s * invR;
        double drag = METEOR_DRAG_BASE + METEOR_DRAG_NEAR * (rScale * rScale);
        double damp = 1.0 - drag * h;
        if (damp < 0.0) {
            damp = 0.0;
        }

        vel = vel * damp + grav * h;
        pos += vel * h;

        double nextR = glm::length(pos - center);
        if (nextR > METEOR_DESPAWN_RADIUS) {
            clearMeteorObjects();
            return;
        }
    }

    meteor.velocity = vec3(vel);
    meteor.posRadius = vec4(vec3(pos), float(METEOR_RADIUS * stretchScale));
}

struct Engine {
    GLuint gridShaderProgram;
    // -- Quad & Texture render -- //
    GLFWwindow* window;
    GLuint quadVAO;
    GLuint texture;
    GLuint shaderProgram;
    GLuint computeProgram = 0;

    // -- HW4 Scheme A: orbit trail rendered by Geometry Shader -- //
    GLuint orbitShaderProgram = 0;
    GLuint orbitVAO = 0;
    GLuint orbitVBO = 0;
    int   orbitPointCount = 512;
    float orbitParticleSize = 2e9f;   // 可自行調整大小
    vector<vec4> orbitTrail;          // xyz = pos, w = spawn time (seconds)

    // -- UBOs -- //
    GLuint cameraUBO = 0;
    GLuint diskUBO = 0;
    GLuint objectsUBO = 0;
    // -- grid mess vars -- //
    GLuint gridVAO = 0;
    GLuint gridVBO = 0;
    GLuint gridEBO = 0;
    int gridIndexCount = 0;

    int WIDTH = 800;  // Window width
    int HEIGHT = 600; // Window height
    int COMPUTE_WIDTH  = 200;   // Compute resolution width
    int COMPUTE_HEIGHT = 150;  // Compute resolution height
    float width = 100000000000.0f; // Width of the viewport in meters
    float height = 75000000000.0f; // Height of the viewport in meters
    
    Engine() {
        if (!glfwInit()) {
            cerr << "GLFW init failed\n";
            exit(EXIT_FAILURE);
        }
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        window = glfwCreateWindow(WIDTH, HEIGHT, "Black Hole", nullptr, nullptr);
        if (!window) {
            cerr << "Failed to create GLFW window\n";
            glfwTerminate();
            exit(EXIT_FAILURE);
        }
        glfwMakeContextCurrent(window);
        glewExperimental = GL_TRUE;
        GLenum glewErr = glewInit();
        if (glewErr != GLEW_OK) {
            cerr << "Failed to initialize GLEW: "
                << (const char*)glewGetErrorString(glewErr)
                << "\n";
            glfwTerminate();
            exit(EXIT_FAILURE);
        }
        cout << "OpenGL " << glGetString(GL_VERSION) << "\n";
        this->shaderProgram = CreateShaderProgram();
        gridShaderProgram = CreateShaderProgram("grid.vert", "grid.frag");

        // (新增) orbit trail shader: VS + GS + FS
        orbitShaderProgram = CreateShaderProgram(
            "orbit_particles.vert",
            "orbit_particles.geom",
            "orbit_particles.frag"
        );
        
        computeProgram = CreateComputeProgram("geodesic.comp");
        glGenBuffers(1, &cameraUBO);
        glBindBuffer(GL_UNIFORM_BUFFER, cameraUBO);
        glBufferData(GL_UNIFORM_BUFFER, 128, nullptr, GL_DYNAMIC_DRAW); // alloc ~128 bytes
        glBindBufferBase(GL_UNIFORM_BUFFER, 1, cameraUBO); // binding = 1 matches shader

        glGenBuffers(1, &diskUBO);
        glBindBuffer(GL_UNIFORM_BUFFER, diskUBO);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(float) * 4, nullptr, GL_DYNAMIC_DRAW); // 3 values + 1 padding
        glBindBufferBase(GL_UNIFORM_BUFFER, 2, diskUBO); // binding = 2 matches compute shader

        glGenBuffers(1, &objectsUBO);
        glBindBuffer(GL_UNIFORM_BUFFER, objectsUBO);
        // allocate space for 16 objects: 
        // sizeof(int) + padding + 16×(vec4 posRadius + vec4 color)
        GLsizeiptr objUBOSize = sizeof(int) + 3 * sizeof(float)
            + 16 * (sizeof(vec4) + sizeof(vec4))
            + 16 * sizeof(float); // 16 floats for mass
        glBufferData(GL_UNIFORM_BUFFER, objUBOSize, nullptr, GL_DYNAMIC_DRAW);
        glBindBufferBase(GL_UNIFORM_BUFFER, 3, objectsUBO);  // binding = 3 matches shader

        auto result = QuadVAO();
        this->quadVAO = result[0];
        this->texture = result[1];

        // (新增) orbit trail VBO/VAO
        orbitTrail.assign(orbitPointCount, vec4(vec3(objects[1].posRadius), 0.0f)); // 先全部填同一點，避免第一幀拉一條超長線
        glGenVertexArrays(1, &orbitVAO);
        glGenBuffers(1, &orbitVBO);

        glBindVertexArray(orbitVAO);
        glBindBuffer(GL_ARRAY_BUFFER, orbitVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(vec4) * orbitTrail.size(), orbitTrail.data(), GL_DYNAMIC_DRAW);
        glEnableVertexAttribArray(0); // layout(location=0) in vec4 aPosTime (xyz,pos; w,time)
        glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(vec4), (void*)0);
        glBindVertexArray(0);
        
    }
    void generateGrid(const vector<ObjectData>& objects) {
        const int gridSize = 25;
        const float spacing = 1e10f;  // tweak this

        vector<vec3> vertices;
        vector<GLuint> indices;

        for (int z = 0; z <= gridSize; ++z) {
            for (int x = 0; x <= gridSize; ++x) {
                float worldX = (x - gridSize / 2) * spacing;
                float worldZ = (z - gridSize / 2) * spacing;

                float y = GRID_BASE_OFFSET;

                // ✅ Warp grid using Schwarzschild geometry
                for (size_t objIndex = 0; objIndex < objects.size(); ++objIndex) {
                    if (isDynamicIndex(static_cast<int>(objIndex))) {
                        continue;
                    }
                    const auto& obj = objects[objIndex];
                    if (obj.mass <= 0.0f) {
                        continue;
                    }
                    vec3 objPos = vec3(obj.posRadius);
                    double mass = obj.mass;

                    double r_s = 2.0 * G * mass / (c * c);
                    double dx = worldX - objPos.x;
                    double dz = worldZ - objPos.z;
                    double dist = sqrt(dx * dx + dz * dz);

                    // prevent sqrt of negative or divide-by-zero (inside or at the black hole center)
                    if (dist > r_s) {
                        double deltaY = 2.0 * sqrt(r_s * (dist - r_s));
                        y += static_cast<float>(deltaY);
                    } else {
                        // ?? For points inside or at r_s: make it dip down sharply
                        y += 2.0f * static_cast<float>(sqrt(r_s * r_s));  // or add a deep pit
                    }
                }

                vertices.emplace_back(worldX, y, worldZ);
            }
        }

        // 🧩 Add indices for GL_LINE rendering
        for (int z = 0; z < gridSize; ++z) {
            for (int x = 0; x < gridSize; ++x) {
                int i = z * (gridSize + 1) + x;
                indices.push_back(i);
                indices.push_back(i + 1);

                indices.push_back(i);
                indices.push_back(i + gridSize + 1);
            }
        }

        // 🔌 Upload to GPU
        if (gridVAO == 0) glGenVertexArrays(1, &gridVAO);
        if (gridVBO == 0) glGenBuffers(1, &gridVBO);
        if (gridEBO == 0) glGenBuffers(1, &gridEBO);

        glBindVertexArray(gridVAO);

        glBindBuffer(GL_ARRAY_BUFFER, gridVBO);
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(vec3), vertices.data(), GL_DYNAMIC_DRAW);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gridEBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), indices.data(), GL_STATIC_DRAW);

        glEnableVertexAttribArray(0); // location = 0
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), (void*)0);

        gridIndexCount = indices.size();

        glBindVertexArray(0);
    }
    void drawGrid(const mat4& viewProj) {
        glUseProgram(gridShaderProgram);
        glUniformMatrix4fv(glGetUniformLocation(gridShaderProgram, "viewProj"),
                        1, GL_FALSE, glm::value_ptr(viewProj));
        glBindVertexArray(gridVAO);

        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        glDrawElements(GL_LINES, gridIndexCount, GL_UNSIGNED_INT, 0);

        glBindVertexArray(0);
        glEnable(GL_DEPTH_TEST);
    }
    void drawFullScreenQuad() {
        glUseProgram(shaderProgram); // fragment + vertex shader
        glBindVertexArray(quadVAO);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texture);
        glUniform1i(glGetUniformLocation(shaderProgram, "screenTexture"), 0);

        glDisable(GL_DEPTH_TEST);  // draw as background
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 6);  // 2 triangles
        glEnable(GL_DEPTH_TEST);
    }
    void pushOrbitTrailSample(const vec3& p, float nowSeconds) {
        // 你的 geom shader 用 age = 1 - t，而 t 由 gl_PrimitiveIDIn 來
        // 所以我們讓「最新的點」放在 index 0，越舊的點往後推
        for (int i = orbitPointCount - 1; i > 0; --i) {
            orbitTrail[i] = orbitTrail[i - 1];
        }
        orbitTrail[0] = vec4(p, nowSeconds);

        glBindBuffer(GL_ARRAY_BUFFER, orbitVBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vec4) * orbitTrail.size(), orbitTrail.data());
    }

    void drawOrbitTrail(const mat4& viewProj, const vec3& camPos, float nowSeconds) {
        glUseProgram(orbitShaderProgram);

        glUniformMatrix4fv(glGetUniformLocation(orbitShaderProgram, "viewProj"),
                        1, GL_FALSE, value_ptr(viewProj));
        glUniform1f(glGetUniformLocation(orbitShaderProgram, "uSizeWorld"), orbitParticleSize);
        glUniform1i(glGetUniformLocation(orbitShaderProgram, "uPointCount"), orbitPointCount);
        glUniform3fv(glGetUniformLocation(orbitShaderProgram, "uCamPos"), 1, value_ptr(camPos));
        glUniform1f(glGetUniformLocation(orbitShaderProgram, "uTime"), nowSeconds);

        glBindVertexArray(orbitVAO);

        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE); // 加色發光

        glDrawArrays(GL_POINTS, 0, orbitPointCount);

        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glBindVertexArray(0);
        glEnable(GL_DEPTH_TEST);
    }
    GLuint CreateShaderProgram(){
        const char* vertexShaderSource = R"(
        #version 330 core
        layout (location = 0) in vec2 aPos;  // Changed to vec2
        layout (location = 1) in vec2 aTexCoord;
        out vec2 TexCoord;
        void main() {
            gl_Position = vec4(aPos, 0.0, 1.0);  // Explicit z=0
            TexCoord = aTexCoord;
        })";

        const char* fragmentShaderSource = R"(
        #version 330 core
        in vec2 TexCoord;
        out vec4 FragColor;
        uniform sampler2D screenTexture;
        void main() {
            FragColor = texture(screenTexture, TexCoord);
        })";

        // vertex shader
        GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertexShader, 1, &vertexShaderSource, nullptr);
        glCompileShader(vertexShader);

        // fragment shader
        GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragmentShader, 1, &fragmentShaderSource, nullptr);
        glCompileShader(fragmentShader);

        GLuint shaderProgram = glCreateProgram();
        glAttachShader(shaderProgram, vertexShader);
        glAttachShader(shaderProgram, fragmentShader);
        glLinkProgram(shaderProgram);

        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);

        return shaderProgram;
    };
    GLuint CreateShaderProgram(const char* vertPath, const char* fragPath) {
        auto loadShader = [](const char* path, GLenum type) -> GLuint {
            std::ifstream in(path);
            if (!in.is_open()) {
                std::cerr << "Failed to open shader: " << path << "\n";
                exit(EXIT_FAILURE);
            }
            std::stringstream ss;
            ss << in.rdbuf();
            std::string srcStr = ss.str();
            const char* src = srcStr.c_str();

            GLuint shader = glCreateShader(type);
            glShaderSource(shader, 1, &src, nullptr);
            glCompileShader(shader);

            GLint success;
            glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
            if (!success) {
                GLint logLen;
                glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &logLen);
                std::vector<char> log(logLen);
                glGetShaderInfoLog(shader, logLen, nullptr, log.data());
                std::cerr << "Shader compile error (" << path << "):\n" << log.data() << "\n";
                exit(EXIT_FAILURE);
            }
            return shader;
        };

        GLuint vertShader = loadShader(vertPath, GL_VERTEX_SHADER);
        GLuint fragShader = loadShader(fragPath, GL_FRAGMENT_SHADER);

        GLuint program = glCreateProgram();
        glAttachShader(program, vertShader);
        glAttachShader(program, fragShader);
        glLinkProgram(program);

        GLint linkSuccess;
        glGetProgramiv(program, GL_LINK_STATUS, &linkSuccess);
        if (!linkSuccess) {
            GLint logLen;
            glGetProgramiv(program, GL_INFO_LOG_LENGTH, &logLen);
            std::vector<char> log(logLen);
            glGetProgramInfoLog(program, logLen, nullptr, log.data());
            std::cerr << "Shader link error:\n" << log.data() << "\n";
            exit(EXIT_FAILURE);
        }

        glDeleteShader(vertShader);
        glDeleteShader(fragShader);

        return program;
    }
    GLuint CreateShaderProgram(const char* vertPath, const char* geomPath, const char* fragPath) {
        auto loadShader = [](const char* path, GLenum type) -> GLuint {
            std::ifstream in(path);
            if (!in.is_open()) {
                std::cerr << "Failed to open shader: " << path << "\n";
                exit(EXIT_FAILURE);
            }
            std::stringstream ss;
            ss << in.rdbuf();
            std::string srcStr = ss.str();
            const char* src = srcStr.c_str();

            GLuint shader = glCreateShader(type);
            glShaderSource(shader, 1, &src, nullptr);
            glCompileShader(shader);

            GLint success;
            glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
            if (!success) {
                GLint logLen;
                glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &logLen);
                std::vector<char> log(logLen);
                glGetShaderInfoLog(shader, logLen, nullptr, log.data());
                std::cerr << "Shader compile error (" << path << "):\n" << log.data() << "\n";
                exit(EXIT_FAILURE);
            }
            return shader;
        };

        GLuint vertShader = loadShader(vertPath, GL_VERTEX_SHADER);
        GLuint geomShader = loadShader(geomPath, GL_GEOMETRY_SHADER);
        GLuint fragShader = loadShader(fragPath, GL_FRAGMENT_SHADER);

        GLuint program = glCreateProgram();
        glAttachShader(program, vertShader);
        glAttachShader(program, geomShader);
        glAttachShader(program, fragShader);
        glLinkProgram(program);

        GLint linkSuccess;
        glGetProgramiv(program, GL_LINK_STATUS, &linkSuccess);
        if (!linkSuccess) {
            GLint logLen;
            glGetProgramiv(program, GL_INFO_LOG_LENGTH, &logLen);
            std::vector<char> log(logLen);
            glGetProgramInfoLog(program, logLen, nullptr, log.data());
            std::cerr << "Shader link error:\n" << log.data() << "\n";
            exit(EXIT_FAILURE);
        }

        glDeleteShader(vertShader);
        glDeleteShader(geomShader);
        glDeleteShader(fragShader);

        return program;
    }
    GLuint CreateComputeProgram(const char* path) {
        // 1) read GLSL source
        std::ifstream in(path);
        if(!in.is_open()) {
            std::cerr << "Failed to open compute shader: " << path << "\n";
            exit(EXIT_FAILURE);
        }
        std::stringstream ss;
        ss << in.rdbuf();
        std::string srcStr = ss.str();
        const char* src = srcStr.c_str();

        // 2) compile
        GLuint cs = glCreateShader(GL_COMPUTE_SHADER);
        glShaderSource(cs, 1, &src, nullptr);
        glCompileShader(cs);
        GLint ok; 
        glGetShaderiv(cs, GL_COMPILE_STATUS, &ok);
        if(!ok) {
            GLint logLen;
            glGetShaderiv(cs, GL_INFO_LOG_LENGTH, &logLen);
            std::vector<char> log(logLen);
            glGetShaderInfoLog(cs, logLen, nullptr, log.data());
            std::cerr << "Compute shader compile error:\n" << log.data() << "\n";
            exit(EXIT_FAILURE);
        }

        // 3) link
        GLuint prog = glCreateProgram();
        glAttachShader(prog, cs);
        glLinkProgram(prog);
        glGetProgramiv(prog, GL_LINK_STATUS, &ok);
        if(!ok) {
            GLint logLen;
            glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &logLen);
            std::vector<char> log(logLen);
            glGetProgramInfoLog(prog, logLen, nullptr, log.data());
            std::cerr << "Compute shader link error:\n" << log.data() << "\n";
            exit(EXIT_FAILURE);
        }

        glDeleteShader(cs);
        return prog;
    }
    void dispatchCompute(const Camera& cam) {
        // determine target compute‐res
        int cw = cam.moving ? COMPUTE_WIDTH  : 200;
        int ch = cam.moving ? COMPUTE_HEIGHT : 150;

        // 1) reallocate the texture if needed
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexImage2D(GL_TEXTURE_2D,
                    0,                // mip
                    GL_RGBA8,         // internal format
                    cw,               // width
                    ch,               // height
                    0, GL_RGBA, 
                    GL_UNSIGNED_BYTE, 
                    nullptr);

        // 2) bind compute program & UBOs
        glUseProgram(computeProgram);
        uploadCameraUBO(cam);
        uploadDiskUBO();
        uploadObjectsUBO(objects);

        // 3) bind it as image unit 0
        glBindImageTexture(0, texture, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA8);

        // 4) dispatch grid
        GLuint groupsX = (GLuint)std::ceil(cw / 16.0f);
        GLuint groupsY = (GLuint)std::ceil(ch / 16.0f);
        glDispatchCompute(groupsX, groupsY, 1);

        // 5) sync
        glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
    }
    void uploadCameraUBO(const Camera& cam) {
        struct UBOData {
            vec3 pos; float _pad0;
            vec3 right; float _pad1;
            vec3 up; float _pad2;
            vec3 forward; float _pad3;
            float tanHalfFov;
            float aspect;
            bool moving;
            int _pad4;
        } data;
        vec3 fwd = normalize(cam.target - cam.position());
        vec3 up = vec3(0, 1, 0); // y axis is up, so disk is in x-z plane
        vec3 right = normalize(cross(fwd, up));
        up = cross(right, fwd);

        data.pos = cam.position();
        data.right = right;
        data.up = up;
        data.forward = fwd;
        data.tanHalfFov = tan(radians(60.0f * 0.5f));
        data.aspect = float(WIDTH) / float(HEIGHT);
        data.moving = cam.dragging || cam.panning;

        glBindBuffer(GL_UNIFORM_BUFFER, cameraUBO);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(UBOData), &data);
    }
    void uploadObjectsUBO(const vector<ObjectData>& objs) {
        struct UBOData {
            int   numObjects;
            float _pad0, _pad1, _pad2;        // <-- pad out to 16 bytes
            vec4  posRadius[16];
            vec4  color[16];
            float  mass[16]; 
        } data;

        size_t count = std::min(objs.size(), size_t(16));
        data.numObjects = static_cast<int>(count);

        for (size_t i = 0; i < count; ++i) {
            data.posRadius[i] = objs[i].posRadius;
            data.color[i] = objs[i].color;
            data.mass[i] = objs[i].mass;
        }

        // Upload
        glBindBuffer(GL_UNIFORM_BUFFER, objectsUBO);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(data), &data);
    }
    void uploadDiskUBO() {
        // disk
        float r1 = SagA.r_s * 2.2f;    // inner radius just outside the event horizon
        float r2 = SagA.r_s * 5.2f;   // outer radius of the disk
        float num = 2.0;               // number of rays
        float thickness = 1e9f;          // padding for std140 alignment
        float diskData[4] = { r1, r2, num, thickness };

        glBindBuffer(GL_UNIFORM_BUFFER, diskUBO);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(diskData), diskData);
    }
    
    vector<GLuint> QuadVAO(){
        float quadVertices[] = {
            // positions   // texCoords
            -1.0f,  1.0f,  0.0f, 1.0f,  // top left
            -1.0f, -1.0f,  0.0f, 0.0f,  // bottom left
            1.0f, -1.0f,  1.0f, 0.0f,  // bottom right

            -1.0f,  1.0f,  0.0f, 1.0f,  // top left
            1.0f, -1.0f,  1.0f, 0.0f,  // bottom right
            1.0f,  1.0f,  1.0f, 1.0f   // top right
        };
        
        GLuint VAO, VBO;
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);

        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), quadVertices, GL_STATIC_DRAW);

        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
        glEnableVertexAttribArray(1);

        GLuint texture;
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexImage2D(GL_TEXTURE_2D,
                    0,             // mip
                    GL_RGBA8,      // internal format
                    COMPUTE_WIDTH,
                    COMPUTE_HEIGHT,
                    0,
                    GL_RGBA,
                    GL_UNSIGNED_BYTE,
                    nullptr);
        vector<GLuint> VAOtexture = {VAO, texture};
        return VAOtexture;
    }
    void renderScene() {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glUseProgram(shaderProgram);
        glBindVertexArray(quadVAO);
        // make sure your fragment shader samples from texture unit 0:
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texture);
        glDrawArrays(GL_TRIANGLES, 0, 6);
        glfwSwapBuffers(window);
        glfwPollEvents();
    };
};
Engine engine;
void setupCameraCallbacks(GLFWwindow* window) {
    glfwSetWindowUserPointer(window, &camera);

    glfwSetMouseButtonCallback(window, [](GLFWwindow* win, int button, int action, int mods) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->processMouseButton(button, action, mods, win);
    });

    glfwSetCursorPosCallback(window, [](GLFWwindow* win, double x, double y) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->processMouseMove(x, y);
    });

    glfwSetScrollCallback(window, [](GLFWwindow* win, double xoffset, double yoffset) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->processScroll(xoffset, yoffset);
    });

    glfwSetKeyCallback(window, [](GLFWwindow* win, int key, int scancode, int action, int mods) {
        Camera* cam = (Camera*)glfwGetWindowUserPointer(win);
        cam->processKey(key, scancode, action, mods);
    });
}


// -- MAIN -- //
int main() {
    setupCameraCallbacks(engine.window);
    vector<unsigned char> pixels(engine.WIDTH * engine.HEIGHT * 3);

    auto t0 = Clock::now();
    lastPrintTime = chrono::duration<double>(t0.time_since_epoch()).count();

    double lastTime = glfwGetTime();
    int   renderW  = 800, renderH = 600, numSteps = 80000;
    while (!glfwWindowShouldClose(engine.window)) {
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);  // optional, but good practice
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        double now   = glfwGetTime();
        double dt    = now - lastTime;   // seconds since last frame
        if (dt > 0.05) {
            dt = 0.05;
        }
        lastTime     = now;

        // Gravity
        for (size_t i = 0; i < objects.size(); ++i) {
            if (isDynamicIndex(static_cast<int>(i))) continue;
            for (size_t j = 0; j < objects.size(); ++j) {
                if (i == j || isDynamicIndex(static_cast<int>(j))) continue; // skip self and meteor/debris
                auto& obj = objects[i];
                auto& obj2 = objects[j];
                float dx  = obj2.posRadius.x - obj.posRadius.x;
                float dy = obj2.posRadius.y - obj.posRadius.y;
                float dz = obj2.posRadius.z - obj.posRadius.z;
                float distance = sqrt(dx * dx + dy * dy + dz * dz);
                if (distance > 0) {
                    vector<double> direction = {dx / distance, dy / distance, dz / distance};
                    //distance *= 1000;
                    double Gforce = (G * obj.mass * obj2.mass) / (distance * distance);

                    double acc1 = Gforce / obj.mass;
                    std::vector<double> acc = {direction[0] * acc1, direction[1] * acc1, direction[2] * acc1};
                    if (Gravity) {
                        obj.velocity.x += acc[0];
                        obj.velocity.y += acc[1];
                        obj.velocity.z += acc[2];

                        obj.posRadius.x += obj.velocity.x;
                        obj.posRadius.y += obj.velocity.y;
                        obj.posRadius.z += obj.velocity.z;
                        cout << "velocity: " <<obj.velocity.x<<", " <<obj.velocity.y<<", " <<obj.velocity.z<<endl;
                    }
                }
            }
        }
        
        // --- simple scripted orbit (animation) ---
        // 讓 objects[1] 沿 XZ 平面做圓周運動：這樣不用碰重力/數值穩定性
        {
        // 只初始化一次，拿目前位置當作半徑與初相位
        static float orbitR = glm::length(glm::vec2(objects[1].posRadius.x, objects[1].posRadius.z));
        static float phase0 = std::atan2(objects[1].posRadius.z, objects[1].posRadius.x);

        const float omega = 0.8f; // 角速度(rad/s)，想更快就調大
        float ang = phase0 + omega * (float)now;

        objects[1].posRadius.x = orbitR * std::cos(ang);
        objects[1].posRadius.z = orbitR * std::sin(ang);
        objects[1].posRadius.y = 0.0f; // 固定在平面上
        }


        updateMeteor(dt);
        engine.pushOrbitTrailSample(vec3(objects[1].posRadius), float(now));

        // ---------- GRID ------------- //
        // 2) rebuild grid mesh on CPU
        engine.generateGrid(objects);
        // 5) overlay the bent grid
        mat4 view = lookAt(camera.position(), camera.target, vec3(0,1,0));
        mat4 proj = perspective(radians(60.0f), float(engine.COMPUTE_WIDTH)/engine.COMPUTE_HEIGHT, 1e9f, 1e14f);
        mat4 viewProj = proj * view;
        engine.drawGrid(viewProj);

        // ---------- RUN RAYTRACER ------------- //
        glViewport(0, 0, engine.WIDTH, engine.HEIGHT);
        engine.dispatchCompute(camera);
        engine.drawFullScreenQuad();

        engine.drawOrbitTrail(viewProj, camera.position(), float(now));

        // 6) present to screen
        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }

    glfwDestroyWindow(engine.window);
    glfwTerminate();
    return 0;
}
