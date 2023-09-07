/*
See LICENSE folder for this sample’s licensing information.

Abstract:
The header that contains the types and enumeration constants that the Metal shaders and the C/Objective-C source share.
*/

#ifndef ShaderTypes_h
#define ShaderTypes_h

#include <simd/simd.h>

#define GEOMETRY_MASK_TRIANGLE_LIGHT    4
#define GEOMETRY_MASK_SPHERE_LIGHT    8
#define GEOMETRY_MASK_VOLUME_CONTAINER_TRIANGLE    16
#define GEOMETRY_MASK_VOLUME_CONTAINER_SPHERE    32
#define GEOMETRY_MASK_TRIANGLE (1 | GEOMETRY_MASK_VOLUME_CONTAINER_TRIANGLE)
#define GEOMETRY_MASK_SPHERE   (2 | GEOMETRY_MASK_VOLUME_CONTAINER_SPHERE)

#define GEOMETRY_MASK_VOLUME_CONTAINER (GEOMETRY_MASK_VOLUME_CONTAINER_TRIANGLE | GEOMETRY_MASK_VOLUME_CONTAINER_SPHERE)
#define GEOMETRY_MASK_LIGHT (GEOMETRY_MASK_TRIANGLE_LIGHT | GEOMETRY_MASK_SPHERE_LIGHT)
#define GEOMETRY_MASK_GEOMETRY (GEOMETRY_MASK_TRIANGLE | GEOMETRY_MASK_SPHERE | GEOMETRY_MASK_VOLUME_CONTAINER)

#define RAY_MASK_PRIMARY   (GEOMETRY_MASK_GEOMETRY | GEOMETRY_MASK_LIGHT)
#define RAY_MASK_SHADOW_ALL   (GEOMETRY_MASK_GEOMETRY | GEOMETRY_MASK_LIGHT)
#define RAY_MASK_SHADOW    GEOMETRY_MASK_GEOMETRY
#define RAY_MASK_SECONDARY GEOMETRY_MASK_GEOMETRY

#ifndef __METAL_VERSION__
struct packed_float3 {
#ifdef __cplusplus
    packed_float3() = default;
    packed_float3(vector_float3 v) : x(v.x), y(v.y), z(v.z) {}
#endif
    float x;
    float y;
    float z;
};
#endif

struct Camera {
    vector_float3 position;
    vector_float3 right;
    vector_float3 up;
    vector_float3 forward;
    float fov = 45.0f;
};

struct AreaLight {
    vector_float3 position;
    vector_float3 forward;
    vector_float3 right;
    vector_float3 up;
    vector_float3 color;
};

struct Uniforms {
    unsigned int width;
    unsigned int height;
    unsigned int frameIndex;
    unsigned int lightCount; // This is light instance count, there might be more than one light.
    unsigned int totalLightCount;
    unsigned int instanceCount;
    Camera camera;
};

struct Material {
    // color of surface
    vector_float3 color;
    
    // material type
    bool is_metal = false;
    bool is_glass = false;
    bool is_phong = false;
    bool is_contain_volume = false;
    
    // phong parameter
    float exponent;
    
    // volume parameters
    float density;
    vector_float3 albedo;
    vector_float3 emission;
    
};

struct Sphere {
    packed_float3 origin;
    float radiusSquared;
    packed_float3 color;
    float radius;
    Material material;
};

struct Triangle {
    vector_float3 normals[3];
    vector_float3 colors[3];
    vector_float3 positions[3];
    Material material;
};

//struct InverseTransformMatrix {
//    vector_float3 invMatrix[3];
//};
#endif
