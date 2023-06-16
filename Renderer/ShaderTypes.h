/*
See LICENSE folder for this sample’s licensing information.

Abstract:
The header that contains the types and enumeration constants that the Metal shaders and the C/Objective-C source share.
*/

#ifndef ShaderTypes_h
#define ShaderTypes_h

#include <simd/simd.h>

#define GEOMETRY_MASK_TRIANGLE 1
#define GEOMETRY_MASK_SPHERE   2
#define GEOMETRY_MASK_LIGHT    4

#define GEOMETRY_MASK_GEOMETRY (GEOMETRY_MASK_TRIANGLE | GEOMETRY_MASK_SPHERE)

#define RAY_MASK_PRIMARY   (GEOMETRY_MASK_GEOMETRY | GEOMETRY_MASK_LIGHT)
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
    unsigned int lightCount;
    Camera camera;
};

struct Material {
    vector_float3 color;
    bool is_metal = false;
    bool is_glass = false;
//    float metallic;
//    float roughness;
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
    Material material;
};


#endif
