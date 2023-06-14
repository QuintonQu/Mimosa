//
//  scene.hpp
//  helloRayTracing
//
//  Created by Ziyuan Qu on 2023/6/13.
//

#ifndef scene_hpp
#define scene_hpp

#include <stdio.h>
#include <string>
#include <vector>
#import <Metal/Metal.hpp>
#import <simd/simd.h>

struct BoundingBox {
    MTL::PackedFloat3 min;
    MTL::PackedFloat3 max;
};

struct Camera {
    simd::float3 position;
    simd::float3 look_at;
    simd::float3 up;
};

struct Sphere {
    simd::float3 origin;
    float radius;
    simd::float3 color;
};

class Geometry : public NS::Object {
public:
    // Initializer.
    static Geometry* alloc();
    
    Geometry* init(MTL::Device* device);

    // Upload the primitives to Metal buffers so the GPU can access them.
    void uploadToBuffers(){};
    
    // Reset the geometry, removing all primitives.
    void clear(){};

    // Get the acceleration structure geometry descriptor for this piece of
    // geometry.
    MTL::AccelerationStructureGeometryDescriptor* geometryDescriptor(){return NULL;};

    // Get the array of Metal resources, such as buffers and textures, to pass
    // to the geometry's intersection function.
    NS::Array* resources(){return NULL;};

    // Name of the intersection function to use for this geometry, or empty
    // for triangles.
    NS::String* intersectionFunctionName(){return NULL;};
    
    // Metal device used to create the acceleration structures.
    MTL::Device* device;
};

class SphereGeometry : public Geometry {
public:
    void addSphereWithOrigin(simd::float3 origin, float radius, simd::float3 color);
    
    static SphereGeometry* alloc();
    SphereGeometry* init(MTL::Device* device);
    void uploadToBuffers();
    void clear();
    MTL::AccelerationStructureGeometryDescriptor* geometryDescriptor();
    NS::Array* resources();
    NS::String* intersectionFunctionName();
    
    std::vector<Sphere> spheres;
    MTL::Buffer* sphere_buffer;
    MTL::Buffer* bounding_box_buffer;
};

class GeometryInstance : public NS::Object {
public:
    void init(Geometry* geometry, simd::float4x4 transform, unsigned int mask);
    Geometry* geometry;
    simd::float4x4 transform;
    unsigned int mask;
};

class Scene : public NS::Object {
public:
    static Scene* alloc();
    
    Scene* init(MTL::Device* device);

    // Add a piece of geometry to the scene.
    void addGeometry(Geometry* mesh);

    // Add an instance of a piece of geometry to the scene.
    void addInstance(GeometryInstance* instance);

    // Add a light to the scene.
//    void addLight(AreaLight light);

    // Remove all geometry, instances, and lights from the scene.
    void clear();

    // Upload all scene data to Metal buffers so the GPU can access the data.
    void uploadToBuffers();
    
    
    MTL::Device* device;
    std::vector<Geometry*>* geometries;
    std::vector<GeometryInstance*>* instances;
    MTL::Buffer light_buffer;
    int light_count;
    Camera camera;
    
};

static Scene* newScene(MTL::Device* device);

#endif /* scene_hpp */
