//
//  scene.cpp
//  helloRayTracing
//
//  Created by Ziyuan Qu on 2023/6/13.
//

#include "scene.hpp"
#import <vector>
using namespace simd;

#pragma mark - Geometry
#pragma region Geometry {

Geometry* Geometry::init(MTL::Device* device){
    this->device = device;
    return this;
}

#pragma endregion Geometry }

#pragma mark - SphereGeometry
#pragma region SphereGeometry {

SphereGeometry* SphereGeometry::init(MTL::Device* device){
    this->device = device;
    return this;
}

void SphereGeometry::uploadToBuffers(){
    MTL::Device* device = this->device;
    this->sphere_buffer = device->newBuffer(this->spheres.size() * sizeof(Sphere), MTL::ResourceStorageModeManaged);
    this->bounding_box_buffer = device->newBuffer(this->spheres.size() * sizeof(BoundingBox), MTL::ResourceStorageModeManaged);
    std::vector<BoundingBox> bounding_boxes;
    for (Sphere sphere : spheres){
        BoundingBox bbox;

        bbox.min.x = sphere.origin.x - sphere.radius;
        bbox.min.y = sphere.origin.y - sphere.radius;
        bbox.min.z = sphere.origin.z - sphere.radius;

        bbox.max.x = sphere.origin.x + sphere.radius;
        bbox.max.y = sphere.origin.y + sphere.radius;
        bbox.max.z = sphere.origin.z + sphere.radius;

        bounding_boxes.push_back(bbox);
    }

    memcpy(this->sphere_buffer, this->spheres.data(), this->sphere_buffer->length());
    memcpy(this->bounding_box_buffer, bounding_boxes.data(), this->bounding_box_buffer->length());
}

void SphereGeometry::clear(){
    this->spheres.clear();
}

MTL::AccelerationStructureGeometryDescriptor* SphereGeometry::geometryDescriptor(){
    MTL::AccelerationStructureBoundingBoxGeometryDescriptor* discriptor = MTL::AccelerationStructureBoundingBoxGeometryDescriptor::alloc()->init();
    discriptor->setBoundingBoxBuffer(this->bounding_box_buffer);
    discriptor->setBoundingBoxCount(this->spheres.size());
    discriptor->setPrimitiveDataBuffer(this->sphere_buffer);
    discriptor->setPrimitiveDataStride(sizeof(Sphere));
    discriptor->setPrimitiveDataElementSize(sizeof(Sphere));
    return discriptor;
}

NS::Array* SphereGeometry::resources(){
    return NS::Array::array(this->sphere_buffer);
}

NS::String* SphereGeometry::intersectionFunctionName(){
    return NS::String::string( "sphereIntersectionFunction ", NS::UTF8StringEncoding );
}

#pragma endregion SphereGeometry }

#pragma mark - GeometryInstance
#pragma region GeometryInstance {

void GeometryInstance::init(Geometry* geometry, float4x4 transform, unsigned int mask){
    this->geometry = geometry;
    this->transform = transform;
    this->mask = mask;
}

#pragma endregion GeometryInstance }

#pragma mark - Scene
#pragma region Scene {

Scene* Scene::init(MTL::Device* device){
    this->device = device;
    this->geometries = new std::vector<Geometry*>();
    this->instances = new std::vector<GeometryInstance*>();
    this->camera.position = {0.0f, 0.0f, -1.0f};
    this->camera.look_at = {0.0f, 0.0f, 0.0f};
    this->camera.up = {0.0f, 1.0f, 0.0f};
    return this;
}

void Scene::clear(){
    this->geometries->clear();
    this->instances->clear();
}

void Scene::addGeometry(Geometry* mesh){
    this->geometries->push_back(mesh);
}

void Scene::addInstance(GeometryInstance* instance){
    this->instances->push_back(instance);
}

void Scene::uploadToBuffers(){
    for (Geometry* geometry : *this->geometries){
        geometry->uploadToBuffers();
    }
    // Lights here
}

#pragma endregion Scene }

#pragma mark - New Scene
static Scene* newScene(MTL::Device* device){
    Scene* scene = Scene::alloc()->init(device);
    scene->camera.position = {0.0f, 1.0f, -10.0f};
    scene->camera.look_at = {0.0f, 1.0f, 0.0f};
    scene->camera.up = {0.0f, 1.0f, 0.0f};
    
    // Add a sphere.
    SphereGeometry* sphere = SphereGeometry::alloc()->init(device);
    sphere->addSphereWithOrigin({0.0f, 1.0f, 10.0f}, 2.0f, {0.725f, 0.71f, 0.68f});
    scene->addGeometry(sphere);
    
    return scene;
}

