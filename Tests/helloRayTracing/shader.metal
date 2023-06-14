//
//  shader.metal
//  helloRayTracing
//
//  Created by Ziyuan Qu on 2023/6/12.
//

#include <metal_stdlib>
using namespace metal;
using namespace raytracing;

#define GEOMETRY_MASK_TRIANGLE 1
#define GEOMETRY_MASK_SPHERE   2
#define GEOMETRY_MASK_LIGHT    4

#define GEOMETRY_MASK_GEOMETRY (GEOMETRY_MASK_TRIANGLE | GEOMETRY_MASK_SPHERE)

#pragma region ShaderTypes {

struct Sphere {
    float3 origin;
    float radius;
    float3 color;
};

struct Camera {
    float3 position;
    float3 look_at;
    float3 up;
    float4x4 transform;
    float fov;
    float focal_distance;
};

struct Uniforms {
    unsigned int width;
    unsigned int height;
    unsigned int frameIndex;
    unsigned int lightCount;
    Camera camera;
};

struct AreaLight {
    vector_float3 position;
    vector_float3 forward;
    vector_float3 right;
    vector_float3 up;
    vector_float3 color;
};

#pragma endregion ShaderTypes }


#pragma region TransformFunc {

float3 transform_point(float3 p, float4x4 transform) {
    return (transform * float4(p.x, p.y, p.z, 1.0f)).xyz;
}

float3 transform_direction(float3 p, float4x4 transform) {
    return (transform * float4(p.x, p.y, p.z, 0.0f)).xyz;
}

#pragma endregion TransformFunc }


#pragma region RandomFunc {

constant unsigned int primes[] = {
    2,   3,  5,  7,
    11, 13, 17, 19,
    23, 29, 31, 37,
    41, 43, 47, 53,
    59, 61, 67, 71,
    73, 79, 83, 89
};

float halton(unsigned int i, unsigned int d) {
    unsigned int b = primes[d];

    float f = 1.0f;
    float invB = 1.0f / b;

    float r = 0;

    while (i > 0) {
        f = f * invB;
        r = r + f * (i % b);
        i = i / b;
    }

    return r;
}

#pragma region RandomFunc }


#pragma region GenerateRay {

ray generate_ray(uint2 tid [[thread_position_in_grid]], constant Uniforms &uniforms [[buffer(0)]], texture2d<unsigned int> randomTex [[texture(0)]]) {
    ray ray;
    float2 pixel = (float2)tid;
    unsigned int offset = randomTex.read(tid).x;
    float2 random_offset = float2(halton(offset + uniforms.frameIndex, 0), halton(offset + uniforms.frameIndex, 1));
    pixel += random_offset;
    float2 uv = (float2)pixel / float2(uniforms.width, uniforms.height);
    uv = uv * 2.0f - 1.0f;
    
    constant Camera & camera = uniforms.camera;
    float3 forward = camera.look_at - camera.position;
    forward = normalize(forward);
    float3 right = cross(camera.up, forward);
    right = normalize(right);
    float3 up = cross(forward, right);
    
    ray.origin = camera.position;
    ray.direction = normalize(uv.x * right + uv.y * up + forward);
    ray.max_distance = INFINITY;
    ray.min_distance = 0.0001;
    return ray;
}

#pragma endregion GenerateRay }


#pragma region Intersection {

struct BoundingBoxIntersection {
    bool accept    [[accept_intersection]]; // Whether to accept or reject the intersection.
    float distance [[distance]];            // Distance from the ray origin to the intersection point.
};

[[intersection(bounding_box, triangle_data, instancing)]]
BoundingBoxIntersection sphereIntersectionFunction(// Ray parameters passed to the ray intersector below
                                                   float3 origin                        [[origin]],
                                                   float3 direction                     [[direction]],
                                                   float minDistance                    [[min_distance]],
                                                   float maxDistance                    [[max_distance]],
                                                   // Information about the primitive.
                                                   unsigned int primitiveIndex          [[primitive_id]],
                                                   unsigned int geometryIndex           [[geometry_intersection_function_table_offset]],
                                                   // Custom resources bound to the intersection function table.
                                                   const device void* perPrimitiveData [[primitive_data]]
                                                   )
{
    Sphere sphere;
    // Look up the resources for this piece of sphere geometry.
    // Per-primitive data points to data from the specified buffer as was configured in the MTLAccelerationStructureBoundingBoxGeometryDescriptor.
    sphere = *(const device Sphere*)perPrimitiveData;

    // Check for intersection between the ray and sphere mathematically.
    float3 oc = origin - sphere.origin;

    float a = dot(direction, direction);
    float b = 2 * dot(oc, direction);
    float c = dot(oc, oc) - sphere.radius * sphere.radius;

    float disc = b * b - 4 * a * c;

    BoundingBoxIntersection ret;

    if (disc <= 0.0f) {
        // If the ray missed the sphere, return false.
        ret.accept = false;
    }
    else {
        // Otherwise, compute the intersection distance.
        ret.distance = (-b - sqrt(disc)) / (2 * a);

        // The intersection function must also check whether the intersection distance is
        // within the acceptable range. Intersection functions do not run in any particular order,
        // so the maximum distance may be different from the one passed into the ray intersector.
        ret.accept = ret.distance >= minDistance && ret.distance <= maxDistance;
    }

    return ret;
}

#pragma endregion Intersection }

// Main ray tracing kernel.
kernel void raytracingKernel(
     uint2                                                  tid                       [[thread_position_in_grid]],
     constant Uniforms &                                    uniforms                  [[buffer(0)]],
     texture2d<unsigned int>                                randomTex                 [[texture(0)]],
     texture2d<float>                                       prevTex                   [[texture(1)]],
     texture2d<float, access::write>                        dstTex                    [[texture(2)]],
     constant MTLAccelerationStructureInstanceDescriptor   *instances                 [[buffer(2)]],
     constant AreaLight                                    *areaLights                [[buffer(3)]],
     instance_acceleration_structure                        accelerationStructure     [[buffer(4)]],
     intersection_function_table<triangle_data, instancing> intersectionFunctionTable [[buffer(5)]]
)
{
    if (tid.x > uniforms.width || tid.y > uniforms.height) return;
    // Generate Ray.
    ray ray = generate_ray(tid, uniforms, randomTex);
    
    // Intersect with the scene.
    float3 accumulatedColor = float3(0.f, 0.f, 0.f);
    intersector<triangle_data, instancing> i;
    typename intersector<triangle_data, instancing>::result_type intersection;
    
    // Recursion
    for(int bounce = 0; bounce < 1; bounce++){
        intersection = i.intersect(ray, accelerationStructure, GEOMETRY_MASK_GEOMETRY, intersectionFunctionTable);
        if (intersection.type == intersection_type::none) break;
        unsigned int instanceIndex = intersection.instance_id;
        // Look up the mask for this instance, which indicates what type of geometry the ray hit.
        unsigned int mask = instances[instanceIndex].mask;
        if (mask == GEOMETRY_MASK_LIGHT) {
            accumulatedColor = float3(1.0f, 1.0f, 1.0f); // Why not area_light.color ?
            break;
        }
        
        // The ray hit something. Look up the transformation matrix for this instance.
        float4x4 objectToWorldSpaceTransform(1.0f);

        for (int column = 0; column < 4; column++)
            for (int row = 0; row < 3; row++)
                objectToWorldSpaceTransform[column][row] = instances[instanceIndex].transformationMatrix[column][row];

        // Compute the intersection point in world space.
        float3 worldSpaceIntersectionPoint = ray.origin + ray.direction * intersection.distance;

//        unsigned primitiveIndex = intersection.primitive_id;
//        unsigned int geometryIndex = instances[instanceIndex].accelerationStructureIndex;
//        float2 barycentric_coords = intersection.triangle_barycentric_coord;

        float3 worldSpaceSurfaceNormal = 0.0f;
        float3 surfaceColor = 0.0f;
        
        if (mask & GEOMETRY_MASK_SPHERE) {
            Sphere sphere;
            sphere = *(const device Sphere*)intersection.primitive_data;
            // Transform the sphere's origin from object space to world space.
            float3 worldSpaceOrigin = transform_point(sphere.origin, objectToWorldSpaceTransform);

            // Compute the surface normal directly in world space.
            worldSpaceSurfaceNormal = normalize(worldSpaceIntersectionPoint - worldSpaceOrigin);

            // The sphere is a uniform color, so you don't need to interpolate the color across the surface.
            surfaceColor = sphere.color;
        }
        accumulatedColor += surfaceColor;
        dstTex.write(float4(accumulatedColor, 1.0f), tid);
    }
}

struct v2f
{
    float4 position [[position]];
    half3 color;
};

v2f vertex vertexMain( uint vertexId [[vertex_id]],
                       device const float3* positions [[buffer(0)]],
                       device const float3* colors [[buffer(1)]] )
{
    v2f o;
    o.position = float4( positions[ vertexId ], 1.0 );
    o.color = half3 ( colors[ vertexId ] );
    return o;
}

half4 fragment fragmentMain( v2f in [[stage_in]] )
{
    return half4( in.color, 1.0 );
}

