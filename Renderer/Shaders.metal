/*
See LICENSE folder for this sample’s licensing information.

Abstract:
The Metal shaders used for this sample.
*/
#include "ShaderTypes.h"

#include <metal_stdlib>
#include <simd/simd.h>

using namespace metal;

using namespace raytracing;

constant unsigned int resourcesStride   [[function_constant(0)]];
constant bool useIntersectionFunctions  [[function_constant(1)]];
constant bool usePerPrimitiveData       [[function_constant(2)]];
//constant bool useResourcesBuffer = !usePerPrimitiveData;

constant unsigned int primes[] = {
    2,   3,  5,  7,
    11, 13, 17, 19,
    23, 29, 31, 37,
    41, 43, 47, 53,
    59, 61, 67, 71,
    73, 79, 83, 89
};

#pragma mark - Generate Random Sample
// Returns the i'th element of the Halton sequence using the d'th prime number as a
// base. The Halton sequence is a low discrepency sequence: the values appear
// random, but are more evenly distributed than a purely random sequence. Each random
// value used to render the image uses a different independent dimension, `d`,
// and each sample (frame) uses a different index `i`. To decorrelate each pixel,
// you can apply a random offset to `i`.
float halton(unsigned int i, unsigned int d) {
    unsigned int b = primes[d % 24];

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

#pragma mark - Inline Utility Interpolate
// Interpolates the vertex attribute of an arbitrary type across the surface of a triangle
// given the barycentric coordinates and triangle index in an intersection structure.
template<typename T, typename IndexType>
inline T interpolateVertexAttribute(device T *attributes,
                                    IndexType i0,
                                    IndexType i1,
                                    IndexType i2,
                                    float2 uv) {
    // Look up value for each vertex.
    const T T0 = attributes[i0];
    const T T1 = attributes[i1];
    const T T2 = attributes[i2];

    // Compute the sum of the vertex attributes weighted by the barycentric coordinates.
    // The barycentric coordinates sum to one.
    return (1.0f - uv.x - uv.y) * T0 + uv.x * T1 + uv.y * T2;
}

template<typename T>
inline T interpolateVertexAttribute(thread T *attributes, float2 uv) {
    // Look up the value for each vertex.
    const T T0 = attributes[0];
    const T T1 = attributes[1];
    const T T2 = attributes[2];

    // Compute the sum of the vertex attributes weighted by the barycentric coordinates.
    // The barycentric coordinates sum to one.
    return (1.0f - uv.x - uv.y) * T0 + uv.x * T1 + uv.y * T2;
}

#pragma mark - Inline Utility Transform

__attribute__((always_inline))
float3 transformPoint(float3 p, float4x4 transform) {
    return (transform * float4(p.x, p.y, p.z, 1.0f)).xyz;
}

__attribute__((always_inline))
float3 transformDirection(float3 p, float4x4 transform) {
    return (transform * float4(p.x, p.y, p.z, 0.0f)).xyz;
}

#pragma mark - Inline Utility Sample
// Uses the inversion method to map two uniformly random numbers to a 3D
// unit hemisphere, where the probability of a given sample is proportional to the cosine
// of the angle between the sample direction and the "up" direction (0, 1, 0).
inline float3 sampleCosineWeightedHemisphere(float2 u) {
    float phi = 2.0f * M_PI_F * u.x;

    float cos_phi;
    float sin_phi = sincos(phi, cos_phi);

    float cos_theta = sqrt(u.y);
    float sin_theta = sqrt(1.0f - cos_theta * cos_theta);

    return float3(sin_theta * cos_phi, cos_theta, sin_theta * sin_phi);
}

// Sample cosine weighted hemisphere with power
inline float3 sampleCosineWeightedHemisphereWithPower(float2 u, float exponent) {
    float phi = 2.0f * M_PI_F * u.x;

    float cos_phi;
    float sin_phi = sincos(phi, cos_phi);

    float cos_theta = pow(u.y, 1/(exponent + 1));
    float sin_theta = sqrt(1.0f - cos_theta * cos_theta);

    return float3(sin_theta * cos_phi, cos_theta, sin_theta * sin_phi);
}

// Sample sphere cap
inline float3 sampleSphereCap(float2 u, float cos_theta_max) {
    float phi = 2.0f * M_PI_F * u.x;

    float cos_phi;
    float sin_phi = sincos(phi, cos_phi);

    float cos_theta = mix(cos_theta_max, 1.0f, u.y);
    float sin_theta = sqrt(1.0f - cos_theta * cos_theta);

    return float3(sin_theta * cos_phi, cos_theta, sin_theta * sin_phi);
}

// Maps two uniformly random numbers to the surface of a 2D area light
// source and returns the direction to this point, the amount of light that travels
// between the intersection point and the sample point on the light source, as well
// as the distance between these two points.

inline void sampleAreaLight(constant AreaLight & light,
                            float2 u,
                            float3 position,
                            thread float3 & lightDirection,
                            thread float3 & lightColor,
                            thread float & lightDistance)
{
    // Map to -1..1
    u = u * 2.0f - 1.0f;

    // Transform into the light's coordinate system.
    float3 samplePosition = light.position +
                            light.right * u.x +
                            light.up * u.y;

    // Compute the vector from sample point on  the light source to intersection point.
    lightDirection = samplePosition - position;

    lightDistance = length(lightDirection);

    float inverseLightDistance = 1.0f / max(lightDistance, 1e-3f);

    // Normalize the light direction.
    lightDirection *= inverseLightDistance;

    // Start with the light's color.
    lightColor = light.color;

    // Light falls off with the inverse square of the distance to the intersection point.
    lightColor *= (inverseLightDistance * inverseLightDistance);

    // Light also falls off with the cosine of the angle between the intersection point
    // and the light source.
    lightColor *= saturate(dot(-lightDirection, light.forward));
}

// Return the type for Scatter Record for different light.
struct EmitterRecord {
    bool accept;
    float3 out_direction;
    float3 emit;
    float pdf;
    float distance;
    bool test = false;
};

// Aligns a direction on the unit hemisphere such that the hemisphere's "up" direction
// (0, 1, 0) maps to the given surface normal direction.
inline float3 alignHemisphereWithNormal(float3 sample, float3 normal) {
    // Set the "up" vector to the normal
    float3 up = normal;

    // Find an arbitrary direction perpendicular to the normal, which becomes the
    // "right" vector.
    float3 right = normalize(cross(normal, float3(0.0072f, 1.0f, 0.0034f)));

    // Find a third vector perpendicular to the previous two, which becomes the
    // "forward" vector.
    float3 forward = cross(right, up);

    // Map the direction on the unit hemisphere to the coordinate system aligned
    // with the normal.
    return sample.x * right + sample.y * up + sample.z * forward;
}

// Return the type for a bounding box intersection function.
struct BoundingBoxIntersection {
    bool accept    [[accept_intersection]]; // Whether to accept or reject the intersection.
    float distance [[distance]];            // Distance from the ray origin to the intersection point.
};

// Sample a sphere light
inline EmitterRecord sampleSphereLight(Sphere sphere, float3 hit_point, float2 random_variable, float4x4 transform){
    EmitterRecord emitter_record;
    float3 origin = transformPoint(sphere.origin, transform);
    float distance = length(hit_point - origin);
    float radius = sphere.radius;
    if(distance < radius){
        emitter_record.accept = false;
        emitter_record.pdf = 0.f;
        return emitter_record;
    }
    float cos_max = sqrt(distance * distance - radius * radius) / distance;
    emitter_record.out_direction = alignHemisphereWithNormal(sampleSphereCap(random_variable, cos_max), normalize(origin - hit_point));
    float area = 2 * M_PI_F * (1-cos_max);
    emitter_record.pdf = abs(1.f / area);
    emitter_record.emit = sphere.material.color;
    emitter_record.distance = distance;
    return emitter_record;
}

inline float spherePdf(Sphere sphere, float3 hit_point, float4x4 transform){
    float3 origin = transformPoint(sphere.origin, transform);
    float distance = length(hit_point - origin);
    float radius = sphere.radius;
    if(distance < radius){
        return 0.f;
    }
    float cos_max = sqrt(distance * distance - radius * radius) / distance;
    float area = 2 * M_PI_F * (1-cos_max);
    return abs(1.f / area);
}

# pragma mark - Materials' Scatter Methods

// Return the type for Scatter Record for different material.
struct ScatterRecord {
    bool accept;
    float3 out_direction;
    float3 attenuation;
    float pdf;
    bool test = false;
};

// Material scatter function -> Metallic
inline ScatterRecord metallicScatter(float3 normal, float3 in_direction) {
    ScatterRecord scatter_record;
    
    scatter_record.out_direction = reflect(in_direction, normal);
    scatter_record.attenuation = float3(1.0f, 1.0f, 1.0f);
    scatter_record.accept = true;
    
    return scatter_record;
}

// Fresnel Dieletric
inline float fresnel_dielectric(float cos_theta_i, float eta_i, float eta_t)
{
    if (eta_t == eta_i)
        return 0.f;

    // Swap the indices of refraction if the interaction starts at the inside of the object
    bool entering = cos_theta_i > 0.0f;
    if (!entering)
    {
        float eta = 0.0f;
        eta = eta_i;
        eta_i = eta_t;
        eta_t = eta;
        cos_theta_i = -cos_theta_i;
    }

    // Using Sahl-Snell's law, calculate the squared sine of the angle between the normal and the transmitted ray
    float eta          = eta_i / eta_t;
    float sin_theta_t2 = eta * eta * (1 - cos_theta_i * cos_theta_i);

    // Total internal reflection!
    if (sin_theta_t2 > 1.0f)
        return 1.0f;

    float cos_theta_t = sqrt(1.0f - sin_theta_t2);

    float Rs = (eta_i * cos_theta_i - eta_t * cos_theta_t) / (eta_i * cos_theta_i + eta_t * cos_theta_t);
    float Rp = (eta_t * cos_theta_i - eta_i * cos_theta_t) / (eta_t * cos_theta_i + eta_i * cos_theta_t);

    return 0.5f * (Rs * Rs + Rp * Rp);
}


// Glossy scatter function -> Glossy
inline ScatterRecord glossyScatter(float3 normal, float3 in_direction, float random_variable) {
    ScatterRecord scatter_record;
    
    float ior = 1.5f;
    float cosine = dot(normal, in_direction);
    float reflected = fresnel_dielectric(cosine, ior, 1.0f);
    if (cosine > 0.0f){
        normal = -normal;
        ior = 1.0f / ior;
    }
    
    if(random_variable < reflected){
        scatter_record.out_direction = reflect(in_direction, normal);
        scatter_record.attenuation = float3(1.0f, 1.0f, 1.0f);
    } else {
        scatter_record.out_direction = refract(in_direction, normal, 1.0f/ior);
        scatter_record.attenuation = float3(1.0f, 1.0f, 1.0f);
    }
    scatter_record.accept = true;
    
    return scatter_record;
}

// Phong Material
inline ScatterRecord phongScatter(float3 normal, float3 in_direction, float2 random_variable, Material material ){
    ScatterRecord scatter_record;
    
    float3 reflection_direction = reflect(in_direction, normal);
    float3 random_sample_direction = sampleCosineWeightedHemisphereWithPower(random_variable, material.exponent);
    scatter_record.out_direction = alignHemisphereWithNormal(random_sample_direction, reflection_direction);
    scatter_record.attenuation = material.color;
    scatter_record.accept = dot(in_direction, normal) < 0.f && dot(scatter_record.out_direction, normal) > 0.f;
    
    return scatter_record;
}

inline float phongPdf(float3 normal, float3 in_direction, float3 out_direction, Material material){
    if(dot(in_direction, normal) > 0.f || dot(out_direction, normal) < 0.f){
        return 0.f;
    }
    float3 mirror_direction = normalize(reflect(in_direction, normal));
    float cosine = saturate(dot(out_direction, mirror_direction));
    return ((material.exponent + 1.f)/(2.f*M_PI_F))*pow(cosine, material.exponent);
}

// Resources for a piece of triangle geometry.
struct TriangleResources {
    device uint16_t *indices;
    device float3 *vertexNormals;
    device float3 *vertexColors;
};

// Resources for a piece of sphere geometry.
struct SphereResources {
    device Sphere *spheres;
};

#pragma mark - Sphere Intersection
/*
 Custom sphere intersection function. The [[intersection]] keyword marks this as an intersection
 function. The [[bounding_box]] keyword means that this intersection function handles intersecting rays
 with bounding box primitives. To create sphere primitives, the sample creates bounding boxes that
 enclose the sphere primitives.

 The [[triangle_data]] and [[instancing]] keywords indicate that the intersector that calls this
 intersection function returns barycentric coordinates for triangle intersections and traverses
 an instance acceleration structure. These keywords must match between the intersection functions,
 intersection function table, intersector, and intersection result to ensure that Metal propagates
 data correctly between stages. Using fewer tags when possible may result in better performance,
 as Metal may need to store less data and pass less data between stages. For example, if you do not
 need barycentric coordinates, omitting [[triangle_data]] means Metal can avoid computing and storing
 them.

 The arguments to the intersection function contain information about the ray, primitive to be
 tested, and so on. The ray intersector provides this datas when it calls the intersection function.
 Metal provides other built-in arguments, but this sample doesn't use them.
 */
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
//                                                   device void *resources               [[buffer(0), function_constant(useResourcesBuffer)]]
//#if SUPPORTS_METAL_3
                                                   const device void* perPrimitiveData [[primitive_data]]
//#endif
                                                   )
{
    Sphere sphere;
    sphere = *(const device Sphere*)perPrimitiveData;
    
//#if SUPPORTS_METAL_3
//    // Look up the resources for this piece of sphere geometry.
//    if (usePerPrimitiveData) {
//        // Per-primitive data points to data from the specified buffer as was configured in the MTLAccelerationStructureBoundingBoxGeometryDescriptor.
//        sphere = *(const device Sphere*)perPrimitiveData;
//    } else
//#endif
//    {
//        device SphereResources& sphereResources = *(device SphereResources *)((device char *)resources + resourcesStride * geometryIndex);
//        // Get the actual sphere enclosed in this bounding box.
//        sphere = sphereResources.spheres[primitiveIndex];
//    }

    // Check for intersection between the ray and sphere mathematically.
    float3 oc = origin - sphere.origin;

    float a = dot(direction, direction);
    float b = 2 * dot(oc, direction);
    float c = dot(oc, oc) - sphere.radiusSquared;

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
        
        // Modify to implement the transparent glass with refrection. Both root are needed to check.
        if (ret.distance > 1e-3f){
            ret.accept = ret.distance >= minDistance && ret.distance <= maxDistance;
        } else {
            ret.distance = (-b + sqrt(disc)) / (2 * a);
            ret.accept = ret.distance >= minDistance && ret.distance <= maxDistance;
        }
    }
    return ret;
}

#pragma mark - Generate Ray
#pragma region GenerateRay {

ray generate_ray(uint2 tid [[thread_position_in_grid]], constant Uniforms &uniforms [[buffer(0)]], texture2d<unsigned int> randomTex [[texture(0)]]) {
    // The ray to cast.
    ray ray;

    // Pixel coordinates for this thread.
    float2 pixel = (float2)tid;

    // Apply a random offset to the random number index to decorrelate pixels.
    unsigned int offset = randomTex.read(tid).x;

    // Add a random offset to the pixel coordinates for antialiasing.
    float2 r = float2(halton(offset + uniforms.frameIndex, 0),
                      halton(offset + uniforms.frameIndex, 1));

    pixel += r;

    // Map pixel coordinates to -1..1.
    float2 uv = (float2)pixel / float2(uniforms.width, uniforms.height);
    uv = uv * 2.0f - 1.0f;

    constant Camera & camera = uniforms.camera;

    // Rays start at the camera position.
    ray.origin = camera.position;

    // Map normalized pixel coordinates into camera's coordinate system.
    ray.direction = normalize(uv.x * camera.right +
                              uv.y * camera.up +
                              camera.forward);

    // Don't limit intersection distance.
    ray.max_distance = INFINITY;
    ray.min_distance = 0.0001;
    
    return ray;
}

#pragma endregion GenerateRay }

#pragma mark - Ray Tracing Kernel - Mats

// Main ray tracing kernel.
kernel void raytracingKernelMats(
     uint2                                                  tid                       [[thread_position_in_grid]],
     constant Uniforms &                                    uniforms                  [[buffer(0)]],
     texture2d<unsigned int>                                randomTex                 [[texture(0)]],
     texture2d<float>                                       prevTex                   [[texture(1)]],
     texture2d<float, access::write>                        dstTex                    [[texture(2)]],
//     device void                                           *resources                 [[buffer(1), function_constant(useResourcesBuffer)]],
     constant MTLAccelerationStructureInstanceDescriptor   *instances                 [[buffer(2)]],
     constant AreaLight                                    *areaLights                [[buffer(3)]],
     instance_acceleration_structure                        accelerationStructure     [[buffer(4)]],
     intersection_function_table<triangle_data, instancing> intersectionFunctionTable [[buffer(5)]]
)
{
    // The sample aligns the thread count to the threadgroup size, which means the thread count
    // may be different than the bounds of the texture. Test to make sure this thread
    // is referencing a pixel within the bounds of the texture.
    if (tid.x > uniforms.width || tid.y > uniforms.height) return;
    
    // Apply a random offset to the random number index to decorrelate pixels.
    unsigned int offset = randomTex.read(tid).x;
    
    // Generate Ray.
    ray ray = generate_ray(tid, uniforms, randomTex);

    // Start with a fully white color. The kernel scales the light each time the
    // ray bounces off of a surface, based on how much of each light component
    // the surface absorbs.
    float3 color = float3(1.0f, 1.0f, 1.0f);

    float3 accumulatedColor = float3(0.0f, 0.0f, 0.0f);
    
    float3 background_color = float3(0.0f, 0.0f, 0.0f);

    // Create an intersector to test for intersection between the ray and the geometry in the scene.
    intersector<triangle_data, instancing> i;

    // If the sample isn't using intersection functions, provide some hints to Metal for
    // better performance.
//    if (!useIntersectionFunctions) {
//        i.assume_geometry_type(geometry_type::triangle);
//        i.force_opacity(forced_opacity::opaque);
//    }

    typename intersector<triangle_data, instancing>::result_type intersection;
    
    ScatterRecord scatter_record;

    // Simulate up to three ray bounces. Each bounce propagates light backward along the
    // ray's path toward the camera.
    for (int bounce = 0; bounce < 2; bounce++) {
        // Get the closest intersection, not the first intersection. This is the default, but
        // the sample adjusts this property below when it casts shadow rays.
        i.accept_any_intersection(false);

        // Check for intersection between the ray and the acceleration structure. If the sample
        // isn't using intersection functions, it doesn't need to include one.
        if (useIntersectionFunctions)
            intersection = i.intersect(ray, accelerationStructure, bounce == 0 ? RAY_MASK_PRIMARY : RAY_MASK_PRIMARY, intersectionFunctionTable);
        else
            intersection = i.intersect(ray, accelerationStructure, bounce == 0 ? RAY_MASK_PRIMARY : RAY_MASK_SECONDARY);

        // Stop if the ray didn't hit anything and has bounced out of the scene.
        if (intersection.type == intersection_type::none) {
            accumulatedColor = background_color * color;
            break;
        }
            
        unsigned int instanceIndex = intersection.instance_id;

        // Look up the mask for this instance, which indicates what type of geometry the ray hit.
        unsigned int mask = instances[instanceIndex].mask;

        // If the ray hit a light source, set the color to white, and stop immediately.
        if(mask & GEOMETRY_MASK_LIGHT){
            if(mask & GEOMETRY_MASK_SPHERE_LIGHT){
                Sphere area_light;
                area_light = *(const device Sphere*)intersection.primitive_data;
                float3 light_color = area_light.material.color;
                accumulatedColor = color * light_color;
                break;
            }
            if(mask & GEOMETRY_MASK_TRIANGLE_LIGHT){
                Triangle area_light;
                area_light = *(const device Triangle*)intersection.primitive_data;
                float3 light_color = area_light.material.color;
                accumulatedColor = color * light_color;
                break;
            }
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
        float2 barycentric_coords = intersection.triangle_barycentric_coord;

        float3 worldSpaceSurfaceNormal = 0.0f;
        float3 surfaceColor = 0.0f;
        Material material;

        if (mask & GEOMETRY_MASK_TRIANGLE) {
            Triangle triangle;
            
            float3 objectSpaceSurfaceNormal;
            
            triangle = *(const device Triangle*)intersection.primitive_data;
            
//#if SUPPORTS_METAL_3
//            if (usePerPrimitiveData) {
//                // Per-primitive data points to data from the specified buffer as was configured in the MTLAccelerationStructureTriangleGeometryDescriptor.
//                triangle = *(const device Triangle*)intersection.primitive_data;
//            } else
//#endif
//            {
//                // The ray hit a triangle. Look up the corresponding geometry's normal and UV buffers.
//                device TriangleResources & triangleResources = *(device TriangleResources *)((device char *)resources + resourcesStride * geometryIndex);
//
//                triangle.normals[0] =  triangleResources.vertexNormals[triangleResources.indices[primitiveIndex * 3 + 0]];
//                triangle.normals[1] =  triangleResources.vertexNormals[triangleResources.indices[primitiveIndex * 3 + 1]];
//                triangle.normals[2] =  triangleResources.vertexNormals[triangleResources.indices[primitiveIndex * 3 + 2]];
//
//                triangle.colors[0] =  triangleResources.vertexColors[triangleResources.indices[primitiveIndex * 3 + 0]];
//                triangle.colors[1] =  triangleResources.vertexColors[triangleResources.indices[primitiveIndex * 3 + 1]];
//                triangle.colors[2] =  triangleResources.vertexColors[triangleResources.indices[primitiveIndex * 3 + 2]];
//            }

            // Interpolate the vertex normal at the intersection point.
            objectSpaceSurfaceNormal = interpolateVertexAttribute(triangle.normals, barycentric_coords);

            // Interpolate the vertex color at the intersection point.
//            surfaceColor = interpolateVertexAttribute(triangle.colors, barycentric_coords);
            surfaceColor = triangle.material.color;
            material = triangle.material;

            // Transform the normal from object to world space.
            worldSpaceSurfaceNormal = normalize(transformDirection(objectSpaceSurfaceNormal, objectToWorldSpaceTransform));
        }
        else if (mask & GEOMETRY_MASK_SPHERE) {
            Sphere sphere;
            
            sphere = *(const device Sphere*)intersection.primitive_data;
            
//#if SUPPORTS_METAL_3
//            if (usePerPrimitiveData) {
//                // Per-primitive data points to data from the specified buffer as was configured in the MTLAccelerationStructureBoundingBoxGeometryDescriptor.
//                sphere = *(const device Sphere*)intersection.primitive_data;
//            } else
//#endif
//            {
//                // The ray hit a sphere. Look up the corresponding sphere buffer.
//                device SphereResources & sphereResources = *(device SphereResources *)((device char *)resources + resourcesStride * geometryIndex);
//                sphere = sphereResources.spheres[primitiveIndex];
//            }

            // Transform the sphere's origin from object space to world space.
            float3 worldSpaceOrigin = transformPoint(sphere.origin, objectToWorldSpaceTransform);

            // Compute the surface normal directly in world space.
            worldSpaceSurfaceNormal = normalize(worldSpaceIntersectionPoint - worldSpaceOrigin);

            // The sphere is a uniform color, so you don't need to interpolate the color across the surface.
            surfaceColor = sphere.color;
            surfaceColor = sphere.material.color;
            material = sphere.material;
        }

        // Choose a random light source to sample.
//        float lightSample = halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 0);
//        unsigned int lightIndex = min((unsigned int)(lightSample * uniforms.lightCount), uniforms.lightCount - 1);

        // Choose a random point to sample on the light source.
        float2 r = float2(halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 1),
                          halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 2));

//        float3 worldSpaceLightDirection;
//        float3 lightColor;
//        float lightDistance;
//
//        // Sample the lighting between the intersection point and the point on the area light.
//        sampleAreaLight(areaLights[lightIndex], r, worldSpaceIntersectionPoint, worldSpaceLightDirection,
//                        lightColor, lightDistance);
//
//        // Scale the light color by the cosine of the angle between the light direction and
//        // surface normal.
//        lightColor *= saturate(dot(worldSpaceSurfaceNormal, worldSpaceLightDirection));
//
//        // Scale the light color by the number of lights to compensate for the fact that
//        // the sample samples only one light source at random.
//        lightColor *= uniforms.lightCount;
//
        // Scale the ray color by the color of the surface to simulate the surface absorbing light.
//        color *= surfaceColor;
//
//        // Compute the shadow ray. The shadow ray checks whether the sample position on the
//        // light source is visible from the current intersection point.
//        // If it is, the kernel adds lighting to the output image.
//        struct ray shadowRay;
//
//        // Add a small offset to the intersection point to avoid intersecting the same
//        // triangle again.
//        shadowRay.origin = worldSpaceIntersectionPoint + worldSpaceSurfaceNormal * 1e-3f;
//
//        // Travel toward the light source.
//        shadowRay.direction = worldSpaceLightDirection;
//
//        // Don't overshoot the light source.
//        shadowRay.max_distance = lightDistance - 1e-3f;
//
//        // Shadow rays check only whether there is an object between the intersection point
//        // and the light source. Tell Metal to return after finding any intersection.
//        i.accept_any_intersection(true);
//
//        if (useIntersectionFunctions)
//            intersection = i.intersect(shadowRay, accelerationStructure, RAY_MASK_SHADOW, intersectionFunctionTable);
//        else
//            intersection = i.intersect(shadowRay, accelerationStructure, RAY_MASK_SHADOW);

//        // If there was no intersection, then the light source is visible from the original
//        // intersection  point. Add the light's contribution to the image.
//        if (intersection.type == intersection_type::none)
//            accumulatedColor += lightColor * color;

        // Choose a random direction to continue the path of the ray. This causes light to
        // bounce between surfaces. An app might evaluate a more complicated equation to
        // calculate the amount of light that reflects between intersection points.  However,
        // all the math in this kernel cancels out because this app assumes a simple diffuse
        // BRDF and samples the rays with a cosine distribution over the hemisphere (importance
        // sampling). This requires that the kernel only multiply the colors together. This
        // sampling strategy also reduces the amount of noise in the output image.
        r = float2(halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 3),
                   halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 4));
        
        ray.origin = worldSpaceIntersectionPoint;
        
        // Deal with scattering.
        if(material.is_metal){
            scatter_record = metallicScatter(worldSpaceSurfaceNormal, ray.direction);
            if(scatter_record.accept){
                ray.direction = scatter_record.out_direction;
            }else{
                break;
            }
        } else if(material.is_glass){
            float random_variable = r.x;
            scatter_record = glossyScatter(worldSpaceSurfaceNormal, ray.direction, random_variable);
            
            if(scatter_record.accept){
                ray.direction = scatter_record.out_direction;
            }else{
                break;
            }
        } else if(material.is_phong){
            float2 random_variable = r;
            scatter_record = phongScatter(worldSpaceSurfaceNormal, ray.direction, random_variable, material);
            if(scatter_record.accept){
                ray.direction = scatter_record.out_direction;
            }else{
                break;
            }
            color *= scatter_record.attenuation;
        } else{
            if(dot(worldSpaceSurfaceNormal, ray.direction) > 0){
                worldSpaceSurfaceNormal = -worldSpaceSurfaceNormal;
            }
            float3 worldSpaceSampleDirection = sampleCosineWeightedHemisphere(r);
            worldSpaceSampleDirection = alignHemisphereWithNormal(worldSpaceSampleDirection, worldSpaceSurfaceNormal);
            ray.direction = worldSpaceSampleDirection;
            color *= material.color;
        }
        
        ray.direction = normalize(ray.direction);
    }

    // Average this frame's sample with all of the previous frames.
    if (uniforms.frameIndex > 0) {
        float3 prevColor = prevTex.read(tid).xyz;
        prevColor *= uniforms.frameIndex;

        accumulatedColor += prevColor;
        accumulatedColor /= (uniforms.frameIndex + 1);
    }

    dstTex.write(float4(accumulatedColor, 1.0f), tid);
}

#pragma mark - Rat Tracing Kernel - NEE

// Main ray tracing kernel.
kernel void raytracingKernelNEE(
     uint2                                                  tid                       [[thread_position_in_grid]],
     constant Uniforms &                                    uniforms                  [[buffer(0)]],
     texture2d<unsigned int>                                randomTex                 [[texture(0)]],
     texture2d<float>                                       prevTex                   [[texture(1)]],
     texture2d<float, access::write>                        dstTex                    [[texture(2)]],
     device void                                           *resources                 [[buffer(1)]],
     constant MTLAccelerationStructureInstanceDescriptor   *instances                 [[buffer(2)]],
     constant AreaLight                                    *areaLights                [[buffer(3)]],
     instance_acceleration_structure                        accelerationStructure     [[buffer(4)]],
     intersection_function_table<triangle_data, instancing> intersectionFunctionTable [[buffer(5)]],
     constant unsigned int                                 *lightIndexs               [[buffer(6)]],
     constant unsigned int                                 *lightCounts               [[buffer(7)]]
)
{
    // The sample aligns the thread count to the threadgroup size, which means the thread count
    // may be different than the bounds of the texture. Test to make sure this thread
    // is referencing a pixel within the bounds of the texture.
    if (tid.x > uniforms.width || tid.y > uniforms.height) return;
    
    // Apply a random offset to the random number index to decorrelate pixels.
    unsigned int offset = randomTex.read(tid).x;
    
    // Generate Ray.
    ray ray = generate_ray(tid, uniforms, randomTex);

    // Start with a fully white color. The kernel scales the light each time the
    // ray bounces off of a surface, based on how much of each light component
    // the surface absorbs.
    float3 color = float3(1.0f, 1.0f, 1.0f);

    float3 accumulatedColor = float3(0.0f, 0.0f, 0.0f);
    
    float3 background_color = float3(0.0f, 0.0f, 0.0f);

    // Create an intersector to test for intersection between the ray and the geometry in the scene.
    intersector<triangle_data, instancing> i;

    typename intersector<triangle_data, instancing>::result_type intersection;
    
    ScatterRecord scatter_record;

    // Simulate up to three ray bounces. Each bounce propagates light backward along the
    // ray's path toward the camera.
    for (int bounce = 0; bounce < 4; bounce++) {
        // Get the closest intersection, not the first intersection. This is the default, but
        // the sample adjusts this property below when it casts shadow rays.
        i.accept_any_intersection(false);

        // Check for intersection between the ray and the acceleration structure. If the sample
        // isn't using intersection functions, it doesn't need to include one.
        intersection = i.intersect(ray, accelerationStructure, bounce == 0 ? RAY_MASK_PRIMARY : RAY_MASK_SECONDARY, intersectionFunctionTable);

        // Stop if the ray didn't hit anything and has bounced out of the scene.
        if (intersection.type == intersection_type::none) {
            accumulatedColor = background_color * color;
            break;
        }
            
        unsigned int instanceIndex = intersection.instance_id;

        // Look up the mask for this instance, which indicates what type of geometry the ray hit.
        unsigned int mask = instances[instanceIndex].mask;

        // If the ray hit a light source, set the color to white, and stop immediately.
        if(mask & GEOMETRY_MASK_LIGHT){
            if(mask & GEOMETRY_MASK_SPHERE_LIGHT){
                Sphere area_light;
                area_light = *(const device Sphere*)intersection.primitive_data;
                float3 light_color = area_light.material.color;
                accumulatedColor = light_color;
                break;
            }
            if(mask & GEOMETRY_MASK_TRIANGLE_LIGHT){
                Triangle area_light;
                area_light = *(const device Triangle*)intersection.primitive_data;
                float3 light_color = area_light.material.color;
                accumulatedColor = light_color;
                break;
            }
        }

        // The ray hit something. Look up the transformation matrix for this instance.
        float4x4 objectToWorldSpaceTransform(1.0f);

        for (int column = 0; column < 4; column++)
            for (int row = 0; row < 3; row++)
                objectToWorldSpaceTransform[column][row] = instances[instanceIndex].transformationMatrix[column][row];

        // Compute the intersection point in world space.
        float3 worldSpaceIntersectionPoint = ray.origin + ray.direction * intersection.distance;

        float2 barycentric_coords = intersection.triangle_barycentric_coord;

        float3 worldSpaceSurfaceNormal = 0.0f;
        float3 surfaceColor = 0.0f;
        Material material;

        if (mask & GEOMETRY_MASK_TRIANGLE) {
            Triangle triangle;
            
            float3 objectSpaceSurfaceNormal;
            
            triangle = *(const device Triangle*)intersection.primitive_data;

            // Interpolate the vertex normal at the intersection point.
            objectSpaceSurfaceNormal = interpolateVertexAttribute(triangle.normals, barycentric_coords);

            surfaceColor = triangle.material.color;
            material = triangle.material;

            // Transform the normal from object to world space.
            worldSpaceSurfaceNormal = normalize(transformDirection(objectSpaceSurfaceNormal, objectToWorldSpaceTransform));
        }
        else if (mask & GEOMETRY_MASK_SPHERE) {
            Sphere sphere;
            
            sphere = *(const device Sphere*)intersection.primitive_data;

            // Transform the sphere's origin from object space to world space.
            float3 worldSpaceOrigin = transformPoint(sphere.origin, objectToWorldSpaceTransform);

            // Compute the surface normal directly in world space.
            worldSpaceSurfaceNormal = normalize(worldSpaceIntersectionPoint - worldSpaceOrigin);

            // The sphere is a uniform color, so you don't need to interpolate the color across the surface.
            surfaceColor = sphere.material.color;
            material = sphere.material;
        }

        // Sample a random light instance
        EmitterRecord emitter_record;
        float light_instance_sample = halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 7);
        unsigned int light_instance_index = min((unsigned int)(light_instance_sample * uniforms.lightCount), uniforms.lightCount - 1);
        unsigned int light_instance = lightIndexs[light_instance_index];
        
        unsigned int geometries_index = instances[light_instance].accelerationStructureIndex;
        unsigned int light_mask = instances[light_instance_index].mask;
        
        float4x4 lightToWorldSpaceTransform(1.0f);

        for (int column = 0; column < 4; column++)
            for (int row = 0; row < 3; row++)
                lightToWorldSpaceTransform[column][row] = instances[instanceIndex].transformationMatrix[column][row];
        
        // Sample a random light geometry from one instance.
        float light_geometry_sample = halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 8);
        unsigned int light_geometry_index = min((unsigned int)(light_geometry_sample * lightCounts[light_instance_index]), lightCounts[light_instance_index] - 1);
        
        // Choose a random point to sample on the light source.
        float2 r = float2(halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 1),
                          halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 2));
        
        if(light_mask & GEOMETRY_MASK_SPHERE_LIGHT){
            device SphereResources & sphereResources = *(device SphereResources *)((device char *)resources + resourcesStride * geometries_index);
            Sphere sphere = sphereResources.spheres[light_geometry_index];
            emitter_record = sampleSphereLight(sphere, worldSpaceIntersectionPoint, r, lightToWorldSpaceTransform);
        }
        

        // Choose a random direction to continue the path of the ray. This causes light to
        // bounce between surfaces. An app might evaluate a more complicated equation to
        // calculate the amount of light that reflects between intersection points.  However,
        // all the math in this kernel cancels out because this app assumes a simple diffuse
        // BRDF and samples the rays with a cosine distribution over the hemisphere (importance
        // sampling). This requires that the kernel only multiply the colors together. This
        // sampling strategy also reduces the amount of noise in the output image.
        r = float2(halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 3),
                   halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 4));
        
        ray.origin = worldSpaceIntersectionPoint;
        float3 in_direction = ray.direction;
        
        // Deal with scattering.
        if(material.is_metal){
            scatter_record = metallicScatter(worldSpaceSurfaceNormal, ray.direction);
            if(scatter_record.accept){
                ray.direction = scatter_record.out_direction;
            }else{
                break;
            }
        } else if(material.is_glass){
            float random_variable = r.x;
            scatter_record = glossyScatter(worldSpaceSurfaceNormal, ray.direction, random_variable);
            
            if(scatter_record.accept){
                ray.direction = scatter_record.out_direction;
            }else{
                break;
            }
        } else if(material.is_phong){
            float2 random_variable = r;
            scatter_record = phongScatter(worldSpaceSurfaceNormal, ray.direction, random_variable, material);
            if(scatter_record.accept){
                ray.direction = scatter_record.out_direction;
            }else{
                break;
            }
        } else{
            if(dot(worldSpaceSurfaceNormal, ray.direction) > 0){
                worldSpaceSurfaceNormal = -worldSpaceSurfaceNormal;
            }
            float3 worldSpaceSampleDirection = sampleCosineWeightedHemisphere(r);
            worldSpaceSampleDirection = alignHemisphereWithNormal(worldSpaceSampleDirection, worldSpaceSurfaceNormal);
            ray.direction = worldSpaceSampleDirection;
        }
        
        // Deal with light sampling
        if(material.is_metal){

        }else if(material.is_glass){

        }else{
            // Sample light for diffuse material
            struct ray shadowRay;
            shadowRay.origin = worldSpaceIntersectionPoint + worldSpaceSurfaceNormal * 1e-3f;;
            shadowRay.direction = emitter_record.out_direction;
            shadowRay.max_distance = emitter_record.distance - 1e-3f;
            i.accept_any_intersection(true);
            intersection = i.intersect(shadowRay, accelerationStructure, RAY_MASK_SHADOW, intersectionFunctionTable);
            if (intersection.type != intersection_type::none){
                break;
            }
            float material_pdf;
            if(material.is_phong){
                material_pdf = phongPdf(worldSpaceSurfaceNormal, in_direction, shadowRay.direction, material);
            }else{
                material_pdf = saturate(dot(worldSpaceSurfaceNormal, shadowRay.direction)) / M_PI_F;
            }
            accumulatedColor += emitter_record.emit * color * material_pdf / emitter_record.pdf;
        }
        
        // Scale the ray color by the color of the surface to simulate the surface absorbing light.
        color *= surfaceColor;
        
        ray.direction = normalize(ray.direction);
    }

    // Average this frame's sample with all of the previous frames.
    if (uniforms.frameIndex > 0) {
        float3 prevColor = prevTex.read(tid).xyz;
        prevColor *= uniforms.frameIndex;

        accumulatedColor += prevColor;
        accumulatedColor /= (uniforms.frameIndex + 1);
    }

    dstTex.write(float4(accumulatedColor, 1.0f), tid);
}

#pragma mark - Rat Tracing Kernel - MIS

// Main ray tracing kernel.
kernel void raytracingKernelMIS(
     uint2                                                  tid                       [[thread_position_in_grid]],
     constant Uniforms &                                    uniforms                  [[buffer(0)]],
     texture2d<unsigned int>                                randomTex                 [[texture(0)]],
     texture2d<float>                                       prevTex                   [[texture(1)]],
     texture2d<float, access::write>                        dstTex                    [[texture(2)]],
     device void                                           *resources                 [[buffer(1)]],
     constant MTLAccelerationStructureInstanceDescriptor   *instances                 [[buffer(2)]],
     constant AreaLight                                    *areaLights                [[buffer(3)]],
     instance_acceleration_structure                        accelerationStructure     [[buffer(4)]],
     intersection_function_table<triangle_data, instancing> intersectionFunctionTable [[buffer(5)]],
     constant unsigned int                                 *lightIndexs               [[buffer(6)]],
     constant unsigned int                                 *lightCounts               [[buffer(7)]]
)
{
    // The sample aligns the thread count to the threadgroup size, which means the thread count
    // may be different than the bounds of the texture. Test to make sure this thread
    // is referencing a pixel within the bounds of the texture.
    if (tid.x > uniforms.width || tid.y > uniforms.height) return;
    
    // Apply a random offset to the random number index to decorrelate pixels.
    unsigned int offset = randomTex.read(tid).x;
    
    // Generate Ray.
    ray ray = generate_ray(tid, uniforms, randomTex);

    // Start with a fully white color. The kernel scales the light each time the
    // ray bounces off of a surface, based on how much of each light component
    // the surface absorbs.
    float3 color = float3(1.0f, 1.0f, 1.0f);

    float3 accumulatedColor = float3(0.0f, 0.0f, 0.0f);
    
    float3 background_color = float3(0.0f, 0.0f, 0.0f);
    
    float Le_weight = 1.f;
    float nee_pdf = 1.f;
    float brdf_pdf = 1.f;
    float mis_power = 3.f;

    // Create an intersector to test for intersection between the ray and the geometry in the scene.
    intersector<triangle_data, instancing> i;

    typename intersector<triangle_data, instancing>::result_type intersection;
    
    ScatterRecord scatter_record;

    // Simulate up to three ray bounces. Each bounce propagates light backward along the
    // ray's path toward the camera.
    for (int bounce = 0; bounce < 2; bounce++) {
        i.accept_any_intersection(false);
        intersection = i.intersect(ray, accelerationStructure, RAY_MASK_PRIMARY, intersectionFunctionTable);

        // Stop if the ray didn't hit anything and has bounced out of the scene.
        if (intersection.type == intersection_type::none) {
            accumulatedColor = background_color * color;
            break;
        }
            
        unsigned int instanceIndex = intersection.instance_id;

        // Look up the mask for this instance, which indicates what type of geometry the ray hit.
        unsigned int mask = instances[instanceIndex].mask;

        // The ray hit something. Look up the transformation matrix for this instance.
        float4x4 objectToWorldSpaceTransform(1.0f);

        for (int column = 0; column < 4; column++)
            for (int row = 0; row < 3; row++)
                objectToWorldSpaceTransform[column][row] = instances[instanceIndex].transformationMatrix[column][row];
        
        // Compute the intersection point in world space.
        float3 worldSpaceIntersectionPoint = ray.origin + ray.direction * intersection.distance;
        
        // If the ray hit a light source, set the color to white, and stop immediately.
        if(mask & GEOMETRY_MASK_LIGHT){
            if(mask & GEOMETRY_MASK_SPHERE_LIGHT){
                Sphere area_light;
                area_light = *(const device Sphere*)intersection.primitive_data;
                float3 light_color = area_light.material.color;
                if(bounce > 0){
                    nee_pdf = spherePdf(area_light, ray.origin, objectToWorldSpaceTransform);
                }
                Le_weight = pow(brdf_pdf, mis_power) / (pow(nee_pdf, mis_power) + pow(brdf_pdf, mis_power));
                accumulatedColor += light_color * Le_weight * color;
            }
            if(mask & GEOMETRY_MASK_TRIANGLE_LIGHT){
                Triangle area_light;
                area_light = *(const device Triangle*)intersection.primitive_data;
                float3 light_color = area_light.material.color;
                if(bounce > 0){
//                    nee_pdf = spherePdf(area_light, ray.origin, objectToWorldSpaceTransform);
                }
                Le_weight = pow(brdf_pdf, mis_power) / (pow(nee_pdf, mis_power) + pow(brdf_pdf, mis_power));
                accumulatedColor += light_color * Le_weight * color;
            }
            break;
        }
        color *= Le_weight;


        float3 worldSpaceSurfaceNormal = 0.0f;
        float3 surfaceColor = 0.0f;
        Material material;

        if (mask & GEOMETRY_MASK_TRIANGLE) {
            Triangle triangle;
            
            float3 objectSpaceSurfaceNormal;
            
            float2 barycentric_coords = intersection.triangle_barycentric_coord;
            
            triangle = *(const device Triangle*)intersection.primitive_data;

            // Interpolate the vertex normal at the intersection point.
            objectSpaceSurfaceNormal = interpolateVertexAttribute(triangle.normals, barycentric_coords);

            surfaceColor = triangle.material.color;
            material = triangle.material;

            // Transform the normal from object to world space.
            worldSpaceSurfaceNormal = normalize(transformDirection(objectSpaceSurfaceNormal, objectToWorldSpaceTransform));
        }
        else if (mask & GEOMETRY_MASK_SPHERE) {
            Sphere sphere;
            
            sphere = *(const device Sphere*)intersection.primitive_data;

            // Transform the sphere's origin from object space to world space.
            float3 worldSpaceOrigin = transformPoint(sphere.origin, objectToWorldSpaceTransform);

            // Compute the surface normal directly in world space.
            worldSpaceSurfaceNormal = normalize(worldSpaceIntersectionPoint - worldSpaceOrigin);

            // The sphere is a uniform color, so you don't need to interpolate the color across the surface.
            surfaceColor = sphere.color;
            surfaceColor = sphere.material.color;
            material = sphere.material;
        }

        // Sample a random light instance
        EmitterRecord emitter_record;
        float light_instance_sample = halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 7);
        unsigned int light_instance_index = min((unsigned int)(light_instance_sample * uniforms.lightCount), uniforms.lightCount - 1);
        unsigned int light_instance = lightIndexs[light_instance_index];
        
        unsigned int geometries_index = instances[light_instance].accelerationStructureIndex;
        unsigned int light_mask = instances[light_instance_index].mask;
        
        float4x4 lightToWorldSpaceTransform(1.0f);

        for (int column = 0; column < 4; column++)
            for (int row = 0; row < 3; row++)
                lightToWorldSpaceTransform[column][row] = instances[instanceIndex].transformationMatrix[column][row];
        
        // Sample a random light geometry from one instance.
        float light_geometry_sample = halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 8);
        unsigned int light_geometry_index = min((unsigned int)(light_geometry_sample * lightCounts[light_instance_index]), lightCounts[light_instance_index] - 1);
        
        // Choose a random point to sample on the light source.
        float2 r = float2(halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 1),
                          halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 2));
        
        if(light_mask & GEOMETRY_MASK_SPHERE_LIGHT){
            device SphereResources & sphereResources = *(device SphereResources *)((device char *)resources + resourcesStride * geometries_index);
            Sphere sphere = sphereResources.spheres[light_geometry_index];
            emitter_record = sampleSphereLight(sphere, worldSpaceIntersectionPoint, r, lightToWorldSpaceTransform);
        }

        // Choose a random direction to continue the path of the ray. This causes light to
        // bounce between surfaces. An app might evaluate a more complicated equation to
        // calculate the amount of light that reflects between intersection points.  However,
        // all the math in this kernel cancels out because this app assumes a simple diffuse
        // BRDF and samples the rays with a cosine distribution over the hemisphere (importance
        // sampling). This requires that the kernel only multiply the colors together. This
        // sampling strategy also reduces the amount of noise in the output image.
        r = float2(halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 3),
                   halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 4));
        
        ray.origin = worldSpaceIntersectionPoint;
        float3 in_direction = ray.direction;
        
        // Deal with scattering.
        if(material.is_metal){
            scatter_record = metallicScatter(worldSpaceSurfaceNormal, ray.direction);
            if(scatter_record.accept){
                ray.direction = scatter_record.out_direction;
            }else{
                break;
            }
        } else if(material.is_glass){
            float random_variable = r.x;
            scatter_record = glossyScatter(worldSpaceSurfaceNormal, ray.direction, random_variable);
            
            if(scatter_record.accept){
                ray.direction = scatter_record.out_direction;
            }else{
                break;
            }
        } else if(material.is_phong){
            float2 random_variable = r;
            scatter_record = phongScatter(worldSpaceSurfaceNormal, ray.direction, random_variable, material);
            if(scatter_record.accept){
                ray.direction = scatter_record.out_direction;
            }else{
                break;
            }
        } else{
            if(dot(worldSpaceSurfaceNormal, ray.direction) > 0){
//                worldSpaceSurfaceNormal = -worldSpaceSurfaceNormal;
            }
            float3 worldSpaceSampleDirection = sampleCosineWeightedHemisphere(r);
            worldSpaceSampleDirection = alignHemisphereWithNormal(worldSpaceSampleDirection, worldSpaceSurfaceNormal);
            ray.direction = worldSpaceSampleDirection;
        }
        
        ray.direction = normalize(ray.direction);
        
        // Deal with light sampling
        if(material.is_metal){

        }else if(material.is_glass){

        }else{
            // Sample light for diffuse material
            struct ray shadowRay;
            shadowRay.origin = worldSpaceIntersectionPoint + worldSpaceSurfaceNormal * 1e-3f;;
            shadowRay.direction = emitter_record.out_direction;
            shadowRay.max_distance = emitter_record.distance - 1e-3f;
            i.accept_any_intersection(true);
            intersection = i.intersect(shadowRay, accelerationStructure, RAY_MASK_SHADOW, intersectionFunctionTable);
            if (intersection.type != intersection_type::none){
                break;
            }
            if(material.is_phong){
                brdf_pdf = phongPdf(worldSpaceSurfaceNormal, in_direction, shadowRay.direction, material);
            }else{
                brdf_pdf = saturate(dot(shadowRay.direction, worldSpaceSurfaceNormal)) / M_PI_F;
            }
            nee_pdf = emitter_record.pdf;
            Le_weight = pow(nee_pdf, mis_power) / (pow(nee_pdf, mis_power) + pow(brdf_pdf, mis_power));
//            Le_weight = nee_pdf / (nee_pdf + brdf_pdf);
            accumulatedColor += Le_weight * emitter_record.emit * color * brdf_pdf / nee_pdf;
        }
        
        // Scale the ray color by the color of the surface to simulate the surface absorbing light.
        color *= surfaceColor;
        
        // Calculate BRDF_PDF for next trace.
        if(material.is_metal){

        }else if(material.is_glass){

        }else if(material.is_phong){
            brdf_pdf = phongPdf(worldSpaceSurfaceNormal, in_direction, ray.direction, material);
        }else{
            brdf_pdf = max(0.f, dot(ray.direction, worldSpaceSurfaceNormal)) / M_PI_F;
        }
    }

    // Average this frame's sample with all of the previous frames.
    if (uniforms.frameIndex > 0) {
        float3 prevColor = prevTex.read(tid).xyz;
        prevColor *= uniforms.frameIndex;

        accumulatedColor += prevColor;
        accumulatedColor /= (uniforms.frameIndex + 1);
    }

    dstTex.write(float4(accumulatedColor, 1.0f), tid);
}

// Screen filling quad in normalized device coordinates.
constant float2 quadVertices[] = {
    float2(-1, -1),
    float2(-1,  1),
    float2( 1,  1),
    float2(-1, -1),
    float2( 1,  1),
    float2( 1, -1)
};

struct CopyVertexOut {
    float4 position [[position]];
    float2 uv;
};

// Simple vertex shader that passes through NDC quad positions.
vertex CopyVertexOut copyVertex(unsigned short vid [[vertex_id]]) {
    float2 position = quadVertices[vid];

    CopyVertexOut out;

    out.position = float4(position, 0, 1);
    out.uv = position * 0.5f + 0.5f;

    return out;
}

// Simple fragment shader that copies a texture and applies a simple tonemapping function.
fragment float4 copyFragment(CopyVertexOut in [[stage_in]],
                             texture2d<float> tex)
{
    constexpr sampler sam(min_filter::nearest, mag_filter::nearest, mip_filter::none);

    float3 color = tex.sample(sam, in.uv).xyz;

    // Apply a simple tonemapping function to reduce the dynamic range of the
    // input image into a range which the screen can display.
    color = color / (1.0f + color);

    return float4(color, 1.0f);
}
