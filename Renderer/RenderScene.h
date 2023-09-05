/*
See LICENSE folder for this sample’s licensing information.

Abstract:
The header for the class that describes objects in a scene.
*/

#ifndef Scene_h
#define Scene_h

#import <Metal/Metal.h>

#import "Transforms.h"
#import "ShaderTypes.h"

#define FACE_MASK_NONE       0
#define FACE_MASK_NEGATIVE_X (1 << 0)
#define FACE_MASK_POSITIVE_X (1 << 1)
#define FACE_MASK_NEGATIVE_Y (1 << 2)
#define FACE_MASK_POSITIVE_Y (1 << 3)
#define FACE_MASK_NEGATIVE_Z (1 << 4)
#define FACE_MASK_POSITIVE_Z (1 << 5)
#define FACE_MASK_ALL        ((1 << 6) - 1)

struct BoundingBox {
    MTLPackedFloat3 min;
    MTLPackedFloat3 max;
};

MTLResourceOptions getManagedBufferStorageMode();

#pragma mark - Geometry

// Represents a piece of geometry in a scene, The sample composes Geometry objects
// from primitives such as triangles or spheres. Each geometry object has its
// own primitive acceleration structure and, optionally, an intersection function.
// The sample creates copies, or "instances" of geometry objects using the GeometryInstance
// class.
@interface Geometry : NSObject

// Metal device used to create the acceleration structures.
@property (nonatomic, readonly) id <MTLDevice> device;

// Name of the intersection function to use for this geometry, or nil
// for triangles.
@property (nonatomic, readonly) NSString *intersectionFunctionName;

// Number of geometries
@property (nonatomic, readonly) NSUInteger geometryCount;

// Initializer.
- (instancetype)initWithDevice:(id <MTLDevice>)device;

// Reset the geometry, removing all primitives.
- (void)clear;

// Upload the primitives to Metal buffers so the GPU can access them.
- (void)uploadToBuffers;

// Get the acceleration structure geometry descriptor for this piece of
// geometry.
- (MTLAccelerationStructureGeometryDescriptor *)geometryDescriptor;

// Get the array of Metal resources, such as buffers and textures, to pass
// to the geometry's intersection function.
- (NSArray <id <MTLResource>> *)resources;

@end

#pragma mark - Specific Geometry

// Represents a piece of geometry made of triangles.
@interface TriangleGeometry : Geometry

// Add a cube to the triangle geometry.
- (void)addCubeWithFaces:(unsigned int)faceMask
                   color:(vector_float3)color
               transform:(matrix_float4x4)transform
           inwardNormals:(bool)inwardNormals
                material:(Material)material;

// Add a XY-Plane
- (void)addXYPlane:(matrix_float4x4)transform
     inwardNormals:(bool)inwardNormals
          material:(Material)material;

// Add a mesh
//- (void)addMesh:(NSString*)file_name
//      transform:(matrix_float4x4)transform
//  inwardNormals:(bool)inwardNormals
//       material:(Material)material;

- (void)addGeometryWithURL:(NSURL *)URL
                 transform:(matrix_float4x4)transform
                  material:(Material)material;
            
@end

// Represents a piece of geometry made of spheres.
@interface SphereGeometry : Geometry

- (void)addSphereWithOrigin:(vector_float3)origin
                     radius:(float)radius
                      color:(vector_float3)color
                   material:(Material)material;

@end

#pragma mark - GeometryInstance

// Represents an instance, or copy, of a piece of geometry in a scene.
// Each instance has its own transformation matrix that determines
// where to place it in the scene.
@interface GeometryInstance : NSObject

// The geometry to use in the instance.
@property (nonatomic, readonly) Geometry *geometry;

// Transformation matrix describing where to place the geometry in the
// scene.
@property (nonatomic, readonly) matrix_float4x4 transform;

// Mask used to filter out intersections between rays and different
// types of geometry.
@property (nonatomic, readonly) unsigned int mask;

// Pointer to the material type of geometry instance
//@property (nonatomic, readonly) Material* material;

// Initializer.
- (instancetype)initWithGeometry:(Geometry *)geometry
                       transform:(matrix_float4x4)transform
                            mask:(unsigned int)mask;

@end

#pragma mark - Material

#pragma mark - Scene

// Represents an entire scene, including different types of geometry,
// instances of that geometry, lights, and a camera.
@interface RenderScene : NSObject

// The device used to create the scene.
@property (nonatomic, readonly) id <MTLDevice> device;

// Array of geometries in the scene.
@property (nonatomic, readonly) NSArray <Geometry *> *geometries;

// Array of geometry instances in the scene.
@property (nonatomic, readonly) NSArray <GeometryInstance *> *instances;

// Array of geometry instances in the scene. (NOT USED!)
@property (nonatomic, readonly) NSArray <GeometryInstance *> *light_instances;

// Buffer containing lights. (NOT USED!)
@property (nonatomic, readonly) id <MTLBuffer> lightBuffer;

// Buffer containing light instance index
@property (nonatomic, readonly) id <MTLBuffer> lightIndexBuffer;

// Buffer containing light counts
@property (nonatomic, readonly) id <MTLBuffer> lightCountBuffer;

// Buffer containing max volume density
@property (nonatomic, readonly) id <MTLBuffer> maxDensityBuffer;

// Number of light instances.
@property (nonatomic, readonly) NSUInteger lightCount;

// Number of lights.
@property (nonatomic, readonly) NSUInteger totalLightCount;

// Number of instances in the light buffer.
@property (nonatomic, readonly) NSUInteger instanceCount;

// Camera position vector.
@property (nonatomic) vector_float3 cameraPosition;

// Camera target vector. The camera faces this point.
@property (nonatomic) vector_float3 cameraTarget;

// Camera up vector.
@property (nonatomic) vector_float3 cameraUp;

// Camera field of view
@property (nonatomic) float cameraFov;

// Environment map texture
@property (nonatomic) id <MTLTexture> envmapTexture;

// Volume 3D density grid
@property (nonatomic) id <MTLTexture> volumeDensityGridTex;

// Initializer
- (instancetype)initWithDevice:(id <MTLDevice>)device;

#pragma mark - Create Scene
#pragma region CreateScene {
// Create scene with instances of a Cornell box. Each box can optionally
// contain a sphere primitive that uses an intersection function.
+ (RenderScene *)newInstancedMultipleCornellBoxSceneWithDevice:(id <MTLDevice>)device
                                useIntersectionFunctions:(BOOL)useIntersectionFunctions;

// Create one Cornell box for test.
+ (RenderScene *)newInstancedCornellBoxSceneWithDevice:(id <MTLDevice>)device
                        useIntersectionFunctions:(BOOL)useIntersectionFunctions;

// Create test scene from darts.
+ (RenderScene *)newTestScene:(id <MTLDevice>)device;

// Create test scene from Obj.
+ (RenderScene *)newTestSceneObj:(id <MTLDevice>)device;

// Create test scene for MIS Implementation
+ (RenderScene *)newTestSceneMIS:(id <MTLDevice>)device;

// Create test scene for EnvMap
+ (RenderScene *)newTestSceneEnv:(id <MTLDevice>)device;

// Create test scene for Homogeneous Media
+ (RenderScene *)newTestSceneVolHomo:(id <MTLDevice>)device;

// Create test scene for Hetreogeneous Media (VDB File)
+ (RenderScene *)newTestSceneVolHetero:(id <MTLDevice>)device;

// Create test scene for Hetreogeneous Media (VDB File)
+ (RenderScene *)newTestSceneVolHeteroBunny:(id <MTLDevice>)device;

// Create test scene for Hetreogeneous Media (VDB File)
+ (RenderScene *)newTestSceneVolHeteroForResearch:(id <MTLDevice>)device;

#pragma endregion CreateScene }

// Add a piece of geometry to the scene.
- (void)addGeometry:(Geometry *)mesh;

// Add an instance of a piece of geometry to the scene.
- (void)addInstance:(GeometryInstance *)instance;

// Add a light to the scene.
- (void)addLight:(AreaLight)light;

// Remove all geometry, instances, and lights from the scene.
- (void)clear;

// Upload all scene data to Metal buffers so the GPU can access the data.
- (void)uploadToBuffers;

@end

#endif
