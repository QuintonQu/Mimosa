/*
See LICENSE folder for this sample’s licensing information.

Abstract:
The implementation of the class that describes objects in a scene.
*/

#import "Scene.h"

#import <vector>
#import <set>
#import <iostream>

#import <ModelIO/ModelIO.h>

using namespace simd;

MTLResourceOptions getManagedBufferStorageMode() {
#if !TARGET_OS_IPHONE
    return MTLResourceStorageModeManaged;
#else
    return MTLResourceStorageModeShared;
#endif
}

#pragma mark - Geometry
@implementation Geometry

- (instancetype)initWithDevice:(id <MTLDevice>)device
{
    self = [super init];

    if (self) {
        _device = device;
    }

    return self;
}

- (void)uploadToBuffers {
}

- (void)clear {
}

- (MTLAccelerationStructureGeometryDescriptor *)geometryDescriptor {
    return nil;
}

- (NSArray <id <MTLResource>> *)resources {
    return @[];
}

- (NSString *)intersectionFunctionName {
    return nil;
}

@end

float3 getTriangleNormal(float3 v0, float3 v1, float3 v2) {
    float3 e1 = normalize(v1 - v0);
    float3 e2 = normalize(v2 - v0);

    return cross(e1, e2);
}

#pragma mark - Specific Geometry Triangle
@implementation TriangleGeometry {
    id <MTLBuffer> _indexBuffer;
    id <MTLBuffer> _vertexPositionBuffer;
    id <MTLBuffer> _vertexNormalBuffer;
    id <MTLBuffer> _vertexColorBuffer;
    id <MTLBuffer> _perPrimitiveDataBuffer;

    std::vector<uint32_t> _indices;
    std::vector<vector_float3> _vertices;
    std::vector<vector_float3> _normals;
    std::vector<vector_float3> _colors;
    std::vector<Triangle> _triangles;
};

- (NSUInteger)geometryCount{
    return _triangles.size();
}

- (void)uploadToBuffers {
    MTLResourceOptions options = getManagedBufferStorageMode();

    id <MTLDevice> device = self.device;

    _indexBuffer = [device newBufferWithLength:_indices.size() * sizeof(uint32_t) options:options];
    _vertexPositionBuffer = [device newBufferWithLength:_vertices.size() * sizeof(vector_float3) options:options];
    _vertexNormalBuffer = [device newBufferWithLength:_normals.size() * sizeof(vector_float3) options:options];
    _vertexColorBuffer = [device newBufferWithLength:_colors.size() * sizeof(vector_float3) options:options];
    _perPrimitiveDataBuffer = [device newBufferWithLength:_triangles.size() * sizeof(Triangle) options:options];

    memcpy(_indexBuffer.contents, _indices.data(), _indexBuffer.length);
    memcpy(_vertexPositionBuffer.contents, _vertices.data(), _vertexPositionBuffer.length);
    memcpy(_vertexNormalBuffer.contents, _normals.data(), _vertexNormalBuffer.length);
    memcpy(_vertexColorBuffer.contents, _colors.data(), _vertexColorBuffer.length);
    memcpy(_perPrimitiveDataBuffer.contents, _triangles.data(), _perPrimitiveDataBuffer.length);

#if !TARGET_OS_IPHONE
    [_indexBuffer didModifyRange:NSMakeRange(0, _indexBuffer.length)];
    [_vertexPositionBuffer didModifyRange:NSMakeRange(0, _vertexPositionBuffer.length)];
    [_vertexNormalBuffer didModifyRange:NSMakeRange(0, _vertexNormalBuffer.length)];
    [_vertexColorBuffer didModifyRange:NSMakeRange(0, _vertexColorBuffer.length)];
    [_perPrimitiveDataBuffer didModifyRange:NSMakeRange(0, _perPrimitiveDataBuffer.length)];
#endif
}

- (void)clear {
    _indices.clear();
    _vertices.clear();
    _normals.clear();
    _colors.clear();
    _triangles.clear();
}

- (void)addCubeFaceWithCubeVertices:(float3 *)cubeVertices
                              color:(float3)color
                                 i0:(uint32_t)i0
                                 i1:(uint32_t)i1
                                 i2:(uint32_t)i2
                                 i3:(uint32_t)i3
                      inwardNormals:(bool)inwardNormals
                           material:(Material)material
{
    const float3 v0 = cubeVertices[i0];
    const float3 v1 = cubeVertices[i1];
    const float3 v2 = cubeVertices[i2];
    const float3 v3 = cubeVertices[i3];

    float3 n0 = getTriangleNormal(v0, v1, v2);
    float3 n1 = getTriangleNormal(v0, v2, v3);

    if (inwardNormals) {
        n0 = -n0;
        n1 = -n1;
    }
    
    const size_t firstIndex = _indices.size();

    const uint32_t baseIndex = (uint32_t)_vertices.size();
    _indices.push_back(baseIndex + 0);
    _indices.push_back(baseIndex + 1);
    _indices.push_back(baseIndex + 2);
    _indices.push_back(baseIndex + 0);
    _indices.push_back(baseIndex + 2);
    _indices.push_back(baseIndex + 3);

    _vertices.push_back(v0);
    _vertices.push_back(v1);
    _vertices.push_back(v2);
    _vertices.push_back(v3);
    
    _normals.push_back(normalize(n0 + n1));
    _normals.push_back(n0);
    _normals.push_back(normalize(n0 + n1));
    _normals.push_back(n1);

    for (int i = 0; i < 4; i++)
        _colors.push_back(color);

    for (size_t triangleIndex = 0; triangleIndex < 2; triangleIndex++)
    {
        Triangle triangle;
        for (size_t i = 0; i < 3; i++)
        {
            const uint32_t index = _indices[firstIndex + triangleIndex * 3 + i];
            triangle.normals[i] = _normals[index];
            triangle.colors[i] = _colors[index];
        }
        triangle.material = material;
        _triangles.push_back(triangle);
    }
}

- (void)addXYPlane:(matrix_float4x4)transform
     inwardNormals:(bool)inwardNormals
          material:(Material)material
{
    float3 v0 = (transform * vector4(0.5f, 0.5f, 0.0f, 1.0f)).xyz;
    float3 v1 = (transform * vector4(-0.5f, 0.5f, 0.0f, 1.0f)).xyz;
    float3 v2 = (transform * vector4(-0.5f, -0.5f, 0.0f, 1.0f)).xyz;
    float3 v3 = (transform * vector4(0.5f, -0.5f, 0.0f, 1.0f)).xyz;
    

    float3 n0 = getTriangleNormal(v0, v1, v2);
    float3 n1 = getTriangleNormal(v0, v2, v3);

    if (inwardNormals) {
        n0 = -n0;
        n1 = -n1;
    }
    
    const size_t firstIndex = _indices.size();

    const uint32_t baseIndex = (uint32_t)_vertices.size();
    _indices.push_back(baseIndex + 0);
    _indices.push_back(baseIndex + 1);
    _indices.push_back(baseIndex + 2);
    _indices.push_back(baseIndex + 0);
    _indices.push_back(baseIndex + 2);
    _indices.push_back(baseIndex + 3);

    _vertices.push_back(v0);
    _vertices.push_back(v1);
    _vertices.push_back(v2);
    _vertices.push_back(v3);
    
    _normals.push_back(normalize(n0 + n1));
    _normals.push_back(n0);
    _normals.push_back(normalize(n0 + n1));
    _normals.push_back(n1);

    for (int i = 0; i < 4; i++)
        _colors.push_back(vector3(0.0f, 0.0f, 0.0f));

    for (size_t triangleIndex = 0; triangleIndex < 2; triangleIndex ++)
    {
        Triangle triangle;
        for (size_t i = 0; i < 3; i++)
        {
            const uint32_t index = _indices[firstIndex + triangleIndex * 3 + i];
            triangle.normals[i] = _normals[index];
            triangle.colors[i] = _colors[index];
        }
        triangle.material = material;
        _triangles.push_back(triangle);
    }
}

- (void)addCubeWithFaces:(unsigned int)faceMask
                   color:(vector_float3)color
               transform:(matrix_float4x4)transform
           inwardNormals:(bool)inwardNormals
                material:(Material)material
{
    float3 cubeVertices[] = {
        vector3(-0.5f, -0.5f, -0.5f),
        vector3( 0.5f, -0.5f, -0.5f),
        vector3(-0.5f,  0.5f, -0.5f),
        vector3( 0.5f,  0.5f, -0.5f),
        vector3(-0.5f, -0.5f,  0.5f),
        vector3( 0.5f, -0.5f,  0.5f),
        vector3(-0.5f,  0.5f,  0.5f),
        vector3( 0.5f,  0.5f,  0.5f),
    };

    for (int i = 0; i < 8; i++) {
        float3 vertex = cubeVertices[i];

        float4 transformedVertex = vector4(vertex.x, vertex.y, vertex.z, 1.0f);
        transformedVertex = transform * transformedVertex;

        cubeVertices[i] = transformedVertex.xyz;
    }

    uint32_t cubeIndices[][4] = {
        { 0, 4, 6, 2 },
        { 1, 3, 7, 5 },
        { 0, 1, 5, 4 },
        { 2, 6, 7, 3 },
        { 0, 2, 3, 1 },
        { 4, 5, 7, 6 }
    };

    for (unsigned face = 0; face < 6; face++) {
        if (faceMask & (1 << face)) {
            [self addCubeFaceWithCubeVertices:cubeVertices
                                        color:color
                                           i0:cubeIndices[face][0]
                                           i1:cubeIndices[face][1]
                                           i2:cubeIndices[face][2]
                                           i3:cubeIndices[face][3]
                                inwardNormals:inwardNormals
                                     material:(Material)material];
        }
    }
}

- (void)addGeometryWithURL:(NSURL *)URL
                 transform:(matrix_float4x4)transform
                  material:(Material)material
{
    MDLAsset *asset = [[MDLAsset alloc] initWithURL:URL];
    
    NSAssert(asset, @"Could not open %@", URL);
    
    MDLMesh *mesh = (MDLMesh *)asset[0];
    
    // Warning: Need to change with or without uv!
    struct MeshVertex {
        float position[3];
        float normal[3];
//        float uv[2];
    };
    
    MeshVertex *meshVertices = (MeshVertex *)mesh.vertexBuffers[0].map.bytes;
    std::set<MeshVertex> mesh_vertices;
    
    const size_t firstIndex = _indices.size();
    const uint32_t baseIndex = (uint32_t)_vertices.size();
    
    for (MDLSubmesh *submesh in mesh.submeshes) {
        uint32_t *indices = (uint32_t *)submesh.indexBuffer.map.bytes;
        
        vector_float3 color = [submesh.material propertyWithSemantic:MDLMaterialSemanticBaseColor].float3Value;
        
//        for (NSUInteger i = 0; i < submesh.indexCount; i++) {
//            uint32_t index = indices[i];
////            _indices.push_back(baseIndex + i);
//
//            MeshVertex vertex = meshVertices[index];
//            mesh_vertices.insert(vertex);
//        }
//
//        std::vector<MeshVertex> v(mesh_vertices.begin(), mesh_vertices.end());
//
//        for(MeshVertex vertex : v){
//            _vertices.push_back((transform * vector4(vertex.position[0], vertex.position[1], vertex.position[2], 1.0f)).xyz);
//            _normals.push_back(normalize((transpose(inverse(transform)) * vector4(vertex.normal[0], vertex.normal[1], vertex.normal[2], 0.0f)).xyz));
//            _colors.push_back(color);
//        }
        
        for (size_t triangle_index = 0; triangle_index < submesh.indexCount / 3; triangle_index++) {
            Triangle triangle;
            for (size_t i = 0; i < 3; i++)
            {
                const uint32_t index = indices[triangle_index * 3 + i];
                MeshVertex & vertex = meshVertices[index];

                _indices.push_back(baseIndex + (uint32_t)triangle_index * 3 + (uint32_t)i);
                _vertices.push_back((transform * vector4(vertex.position[0], vertex.position[1], vertex.position[2], 1.0f)).xyz);
                _normals.push_back(normalize((transpose(inverse(transform)) * vector4(vertex.normal[0], vertex.normal[1], vertex.normal[2], 0.0f)).xyz));
                _colors.push_back(color);

                triangle.normals[i] = _normals[_indices[firstIndex + triangle_index * 3 + i]];
                triangle.colors[i] = _colors[_indices[firstIndex + triangle_index * 3 + i]];
            }
            triangle.material = material;
            _triangles.push_back(triangle);
        }
    }
}

- (MTLAccelerationStructureGeometryDescriptor *)geometryDescriptor {
    // Metal represents each piece of piece of geometry in an acceleration structure using
    // a geometry descriptor. The sample uses a triangle geometry descriptor to represent
    // triangle geometry. Each triangle geometry descriptor can have its own
    // vertex buffer, index buffer, and triangle count. The sample use a single geometry
    // descriptor because it already packed all of the vertex data into a single buffer.
    MTLAccelerationStructureTriangleGeometryDescriptor *descriptor = [MTLAccelerationStructureTriangleGeometryDescriptor descriptor];

    descriptor.indexBuffer = _indexBuffer;
    descriptor.indexType = MTLIndexTypeUInt32;

    descriptor.vertexBuffer = _vertexPositionBuffer;
    descriptor.vertexStride = sizeof(float3);
    descriptor.triangleCount = _indices.size() / 3;

#if SUPPORTS_METAL_3
    if (@available(iOS 16, macOS 13, *)) {
        descriptor.primitiveDataBuffer = _perPrimitiveDataBuffer;
        descriptor.primitiveDataStride = sizeof(Triangle);
        descriptor.primitiveDataElementSize = sizeof(Triangle);
    }
#endif

    return descriptor;
}

- (NSArray <id <MTLResource>> *)resources {
    // The colors and normals for the vertices.
    return @[ _indexBuffer, _vertexNormalBuffer, _vertexColorBuffer ];
}

@end

#pragma mark - Specific Geometry Sphere
@implementation SphereGeometry {
    id <MTLBuffer> _sphereBuffer;
    id <MTLBuffer> _boundingBoxBuffer;
    id <MTLBuffer> _perPrimitiveDataBuffer;

    std::vector<Sphere> _spheres;
};

- (NSUInteger)geometryCount{
    return _spheres.size();
}

- (void)uploadToBuffers {
    MTLResourceOptions options = getManagedBufferStorageMode();

    id <MTLDevice> device = self.device;

    _sphereBuffer = [device newBufferWithLength:_spheres.size() * sizeof(Sphere) options:options];
    _boundingBoxBuffer = [device newBufferWithLength:_spheres.size() * sizeof(BoundingBox) options:options];

    std::vector<BoundingBox> boundingBoxes;

    // Geometry types that use custom intersection functions provide bounding boxes that enclose
    // each primitive. Metal invokes the intersection function whenever a ray potentially intersects
    // one of these bounding boxes.
    for (Sphere & sphere : _spheres) {
        BoundingBox bounds;

        bounds.min.x = sphere.origin.x - sphere.radius;
        bounds.min.y = sphere.origin.y - sphere.radius;
        bounds.min.z = sphere.origin.z - sphere.radius;

        bounds.max.x = sphere.origin.x + sphere.radius;
        bounds.max.y = sphere.origin.y + sphere.radius;
        bounds.max.z = sphere.origin.z + sphere.radius;

        boundingBoxes.push_back(bounds);
    }

    memcpy(_sphereBuffer.contents, _spheres.data(), _sphereBuffer.length);
    memcpy(_boundingBoxBuffer.contents, boundingBoxes.data(), _boundingBoxBuffer.length);

#if !TARGET_OS_IPHONE
    [_sphereBuffer didModifyRange:NSMakeRange(0, _sphereBuffer.length)];
    [_boundingBoxBuffer didModifyRange:NSMakeRange(0, _boundingBoxBuffer.length)];
#endif
}

- (void)clear {
    _spheres.clear();
}

- (void)addSphereWithOrigin:(vector_float3)origin
                     radius:(float)radius
                      color:(vector_float3)color
                    material:(Material)material
{
    Sphere sphere;

    sphere.origin = origin;
    sphere.radiusSquared = radius * radius;
    sphere.color = color;
    sphere.radius = radius;
    sphere.material = material;

    _spheres.push_back(sphere);
}

- (MTLAccelerationStructureGeometryDescriptor *)geometryDescriptor {
    // Metal represents each piece of geometry in an acceleration structure using
    // a geometry descriptor. The sample uses a bounding box geometry descriptor to
    // represent a custom primitive type.

    MTLAccelerationStructureBoundingBoxGeometryDescriptor *descriptor = [MTLAccelerationStructureBoundingBoxGeometryDescriptor descriptor];

    descriptor.boundingBoxBuffer = _boundingBoxBuffer;
    descriptor.boundingBoxCount = _spheres.size();

#if SUPPORTS_METAL_3
    if (@available(iOS 16, macOS 13, *)) {
        descriptor.primitiveDataBuffer = _sphereBuffer;
        descriptor.primitiveDataStride = sizeof(Sphere);
        descriptor.primitiveDataElementSize = sizeof(Sphere);
    }
#endif

    return descriptor;
}

- (NSArray <id <MTLResource>> *)resources {
    // The sphere intersection function uses the sphere origins and radii to check for
    // intersection with rays.
    return @[ _sphereBuffer ];
}

- (NSString *)intersectionFunctionName {
    return @"sphereIntersectionFunction";
}

@end

#pragma mark - GeometryInstance
@implementation GeometryInstance

- (instancetype)initWithGeometry:(Geometry *)geometry
                       transform:(matrix_float4x4)transform
                            mask:(unsigned int)mask
{
    self = [super init];

    if (self) {
        _geometry = geometry;
        _transform = transform;
        _mask = mask;
    }

    return self;
}

@end

#pragma mark - Scene
@implementation Scene {
    NSMutableArray <Geometry *> *_geometries;
    NSMutableArray <GeometryInstance *> *_instances;

    std::vector<AreaLight> _lights;
    std::vector<unsigned int> _light_indexs;
    std::vector<unsigned int> _light_counts; //_light_indexs.size() = _light_counts.size()
}

- (NSArray <Geometry *> *)geometries {
    return _geometries;
}

- (NSUInteger)lightCount {
    return (NSUInteger)_light_indexs.size();

}

- (NSUInteger)instanceCount {
    return _instances.count;
}

- (instancetype)initWithDevice:(id<MTLDevice>)device {
    self = [super init];

    if (self) {
        _device = device;

        _geometries = [[NSMutableArray alloc] init];
        _instances = [[NSMutableArray alloc] init];

        _cameraPosition = vector3(0.0f, 0.0f, -1.0f);
        _cameraTarget = vector3(0.0f, 0.0f, 0.0f);
        _cameraUp = vector3(0.0f, 1.0f, 0.0f);
    }

    return self;
}

- (void)clear {
    [_geometries removeAllObjects];
    [_instances removeAllObjects];

    _lights.clear();
    _light_indexs.clear();
    _light_counts.clear();
}

- (void)addGeometry:(Geometry *)mesh {
    [_geometries addObject:mesh];
}

- (void)addInstance:(GeometryInstance *)instance {
    [_instances addObject:instance];
    if (instance.mask & GEOMETRY_MASK_LIGHT){
        _light_indexs.push_back((unsigned int)(_instances.count-1));
        Geometry* geometry = instance.geometry;
        _light_counts.push_back((unsigned int)geometry.geometryCount);
    }
}

- (void)addLight:(AreaLight)light {
    _lights.push_back(light);
}

- (void)uploadToBuffers {
    for (Geometry *geometry in _geometries)
        [geometry uploadToBuffers];

    MTLResourceOptions options = getManagedBufferStorageMode();

    _lightBuffer = [_device newBufferWithLength:_lights.size() * sizeof(AreaLight) options:options];
    
    _lightIndexBuffer = [_device newBufferWithLength:_light_indexs.size() * sizeof(unsigned int) options:options];
    
    _lightCountBuffer = [_device newBufferWithLength:_light_counts.size() * sizeof(unsigned int) options:options];

    memcpy(_lightBuffer.contents, &_lights[0], _lightBuffer.length);
    memcpy(_lightIndexBuffer.contents, &_light_indexs[0], _lightIndexBuffer.length);
    memcpy(_lightCountBuffer.contents, &_light_counts[0], _lightCountBuffer.length);

#if !TARGET_OS_IPHONE
    [_lightBuffer didModifyRange:NSMakeRange(0, _lightBuffer.length)];
    [_lightIndexBuffer didModifyRange:NSMakeRange(0, _lightIndexBuffer.length)];
    [_lightCountBuffer didModifyRange:NSMakeRange(0, _lightCountBuffer.length)];
#endif
}

#pragma mark - Create Scene
+ (Scene *)newInstancedMultipleCornellBoxSceneWithDevice:(id <MTLDevice>)device
                        useIntersectionFunctions:(BOOL)useIntersectionFunctions
{
    Scene *scene = [[Scene alloc] initWithDevice:device];

    // Set up the camera.
    scene.cameraPosition = vector3(0.0f, 1.0f, 10.0f);
    scene.cameraTarget = vector3(0.0f, 1.0f, 0.0f);
    scene.cameraUp = vector3(0.0f, 1.0f, 0.0f);
    
    // Add a default material
    Material* default_material = new Material();
    default_material->color = vector3(0.8f, 0.8f, 0.8f);

    // Create a piece of triangle geometry for the light source.
    TriangleGeometry *lightMesh = [[TriangleGeometry alloc] initWithDevice:device];

    [scene addGeometry:lightMesh];

    matrix_float4x4 transform = matrix4x4_translation(0.0f, 1.0f, 0.0f) * matrix4x4_scale(0.5f, 1.98f, 0.5f);

    // Add the light source.
    [lightMesh addCubeWithFaces:FACE_MASK_POSITIVE_Y
                          color:vector3(1.0f, 1.0f, 1.0f)
                      transform:transform
                  inwardNormals:true
                    material:*default_material];

    // Create a piece of triangle geometry for the Cornell box.
    TriangleGeometry *geometryMesh = [[TriangleGeometry alloc] initWithDevice:device];

    [scene addGeometry:geometryMesh];

    transform = matrix4x4_translation(0.0f, 1.0f, 0.0f) * matrix4x4_scale(2.0f, 2.0f, 2.0f);

    // Add the top, bottom, and back walls.
    [geometryMesh addCubeWithFaces:FACE_MASK_NEGATIVE_Y | FACE_MASK_POSITIVE_Y | FACE_MASK_NEGATIVE_Z
                             color:vector3(0.725f, 0.71f, 0.68f)
                         transform:transform
                     inwardNormals:true
                          material:*default_material];

    // Add the left wall.
    [geometryMesh addCubeWithFaces:FACE_MASK_NEGATIVE_X
                             color:vector3(0.63f, 0.065f, 0.05f)
                         transform:transform
                     inwardNormals:true
                          material:*default_material];

    // Add the right wall.
    [geometryMesh addCubeWithFaces:FACE_MASK_POSITIVE_X
                             color:vector3(0.14f, 0.45f, 0.091f)
                         transform:transform
                     inwardNormals:true
                          material:*default_material];

    transform = matrix4x4_translation(-0.335f, 0.6f, -0.29f) *
                matrix4x4_rotation(0.3f, vector3(0.0f, 1.0f, 0.0f)) *
                matrix4x4_scale(0.6f, 1.2f, 0.6f);

    // Add the tall box.
    [geometryMesh addCubeWithFaces:FACE_MASK_ALL
                             color:vector3(0.725f, 0.71f, 0.68f)
                         transform:transform
                     inwardNormals:false
                          material:*default_material];

    SphereGeometry *sphereGeometry = nil;

    if (!useIntersectionFunctions) {
        transform = matrix4x4_translation(0.3275f, 0.3f, 0.3725f) *
                    matrix4x4_rotation(-0.3f, vector3(0.0f, 1.0f, 0.0f)) *
                    matrix4x4_scale(0.6f, 0.6f, 0.6f);

        // If the sample isn't using intersection functions, add the short box.
        [geometryMesh addCubeWithFaces:FACE_MASK_ALL
                                 color:vector3(0.725f, 0.71f, 0.68f)
                             transform:transform
                         inwardNormals:false
                              material:*default_material];
    }
    else {
        // Otherwise, create a piece of sphere geometry.
        sphereGeometry = [[SphereGeometry alloc] initWithDevice:device];

        [scene addGeometry:sphereGeometry];

        [sphereGeometry addSphereWithOrigin:vector3(0.3275f, 0.3f, 0.3725f)
                                     radius:0.3f
                                      color:vector3(0.725f, 0.71f, 0.68f)
                                   material:*default_material];
    }

    // Create nine instances of the scene.
    for (int y = -1; y <= 1; y++) {
        for (int x = -1; x <= 1; x++) {
            matrix_float4x4 transform = matrix4x4_translation(x * 2.5f, y * 2.5f, 0.0f);

            // Create an instance of the light.
            GeometryInstance *lightMeshInstance = [[GeometryInstance alloc] initWithGeometry:lightMesh
                                                                                   transform:transform
                                                                                        mask:GEOMETRY_MASK_LIGHT];

            [scene addInstance:lightMeshInstance];

            // Create an instance of the Cornell box.
            GeometryInstance *geometryMeshInstance = [[GeometryInstance alloc] initWithGeometry:geometryMesh
                                                                                      transform:transform
                                                                                           mask:GEOMETRY_MASK_TRIANGLE];

            [scene addInstance:geometryMeshInstance];

            // Create an instance of the sphere.
            if (useIntersectionFunctions) {
                GeometryInstance *sphereGeometryInstance = [[GeometryInstance alloc] initWithGeometry:sphereGeometry
                                                                                            transform:transform
                                                                                                 mask:GEOMETRY_MASK_SPHERE];

                [scene addInstance:sphereGeometryInstance];
            }

            // Add a light for each box.
            AreaLight light;

            light.position = vector3(x * 2.5f, y * 2.5f + 1.98f, 0.0f);
            light.forward = vector3(0.0f, -1.0f, 0.0f);
            light.right = vector3(0.25f, 0.0f, 0.0f);
            light.up = vector3(0.0f, 0.0f, 0.25f);

            float r = (float)rand() / (float)RAND_MAX;
            float g = (float)rand() / (float)RAND_MAX;
            float b = (float)rand() / (float)RAND_MAX;

            light.color = vector3(r * 4.0f, g * 4.0f, b * 4.0f);
            light.color = vector3(1.0f, 1.0f, 1.0f);

            [scene addLight:light];
        }
    }

    return scene;
}

+ (Scene *)newInstancedCornellBoxSceneWithDevice:(id <MTLDevice>)device
                        useIntersectionFunctions:(BOOL)useIntersectionFunctions
{
    Scene *scene = [[Scene alloc] initWithDevice:device];

    // Set up the camera.
    scene.cameraPosition = vector3(0.0f, 1.0f, 3.5f);
    scene.cameraTarget = vector3(0.0f, 1.0f, 0.0f);
    scene.cameraUp = vector3(0.0f, 1.0f, 0.0f);
    
    // Add a default material
    Material* default_material = new Material();
    default_material->color = vector3(0.8f, 0.8f, 0.8f);
    
    // Add multiple materials
    Material* left_wall_material = new Material();
    Material* right_wall_material = new Material();
    Material* tall_box_material = new Material();
    Material* sphere_material = new Material();
    Material* back_wall_material  = new Material();
    
    left_wall_material->color = vector3(0.63f, 0.065f, 0.05f);
    right_wall_material->color = vector3(0.14f, 0.45f, 0.091f);
    tall_box_material->color = vector3(1.0f, 1.0f, 1.0f);
    tall_box_material->is_metal = true;
    sphere_material->color = vector3(1.0f, 1.0f, 1.0f);
    sphere_material->is_glass = true;
    back_wall_material->color = vector3(0.4f, 0.3f, 0.68f);

    // Create a piece of triangle geometry for the light source.
    TriangleGeometry *lightMesh = [[TriangleGeometry alloc] initWithDevice:device];

    [scene addGeometry:lightMesh];

    matrix_float4x4 transform = matrix4x4_translation(0.0f, 1.0f, 0.0f) * matrix4x4_scale(0.5f, 1.98f, 0.5f);

    // Add the light source.
    [lightMesh addCubeWithFaces:FACE_MASK_POSITIVE_Y
                          color:vector3(1.0f, 1.0f, 1.0f)
                      transform:transform
                  inwardNormals:true
                    material:*default_material];

    // Create a piece of triangle geometry for the Cornell box.
    TriangleGeometry *geometryMesh = [[TriangleGeometry alloc] initWithDevice:device];

    [scene addGeometry:geometryMesh];

    transform = matrix4x4_translation(0.0f, 1.0f, 0.0f) * matrix4x4_scale(2.0f, 2.0f, 2.0f);

    // Add the top, bottom walls.
    [geometryMesh addCubeWithFaces:FACE_MASK_NEGATIVE_Y | FACE_MASK_POSITIVE_Y
                             color:vector3(0.725f, 0.71f, 0.68f)
                         transform:transform
                     inwardNormals:true
                          material:*default_material];
    
    // Add the back wall.
    [geometryMesh addCubeWithFaces:FACE_MASK_NEGATIVE_Z
                             color:vector3(0.725f, 0.71f, 0.68f)
                         transform:transform
                     inwardNormals:true
                          material:*back_wall_material];

    // Add the left wall.
    [geometryMesh addCubeWithFaces:FACE_MASK_NEGATIVE_X
                             color:vector3(0.63f, 0.065f, 0.05f)
                         transform:transform
                     inwardNormals:true
                          material:*left_wall_material];

    // Add the right wall.
    [geometryMesh addCubeWithFaces:FACE_MASK_POSITIVE_X
                             color:vector3(0.14f, 0.45f, 0.091f)
                         transform:transform
                     inwardNormals:true
                          material:*right_wall_material];

    transform = matrix4x4_translation(-0.335f, 0.6f, -0.29f) *
                matrix4x4_rotation(0.3f, vector3(0.0f, 1.0f, 0.0f)) *
                matrix4x4_scale(0.6f, 1.2f, 0.6f);

    // Add the tall box.
    [geometryMesh addCubeWithFaces:FACE_MASK_ALL
                             color:vector3(0.725f, 0.71f, 0.68f)
                         transform:transform
                     inwardNormals:false
                          material:*tall_box_material];

    SphereGeometry *sphereGeometry = nil;

    if (!useIntersectionFunctions) {
//        transform = matrix4x4_translation(0.3275f, 0.3f, 0.3725f) *
//                    matrix4x4_rotation(-0.3f, vector3(0.0f, 1.0f, 0.0f)) *
//                    matrix4x4_scale(0.6f, 0.6f, 0.6f);
//
//        // If the sample isn't using intersection functions, add the short box.
//        [geometryMesh addCubeWithFaces:FACE_MASK_ALL
//                                 color:vector3(0.725f, 0.71f, 0.68f)
//                             transform:transform
//                         inwardNormals:false
//                              material:*sphere_material];
    }
    else {
        // Otherwise, create a piece of sphere geometry.
        sphereGeometry = [[SphereGeometry alloc] initWithDevice:device];

        [scene addGeometry:sphereGeometry];

        [sphereGeometry addSphereWithOrigin:vector3(0.3275f, 0.3f, 0.3725f)
                                     radius:0.3f
                                      color:vector3(0.725f, 0.71f, 0.68f)
                                   material:*sphere_material];
    }

    transform = matrix4x4_translation(0.0f, 0.0f, 0.0f);

    // Create an instance of the light.
    GeometryInstance *lightMeshInstance = [[GeometryInstance alloc] initWithGeometry:lightMesh
                                                                           transform:transform
                                                                                mask:GEOMETRY_MASK_LIGHT];

    [scene addInstance:lightMeshInstance];

    // Create an instance of the Cornell box.
    GeometryInstance *geometryMeshInstance = [[GeometryInstance alloc] initWithGeometry:geometryMesh
                                                                              transform:transform
                                                                                   mask:GEOMETRY_MASK_TRIANGLE];

    [scene addInstance:geometryMeshInstance];

    // Create an instance of the sphere.
    if (useIntersectionFunctions) {
        GeometryInstance *sphereGeometryInstance = [[GeometryInstance alloc] initWithGeometry:sphereGeometry
                                                                                    transform:transform
                                                                                         mask:GEOMETRY_MASK_SPHERE];

        [scene addInstance:sphereGeometryInstance];
    }

    // Add a light for each box.
    AreaLight light;

    light.position = vector3(0.0f, 1.98f, 0.0f);
    light.forward = vector3(0.0f, -1.0f, 0.0f);
    light.right = vector3(0.25f, 0.0f, 0.0f);
    light.up = vector3(0.0f, 0.0f, 0.25f);

//    float r = (float)rand() / (float)RAND_MAX;
//    float g = (float)rand() / (float)RAND_MAX;
//    float b = (float)rand() / (float)RAND_MAX;
//
//    light.color = vector3(r * 4.0f, g * 4.0f, b * 4.0f);
    light.color = vector3(1.0f, 1.0f, 1.0f);

    [scene addLight:light];

    return scene;
}

+ (Scene *)newTestScene:(id<MTLDevice>)device
{
    Scene *scene = [[Scene alloc] initWithDevice:device];

    // Set up the camera.
    scene.cameraPosition = vector3(0.0f, 1.75f, 6.0f);
    scene.cameraTarget = vector3(0.0f, 1.0f, 0.0f);
    scene.cameraUp = vector3(0.0f, 1.0f, 0.0f);
    
    // Add a default material
    Material* default_material = new Material();
    default_material->color = vector3(0.8f, 0.8f, 0.8f);
    
    Material* m_light = new Material();
    m_light->color = vector3(20.f, 20.f, 20.f);
    
    Material* m_white = new Material();
    m_white->color = vector3(0.6f, 0.6f, 0.6f);
    
    Material* m_metal = new Material();
    m_metal->is_metal = true;
    m_metal->color = vector3(0.6f, 0.6f, 0.6f);
    
    Material* m_green = new Material();
    m_green->color = vector3(0.4f, 0.6f, 0.4f);
    
    Material* m_glass = new Material();
    m_glass->is_glass = true;
    m_glass->color = vector3(1.0f, 1.0f, 1.0f);

    // Create a piece of triangle geometry for the light source.
    TriangleGeometry *lightMesh = [[TriangleGeometry alloc] initWithDevice:device];

    [scene addGeometry:lightMesh];

    // Add the light source.
    
    [lightMesh addXYPlane:matrix4x4_rotation(65.f/180.f*M_PI, vector3(0.0f, 1.0f, 0.0f)) * matrix4x4_translation(4.0f, 3.5f, 0.0f) * matrix4x4_rotation(50.f/180.f*M_PI, vector3(0.0f, 0.0f, 1.0f)) * matrix4x4_rotation(-90.f/180.f*M_PI, vector3(0.0f, 1.0f, 0.0f))
               inwardNormals:true
                    material:*m_light];

    // Create a piece of triangle geometry for the Cornell box.
    TriangleGeometry *geometryMesh = [[TriangleGeometry alloc] initWithDevice:device];

    [scene addGeometry:geometryMesh];

    // Add quad
    [geometryMesh addXYPlane:matrix4x4_rotation(-90.f/180.f*M_PI, vector3(1.0f, 0.0f, 0.0f)) * matrix4x4_scale(100.0f, 10.0f, 0.0f)
               inwardNormals:false
                    material:*m_white];


    SphereGeometry *sphereGeometry = nil;


    // Otherwise, create a piece of sphere geometry.
    sphereGeometry = [[SphereGeometry alloc] initWithDevice:device];

    [scene addGeometry:sphereGeometry];

    [sphereGeometry addSphereWithOrigin:vector3(0.5f, 1.0f, 0.0f)
                                 radius:1.0f
                                  color:vector3(0.725f, 0.71f, 0.68f)
                               material:*m_glass];
    
    [sphereGeometry addSphereWithOrigin:vector3(2.2f, 1.0f, -4.0f)
                                 radius:1.0f
                                  color:vector3(0.725f, 0.71f, 0.68f)
                               material:*m_green];
    
    [sphereGeometry addSphereWithOrigin:vector3(-2.5f, 1.5f, -3.5f)
                                 radius:1.5f
                                  color:vector3(0.725f, 0.71f, 0.68f)
                               material:*m_metal];

    
    // Create instances
    matrix_float4x4 transform = matrix4x4_translation(0.0f, 0.0f, 0.0f);

    // Create an instance of the light.
    GeometryInstance *lightMeshInstance = [[GeometryInstance alloc] initWithGeometry:lightMesh
                                                                           transform:transform
                                                                                mask:GEOMETRY_MASK_LIGHT];

    [scene addInstance:lightMeshInstance];

    // Create an instance of the Cornell box.
    GeometryInstance *geometryMeshInstance = [[GeometryInstance alloc] initWithGeometry:geometryMesh
                                                                              transform:transform
                                                                                   mask:GEOMETRY_MASK_TRIANGLE];

    [scene addInstance:geometryMeshInstance];

    // Create an instance of the sphere.
    GeometryInstance *sphereGeometryInstance = [[GeometryInstance alloc] initWithGeometry:sphereGeometry
                                                                                transform:transform
                                                                                     mask:GEOMETRY_MASK_SPHERE];

    [scene addInstance:sphereGeometryInstance];


    // Add a light for each box.
    AreaLight light;

//    light.position = vector3(0.0f, 1.98f, 0.0f);
//    light.forward = vector3(0.0f, -1.0f, 0.0f);
//    light.right = vector3(0.25f, 0.0f, 0.0f);
//    light.up = vector3(0.0f, 0.0f, 0.25f);

//    float r = (float)rand() / (float)RAND_MAX;
//    float g = (float)rand() / (float)RAND_MAX;
//    float b = (float)rand() / (float)RAND_MAX;
//
//    light.color = vector3(r * 4.0f, g * 4.0f, b * 4.0f);
//    light.color = vector3(20.0f, 20.0f, 20.0f);

    [scene addLight:light];

    return scene;
}

+ (Scene *)newTestSceneObj:(id <MTLDevice>)device
{
    Scene *scene = [[Scene alloc] initWithDevice:device];

    // Set up the camera.
    scene.cameraPosition = vector3(0.0f, 1.0f, 3.5f);
    scene.cameraTarget = vector3(0.0f, 1.0f, 0.0f);
    scene.cameraUp = vector3(0.0f, 1.0f, 0.0f);
    scene.cameraFov = 45.0f;
    
    // Add a default material
    Material* default_material = new Material();
    default_material->color = vector3(0.8f, 0.8f, 0.8f);
    
    // Add multiple materials
    Material* left_wall_material = new Material();
    Material* right_wall_material = new Material();
    Material* tall_box_material = new Material();
    Material* sphere_material = new Material();
    Material* back_wall_material  = new Material();
    Material* light_material = new Material();
    Material* glass = new Material();
    
    left_wall_material->color = vector3(0.63f, 0.065f, 0.05f);
    right_wall_material->color = vector3(0.14f, 0.45f, 0.091f);
    tall_box_material->color = vector3(1.0f, 1.0f, 1.0f);
    tall_box_material->is_metal = true;
    sphere_material->color = vector3(1.0f, 1.0f, 1.0f);
    sphere_material->is_glass = true;
    back_wall_material->color = vector3(0.4f, 0.3f, 0.68f);
    light_material->color = vector3(20.f, 20.f, 20.f);
    glass->is_glass = true;
    glass->color = vector3(1.0f, 1.0f, 1.0f);

    // Create a piece of triangle geometry for the light source.
    TriangleGeometry *lightMesh = [[TriangleGeometry alloc] initWithDevice:device];

    [scene addGeometry:lightMesh];

    matrix_float4x4 transform = matrix4x4_translation(0.0f, 1.0f, 0.0f) * matrix4x4_scale(0.5f, 1.98f, 0.5f);

    // Add the light source.
    [lightMesh addCubeWithFaces:FACE_MASK_POSITIVE_Y
                          color:vector3(1.0f, 1.0f, 1.0f)
                      transform:transform
                  inwardNormals:true
                    material:*light_material];

    // Create a piece of triangle geometry for the Cornell box.
    TriangleGeometry *geometryMesh = [[TriangleGeometry alloc] initWithDevice:device];

    [scene addGeometry:geometryMesh];

    transform = matrix4x4_translation(0.0f, 1.0f, 0.0f) * matrix4x4_scale(2.0f, 2.0f, 2.0f);

    // Add the top, bottom walls.
    [geometryMesh addCubeWithFaces:FACE_MASK_NEGATIVE_Y | FACE_MASK_POSITIVE_Y
                             color:vector3(0.725f, 0.71f, 0.68f)
                         transform:transform
                     inwardNormals:true
                          material:*default_material];
    
    // Add the back wall.
    [geometryMesh addCubeWithFaces:FACE_MASK_NEGATIVE_Z
                             color:vector3(0.725f, 0.71f, 0.68f)
                         transform:transform
                     inwardNormals:true
                          material:*back_wall_material];

    // Add the left wall.
    [geometryMesh addCubeWithFaces:FACE_MASK_NEGATIVE_X
                             color:vector3(0.63f, 0.065f, 0.05f)
                         transform:transform
                     inwardNormals:true
                          material:*left_wall_material];

    // Add the right wall.
    [geometryMesh addCubeWithFaces:FACE_MASK_POSITIVE_X
                             color:vector3(0.14f, 0.45f, 0.091f)
                         transform:transform
                     inwardNormals:true
                          material:*right_wall_material];
    
    // Import the OBJ File
    NSURL *URL = [[NSBundle mainBundle] URLForResource:@"bunny-fine" withExtension:@"obj"];
    [geometryMesh addGeometryWithURL:URL
                           transform:transform * matrix4x4_scale(0.5f, 0.5f, 0.5f) * matrix4x4_translation(0.0f, -1.0f, 0.0f)
                            material:*glass];

//    transform = matrix4x4_translation(-0.335f, 0.6f, -0.29f) *
//                matrix4x4_rotation(0.3f, vector3(0.0f, 1.0f, 0.0f)) *
//                matrix4x4_scale(0.6f, 1.2f, 0.6f);
//
//    // Add the tall box.
//    [geometryMesh addCubeWithFaces:FACE_MASK_ALL
//                             color:vector3(0.725f, 0.71f, 0.68f)
//                         transform:transform
//                     inwardNormals:false
//                          material:*tall_box_material];

    SphereGeometry *sphereGeometry = nil;

    // Otherwise, create a piece of sphere geometry.
    sphereGeometry = [[SphereGeometry alloc] initWithDevice:device];

    [scene addGeometry:sphereGeometry];

    [sphereGeometry addSphereWithOrigin:vector3(0.3275f, 0.3f, 0.3725f)
                                 radius:0.3f
                                  color:vector3(0.725f, 0.71f, 0.68f)
                               material:*sphere_material];

    transform = matrix4x4_translation(0.0f, 0.0f, 0.0f);

    // Create an instance of the light.
    GeometryInstance *lightMeshInstance = [[GeometryInstance alloc] initWithGeometry:lightMesh
                                                                           transform:transform
                                                                                mask:GEOMETRY_MASK_LIGHT];

    [scene addInstance:lightMeshInstance];

    // Create an instance of the Cornell box.
    GeometryInstance *geometryMeshInstance = [[GeometryInstance alloc] initWithGeometry:geometryMesh
                                                                              transform:transform
                                                                                   mask:GEOMETRY_MASK_TRIANGLE];

    [scene addInstance:geometryMeshInstance];

//    // Create an instance of the sphere.
//    GeometryInstance *sphereGeometryInstance = [[GeometryInstance alloc] initWithGeometry:sphereGeometry
//                                                                                    transform:transform
//                                                                                         mask:GEOMETRY_MASK_SPHERE];
//
//    [scene addInstance:sphereGeometryInstance];

    // Add a light for each box.
    AreaLight light;
//
//    light.position = vector3(0.0f, 1.98f, 0.0f);
//    light.forward = vector3(0.0f, -1.0f, 0.0f);
//    light.right = vector3(0.25f, 0.0f, 0.0f);
//    light.up = vector3(0.0f, 0.0f, 0.25f);
//
//    float r = (float)rand() / (float)RAND_MAX;
//    float g = (float)rand() / (float)RAND_MAX;
//    float b = (float)rand() / (float)RAND_MAX;
//
//    light.color = vector3(r * 4.0f, g * 4.0f, b * 4.0f);
//    light.color = vector3(1.0f, 1.0f, 1.0f);
//
    [scene addLight:light];

    return scene;
}

+ (Scene *)newTestSceneMIS:(id<MTLDevice>)device
{
    Scene *scene = [[Scene alloc] initWithDevice:device];
    
    // Set up the camera
    scene.cameraPosition = vector3(0.f, 6.f, 27.5f);
    scene.cameraTarget = vector3(0.f, -1.5f, 2.5f);
    scene.cameraUp = vector3(0.0f, 1.0f, 0.0f);
    scene.cameraFov = 16.0f;
    
    // Add materials
    Material* diffuse_light_1 = new Material;
    diffuse_light_1->color = vector3(901.f, 901.f, 901.f);
    
    Material* diffuse_light_2 = new Material;
    diffuse_light_2->color = vector3(100.f, 100.f, 100.f);
    
    Material* diffuse_light_3 = new Material;
    diffuse_light_3->color = vector3(11.111f, 11.111f, 11.111f);
    
    Material* diffuse_light_4 = new Material;
    diffuse_light_4->color = vector3(4.23457f, 4.23457f, 4.23457f);
    
    Material* plate_1 = new Material;
    plate_1->color = vector3(0.35f, 0.35f, 0.35f);
    plate_1->exponent = 100000.f;
    plate_1->is_phong = true;
    
    Material* plate_2 = new Material;
    plate_2->color = vector3(0.25f, 0.25f, 0.25f);
    plate_2->exponent = 5000.0f;
    plate_2->is_phong = true;
    
    Material* plate_3 = new Material;
    plate_3->color = vector3(0.2f, 0.2f, 0.2f);
    plate_3->exponent = 400.0f;
    plate_3->is_phong = true;
    
    Material* plate_4 = new Material;
    plate_4->color = vector3(0.2f, 0.2f, 0.2f);
    plate_4->exponent = 100.0f;
    plate_4->is_phong = true;
    
    Material* diffuse = new Material;
    diffuse->color = vector3(0.2f, 0.2f, 0.2f);
    
    matrix_float4x4 non_transform = matrix4x4_translation(0.0f, 0.0f, 0.0f);
    
    SphereGeometry *lightMesh = [[SphereGeometry alloc] initWithDevice:device];

    [scene addGeometry:lightMesh];

    // Add the light sources.
    [lightMesh addSphereWithOrigin:vector3(3.75f, 0.0f, 0.0f)
                            radius:0.0333f
                             color:vector3(0.0f, 0.0f, 0.0f)
                          material:*diffuse_light_1];
    
    [lightMesh addSphereWithOrigin:vector3(1.25f, 0.0f, 0.0f)
                            radius:0.1f
                             color:vector3(0.0f, 0.0f, 0.0f)
                          material:*diffuse_light_2];
    
    [lightMesh addSphereWithOrigin:vector3(-1.25f, 0.0f, 0.0f)
                            radius:0.3f
                             color:vector3(0.0f, 0.0f, 0.0f)
                          material:*diffuse_light_3];
    
    [lightMesh addSphereWithOrigin:vector3(-3.75f, 0.0f, 0.0f)
                            radius:0.9f
                             color:vector3(0.0f, 0.0f, 0.0f)
                          material:*diffuse_light_4];
    
//    [lightMesh addSphereWithOrigin:vector3(-3.75f, 0.0f, 0.0f)
//                            radius:0.89f
//                             color:vector3(0.0f, 0.0f, 0.0f)
//                          material:*diffuse_light];
    
    // Add plate and floor.
    TriangleGeometry *geometryMesh = [[TriangleGeometry alloc] initWithDevice:device];

    [scene addGeometry:geometryMesh];
    
    NSURL* URL;
    
    URL = [[NSBundle mainBundle] URLForResource:@"plate1" withExtension:@"obj"];
    [geometryMesh addGeometryWithURL:URL
                           transform:non_transform
                            material:*plate_1];
    
    URL = [[NSBundle mainBundle] URLForResource:@"plate2" withExtension:@"obj"];
    [geometryMesh addGeometryWithURL:URL
                           transform:non_transform
                            material:*plate_2];

    URL = [[NSBundle mainBundle] URLForResource:@"plate3" withExtension:@"obj"];
    [geometryMesh addGeometryWithURL:URL
                           transform:non_transform
                            material:*plate_3];

    URL = [[NSBundle mainBundle] URLForResource:@"plate4" withExtension:@"obj"];
    [geometryMesh addGeometryWithURL:URL
                           transform:non_transform
                            material:*plate_4];

    URL = [[NSBundle mainBundle] URLForResource:@"floor" withExtension:@"obj"];
    [geometryMesh addGeometryWithURL:URL
                           transform:non_transform
                            material:*diffuse];
    
//    SphereGeometry *sphereGeometry = nil;

    // Otherwise, create a piece of sphere geometry.
//    sphereGeometry = [[SphereGeometry alloc] initWithDevice:device];
//
//    [scene addGeometry:sphereGeometry];
//
//    [sphereGeometry addSphereWithOrigin:vector3(0.3275f, 0.3f, 0.3725f)
//                                 radius:0.3f
//                                  color:vector3(0.725f, 0.71f, 0.68f)
//                               material:*diffuse];

    // Create an instance of the light.
    GeometryInstance *lightMeshInstance = [[GeometryInstance alloc] initWithGeometry:lightMesh
                                                                           transform:non_transform
                                                                                mask:GEOMETRY_MASK_SPHERE_LIGHT];

    [scene addInstance:lightMeshInstance];

    // Create an instance of the Cornell box.
    GeometryInstance *geometryMeshInstance = [[GeometryInstance alloc] initWithGeometry:geometryMesh
                                                                              transform:non_transform
                                                                                   mask:GEOMETRY_MASK_TRIANGLE];

    [scene addInstance:geometryMeshInstance];

    // Create an instance of the sphere.
//    GeometryInstance *sphereGeometryInstance = [[GeometryInstance alloc] initWithGeometry:sphereGeometry
//                                                                                    transform:non_transform
//                                                                                         mask:GEOMETRY_MASK_SPHERE];
//
//    [scene addInstance:sphereGeometryInstance];

    // Add a light for each box.
    AreaLight light;
//
//    light.position = vector3(0.0f, 1.98f, 0.0f);
//    light.forward = vector3(0.0f, -1.0f, 0.0f);
//    light.right = vector3(0.25f, 0.0f, 0.0f);
//    light.up = vector3(0.0f, 0.0f, 0.25f);
//
//    float r = (float)rand() / (float)RAND_MAX;
//    float g = (float)rand() / (float)RAND_MAX;
//    float b = (float)rand() / (float)RAND_MAX;
//
//    light.color = vector3(r * 4.0f, g * 4.0f, b * 4.0f);
//    light.color = vector3(1.0f, 1.0f, 1.0f);
//
    [scene addLight:light];

    return scene;
    
    return scene;
}

@end

