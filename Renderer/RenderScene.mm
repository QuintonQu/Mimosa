/*
See LICENSE folder for this sample’s licensing information.

Abstract:
The implementation of the class that describes objects in a scene.
*/

#import "RenderScene.h"
#import <vector>
#import <set>
#import <iostream>
#import <ModelIO/ModelIO.h>
#import <CoreGraphics/CoreGraphics.h>
#import <ImageIO/ImageIO.h>
#include <random>
#include <stdint.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/ValueTransformer.h>
#import <Foundation/Foundation.h>

using namespace simd;

MTLResourceOptions getManagedBufferStorageMode() {
#if !TARGET_OS_IPHONE
    return MTLResourceStorageModeManaged;
#else
    return MTLResourceStorageModeShared;
#endif
}

#pragma mark - Environment Map Help Functions

static CGImageRef createCGImageFromFile (NSString* path)
{
    // Get the URL for the pathname to pass it to
    // `CGImageSourceCreateWithURL`.
    NSURL *url = [NSURL fileURLWithPath:path];
    CGImageRef        myImage = NULL;
    CGImageSourceRef  myImageSource;
    CFDictionaryRef   myOptions = NULL;
    CFStringRef       myKeys[2];
    CFTypeRef         myValues[2];

    // Set up options if you want them. The options here are for
    // caching the image in a decoded form and for using floating-point
    // values if the image format supports them.
    myKeys[0] = kCGImageSourceShouldCache;
    myValues[0] = (CFTypeRef)kCFBooleanFalse;

    myKeys[1] = kCGImageSourceShouldAllowFloat;
    myValues[1] = (CFTypeRef)kCFBooleanTrue;

    // Create the dictionary.
    myOptions = CFDictionaryCreate(NULL,
                                   (const void **) myKeys,
                                   (const void **) myValues, 2,
                                   &kCFTypeDictionaryKeyCallBacks,
                                   & kCFTypeDictionaryValueCallBacks);

    // Create an image source from the URL.
    myImageSource = CGImageSourceCreateWithURL((CFURLRef)url, myOptions);
    CFRelease(myOptions);

    // Make sure the image source exists before continuing.
    if (myImageSource == NULL)
    {
        fprintf(stderr, "Image source is NULL.");
        return  NULL;
    }

    // Create an image from the first item in the image source.
    myImage = CGImageSourceCreateImageAtIndex(myImageSource, 0, NULL);
    CFRelease(myImageSource);

    // Make sure the image exists before continuing.
    if (myImage == NULL)
    {
         fprintf(stderr, "Image not created from image source.");
         return NULL;
    }

    return myImage;
}


id<MTLTexture> texture_from_radiance_file(NSString * fileName, id<MTLDevice> device, NSError ** error)
{
    // Validate the function inputs.

    if (![fileName containsString:@"."])
    {
        if (error != NULL)
        {
            *error = [[NSError alloc] initWithDomain:@"File load failure."
                                                code:0xdeadbeef
                                            userInfo:@{NSLocalizedDescriptionKey : @"No file extension provided."}];
        }
        return nil;
    }

    NSArray * subStrings = [fileName componentsSeparatedByString:@"."];

    if ([subStrings[1] compare:@"hdr"] != NSOrderedSame)
    {
        if (error != NULL)
        {
            *error = [[NSError alloc] initWithDomain:@"File load failure."
                                                code:0xdeadbeef
                                            userInfo:@{NSLocalizedDescriptionKey : @"Only (.hdr) files are supported."}];
        }
        return nil;
    }

    // Load and validate the image.

    NSString* filePath = [[NSBundle mainBundle] pathForResource:subStrings[0] ofType:subStrings[1]];
    CGImageRef loadedImage = createCGImageFromFile(filePath);

    if (loadedImage == NULL)
    {
        if (error != NULL)
        {
            *error = [[NSError alloc] initWithDomain:@"File load failure."
                                                code:0xdeadbeef
                                            userInfo:@{NSLocalizedDescriptionKey : @"Unable to create CGImage."}];
        }

        return nil;
    }

    size_t bpp = CGImageGetBitsPerPixel(loadedImage);

    size_t kSrcChannelCount;
    if (@available(macOS 14, *)) {
        kSrcChannelCount = 4;
    }else{
        kSrcChannelCount = 3;
    }
    
    const size_t kBitsPerByte = 8;
    const size_t kExpectedBitsPerPixel = sizeof(uint16_t) * kSrcChannelCount * kBitsPerByte;

    if (bpp != kExpectedBitsPerPixel)
    {
        if (error != NULL)
        {
            *error = [[NSError alloc] initWithDomain:@"File load failure."
                                                code:0xdeadbeef
                                            userInfo:@{NSLocalizedDescriptionKey : [NSString stringWithFormat:@"Expected %zu bits per pixel, but file returns %zu", kExpectedBitsPerPixel, bpp]}];
        }
        CFRelease(loadedImage);
        return nil;
    }

    // Copy the image into a tempory buffer.

    size_t width = CGImageGetWidth(loadedImage);
    size_t height = CGImageGetHeight(loadedImage);

    // Make the CG image data accessible.
    CFDataRef cgImageData = CGDataProviderCopyData(CGImageGetDataProvider(loadedImage));

    // Get a pointer to the data.
    const uint16_t * srcData = (const uint16_t * )CFDataGetBytePtr(cgImageData);

    // Metal exposes an RGBA16Float format, but the source data is RGB F16, so
    // this function adds an extra channel of padding.
    const size_t kPixelCount = width * height;
    const size_t kDstChannelCount = 4;
    const size_t kDstSize = kPixelCount * sizeof(uint16_t) * kDstChannelCount;

    uint16_t * dstData = (uint16_t *)malloc(kDstSize);

    for (size_t texIdx = 0; texIdx < kPixelCount; ++texIdx)
    {
        const uint16_t * currSrc = srcData + (texIdx * kSrcChannelCount);
        uint16_t * currDst = dstData + (texIdx * kDstChannelCount);

        currDst[0] = currSrc[0];
        currDst[1] = currSrc[1];
        currDst[2] = currSrc[2];
        currDst[3] = (uint16_t) 1.f;
    }

    // Create an `MTLTexture`.

    MTLTextureDescriptor * texDesc = [MTLTextureDescriptor new];

    texDesc.pixelFormat = MTLPixelFormatRGBA16Float;
    texDesc.width = width;
    texDesc.height = height;

    id<MTLTexture> texture = [device newTextureWithDescriptor:texDesc];

    const NSUInteger kBytesPerRow = sizeof(uint16_t) * kDstChannelCount * width;

    MTLRegion region = { {0,0,0}, {width, height, 1} };

    [texture replaceRegion:region mipmapLevel:0 withBytes:dstData bytesPerRow:kBytesPerRow];

    // Remember to clean things up.
    free(dstData);
    CFRelease(cgImageData);
    CFRelease(loadedImage);

    return texture;
}

#pragma mark - VDB density grid Help Functions
uint16_t floatToHalf(float value)
{
    uint32_t bits = *(uint32_t *)&value;
    uint32_t sign = (bits >> 31) & 0x1;
    int32_t exponent = ((bits >> 23) & 0xFF) - 127;
    uint32_t mantissa = bits & 0x7FFFFF;

    if(exponent == 128)
    {
        // NaN or Infinity
        exponent = 16;
        mantissa >>= 13;
    }
    else if(exponent > 15)
    {
        // Overflow
        return sign << 15 | 0x7C00;
    }
    else if(exponent > -15)
    {
        // Normalized number
        exponent += 15;
        mantissa >>= 13;
    }
    else if(exponent > -25)
    {
        // Subnormal number
        mantissa |= 0x800000;
        mantissa >>= -14 - exponent;
        exponent = -15 + 1;
    }
    else
    {
        // Underflow
        return sign << 15;
    }

    return sign << 15 | exponent << 10 | mantissa;
}

id<MTLTexture> createVolume(NSString * fileName, id<MTLDevice> device, float &_maxDensity, float3 &size, float scale, bool is_grad){
    openvdb::initialize();
    openvdb::io::File file([fileName UTF8String]);
    file.open();
    
    openvdb::GridBase::Ptr base_grid;
    std::string gridname = "";
    
    for (openvdb::io::File::NameIterator name_iter = file.beginName();
        name_iter != file.endName(); ++name_iter)
    {
        // Read in only the grid we are interested in.
        if (gridname == "" || name_iter.gridName() == gridname) {
            std::cout << "reading grid " << name_iter.gridName() << std::endl;
            base_grid = file.readGrid(name_iter.gridName());
            if (gridname == "")
                break;
        } else {
            std::cout << "skipping grid " << name_iter.gridName() << std::endl;
        }
    }
    std::cout << "vdb file reading done!" << std::endl;
    
    file.close();
    openvdb::FloatGrid::Ptr grid = openvdb::gridPtrCast<openvdb::FloatGrid>(base_grid);
    auto bbox = grid->evalActiveVoxelBoundingBox();
    openvdb::FloatGrid::Accessor accessor = grid->getAccessor();

    int width = bbox.max().x() - bbox.min().x();
    int height = bbox.max().y() - bbox.min().y();
    int depth = bbox.max().z() - bbox.min().z();
    
    size[0] = width;
    size[1] = height;
    size[2] = depth;

    uint16_t* values = (uint16_t*)malloc(sizeof(uint16_t) * width * height * depth * 4);
    int value_index = 0;
    float max_density = 0.f;
    for (int k = bbox.min().z(); k < bbox.max().z(); ++k) {
        for (int j = bbox.min().y(); j < bbox.max().y(); ++j) {
            for (int i = bbox.min().x(); i < bbox.max().x(); ++i) {
                // albedo
                if(is_grad){
                    values[value_index] = floatToHalf((float((k - bbox.min().z()))/(bbox.max().z()-bbox.min().z()) / 2 + 0.5));
                    values[value_index + 1] = floatToHalf((float((j - bbox.min().y()))/(bbox.max().y()-bbox.min().y()) / 4 + 0.75));
                    values[value_index + 2] = floatToHalf((float((i - bbox.min().x()))/(bbox.max().x()-bbox.min().x()) / 3 + 0.66));
                }else{
                    values[value_index] = floatToHalf(1.f);
                    values[value_index + 1] = floatToHalf(1.f);
                    values[value_index + 2] = floatToHalf(1.f);
                }
                
                // density
                float density = accessor.getValue(openvdb::Coord(i, j, k)) * scale;
                if(density > max_density) max_density = density;
                values[value_index + 3] = floatToHalf(density);
                
                // next voxel
                value_index += 4;
            }
        }
    }
    _maxDensity = max_density;
    std::cout << "vdb data extraction done!" << std::endl;
    
    MTLTextureDescriptor* textureDescriptor = [MTLTextureDescriptor new];
    textureDescriptor.pixelFormat = MTLPixelFormatRGBA16Float;
    textureDescriptor.textureType = MTLTextureType3D;
    textureDescriptor.width = width;
    textureDescriptor.height = height;
    textureDescriptor.depth = depth;
    
    id<MTLTexture> texture = [device newTextureWithDescriptor:textureDescriptor];
    [texture replaceRegion:MTLRegionMake3D(0, 0, 0, width, height, depth)
               mipmapLevel:    0
                     slice: 0
                 withBytes:values
               bytesPerRow:sizeof(uint16_t) * width * 4
             bytesPerImage:sizeof(uint16_t) * width * height * 4];
    
    free(values);
    return texture;
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
            triangle.positions[i] = _vertices[index];
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
            triangle.positions[i] = _vertices[index];
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
                triangle.positions[i] = _vertices[_indices[firstIndex + triangle_index * 3 + i]];
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
//    return @[ _indexBuffer, _vertexNormalBuffer, _vertexColorBuffer, _vertexPositionBuffer ];
    return @[ _perPrimitiveDataBuffer ];
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
@implementation RenderScene {
    NSMutableArray <Geometry *> *_geometries;
    NSMutableArray <GeometryInstance *> *_instances;

    std::vector<AreaLight> _lights;
    std::vector<unsigned int> _light_indexs;
    std::vector<unsigned int> _light_counts; //_light_indexs.size() = _light_counts.size()
    
    int _totalLightCount;
    float _maxDensity;
}

- (NSArray <Geometry *> *)geometries {
    return _geometries;
}

- (NSUInteger)lightCount {
    if(_light_indexs[0] == -1){
        return 0;
    }
    return (NSUInteger)_light_indexs.size();
}

- (NSUInteger)totalLightCount {
//    std::cout <<_totalLightCount << std::endl;
    return (NSUInteger)_totalLightCount;

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
        
        _totalLightCount = 0;
        _maxDensity = 0.f;
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
        _totalLightCount += (unsigned int)geometry.geometryCount;
    }
//    std::cout << _light_indexs[0] << ", " << _light_counts[0] << std::endl;
}

- (void)addLight:(AreaLight)light {
    _lights.push_back(light);
}

- (void)uploadToBuffers {
    for (Geometry *geometry in _geometries)
        [geometry uploadToBuffers];

    MTLResourceOptions options = getManagedBufferStorageMode();
    
    if(_lights.size() == 0) {
        AreaLight light;
        _lights.push_back(light);
    }
    
    if(_light_indexs.size() == 0){
        _light_indexs.push_back(-1);
        _light_counts.push_back(-1);
    }
    
    _lightBuffer = [_device newBufferWithLength:_lights.size() * sizeof(AreaLight) options:options];
    
    _lightIndexBuffer = [_device newBufferWithLength:_light_indexs.size() * sizeof(unsigned int) options:options];
    
    _lightCountBuffer = [_device newBufferWithLength:_light_counts.size() * sizeof(unsigned int) options:options];

    memcpy(_lightBuffer.contents, &_lights[0], _lightBuffer.length);
    memcpy(_lightIndexBuffer.contents, &_light_indexs[0], _lightIndexBuffer.length);
    memcpy(_lightCountBuffer.contents, &_light_counts[0], _lightCountBuffer.length);
    
    _maxDensityBuffer = [_device newBufferWithLength:sizeof(float) options:options];
    memcpy(_maxDensityBuffer.contents, &_maxDensity, _maxDensityBuffer.length);

#if !TARGET_OS_IPHONE
    [_lightBuffer didModifyRange:NSMakeRange(0, _lightBuffer.length)];
    [_lightIndexBuffer didModifyRange:NSMakeRange(0, _lightIndexBuffer.length)];
    [_lightCountBuffer didModifyRange:NSMakeRange(0, _lightCountBuffer.length)];
    [_maxDensityBuffer didModifyRange:NSMakeRange(0, _maxDensityBuffer.length)];
#endif
}

- (void)setConstantEnvmapTexture:(float3)background_color {
    MTLTextureDescriptor *textureDescriptor = [MTLTextureDescriptor texture2DDescriptorWithPixelFormat:MTLPixelFormatRGBA32Float width:512 height:512 mipmapped:NO];
    
    id<MTLTexture> texture = [_device newTextureWithDescriptor:textureDescriptor];
    
    float color[4] = {background_color.r, background_color.g, background_color.b, 1.0f};
    float *data = (float *)malloc(512 * 512 * 4 * sizeof(float));
    for (int i = 0; i < 512 * 512; i++)
    {
        memcpy(data + i * 4, color, 4 * sizeof(float));
    }

    MTLRegion region = MTLRegionMake2D(0, 0, 512, 512);
    [texture replaceRegion:region mipmapLevel:0 withBytes:data bytesPerRow:512 * 4 * sizeof(float)];

    [self setEnvmapTexture:texture];
    free(data);
}

#pragma mark - Create Scene (Default)
+ (RenderScene *)newInstancedMultipleCornellBoxSceneWithDevice:(id <MTLDevice>)device
                        useIntersectionFunctions:(BOOL)useIntersectionFunctions
{
    RenderScene *scene = [[RenderScene alloc] initWithDevice:device];

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

#pragma mark - Create Scene Cornell Box
+ (RenderScene *)newInstancedCornellBoxSceneWithDevice:(id <MTLDevice>)device
                        useIntersectionFunctions:(BOOL)useIntersectionFunctions
{
    RenderScene *scene = [[RenderScene alloc] initWithDevice:device];

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

#pragma mark - Create Scene Woj-example
+ (RenderScene *)newTestScene:(id<MTLDevice>)device
{
    RenderScene *scene = [[RenderScene alloc] initWithDevice:device];

    // Set up the camera.
    scene.cameraPosition = vector3(0.0f, 1.75f, 6.0f);
    scene.cameraTarget = vector3(0.0f, 1.0f, 0.0f);
    scene.cameraUp = vector3(0.0f, 1.0f, 0.0f);
    scene.cameraFov = 45.0f;
    
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
                                                                                mask:GEOMETRY_MASK_TRIANGLE_LIGHT];

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

    float3 background_color = vector3(0.0f, 0.0f, 0.0f);
    [scene setConstantEnvmapTexture:background_color];

    return scene;
}

#pragma mark - Create Scene Rabbit
+ (RenderScene *)newTestSceneObj:(id <MTLDevice>)device
{
    RenderScene *scene = [[RenderScene alloc] initWithDevice:device];

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
                                                                                mask:GEOMETRY_MASK_TRIANGLE_LIGHT];

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

    float3 background_color = vector3(0.0f, 0.0f, 0.0f);
    [scene setConstantEnvmapTexture:background_color];

    return scene;
}

#pragma mark - Create Scene MIS-Woj-example
+ (RenderScene *)newTestSceneMIS:(id<MTLDevice>)device
{
    RenderScene *scene = [[RenderScene alloc] initWithDevice:device];
    
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

    float3 background_color = vector3(0.0f, 0.0f, 0.0f);
    [scene setConstantEnvmapTexture:background_color];

    return scene;
}

#pragma mark - Create Scene Environment Map
+ (RenderScene *)newTestSceneEnv:(id<MTLDevice>)device
{
    RenderScene *scene = [[RenderScene alloc] initWithDevice:device];
    
    // Set up the camera
    scene.cameraPosition = vector3(13.f, 2.f, 3.f);
    scene.cameraTarget = vector3(0.f, 0.f, 0.f);
    scene.cameraUp = vector3(0.0f, 1.0f, 0.0f);
    scene.cameraFov = 20.0f;
    
    // Add materials
    Material* ground_material = new Material;
    ground_material->color = vector3(0.5f, 0.5f, 0.5f);
    
    Material* diffuse_light_4 = new Material;
    diffuse_light_4->color = vector3(4.23457f, 4.23457f, 4.23457f);
    
    matrix_float4x4 non_transform = matrix4x4_translation(0.0f, 0.0f, 0.0f);
    
    // Add floor.
    TriangleGeometry *geometryMesh = [[TriangleGeometry alloc] initWithDevice:device];

    [scene addGeometry:geometryMesh];
    
    [geometryMesh addXYPlane:matrix4x4_rotation(-90.f/180.f*M_PI, vector3(1.0f, 0.0f, 0.0f)) * matrix4x4_scale(100.0f, 100.0f, 0.0f)
               inwardNormals:false
                    material:*ground_material];
    
    // Add spheres    
    SphereGeometry *sphereGeometry = [[SphereGeometry alloc] initWithDevice:device];
    
    [scene addGeometry:sphereGeometry];
        
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);

    for (int a = -11; a < 11; a++)
    {
        for (int b = -11; b < 11; b++)
        {
            float choose_mat = (float) dis(gen);
            float r1 = (float) dis(gen);
            float r2 = (float) dis(gen);
            float3 origin = vector3(a+0.9f*r1, 0.2f, b+0.9f*r2);
            if(length(origin-vector3(4.0f, 0.2f, 0.0f)) > 0.9f){
                if(choose_mat < 0.8f){
                    // diffuse
                    float r1 = (float) dis(gen);
                    float r2 = (float) dis(gen);
                    float r3 = (float) dis(gen);
                    float r4 = (float) dis(gen);
                    float r5 = (float) dis(gen);
                    float r6 = (float) dis(gen);
                    Material* diffuse = new Material;
                    diffuse->color = vector3(r1 * r2, r3 * r4, r5 * r6);
                    [sphereGeometry addSphereWithOrigin:origin
                                                 radius:0.2f
                                                  color:vector3(0.725f, 0.71f, 0.68f)
                                               material:*diffuse];
                }else if(choose_mat < 0.95f){
                    float r1 = (float) dis(gen);
                    float r2 = (float) dis(gen);
                    float r3 = (float) dis(gen);
//                    float r4 = (float) dis(gen);
                    Material* metal = new Material;
                    metal->color = vector3(0.5f * (1 + r1), 0.5f * (1.0f + r2), 0.5f * (1.0f + r3));
                    metal->is_metal = true;
                    [sphereGeometry addSphereWithOrigin:origin
                                                 radius:0.2f
                                                  color:vector3(0.725f, 0.71f, 0.68f)
                                               material:*metal];
                }else{
                    Material* glass = new Material;
                    glass->color = vector3(1.f, 1.f, 1.f);
                    glass->is_glass = true;
                    [sphereGeometry addSphereWithOrigin:origin
                                                 radius:0.2f
                                                  color:vector3(0.725f, 0.71f, 0.68f)
                                               material:*glass];
                }
            }
        }
    }
    
    Material* big_glass = new Material;
    big_glass->color = vector3(1.f, 1.f, 1.f);
    big_glass->is_glass = true;
    
    Material* big_diffuse = new Material;
    big_diffuse->color = vector3(0.4f, 0.2f, 0.1f);
    
    Material* big_metal = new Material;
    big_metal->color = vector3(0.7f, 0.6f, 0.5f);
    big_metal->is_phong = true;
    big_metal->exponent = 1000.f;
    
    [sphereGeometry addSphereWithOrigin:vector3(0.0f, 1.0f, 0.0f)
                                 radius:1.0f
                                  color:vector3(0.725f, 0.71f, 0.68f)
                               material:*big_glass];
    
    [sphereGeometry addSphereWithOrigin:vector3(-4.0f, 1.0f, 0.0f)
                                 radius:1.0f
                                  color:vector3(0.725f, 0.71f, 0.68f)
                               material:*big_diffuse];
    
    [sphereGeometry addSphereWithOrigin:vector3(4.0f, 1.0f, 0.0f)
                                 radius:1.0f
                                  color:vector3(0.725f, 0.71f, 0.68f)
                               material:*big_metal];

    // Create an instance of the Cornell box.
    GeometryInstance *geometryMeshInstance = [[GeometryInstance alloc] initWithGeometry:geometryMesh
                                                                              transform:non_transform
                                                                                   mask:GEOMETRY_MASK_TRIANGLE];

    [scene addInstance:geometryMeshInstance];

    // Create an instance of the sphere.
    GeometryInstance *sphereGeometryInstance = [[GeometryInstance alloc] initWithGeometry:sphereGeometry
                                                                                    transform:non_transform
                                                                                         mask:GEOMETRY_MASK_SPHERE];

    [scene addInstance:sphereGeometryInstance];
    
    // Create environment map
    NSError *error;
    id<MTLTexture> env_map = texture_from_radiance_file( @"kloppenheim_06_4k.hdr", scene.device, &error );
    NSAssert( env_map, @"Could not load sky texture: %@", error );
    [scene setEnvmapTexture:env_map];

    return scene;
}

#pragma mark - Create Scene Homo-Vol
+ (RenderScene *)newTestSceneVolHomo:(id <MTLDevice>)device
{
    RenderScene *scene = [[RenderScene alloc] initWithDevice:device];

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
    back_wall_material->color = vector3(0.38f, 0.68f, 0.98f);
    
    back_wall_material->is_phong = true;
    back_wall_material->exponent = 500;
//    right_wall_material->is_metal = true;
//    left_wall_material->is_metal = true;
//    default_material->is_metal = true;
    
    light_material->color = vector3(20.f, 20.f, 20.f);
    glass->is_glass = true;
    glass->color = vector3(1.0f, 1.0f, 1.0f);
    glass->is_contain_volume = true;
    glass->density = 2.7f;
    glass->albedo = vector3(0.1f, 0.75f, 0.9f);
//    glass->emission = vector3(0.7f, 0.8f, 0.7f);

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
    
    TriangleGeometry *volumeMesh = [[TriangleGeometry alloc] initWithDevice:device];

    [scene addGeometry:volumeMesh];
    
//     Import the OBJ File
    NSURL *URL = [[NSBundle mainBundle] URLForResource:@"bunny-fine" withExtension:@"obj"];
    [volumeMesh addGeometryWithURL:URL
                           transform:transform * matrix4x4_scale(0.5f, 0.5f, 0.5f) * matrix4x4_translation(0.0f, -1.0f, 0.0f)
                            material:*glass];
    
    // Add a cube in the center
//    [volumeMesh addCubeWithFaces:FACE_MASK_ALL
//                             color:vector3(0.14f, 0.45f, 0.091f)
//                         transform:transform * matrix4x4_scale(0.3f, 0.3f, 0.3f)
//                     inwardNormals:true
//                          material:*glass];

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

    [sphereGeometry addSphereWithOrigin:vector3(0.0f, 1.0f, 0.0f)
                                 radius:0.4f
                                  color:vector3(0.725f, 0.71f, 0.68f)
                               material:*glass];

    transform = matrix4x4_translation(0.0f, 0.0f, 0.0f);

    // Create an instance of the light.
    GeometryInstance *lightMeshInstance = [[GeometryInstance alloc] initWithGeometry:lightMesh
                                                                           transform:transform
                                                                                mask:GEOMETRY_MASK_TRIANGLE_LIGHT];

    [scene addInstance:lightMeshInstance];

    // Create an instance of the Cornell box.
    GeometryInstance *geometryMeshInstance = [[GeometryInstance alloc] initWithGeometry:geometryMesh
                                                                              transform:transform
                                                                                   mask:GEOMETRY_MASK_TRIANGLE];

    [scene addInstance:geometryMeshInstance];
    
    // Create volume
    GeometryInstance *volumeMeshInstance = [[GeometryInstance alloc] initWithGeometry:volumeMesh
                                                                              transform:transform
                                                                                   mask:GEOMETRY_MASK_VOLUME_CONTAINER_TRIANGLE];

    [scene addInstance:volumeMeshInstance];

    // Create an instance of the sphere.
//    GeometryInstance *sphereGeometryInstance = [[GeometryInstance alloc] initWithGeometry:sphereGeometry
//                                                                                    transform:transform
//                                                                                         mask:GEOMETRY_MASK_VOLUME_CONTAINER_SPHERE];
//
//    [scene addInstance:sphereGeometryInstance];

    float3 background_color = vector3(0.0f, 0.0f, 0.0f);
    [scene setConstantEnvmapTexture:background_color];

    return scene;
}

#pragma mark - Create Scene Hetero-Vol
+ (RenderScene *)newTestSceneVolHetero:(id <MTLDevice>)device
{
    RenderScene *scene = [[RenderScene alloc] initWithDevice:device];

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
    back_wall_material->color = vector3(0.2f, 0.1f, 0.3f);
    light_material->color = vector3(20.f, 20.f, 20.f);
    glass->is_glass = true;
    glass->color = vector3(1.0f, 1.0f, 1.0f);
    glass->is_contain_volume = true;

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
    
    TriangleGeometry *volumeMesh = [[TriangleGeometry alloc] initWithDevice:device];

    [scene addGeometry:volumeMesh];
    
//     Import the OBJ File
//    NSURL *URL = [[NSBundle mainBundle] URLForResource:@"bunny-fine" withExtension:@"obj"];
//    [volumeMesh addGeometryWithURL:URL
//                           transform:transform * matrix4x4_scale(0.5f, 0.5f, 0.5f) * matrix4x4_translation(0.0f, -1.0f, 0.0f)
//                            material:*glass];
    
    // Add a cube in the center
    // ASSUME THE VOL CONTAINER IS ALWAYS A STANDARD CUBE WITHOUT TRANSFORM, ALL THE TRANSFORM NEED TO BE DONE ON THE INSTANCE SIDE
    [volumeMesh addCubeWithFaces:FACE_MASK_ALL
                             color:vector3(0.14f, 0.45f, 0.091f)
                         transform:matrix4x4_translation(0.0f, 0.0f, 0.0f)
                     inwardNormals:true
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

    [sphereGeometry addSphereWithOrigin:vector3(0.0f, 1.0f, 0.0f)
                                 radius:0.4f
                                  color:vector3(0.725f, 0.71f, 0.68f)
                               material:*glass];

    transform = matrix4x4_translation(0.0f, 0.0f, 0.0f);

    // Create an instance of the light.
    GeometryInstance *lightMeshInstance = [[GeometryInstance alloc] initWithGeometry:lightMesh
                                                                           transform:transform
                                                                                mask:GEOMETRY_MASK_TRIANGLE_LIGHT];

    [scene addInstance:lightMeshInstance];

    // Create an instance of the Cornell box.
    GeometryInstance *geometryMeshInstance = [[GeometryInstance alloc] initWithGeometry:geometryMesh
                                                                              transform:transform
                                                                                   mask:GEOMETRY_MASK_TRIANGLE];

    [scene addInstance:geometryMeshInstance];
    
    float3 background_color = vector3(0.0f, 0.0f, 0.0f);
    [scene setConstantEnvmapTexture:background_color];
    
    // create volume from VDB file
    float3 size;
    NSString *path = [[NSBundle mainBundle] pathForResource:@"cloud_01" ofType:@"vdb"];
    [scene setVolumeDensityGridTex:createVolume(path, device, scene->_maxDensity, size, 20.f, false)];
    std::cout << "Max density: " << scene->_maxDensity << std::endl;
    size = size / MIN(MIN(size.x, size.y), size.z);
    
    transform = matrix4x4_translation(0.0f, 0.99f, 0.0f) * matrix4x4_rotation(0.f* M_PI, vector3(0.f, 1.f, 0.f)) * matrix4x4_scale(0.618, 0.618, 0.618);
    // Create volume
    GeometryInstance *volumeMeshInstance = [[GeometryInstance alloc] initWithGeometry:volumeMesh
                                                                            transform:transform * matrix4x4_scale(size.x, size.y, size.z)
                                                                                   mask:GEOMETRY_MASK_VOLUME_CONTAINER_TRIANGLE];

    [scene addInstance:volumeMeshInstance];
    
    return scene;
}

#pragma mark - Create Scene Hetero-Vol
+ (RenderScene *)newTestSceneVolHeteroBunny:(id <MTLDevice>)device
{
    RenderScene *scene = [[RenderScene alloc] initWithDevice:device];

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
    back_wall_material->color = vector3(0.9f, 0.6f, 0.6f);
    light_material->color = vector3(20.f, 20.f, 20.f);
    glass->is_glass = false;
    glass->color = vector3(1.0f, 1.0f, 1.0f);
    glass->is_contain_volume = true;

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
    
    TriangleGeometry *volumeMesh = [[TriangleGeometry alloc] initWithDevice:device];

    [scene addGeometry:volumeMesh];
    
//     Import the OBJ File
//    NSURL *URL = [[NSBundle mainBundle] URLForResource:@"bunny-fine" withExtension:@"obj"];
//    [volumeMesh addGeometryWithURL:URL
//                           transform:transform * matrix4x4_scale(0.5f, 0.5f, 0.5f) * matrix4x4_translation(0.0f, -1.0f, 0.0f)
//                            material:*glass];
    
    // Add a cube in the center
    // ASSUME THE VOL CONTAINER IS ALWAYS A STANDARD CUBE WITHOUT TRANSFORM, ALL THE TRANSFORM NEED TO BE DONE ON THE INSTANCE SIDE
    [volumeMesh addCubeWithFaces:FACE_MASK_ALL
                             color:vector3(0.14f, 0.45f, 0.091f)
                         transform:matrix4x4_translation(0.0f, 0.0f, 0.0f)
                     inwardNormals:true
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

    [sphereGeometry addSphereWithOrigin:vector3(0.0f, 1.0f, 0.0f)
                                 radius:0.4f
                                  color:vector3(0.725f, 0.71f, 0.68f)
                               material:*glass];

    transform = matrix4x4_translation(0.0f, 0.0f, 0.0f);

    // Create an instance of the light.
    GeometryInstance *lightMeshInstance = [[GeometryInstance alloc] initWithGeometry:lightMesh
                                                                           transform:transform
                                                                                mask:GEOMETRY_MASK_TRIANGLE_LIGHT];

    [scene addInstance:lightMeshInstance];

    // Create an instance of the Cornell box.
    GeometryInstance *geometryMeshInstance = [[GeometryInstance alloc] initWithGeometry:geometryMesh
                                                                              transform:transform
                                                                                   mask:GEOMETRY_MASK_TRIANGLE];

    [scene addInstance:geometryMeshInstance];
    
    float3 background_color = vector3(0.0f, 0.0f, 0.0f);
    [scene setConstantEnvmapTexture:background_color];
    
    // create volume from VDB file
    float3 size;
    NSString *path = [[NSBundle mainBundle] pathForResource:@"bunny_cloud" ofType:@"vdb"];
    [scene setVolumeDensityGridTex:createVolume(path, device, scene->_maxDensity, size, 30.f, true)];
    std::cout << "Max density: " << scene->_maxDensity << std::endl;
    size = size / MIN(MIN(size.x, size.y), size.z);
    
    transform = matrix4x4_translation(0.0f, 0.6f, 0.2f) * matrix4x4_rotation(0.f* M_PI, vector3(0.f, 1.f, 0.f)) * matrix4x4_scale(0.618, 0.618, 0.618) * matrix4x4_scale(1.3f, 1.3f, 1.3f);
    // Create volume
    GeometryInstance *volumeMeshInstance = [[GeometryInstance alloc] initWithGeometry:volumeMesh
                                                                            transform:transform * matrix4x4_scale(size.x, size.y, size.z)
                                                                                   mask:GEOMETRY_MASK_VOLUME_CONTAINER_TRIANGLE];

    [scene addInstance:volumeMeshInstance];
    
    return scene;
}

@end

