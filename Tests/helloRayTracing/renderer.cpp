//
//  Renderer.cpp
//  helloRayTracing
//
//  Created by Ziyuan Qu on 2023/6/12.
//

#include "renderer.hpp"

using namespace simd;

struct Uniforms {
    unsigned int width;
    unsigned int height;
    unsigned int frameIndex;
    unsigned int lightCount;
    Camera camera;
};

static const NS::UInteger maxFramesInFlight = 3;
static const size_t alignedUniformsSize = (sizeof(Uniforms) + 255) & ~255;

Renderer::Renderer( MTL::Device* pDevice ) : _pDevice( pDevice->retain() )
{
//    _scene = Scene::newScene(_pDevice);
    _pCommandQueue = _pDevice->newCommandQueue();
    buildShaders();
    buildBuffers();
}

Renderer::~Renderer()
{
    _pVertexPositionsBuffer->release();
    _pVertexColorsBuffer->release();
    _pPSO->release();
    _pCommandQueue->release();
    _pDevice->release();
}

void Renderer::loadMetal(){
    NS::Error* pError = nullptr;
    _pLibrary = _pDevice->newDefaultLibrary();
    if ( !_pLibrary )
    {
        __builtin_printf( "%s", pError->localizedDescription()->utf8String() );
        assert( false );
    }
    _pCommandQueue = _pDevice->newCommandQueue();
}

MTL::ComputePipelineState* Renderer::newComputePipelineStateWithFunction(const MTL::Function* function, NS::Array* linked_functions){
    MTL::LinkedFunctions* mtl_linked_functions;
    if(linked_functions){
        mtl_linked_functions = MTL::LinkedFunctions::alloc()->init();
        mtl_linked_functions->setFunctions(linked_functions);
    }
    MTL::ComputePipelineDescriptor* descriptor = MTL::ComputePipelineDescriptor::alloc()->init();
    descriptor->setComputeFunction(function);
    descriptor->setLinkedFunctions(mtl_linked_functions);
    descriptor->setThreadGroupSizeIsMultipleOfThreadExecutionWidth(TRUE);
    
    NS::Error* pError = nullptr;
    MTL::ComputePipelineState* pipeline = _pDevice->newComputePipelineState(descriptor, 0, NULL, &pError);
    if ( !pipeline )
    {
        __builtin_printf( "%s", pError->localizedDescription()->utf8String() );
        assert( false );
    }
    return pipeline;
}

void Renderer::createPipelines(){
    using NS::StringEncoding::UTF8StringEncoding;
    
    MTL::Function* sphere_intersection_function = _pLibrary->newFunction(SphereGeometry::intersectionFunctionName());
    NS::Array* intersection_functions = NS::Array::array(sphere_intersection_function);
    MTL::Function* ray_tracing_function = _pLibrary->newFunction(NS::String::string("raytracingKernel", NS::UTF8StringEncoding));
    _raytracingPipeline = newComputePipelineStateWithFunction(ray_tracing_function, intersection_functions);
    MTL::IntersectionFunctionTableDescriptor* intersection_function_table_descriptor = MTL::IntersectionFunctionTableDescriptor::alloc()->init();
    intersection_function_table_descriptor->setFunctionCount(_scene->geometries->size());
    _intersectionFunctionTable = _raytracingPipeline->newIntersectionFunctionTable(intersection_function_table_descriptor);
    for(NS::UInteger geometry_index=0; geometry_index < _scene->geometries->size(); geometry_index++){
//        Geometry* geometry = _scene->geometries->at(geometry_index);
        MTL::Function* intersect_function = sphere_intersection_function; //Only Sphere for now!
        MTL::FunctionHandle* handle = _raytracingPipeline->functionHandle(intersect_function);
        _intersectionFunctionTable->setFunction(handle, geometry_index);
    }
    
    NS::Error* pError = nullptr;
    MTL::Function* pVertexFn = _pLibrary->newFunction( NS::String::string("copyVertex", UTF8StringEncoding) );
    MTL::Function* pFragFn = _pLibrary->newFunction( NS::String::string("copyFragment", UTF8StringEncoding) );
    MTL::RenderPipelineDescriptor* pDesc = MTL::RenderPipelineDescriptor::alloc()->init();
    pDesc->setVertexFunction( pVertexFn );
    pDesc->setFragmentFunction( pFragFn );
    pDesc->colorAttachments()->object(0)->setPixelFormat( MTL::PixelFormat::PixelFormatBGRA8Unorm_sRGB );

    _copyPipeline = _pDevice->newRenderPipelineState( pDesc, &pError );
    if ( !_copyPipeline )
    {
        __builtin_printf( "%s", pError->localizedDescription()->utf8String() );
        assert( false );
    }
}

void Renderer::createBuffers(){
    NS::UInteger uniform_buffer_size = alignedUniformsSize * maxFramesInFlight;
    _uniformBuffer = _pDevice->newBuffer(uniform_buffer_size, MTL::ResourceStorageModeManaged);
    _scene->uploadToBuffers();
    for (Geometry* geometry : *_scene->geometries){
        _resourcesStride = geometry->resources()->count() * sizeof(uint64_t);
    }
}

MTL::AccelerationStructure* Renderer::newAccelerationStructureWithDescriptor(MTL::AccelerationStructureDescriptor* descriptor){
    MTL::AccelerationStructureSizes accel_sizes = _pDevice->accelerationStructureSizes(descriptor);
    MTL::AccelerationStructure* acceleration_structure = _pDevice->newAccelerationStructure(accel_sizes.accelerationStructureSize);
    MTL::Buffer* scratch_buffer = _pDevice->newBuffer(accel_sizes.accelerationStructureSize, MTL::ResourceStorageModePrivate);
    MTL::CommandBuffer* command_buffer = _pCommandQueue->commandBuffer();
    MTL::AccelerationStructureCommandEncoder* command_encoder = command_buffer->accelerationStructureCommandEncoder();
    MTL::Buffer* compactedSizeBuffer = _pDevice->newBuffer(sizeof(uint32_t), MTL::ResourceStorageModeShared);
    command_encoder->buildAccelerationStructure(acceleration_structure, descriptor, scratch_buffer, 0);
    command_encoder->writeCompactedAccelerationStructureSize(acceleration_structure, compactedSizeBuffer, 0);
    command_encoder->endEncoding();
    command_buffer->commit();
    command_buffer->waitUntilCompleted();
    uint32_t compacted_size = *(uint32_t *)compactedSizeBuffer->contents();
    MTL::AccelerationStructure* compacted_acceleration_structure = _pDevice->newAccelerationStructure(compacted_size);
    command_buffer = _pCommandQueue->commandBuffer();
    command_encoder = command_buffer->accelerationStructureCommandEncoder();
    command_encoder->endEncoding();
    command_buffer->commit();
    return compacted_acceleration_structure;
}

template <typename T>
int findIndex(const std::vector<T>& vec, const T& value) {
    auto it = std::find(vec.begin(), vec.end(), value);
    if (it != vec.end()) {
        return std::distance(vec.begin(), it);
    } else {
        return -1;
    }
}

void Renderer::createAccelerationStructure(){
    for(NS::UInteger geometry_index=0; geometry_index < _scene->geometries->size(); geometry_index++){
        Geometry* geometry = _scene->geometries->at(geometry_index);
        MTL::AccelerationStructureGeometryDescriptor* geometry_descriptor = geometry->geometryDescriptor();
        geometry_descriptor->setIntersectionFunctionTableOffset(geometry_index);
        MTL::PrimitiveAccelerationStructureDescriptor* accele_descriptor = MTL::PrimitiveAccelerationStructureDescriptor::descriptor();
        accele_descriptor->setGeometryDescriptors(NS::Array::array(geometry_descriptor));
        MTL::AccelerationStructure* acceleration_structure = newAccelerationStructureWithDescriptor(accele_descriptor);
        _primitiveAccelerationStructures->push_back(acceleration_structure);
    }
    _instanceBuffer = _pDevice->newBuffer(sizeof(MTL::AccelerationStructureInstanceDescriptor)*_scene->instances->size(), MTL::ResourceStorageModeManaged);
    MTL::AccelerationStructureInstanceDescriptor* instance_descriptors = (MTL::AccelerationStructureInstanceDescriptor*) _instanceBuffer->contents();
    for(NS::UInteger instance_index=0; instance_index < _scene->instances->size(); instance_index++){
        GeometryInstance* instance = _scene->instances->at(instance_index);
        NS::UInteger geometry_index = findIndex(*_scene->geometries, instance->geometry);
        instance_descriptors[instance_index].accelerationStructureIndex = (uint32_t)geometry_index;
        instance_descriptors[instance_index].options = instance->geometry->intersectionFunctionName() == NULL ? MTL::AccelerationStructureInstanceOptionOpaque : 0;
        instance_descriptors[instance_index].intersectionFunctionTableOffset = 0;
        instance_descriptors[instance_index].mask = (uint32_t)instance->mask;
        for (int column = 0; column < 4; column++)
            for (int row = 0; row < 3; row++)
                instance_descriptors[instance_index].transformationMatrix.columns[column][row] = instance->transform.columns[column][row];
    }
    _instanceBuffer->didModifyRange(NS::Range(0, _instanceBuffer->length()));
    
    MTL::InstanceAccelerationStructureDescriptor* accel_descriptor = MTL::InstanceAccelerationStructureDescriptor::descriptor();
    
    NS::Array* ns_array = NS::Array::alloc();
    for(NS::UInteger i = 0; i < _primitiveAccelerationStructures->size(); i++){
        MTL::AccelerationStructure* mtlas= _primitiveAccelerationStructures->at(i);
//        NS::Value* value = NS::Value::alloc()->init(mtlas, "MTL::AccelerationStructure*");
        NS::Array* tmp_ns_array = NS::Array::array(_primitiveAccelerationStructures);
        ns_array->object(i)->
    }
    NS::Array::alloc()->init(_primitiveAccelerationStructures->data(), _primitiveAccelerationStructures->size());
    accel_descriptor->setInstancedAccelerationStructures(ns_array);
    accel_descriptor->setInstanceCount(_scene->instances->size());
    accel_descriptor->setInstanceDescriptorBuffer(_instanceBuffer);
    _instanceAccelerationStructure = newAccelerationStructureWithDescriptor(accel_descriptor);
    
}

void Renderer::buildShaders()
{
    using NS::StringEncoding::UTF8StringEncoding;

    NS::Error* pError = nullptr;
    _pLibrary = _pDevice->newDefaultLibrary();
    if ( !_pLibrary )
    {
        __builtin_printf( "%s", pError->localizedDescription()->utf8String() );
        assert( false );
    }

    MTL::Function* pVertexFn = _pLibrary->newFunction( NS::String::string("vertexMain", UTF8StringEncoding) );
    MTL::Function* pFragFn = _pLibrary->newFunction( NS::String::string("fragmentMain", UTF8StringEncoding) );

    MTL::RenderPipelineDescriptor* pDesc = MTL::RenderPipelineDescriptor::alloc()->init();
    pDesc->setVertexFunction( pVertexFn );
    pDesc->setFragmentFunction( pFragFn );
    pDesc->colorAttachments()->object(0)->setPixelFormat( MTL::PixelFormat::PixelFormatBGRA8Unorm_sRGB );

    _pPSO = _pDevice->newRenderPipelineState( pDesc, &pError );
    if ( !_pPSO )
    {
        __builtin_printf( "%s", pError->localizedDescription()->utf8String() );
        assert( false );
    }

    pVertexFn->release();
    pFragFn->release();
    pDesc->release();
    _pLibrary->release();
}

void Renderer::buildBuffers()
{
    const size_t NumVertices = 6;

    simd::float3 positions[NumVertices] =
    {
        { -0.8f,  0.8f, 0.0f },
        {  0.0f, -0.8f, 0.0f },
        { +0.8f,  0.8f, 0.0f },
        { -0.5f,  0.5f, 0.0f },
        {  0.0f, -0.5f, 0.0f },
        { +0.5f,  0.5f, 0.0f }
    };

    simd::float3 colors[NumVertices] =
    {
        {  1.0, 0.3f, 0.2f },
        {  0.8f, 1.0, 0.0f },
        {  0.8f, 0.0f, 1.0 },
        {  1.0, 0.3f, 0.2f },
        {  0.8f, 1.0, 0.0f },
        {  0.8f, 0.0f, 1.0 }
    };

    const size_t positionsDataSize = NumVertices * sizeof( simd::float3 );
    const size_t colorDataSize = NumVertices * sizeof( simd::float3 );

    MTL::Buffer* pVertexPositionsBuffer = _pDevice->newBuffer( positionsDataSize, MTL::ResourceStorageModeManaged );
    MTL::Buffer* pVertexColorsBuffer = _pDevice->newBuffer( colorDataSize, MTL::ResourceStorageModeManaged );

    _pVertexPositionsBuffer = pVertexPositionsBuffer;
    _pVertexColorsBuffer = pVertexColorsBuffer;

    memcpy( _pVertexPositionsBuffer->contents(), positions, positionsDataSize );
    memcpy( _pVertexColorsBuffer->contents(), colors, colorDataSize );

    _pVertexPositionsBuffer->didModifyRange( NS::Range::Make( 0, _pVertexPositionsBuffer->length() ) );
    _pVertexColorsBuffer->didModifyRange( NS::Range::Make( 0, _pVertexColorsBuffer->length() ) );
}

void Renderer::draw( MTK::View* pView )
{
    NS::AutoreleasePool* pPool = NS::AutoreleasePool::alloc()->init();

    MTL::CommandBuffer* pCmd = _pCommandQueue->commandBuffer();
    MTL::RenderPassDescriptor* pRpd = pView->currentRenderPassDescriptor();
    MTL::RenderCommandEncoder* pEnc = pCmd->renderCommandEncoder( pRpd );

    pEnc->setRenderPipelineState( _pPSO );
    pEnc->setVertexBuffer( _pVertexPositionsBuffer, 0, 0 );
    pEnc->setVertexBuffer( _pVertexColorsBuffer, 0, 1 );
    pEnc->drawPrimitives( MTL::PrimitiveType::PrimitiveTypeTriangle, NS::UInteger(0), NS::UInteger(3) );
    pEnc->drawPrimitives( MTL::PrimitiveType::PrimitiveTypeTriangle, NS::UInteger(3), NS::UInteger(6) );

    pEnc->endEncoding();
    pCmd->presentDrawable( pView->currentDrawable() );
    pCmd->commit();

    pPool->release();
}


