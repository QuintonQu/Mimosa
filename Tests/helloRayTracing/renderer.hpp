//
//  Renderer.hpp
//  helloRayTracing
//
//  Created by Ziyuan Qu on 2023/6/12.
//

#ifndef renderer_hpp
#define renderer_hpp

#define NS_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#define MTK_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#include <Metal/Metal.hpp>
#include <AppKit/AppKit.hpp>
#include <MetalKit/MetalKit.hpp>
#include <simd/simd.h>
#include <cassert>
#import "scene.hpp"

class Renderer
{
public:
    Renderer( MTL::Device* pDevice );
    ~Renderer();
    void buildShaders();
    void buildBuffers();
    void draw( MTK::View* pView );
    void loadMetal();
    MTL::ComputePipelineState* newComputePipelineStateWithFunction(const MTL::Function* function, NS::Array* linkedFunctions);
    // MTL::Function* specializedFunctionWithName(NS::String name);
    void createPipelines();
    void createBuffers();
    MTL::AccelerationStructure* newAccelerationStructureWithDescriptor(MTL::AccelerationStructureDescriptor* descriptor);
    void createAccelerationStructure();

private:
    MTL::Device* _pDevice;
    MTL::CommandQueue* _pCommandQueue;
    MTL::Library* _pLibrary;
    
    MTL::Buffer* _uniformBuffer;
    
    MTL::AccelerationStructure* _instanceAccelerationStructure;
    std::vector<MTL::AccelerationStructure*>* _primitiveAccelerationStructures;
    
    MTL::ComputePipelineState* _raytracingPipeline;
    MTL::RenderPipelineState* _copyPipeline;
    MTL::RenderPipelineState* _pPSO;
    
    MTL::Texture* _accumulationTargets[2];
    MTL::Texture* _randomTexture;
    
    MTL::Buffer* _resourceBuffer;
    MTL::Buffer* _instanceBuffer;
    
    MTL::IntersectionFunctionTable* _intersectionFunctionTable;
    
    dispatch_semaphore_t _semaphore;
    CGSize _size;
    NS::UInteger _uniformBufferOffset;
    NS::UInteger _uniformBufferIndex;
    
    unsigned int _frameIndex;
    
    Scene* _scene;
    
    NS::UInteger _resourcesStride;
    
    MTL::Buffer* _pVertexPositionsBuffer;
    MTL::Buffer* _pVertexColorsBuffer;
};


#endif /* Renderer_hpp */
