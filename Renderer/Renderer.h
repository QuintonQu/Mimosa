/*
See LICENSE folder for this sample’s licensing information.

Abstract:
The header for the renderer class that performs Metal setup and per-frame rendering.
*/

#import <Metal/Metal.h>
#import <MetalKit/MetalKit.h>

#import "Scene.h"

typedef NS_ENUM( uint8_t, RenderMode )
{
    Mats = 0,
    NEE = 1,
    MIS = 2
};


@interface Renderer : NSObject <MTKViewDelegate>

- (instancetype)initWithDevice:(id<MTLDevice>)device
                         scene:(Scene *)scene;
- (void)setRenderMode:(RenderMode)renderMode;

@end
