/*
See LICENSE folder for this sample’s licensing information.

Abstract:
The implementation of the cross-platform view controller.
*/
#import "ViewController.h"
#import "Renderer.h"
#import <iostream>

@implementation ViewController
{
    MTKView *_view;

    Renderer *_renderer;
    
    NSString *_sceneName;
}

- (void)viewDidLoad
{
    [super viewDidLoad];

    _view = (MTKView *)self.view;

#if TARGET_OS_IPHONE
    _view.device = MTLCreateSystemDefaultDevice();
#else
    NSArray<id<MTLDevice>> *devices = MTLCopyAllDevices();

    id<MTLDevice> selectedDevice;

    for(id<MTLDevice> device in devices)
    {
        if(device.supportsRaytracing)
        {
            if(!selectedDevice || !device.isLowPower)
            {
                selectedDevice = device;
            }
        }
    }
    _view.device = selectedDevice;
    
    CAMetalLayer *metalLayer = (CAMetalLayer *)_view.layer;
    metalLayer.displaySyncEnabled = YES;

    NSLog(@"Selected Device: %@", _view.device.name);
#endif

    // Device must support Metal and ray tracing.
    NSAssert(_view.device && _view.device.supportsRaytracing,
             @"Ray tracing isn't supported on this device");

#if TARGET_OS_IPHONE
    _view.backgroundColor = UIColor.blackColor;
#endif
    _view.colorPixelFormat = MTLPixelFormatRGBA16Float;

//    Scene *scene = [Scene newInstancedCornellBoxSceneWithDevice:_view.device
//                                       useIntersectionFunctions:YES];
    
#ifdef HAS_SWIFT_UI
    RenderScene *scene;
    if([_sceneName  isEqual: @"Scene 1"])
        scene = [RenderScene newTestScene:_view.device];
    if([_sceneName  isEqual: @"Scene 2"])
        scene = [RenderScene newTestSceneObj:_view.device];
    if([_sceneName  isEqual: @"Scene 3"])
        scene = [RenderScene newTestSceneMIS:_view.device];
    if([_sceneName  isEqual: @"Scene 4"])
        scene = [RenderScene newTestSceneEnv:_view.device];
    if([_sceneName  isEqual: @"Scene 5"])
        scene = [RenderScene newTestSceneVolHomo:_view.device];
    if([_sceneName  isEqual: @"Scene 6"])
        scene = [RenderScene newTestSceneVolHetero:_view.device];
    if([_sceneName  isEqual: @"Scene 7"])
        scene = [RenderScene newTestSceneVolHeteroBunny:_view.device];
#else
    RenderScene *scene = [RenderScene newTestSceneVolHeteroForResearch:_view.device];
#endif
    
    _renderer = [[Renderer alloc] initWithDevice:_view.device
                                           scene:scene];

    [_renderer mtkView:_view drawableSizeWillChange:_view.bounds.size];

    _view.delegate = _renderer;
    
    std::cout << "I initialize" << std::endl;
}

- (void)setRenderMode:(int)renderModeIndex {
    _renderer.renderMode = (RenderMode)renderModeIndex;
    std::cout << "set render mode to " << renderModeIndex << std::endl;
}

- (void)setSceneName:(NSString *)sceneName {
    _sceneName = sceneName;
    std::cout << "load scene: " << [sceneName UTF8String] << std::endl;
}

- (void)setRenderVolume:(bool)renderVolume {
    
}

@end
