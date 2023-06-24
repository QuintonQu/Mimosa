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
}

- (void)viewDidLoad
{
    [super viewDidLoad];
    
    [self _configureBackdrop:_configBackdrop];

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
    
    Scene *scene = [Scene newTestSceneMIS:_view.device];

    _renderer = [[Renderer alloc] initWithDevice:_view.device
                                           scene:scene];

    [_renderer mtkView:_view drawableSizeWillChange:_view.bounds.size];

    _view.delegate = _renderer;
}

- (void)_configureBackdrop:(NSView *)view
{
    view.wantsLayer = YES;
    view.layer.borderWidth = 1.0f;
    view.layer.cornerRadius = 8.0f;
    
    CGFloat components[] = {1.0, 1.0, 1.0, 0.8}; // RGBA values for white color
    CGColorSpaceRef colorSpace = CGColorSpaceCreateDeviceRGB();
    CGColorRef whiteColor = CGColorCreate(colorSpace, components);
    view.layer.backgroundColor = whiteColor;
}

- (IBAction)renderModeControl:(id)sender {
    if ( sender == _renderModeControl )
    {
        _renderer.renderMode = (RenderMode)_renderModeControl.indexOfSelectedItem;
//        std::cout << _renderModeControl.indexOfSelectedItem << std::endl;
    }
}

@end
