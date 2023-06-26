/*
See LICENSE folder for this sample’s licensing information.

Abstract:
The header for the cross-platform view controller.
*/
#import <Metal/Metal.h>
#import <MetalKit/MetalKit.h>
#import <TargetConditionals.h>

#if !TARGET_OS_IPHONE
#import <Cocoa/Cocoa.h>
#else
#import <UIKit/UIKit.h>
#endif

#if !TARGET_OS_IPHONE
@interface ViewController : NSViewController

@property (nonatomic, strong) NSString *sceneName;

- (void)setRenderMode:(int)renderModeIndex;

#else
@interface ViewController : UIViewController
#endif

@end
