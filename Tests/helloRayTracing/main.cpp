//
//  main.cpp
//  helloRayTracing
//
//  Created by Ziyuan Qu on 2023/6/12.
//

#include "renderer.cpp"

#pragma region Declarations {

class MyMTKViewDelegate : public MTK::ViewDelegate
{
    public:
        MyMTKViewDelegate( MTL::Device* pDevice );
        virtual ~MyMTKViewDelegate() override;
        virtual void drawInMTKView( MTK::View* pView ) override;

    private:
        Renderer* _pRenderer;
};

class MyAppDelegate : public NS::ApplicationDelegate
{
    public:
        ~MyAppDelegate();

        NS::Menu* createMenuBar();

        virtual void applicationWillFinishLaunching( NS::Notification* pNotification ) override;
        virtual void applicationDidFinishLaunching( NS::Notification* pNotification ) override;
        virtual bool applicationShouldTerminateAfterLastWindowClosed( NS::Application* pSender ) override;

    private:
        NS::Window* _pWindow;
        MTK::View* _pMtkView;
        MTL::Device* _pDevice;
        MyMTKViewDelegate* _pViewDelegate = nullptr;
};

#pragma endregion Declarations }


int main( int argc, char* argv[] )
{
    NS::AutoreleasePool* pAutoreleasePool = NS::AutoreleasePool::alloc()->init();

    MyAppDelegate del;

    NS::Application* pSharedApplication = NS::Application::sharedApplication();
    pSharedApplication->setDelegate( &del );
    pSharedApplication->run();

    pAutoreleasePool->release();

    return 0;
}


#pragma mark - AppDelegate
#pragma region AppDelegate {

MyAppDelegate::~MyAppDelegate()
{
    _pMtkView->release();
    _pWindow->release();
    _pDevice->release();
    delete _pViewDelegate;
}

NS::Menu* MyAppDelegate::createMenuBar()
{
    using NS::StringEncoding::UTF8StringEncoding;

    NS::Menu* pMainMenu = NS::Menu::alloc()->init();
    NS::MenuItem* pAppMenuItem = NS::MenuItem::alloc()->init();
    NS::Menu* pAppMenu = NS::Menu::alloc()->init( NS::String::string( "Appname", UTF8StringEncoding ) );

    NS::String* appName = NS::RunningApplication::currentApplication()->localizedName();
    NS::String* quitItemName = NS::String::string( "Quit ", UTF8StringEncoding )->stringByAppendingString( appName );
    SEL quitCb = NS::MenuItem::registerActionCallback( "appQuit", [](void*,SEL,const NS::Object* pSender){
        auto pApp = NS::Application::sharedApplication();
        pApp->terminate( pSender );
    } );

    NS::MenuItem* pAppQuitItem = pAppMenu->addItem( quitItemName, quitCb, NS::String::string( "q", UTF8StringEncoding ) );
    pAppQuitItem->setKeyEquivalentModifierMask( NS::EventModifierFlagCommand );
    pAppMenuItem->setSubmenu( pAppMenu );

    NS::MenuItem* pWindowMenuItem = NS::MenuItem::alloc()->init();
    NS::Menu* pWindowMenu = NS::Menu::alloc()->init( NS::String::string( "Window", UTF8StringEncoding ) );

    SEL closeWindowCb = NS::MenuItem::registerActionCallback( "windowClose", [](void*, SEL, const NS::Object*){
        auto pApp = NS::Application::sharedApplication();
            pApp->windows()->object< NS::Window >(0)->close();
    } );
    NS::MenuItem* pCloseWindowItem = pWindowMenu->addItem( NS::String::string( "Close Window", UTF8StringEncoding ), closeWindowCb, NS::String::string( "w", UTF8StringEncoding ) );
    pCloseWindowItem->setKeyEquivalentModifierMask( NS::EventModifierFlagCommand );

    pWindowMenuItem->setSubmenu( pWindowMenu );

    pMainMenu->addItem( pAppMenuItem );
    pMainMenu->addItem( pWindowMenuItem );

    pAppMenuItem->release();
    pWindowMenuItem->release();
    pAppMenu->release();
    pWindowMenu->release();

    return pMainMenu->autorelease();
}

void MyAppDelegate::applicationWillFinishLaunching( NS::Notification* pNotification )
{
    NS::Menu* pMenu = createMenuBar();
    NS::Application* pApp = reinterpret_cast< NS::Application* >( pNotification->object() );
    pApp->setMainMenu( pMenu );
    pApp->setActivationPolicy( NS::ActivationPolicy::ActivationPolicyRegular );
}

void MyAppDelegate::applicationDidFinishLaunching( NS::Notification* pNotification )
{
    CGRect frame = (CGRect){ {100.0, 100.0}, {512.0, 512.0} };

    _pWindow = NS::Window::alloc()->init(
        frame,
        NS::WindowStyleMaskClosable|NS::WindowStyleMaskTitled,
        NS::BackingStoreBuffered,
        false );

    _pDevice = MTL::CreateSystemDefaultDevice();

    _pMtkView = MTK::View::alloc()->init( frame, _pDevice );
    _pMtkView->setColorPixelFormat( MTL::PixelFormat::PixelFormatBGRA8Unorm_sRGB );
    _pMtkView->setClearColor( MTL::ClearColor::Make( 1.0, 0.0, 0.0, 1.0 ) );

    _pViewDelegate = new MyMTKViewDelegate( _pDevice );
    _pMtkView->setDelegate( _pViewDelegate );

    _pWindow->setContentView( _pMtkView );
    _pWindow->setTitle( NS::String::string( "01 - Primitive", NS::StringEncoding::UTF8StringEncoding ) );

    _pWindow->makeKeyAndOrderFront( nullptr );

    NS::Application* pApp = reinterpret_cast< NS::Application* >( pNotification->object() );
    pApp->activateIgnoringOtherApps( true );
}

bool MyAppDelegate::applicationShouldTerminateAfterLastWindowClosed( NS::Application* pSender )
{
    return true;
}

#pragma endregion AppDelegate }


#pragma mark - ViewDelegate
#pragma region ViewDelegate {

MyMTKViewDelegate::MyMTKViewDelegate( MTL::Device* pDevice )
: MTK::ViewDelegate()
, _pRenderer( new Renderer( pDevice ) )
{
}

MyMTKViewDelegate::~MyMTKViewDelegate()
{
    delete _pRenderer;
}

void MyMTKViewDelegate::drawInMTKView( MTK::View* pView )
{
    _pRenderer->draw( pView );
}

#pragma endregion ViewDelegate }


