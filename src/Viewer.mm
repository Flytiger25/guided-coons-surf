#import <Cocoa/Cocoa.h>

#include "Viewer.h"

#include <AIS_Shape.hxx>
#include <AIS_DisplayMode.hxx>
#include <AIS_InteractiveContext.hxx>
#include <V3d_Viewer.hxx>
#include <V3d_View.hxx>
#include <Aspect_DisplayConnection.hxx>
#include <OpenGl_GraphicDriver.hxx>
#include <Cocoa_Window.hxx>
#include <Cocoa_LocalPool.hxx>
#include <Graphic3d_NameOfMaterial.hxx>

#include <STEPControl_Reader.hxx>
#include <TopExp_Explorer.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>

#include <iostream>

// ============================================================================
// OCC 交互 NSView：只负责鼠标事件和 OpenGL 渲染
// ============================================================================
@interface OCCInteractiveView : NSView
{
    NSPoint        myLastPoint;
    Standard_Real  myZoomDelta;
    Standard_Real  myZoomThreshold;
}
@property (nonatomic) Standard_Boolean isRotating;
@property (nonatomic) Standard_Boolean isPanning;
@end

@implementation OCCInteractiveView

static Handle(V3d_View)                 sOCCView;
static Handle(AIS_InteractiveContext)   sOCCContext;
static Handle(AIS_Shape)                sAISShape;

+ (void)setOCCView:(const Handle(V3d_View)&)v                { sOCCView   = v; }
+ (void)setOCCContext:(const Handle(AIS_InteractiveContext)&)c { sOCCContext = c; }
+ (void)setAISShape:(const Handle(AIS_Shape)&)s               { sAISShape  = s; }
+ (Handle(V3d_View))occView               { return sOCCView; }
+ (Handle(AIS_InteractiveContext))occContext { return sOCCContext; }
+ (Handle(AIS_Shape))aisShape             { return sAISShape; }

- (instancetype)initWithFrame:(NSRect)frame
{
    self = [super initWithFrame:frame];
    if (self)
    {
        _isRotating = Standard_False;
        _isPanning  = Standard_False;
        myZoomDelta     = 0.0;
        myZoomThreshold = 3.0;
    }
    return self;
}

- (BOOL)acceptsFirstResponder { return YES; }
- (BOOL)isFlipped { return YES; }

// ── 鼠标 ─────────────────────────────────────────────────
- (void)mouseDown:(NSEvent*)e
{
    myLastPoint = [self convertPoint:[e locationInWindow] fromView:nil];
    self.isRotating = Standard_True;
    sOCCView->StartRotation(Standard_Integer(myLastPoint.x), Standard_Integer(myLastPoint.y));
}
- (void)rightMouseDown:(NSEvent*)e
{
    myLastPoint = [self convertPoint:[e locationInWindow] fromView:nil];
    self.isPanning = Standard_True;
}
- (void)otherMouseDown:(NSEvent*)e
{
    myLastPoint = [self convertPoint:[e locationInWindow] fromView:nil];
    self.isPanning = Standard_True;
}
- (void)mouseDragged:(NSEvent*)e
{
    NSPoint aPt = [self convertPoint:[e locationInWindow] fromView:nil];
    if (self.isRotating)
        sOCCView->Rotation(Standard_Integer(aPt.x), Standard_Integer(aPt.y));
    if (self.isPanning)
    {
        sOCCView->Pan(Standard_Integer(aPt.x - myLastPoint.x),
                      Standard_Integer(myLastPoint.y - aPt.y));
    }
    myLastPoint = aPt;
}
- (void)rightMouseDragged:(NSEvent*)e
{
    NSPoint aPt = [self convertPoint:[e locationInWindow] fromView:nil];
    if (self.isPanning)
    {
        sOCCView->Pan(Standard_Integer(aPt.x - myLastPoint.x),
                      Standard_Integer(myLastPoint.y - aPt.y));
    }
    myLastPoint = aPt;
}
- (void)otherMouseDragged:(NSEvent*)e { [self rightMouseDragged:e]; }
- (void)mouseUp:(NSEvent*)e    { self.isRotating = Standard_False; }
- (void)rightMouseUp:(NSEvent*)e  { self.isPanning = Standard_False; }
- (void)otherMouseUp:(NSEvent*)e  { self.isPanning = Standard_False; }

- (void)scrollWheel:(NSEvent*)e
{
    myZoomDelta += [e deltaY];
    if (std::abs(myZoomDelta) < myZoomThreshold && ![e hasPreciseScrollingDeltas])
        return;
    NSPoint aPt = [self convertPoint:[e locationInWindow] fromView:nil];
    sOCCView->StartZoomAtPoint(Standard_Integer(aPt.x), Standard_Integer(aPt.y));
    sOCCView->ZoomAtPoint(0, 0, Standard_Integer(myZoomDelta), 0);
    myZoomDelta = 0.0;
}
- (void)magnifyWithEvent:(NSEvent*)e
{
    NSPoint aPt = [self convertPoint:[e locationInWindow] fromView:nil];
    sOCCView->StartZoomAtPoint(Standard_Integer(aPt.x), Standard_Integer(aPt.y));
    sOCCView->ZoomAtPoint(0, 0, Standard_Integer([e magnification] * 50.0), 0);
}
- (void)setFrameSize:(NSSize)newSize
{
    [super setFrameSize:newSize];
    if (!sOCCView.IsNull()) { sOCCView->MustBeResized(); sOCCView->Redraw(); }
}
@end

// ============================================================================
// 控制栏按钮的点击回调
// ============================================================================
@interface ViewerControlPanel : NSView
{
    NSButton* mOpenBtn;
    NSButton* mWireBtn;
    NSButton* mShadeBtn;
    NSButton* mTransBtn;
    NSButton* mFitBtn;
}
@end

@implementation ViewerControlPanel

- (instancetype)initWithFrame:(NSRect)frame
{
    self = [super initWithFrame:frame];
    if (self)
    {
        self.wantsLayer = YES;
        self.layer.backgroundColor = [[NSColor colorWithWhite:0.15 alpha:1.0] CGColor];

        auto makeBtn = ^NSButton*(NSString* title, NSRect r, SEL action) {
            NSButton* b = [[NSButton alloc] initWithFrame:r];
            b.title = title;
            b.bezelStyle = NSBezelStyleRounded;
            [b setButtonType:NSButtonTypeMomentaryPushIn];
            [b setTarget:self];
            [b setAction:action];
            return b;
        };

        CGFloat x = 12, y = 4, w = 72, h = 28, gap = 8;
        mOpenBtn  = makeBtn(@"打开",  NSMakeRect(x, y, w, h), @selector(onOpen:));
        x += w + gap;
        mWireBtn  = makeBtn(@"线框",  NSMakeRect(x, y, w, h), @selector(onWireframe:));
        x += w + gap;
        mShadeBtn = makeBtn(@"着色",  NSMakeRect(x, y, w, h), @selector(onShaded:));
        x += w + gap;
        mTransBtn = makeBtn(@"透明度", NSMakeRect(x, y, w, h), @selector(onTransparency:));
        x += w + gap;
        mFitBtn   = makeBtn(@"适配",  NSMakeRect(x, y, w, h), @selector(onFitAll:));

        [self addSubview:mOpenBtn];
        [self addSubview:mWireBtn];
        [self addSubview:mShadeBtn];
        [self addSubview:mTransBtn];
        [self addSubview:mFitBtn];
    }
    return self;
}

- (void)onOpen:(id)sender
{
    NSOpenPanel* panel = [NSOpenPanel openPanel];
    [panel setTitle:@"选择 STEP 文件"];
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
    [panel setAllowedFileTypes:@[@"step", @"stp"]];
#pragma clang diagnostic pop
    [panel setCanChooseFiles:YES];
    [panel setCanChooseDirectories:NO];
    [panel setAllowsMultipleSelection:NO];

    if ([panel runModal] != NSModalResponseOK) return;

    NSURL* fileURL = [panel URL];
    if (!fileURL) return;

    std::string filePath = [[fileURL path] UTF8String];

    // 读取 STEP 文件
    STEPControl_Reader reader;
    IFSelect_ReturnStatus status = reader.ReadFile(filePath.c_str());
    if (status != IFSelect_RetDone)
    {
        NSAlert* alert = [[NSAlert alloc] init];
        [alert setMessageText:@"无法读取文件"];
        [alert setInformativeText:[NSString stringWithFormat:@"STEP 文件读取失败: %s", filePath.c_str()]];
        [alert runModal];
        return;
    }

    reader.TransferRoots();
    TopoDS_Shape shape = reader.OneShape();

    if (shape.IsNull())
    {
        NSAlert* alert = [[NSAlert alloc] init];
        [alert setMessageText:@"文件为空"];
        [alert setInformativeText:@"STEP 文件中没有有效的几何数据。"];
        [alert runModal];
        return;
    }

    // 在 OCC 视图中追加显示
    auto ctx = [OCCInteractiveView occContext];
    auto view = [OCCInteractiveView occView];
    if (ctx.IsNull() || view.IsNull()) return;

    Handle(AIS_Shape) aisShape = new AIS_Shape(shape);
    aisShape->SetMaterial(Graphic3d_NOM_ALUMINIUM);
    aisShape->SetTransparency(0.1);
    aisShape->SetColor(Quantity_NOC_STEELBLUE3);
    ctx->Display(aisShape, Standard_True);
    ctx->SetDisplayMode(aisShape, AIS_Shaded, Standard_False);
    [OCCInteractiveView setAISShape:aisShape];

    view->FitAll();
    view->Redraw();

    std::cout << "已加载: " << filePath << std::endl;
}

- (void)onWireframe:(id)sender
{
    auto ctx   = [OCCInteractiveView occContext];
    auto shape = [OCCInteractiveView aisShape];
    ctx->SetDisplayMode(shape, AIS_WireFrame, Standard_True);
}
- (void)onShaded:(id)sender
{
    auto ctx   = [OCCInteractiveView occContext];
    auto shape = [OCCInteractiveView aisShape];
    ctx->SetDisplayMode(shape, AIS_Shaded, Standard_True);
}
- (void)onTransparency:(id)sender
{
    Handle(AIS_Shape) aShape = [OCCInteractiveView aisShape];
    if (aShape.IsNull()) return;
    Standard_Real nt = (aShape->Transparency() > 0.5) ? 0.1 : 0.8;
    aShape->SetTransparency(nt);
    [OCCInteractiveView occContext]->Redisplay(aShape, Standard_True);
    [OCCInteractiveView occView]->Redraw();
}
- (void)onFitAll:(id)sender
{
    [OCCInteractiveView occView]->FitAll();
    [OCCInteractiveView occView]->Redraw();
}
@end

// ============================================================================
// Cocoa 代理
// ============================================================================
@interface OccViewerDelegate : NSObject <NSApplicationDelegate>
@end
@implementation OccViewerDelegate
- (BOOL)applicationShouldTerminateAfterLastWindowClosed:(NSApplication*)s { return YES; }
@end

// ============================================================================
// 公开接口
// ============================================================================
void DisplayShape(const TopoDS_Shape& theShape, const char* theTitle)
{
    if (NSApp == nil)
    {
        [NSApplication sharedApplication];
        [NSApp setActivationPolicy:NSApplicationActivationPolicyRegular];
        [NSApp setDelegate:[[OccViewerDelegate alloc] init]];
        [NSApp finishLaunching];
    }

    Standard_Integer winW = 900, winH = 636; // 600 + 36 控制栏
    Standard_Integer ctrlH = 36;

    // ── 窗口容器 ──
    NSRect aFrame = NSMakeRect(100, 100, winW, winH);
    NSView* container = [[NSView alloc] initWithFrame:aFrame];

    // ── 控制栏（顶部） ──
    NSRect ctrlFrame = NSMakeRect(0, winH - ctrlH, winW, ctrlH);
    ViewerControlPanel* panel = [[ViewerControlPanel alloc] initWithFrame:ctrlFrame];
    panel.autoresizingMask = NSViewWidthSizable | NSViewMinYMargin;
    [container addSubview:panel];

    // ── OCC GL 视图（下方填满） ──
    NSRect glFrame = NSMakeRect(0, 0, winW, winH - ctrlH);
    OCCInteractiveView* glView = [[OCCInteractiveView alloc] initWithFrame:glFrame];
    glView.autoresizingMask = NSViewWidthSizable | NSViewHeightSizable;
    [container addSubview:glView];

    // ── NSWindow ──
    NSWindow* aNSWindow = [[NSWindow alloc]
        initWithContentRect:aFrame
                  styleMask:NSWindowStyleMaskTitled
                              | NSWindowStyleMaskClosable
                              | NSWindowStyleMaskResizable
                              | NSWindowStyleMaskMiniaturizable
                    backing:NSBackingStoreBuffered
                      defer:NO];
    [aNSWindow setTitle:[NSString stringWithUTF8String:theTitle]];
    [aNSWindow setContentView:container];
    [aNSWindow setAcceptsMouseMovedEvents:YES];

    // ── OCC 初始化 ──
    Handle(Aspect_DisplayConnection) aDisplayConn = new Aspect_DisplayConnection();
    Handle(OpenGl_GraphicDriver)    aDriver      = new OpenGl_GraphicDriver(aDisplayConn);

    Handle(V3d_Viewer) aViewer = new V3d_Viewer(aDriver);
    aViewer->SetDefaultLights();
    aViewer->SetLightOn();

    Handle(AIS_InteractiveContext) aContext = new AIS_InteractiveContext(aViewer);

    Handle(Cocoa_Window) aWindow = new Cocoa_Window((__bridge NSView*)glView);
    aWindow->Map();
    aWindow->DoMapping();

    Handle(V3d_View) aView = aViewer->CreateView();
    aView->SetWindow(aWindow);
    [OCCInteractiveView setOCCView:aView];
    [OCCInteractiveView setOCCContext:aContext];

    aView->SetBackgroundColor(Quantity_NOC_GRAY30);
    aView->SetBgGradientColors(Quantity_NOC_GRAY20, Quantity_NOC_GRAY60,
                               Aspect_GFM_VER, Standard_False);

    // ── 显示 shape（默认着色模式） ──
    Handle(AIS_Shape) anAisShape = new AIS_Shape(theShape);
    anAisShape->SetMaterial(Graphic3d_NOM_ALUMINIUM);
    anAisShape->SetTransparency(0.1);
    anAisShape->SetColor(Quantity_NOC_STEELBLUE3);
    aContext->Display(anAisShape, Standard_True);
    aContext->SetDisplayMode(anAisShape, AIS_Shaded, Standard_False);
    [OCCInteractiveView setAISShape:anAisShape];
    aView->FitAll();
    aView->Redraw();

    [aNSWindow center];
    [aNSWindow makeKeyAndOrderFront:nil];
    [NSApp activateIgnoringOtherApps:YES];
    [NSApp run];
}
