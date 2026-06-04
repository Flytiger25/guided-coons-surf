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
#include <vector>
#include <string>

// ============================================================================
// 全局模型列表
// ============================================================================
struct ModelEntry
{
    Handle(AIS_Shape)     aisShape;
    std::string           name;
};
static std::vector<ModelEntry>  sModels;
static NSMutableArray<NSString*>* sModelNames = nil;

// ============================================================================
// OCC 交互 NSView（鼠标 + OpenGL 渲染）
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

typedef Handle(V3d_View)               HView;
typedef Handle(AIS_InteractiveContext) HCtx;

+ (void)setOCCView:(HView)v    { sOCCView   = v; }
+ (void)setOCCContext:(HCtx)c  { sOCCContext = c; }
+ (Handle(V3d_View))occView              { return sOCCView; }
+ (Handle(AIS_InteractiveContext))occContext { return sOCCContext; }

// ── 选中模型（高亮）──
+ (void)selectModelAtIndex:(NSInteger)idx
{
    if (sOCCContext.IsNull()) return;
    sOCCContext->ClearSelected(Standard_False);
    if (idx >= 0 && idx < (NSInteger)sModels.size())
    {
        Handle(AIS_Shape) sh = sModels[idx].aisShape;
        if (!sh.IsNull())
        {
            sOCCContext->AddOrRemoveSelected(sh, Standard_True);
        }
    }
}

+ (NSInteger)modelCount { return (NSInteger)sModels.size(); }

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
- (BOOL)isFlipped            { return YES; }

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
        sOCCView->Pan(Standard_Integer(aPt.x - myLastPoint.x),
                      Standard_Integer(myLastPoint.y - aPt.y));
    myLastPoint = aPt;
}
- (void)rightMouseDragged:(NSEvent*)e
{
    NSPoint aPt = [self convertPoint:[e locationInWindow] fromView:nil];
    if (self.isPanning)
        sOCCView->Pan(Standard_Integer(aPt.x - myLastPoint.x),
                      Standard_Integer(myLastPoint.y - aPt.y));
    myLastPoint = aPt;
}
- (void)otherMouseDragged:(NSEvent*)e  { [self rightMouseDragged:e]; }
- (void)mouseUp:(NSEvent*)e           { self.isRotating = Standard_False; }
- (void)rightMouseUp:(NSEvent*)e       { self.isPanning  = Standard_False; }
- (void)otherMouseUp:(NSEvent*)e       { self.isPanning  = Standard_False; }

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
// 模型列表面板（NSTableView）
// ============================================================================
@interface ModelListController : NSObject <NSTableViewDataSource, NSTableViewDelegate>
{
    NSTableView* tableView;
}
- (instancetype)initWithTableView:(NSTableView*)tv;
- (void)reload;
@end

@implementation ModelListController

- (instancetype)initWithTableView:(NSTableView*)tv
{
    self = [super init];
    if (self)
    {
        tableView = tv;
        tv.dataSource = self;
        tv.delegate   = self;
        [tv setTarget:self];
        [tv setAction:@selector(onTableClick:)];
    }
    return self;
}

- (void)reload { [tableView reloadData]; }

- (NSInteger)numberOfRowsInTableView:(NSTableView*)tv
{
    return [OCCInteractiveView modelCount];
}

- (id)tableView:(NSTableView*)tv
    objectValueForTableColumn:(NSTableColumn*)col row:(NSInteger)row
{
    if (row < 0 || row >= [OCCInteractiveView modelCount]) return @"";
    return [NSString stringWithUTF8String:sModels[row].name.c_str()];
}

- (void)onTableClick:(id)sender
{
    NSInteger row = [tableView clickedRow];
    if (row >= 0)
        [OCCInteractiveView selectModelAtIndex:row];
}

- (void)tableViewSelectionDidChange:(NSNotification*)note
{
    NSInteger row = [tableView selectedRow];
    if (row >= 0)
        [OCCInteractiveView selectModelAtIndex:row];
}
@end

// ============================================================================
// 列表管理器（全局引用，供按钮回调刷新列表）
// ============================================================================
static ModelListController* gListCtrl = nil;

// ============================================================================
// 顶部控制栏
// ============================================================================
@interface ViewerControlPanel : NSView
@end

@implementation ViewerControlPanel

- (instancetype)initWithFrame:(NSRect)frame
{
    self = [super initWithFrame:frame];
    if (self)
    {
        self.wantsLayer = YES;
        self.layer.backgroundColor = [[NSColor colorWithWhite:0.15 alpha:1.0] CGColor];

        auto mk = ^NSButton*(NSString* t, NSRect r, SEL a) {
            NSButton* b = [[NSButton alloc] initWithFrame:r];
            b.title = t; b.bezelStyle = NSBezelStyleRounded;
            [b setButtonType:NSButtonTypeMomentaryPushIn];
            [b setTarget:self]; [b setAction:a];
            [self addSubview:b];
            return b;
        };

        CGFloat x=12, y=4, w=72, h=28, g=8;
        mk(@"打开",   NSMakeRect(x, y, w, h), @selector(onOpen:));  x+=w+g;
        mk(@"线框",   NSMakeRect(x, y, w, h), @selector(onWire:));  x+=w+g;
        mk(@"着色",   NSMakeRect(x, y, w, h), @selector(onShade:)); x+=w+g;
        mk(@"透明度", NSMakeRect(x, y, w, h), @selector(onTrans:)); x+=w+g;
        mk(@"适配",   NSMakeRect(x, y, w, h), @selector(onFit:));
    }
    return self;
}

// ── 打开 STEP ──
- (void)onOpen:(id)sender
{
    NSOpenPanel* p = [NSOpenPanel openPanel];
    p.title = @"选择 STEP 文件";
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
    [p setAllowedFileTypes:@[@"step", @"stp"]];
#pragma clang diagnostic pop
    p.canChooseFiles = YES;
    p.canChooseDirectories = NO;
    p.allowsMultipleSelection = NO;
    if ([p runModal] != NSModalResponseOK) return;

    std::string filePath = [[[p URL] path] UTF8String];
    STEPControl_Reader reader;
    if (reader.ReadFile(filePath.c_str()) != IFSelect_RetDone) return;
    reader.TransferRoots();
    TopoDS_Shape shape = reader.OneShape();
    if (shape.IsNull()) return;

    auto ctx  = [OCCInteractiveView occContext];
    auto view = [OCCInteractiveView occView];
    if (ctx.IsNull() || view.IsNull()) return;

    Handle(AIS_Shape) ais = new AIS_Shape(shape);
    ais->SetMaterial(Graphic3d_NOM_ALUMINIUM);
    ais->SetTransparency(0.1);
    ais->SetColor(Quantity_NOC_STEELBLUE3);
    ctx->Display(ais, Standard_True);
    ctx->SetDisplayMode(ais, AIS_Shaded, Standard_False);

    // 记入模型列表
    std::string name = [[[p URL] lastPathComponent] UTF8String];
    sModels.push_back({ais, name});

    // 更新模型列表
    if (gListCtrl) [gListCtrl reload];

    view->FitAll();
    view->Redraw();
    std::cout << "+ " << name << std::endl;
}

- (void)onWire:(id)sender
{
    auto ctx = [OCCInteractiveView occContext];
    for (auto& m : sModels)
        if (!m.aisShape.IsNull())
            ctx->SetDisplayMode(m.aisShape, AIS_WireFrame, Standard_False);
    ctx->UpdateCurrentViewer();
    [OCCInteractiveView occView]->Redraw();
}
- (void)onShade:(id)sender
{
    auto ctx = [OCCInteractiveView occContext];
    for (auto& m : sModels)
        if (!m.aisShape.IsNull())
            ctx->SetDisplayMode(m.aisShape, AIS_Shaded, Standard_False);
    ctx->UpdateCurrentViewer();
    [OCCInteractiveView occView]->Redraw();
}
- (void)onTrans:(id)sender
{
    for (auto& m : sModels)
    {
        if (m.aisShape.IsNull()) continue;
        Standard_Real nt = (m.aisShape->Transparency() > 0.5) ? 0.1 : 0.8;
        m.aisShape->SetTransparency(nt);
        [OCCInteractiveView occContext]->Redisplay(m.aisShape, Standard_False);
    }
    [OCCInteractiveView occContext]->UpdateCurrentViewer();
    [OCCInteractiveView occView]->Redraw();
}
- (void)onFit:(id)sender
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

    Standard_Integer winW = 1100, winH = 636, ctrlH = 36, listW = 180;

    // ── 容器 ──
    NSRect aFrame = NSMakeRect(100, 100, winW, winH);
    NSView* container = [[NSView alloc] initWithFrame:aFrame];

    // ── 控制栏 ──
    NSRect ctrlFrame = NSMakeRect(0, winH - ctrlH, winW, ctrlH);
    ViewerControlPanel* panel = [[ViewerControlPanel alloc] initWithFrame:ctrlFrame];
    panel.autoresizingMask = NSViewWidthSizable | NSViewMinYMargin;
    [container addSubview:panel];

    // ── 左侧模型列表 ──
    NSRect listFrame = NSMakeRect(0, 0, listW, winH - ctrlH);
    NSScrollView* scrollView = [[NSScrollView alloc] initWithFrame:listFrame];
    scrollView.autoresizingMask = NSViewHeightSizable | NSViewMaxXMargin;
    scrollView.hasVerticalScroller = YES;
    scrollView.borderType = NSBezelBorder;

    NSTableView* tableView = [[NSTableView alloc] initWithFrame:scrollView.bounds];
    NSTableColumn* col = [[NSTableColumn alloc] initWithIdentifier:@"name"];
    col.title = @"模型";
    col.width = listW - 20;
    [tableView addTableColumn:col];
    tableView.headerView = [[NSTableHeaderView alloc] init];
    tableView.allowsMultipleSelection = NO;
    scrollView.documentView = tableView;

    ModelListController* listCtrl = [[ModelListController alloc] initWithTableView:tableView];
    gListCtrl = listCtrl;

    [container addSubview:scrollView];

    // ── OCC GL 视图 ──
    NSRect glFrame = NSMakeRect(listW, 0, winW - listW, winH - ctrlH);
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
    Handle(Aspect_DisplayConnection) dc = new Aspect_DisplayConnection();
    Handle(OpenGl_GraphicDriver)     dr = new OpenGl_GraphicDriver(dc);

    Handle(V3d_Viewer) vr = new V3d_Viewer(dr);
    vr->SetDefaultLights();
    vr->SetLightOn();
    Handle(AIS_InteractiveContext) ctx = new AIS_InteractiveContext(vr);

    Handle(Cocoa_Window) cw = new Cocoa_Window((__bridge NSView*)glView);
    cw->Map();
    cw->DoMapping();

    Handle(V3d_View) view = vr->CreateView();
    view->SetWindow(cw);
    [OCCInteractiveView setOCCView:view];
    [OCCInteractiveView setOCCContext:ctx];

    view->SetBackgroundColor(Quantity_NOC_GRAY30);
    view->SetBgGradientColors(Quantity_NOC_GRAY20, Quantity_NOC_GRAY60,
                              Aspect_GFM_VER, Standard_False);

    view->FitAll();
    view->Redraw();

    // ── 初始 shape ──
    if (!theShape.IsNull())
    {
        Handle(AIS_Shape) ais = new AIS_Shape(theShape);
        ais->SetMaterial(Graphic3d_NOM_ALUMINIUM);
        ais->SetTransparency(0.1);
        ais->SetColor(Quantity_NOC_STEELBLUE3);
        ctx->Display(ais, Standard_True);
        ctx->SetDisplayMode(ais, AIS_Shaded, Standard_False);
        sModels.push_back({ais, theTitle ? theTitle : "model"});
        [listCtrl reload];
        view->FitAll();
        view->Redraw();
    }

    [aNSWindow center];
    [aNSWindow makeKeyAndOrderFront:nil];
    [NSApp activateIgnoringOtherApps:YES];
    [NSApp run];
}
