//
//  ContentView.swift
//  Mimosa-Swift
//
//  Created by Ziyuan Qu on 2023/6/26.
//

import SwiftUI

struct ContentView: View {
    @State private var searchText = ""
    @State private var columnVisibility = NavigationSplitViewVisibility.detailOnly
    
    var body: some View {
//        EmptyView()
        NavigationSplitView(columnVisibility: $columnVisibility) {
            Sidebar()
        } detail: {
            ScenesUI()
        }
        .navigationTitle("All Boards")
//        .searchable(text: $searchText)
    }
}

struct ViewControllerView: NSViewControllerRepresentable {
    let sceneName: String
    let selectedIndex: Int
    
    func makeNSViewController(context: Context) -> NSViewController {
        let storyboard = NSStoryboard(name: "RenderScene", bundle: nil)
        if let viewController = storyboard.instantiateInitialController() as? ViewController {
            // pass sceneName to viewController
            viewController.sceneName = sceneName
            viewController.setRenderMode(Int32(selectedIndex))
            return viewController
        } else {
            fatalError("Failed to instantiate initial view controller from storyboard")
        }
    }
    
    func updateNSViewController(_ nsViewController: NSViewController, context: Context) {
        // update view controller
        if let viewController = nsViewController as? ViewController {
            viewController.setRenderMode(Int32(selectedIndex))
        }
    }
}

//struct BoardView: View {
//    let sceneName: String
//    @State private var isInspectorPresented = true
//    @State private var selectedIndex = 2
//    
//    var body: some View {
//        if #available(macOS 14.0, *) {
//            ViewControllerView(sceneName: sceneName, selectedIndex: selectedIndex)
//            .navigationTitle("Rendering")
//            .inspector(isPresented: $isInspectorPresented) {
//            InspectorView(selectedIndex: $selectedIndex)
//                .inspectorColumnWidth(min: 200, ideal: 300, max: 400)
//                .toolbar {
//                    Button(action: { isInspectorPresented.toggle() }) {
//                        Label("Toggle Inspector", systemImage: "sidebar.right")
//                    }
//                }
//            }
//        } else {
//            ZStack(alignment: .trailing) {
//                ViewControllerView(sceneName: sceneName, selectedIndex: selectedIndex)
//                    .navigationTitle("Rendering")
//                    .toolbar {
//                        Button(action: { isInspectorPresented.toggle() }) {
//                            Label("Toggle Inspector", systemImage: "sidebar.right")
//                        }
//                    }
//                
//                if isInspectorPresented {
//                    InspectorView(selectedIndex: $selectedIndex)
//                        .frame(width: 300)
//                        .transition(.move(edge: .trailing))
//                }
//            }
//        }
//    }
//}

struct RenderView: View {
    let sceneName: String
    @State private var isInspectorPresented = false
    @State private var selectedIndex = 2
    
    var body: some View {
        ZStack(alignment: .trailing) {
            ViewControllerView(sceneName: sceneName, selectedIndex: selectedIndex)
                .navigationTitle("Rendering")
                .toolbar {
                    Button(action: { isInspectorPresented.toggle() }) {
                        Label("Toggle Inspector", systemImage: "sidebar.right")
                    }
                }
            
            if isInspectorPresented {
                InspectorView(selectedIndex: $selectedIndex)
                    .frame(width: 300)
//                    .transition(.move(edge: .trailing))
            }
        }
    }
}


struct InspectorView: View {
    @Binding var selectedIndex: Int
//    @State private var selection: Int = 0
    
    var body: some View {
        List {
            Text("Integrators")
                .font(.headline)
            Picker(selection: $selectedIndex, label: EmptyView()) {
                Text("BRDF").tag(0)
                Text("NEE").tag(1)
                Text("MIS").tag(2)
            }
            .pickerStyle(SegmentedPickerStyle())
//            .onChange(of: selectedIndex){
//                print("change selection")
//            }
            Divider()
//            Picker(selection: $selection, label: Text("Picker Label")) {
//                Text("Option 1").tag(0)
//                Text("Option 2").tag(1)
//                Text("Option 3").tag(2)
//            }
        }
        .listStyle(SidebarListStyle())
//        .frame(minWidth: 100)
    }
}

//struct ContentView_Previews: PreviewProvider {
//    static var previews: some View {
//        ContentView()
//    }
//}

struct Sidebar: View {
    var body: some View {
        List {
//            NavigationLink(destination: Text("Item 1")) {
//                Label("Item 1", systemImage: "star.fill")
//            }
//            NavigationLink(destination: Text("Item 2")) {
//                Label("Item 2", systemImage: "star.fill")
//            }
//            NavigationLink(destination: Text("Item 3")) {
//                Label("Item 3", systemImage: "star.fill")
//            }
            EmptyView()
        }
        .listStyle(SidebarListStyle())
        .frame(minWidth: 100)
    }
}

struct ScenesUI: View {
    @State private var isCardView = true
    
    let scenes = [
        ("Scene 1", "1 hour ago", "01"),
        ("Scene 2", "2 hours ago", "02"),
        ("Scene 3", "3 hours ago", "03")
    ]
    
    var body: some View {
        NavigationStack {
            if isCardView {
                ScrollView {
                    LazyVGrid(columns: [GridItem(.adaptive(minimum: 210))]) {
                        ForEach(scenes, id: \.0) { scene in
                            NavigationStack() {
                                NavigationLink(destination:RenderView(sceneName: scene.0)){
                                    ZStack(alignment: .bottom) {
                                        if let nsImage = NSImage(named: scene.2) {
                                            Image(nsImage: nsImage)
                                            .resizable()
                                            .scaledToFit()
//                                            .aspectRatio(contentMode: .fit)
//                                            .frame(width: 200, height: 200)
                                        } else {
                                            Image(systemName: "photo")
                                           .resizable()
                                           .frame(width: 200, height: 200)
                                        }
                                        
                                        Rectangle()
                                            .fill(Color.gray.opacity(0.7))
                                            .frame(height: 50)
                                        
                                        VStack(alignment: .leading) {
                                            Text(scene.0)
                                                .font(.headline)
                                                .foregroundColor(.white)
                                            Text("Last opened: \(scene.1)")
                                                .font(.subheadline)
                                                .foregroundColor(.white)
                                        }
                                        .padding()
                                    }
                                    .cornerRadius(10)
                                    .padding()
//                                    .frame(width: 200, height: 200)
                                }
                                .buttonStyle(PlainButtonStyle())
                            }
                        }
                    }
                }
                .background(Color.white)
            } else {
                List {
                    ForEach(scenes, id: \.0) { scene in
                        NavigationStack() {
                            NavigationLink(destination:RenderView(sceneName: scene.0)) {
                                if let nsImage = NSImage(named: scene.2) {
                                    Image(nsImage: nsImage)
                                    .resizable()
                                    .scaledToFit()
                                    .frame(width: 50, height: 50)
                                    .cornerRadius(5)
                                } else {
                                    Image(systemName: "photo")
                                   .resizable()
                                   .frame(width: 50, height: 50)
                                   .cornerRadius(5)
                                }
                                
                                VStack(alignment: .leading) {
                                    Text(scene.0)
                                    Text("Last opened: \(scene.1)")
                                        .font(.subheadline)
                                        .foregroundColor(.secondary)
                                }
                            }
                        }
                    }
                }
            }
        }
        .navigationTitle("All Scenes")
        .toolbar {
            Button(action: { isCardView.toggle() }) {
                Image(systemName: isCardView ? "list.bullet" : "square.grid.2x2")
            }
        }
    }
}

//
//#Preview {
//    ContentView()
//}
