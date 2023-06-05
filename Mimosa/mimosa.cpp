//
//  main.cpp
//  Mimosa
//
//  Created by Ziyuan Qu on 2023/6/5.
//

#define NS_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>
#include <QuartzCore/QuartzCore.hpp>
#include <iostream>

int main(int argc, const char * argv[]) {
    // insert code here...
    MTL::Device* device = MTL::CreateSystemDefaultDevice();
    std::cout << "Hello, Mimosa!" << std::endl;
    device -> release();
    return 0;
}
