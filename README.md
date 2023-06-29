# Mimosa
Quinton's Ray Tracer with Metal

## Introduction
This is an exercise project for writing ray tracing in Metal shader. The project is based on sample code [Accelerating ray tracing using Metal](https://developer.apple.com/documentation/metal/accelerating_ray_tracing_using_metal) from Apple. The sample code is written in Objective-C and Metal Shading Language. Hopefully, I will add Metal-CPP support to this project in the future.

## Features
### Geometry
- [ ] OpenVDB
- [x] Obj file loader
- [x] Sphere
- [x] (One side) Triangle
- [x] (One side) Cube / Plane with triangle

### Material
- [ ] Blinn-phong
- [x] Phong
- [x] Diffuse
- [x] Metallic
- [x] Dielectric

### Light
- [ ] Point light
- [ ] Important Sampling Environment map
- [x] Sky box / Environment map
- [x] Area light (Mesh light)

### Texture
- [ ] Checker texture
- [ ] Perlin noise
- [ ] Texture mapping
- [ ] Normal mapping
- [ ] Displacement mapping

### Path Tracer
- [ ] Volumetric rendering
- [x] Multiple importance sampling
- [x] Naive Path tracer
- [x] Next event estimation

### Other Functions
- [ ] Moving Camera
- [ ] IOS Deployment
- [x] Choosing different scene

## References
- [Accelerating ray tracing using Metal](https://developer.apple.com/documentation/metal/accelerating_ray_tracing_using_metal)
- [Darts (Wojciech Jarosz)](https://cs87-dartmouth.github.io/Fall2021/darts-overview.html)
