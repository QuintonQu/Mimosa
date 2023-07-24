# Mimosa
Quinton's Ray Tracer with Metal

## Introduction
This is an exercise project for writing ray tracing in Metal shader. The project is based on sample code [Accelerating ray tracing using Metal](https://developer.apple.com/documentation/metal/accelerating_ray_tracing_using_metal) from Apple. The code is written in Objective-C++ and Metal Shading Language. Hopefully, I will add Metal-CPP support to this project in the future.

## App Viewer
The app viewer is written in Swift. It is a simple app that can load a scene and render it. You can simply choose different scenes and path tracing algorithms.
<p float="left">
  <img src="./Results/App View/01.png" height="200" />
  <img src="./Results/App View/06.png" height="200" />
  <img src="./Results/App View/02.png" height="200" />
  <img src="./Results/App View/03.png" height="200" />
  <img src="./Results/App View/04.png" height="200" />
  <img src="./Results/App View/07.png" height="200" />
</p>
Here's a little video to help you learn how to operate the app.
<p float="left">
  <img src="./Results/App View/08.gif" width="400" />
  <img src="./Results/App View/09.gif" width="400" />
  <img src="./Results/App View/10.gif" width="400" />
</p>

## Features
### Geometry
- [x] OpenVDB
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
- [x] Volumetric rendering
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
