/*
See LICENSE folder for this sample’s licensing information.

Abstract:
The Metal shaders used for this sample.
*/
#include "ShaderTypes.h"

#include <metal_stdlib>
#include <simd/simd.h>

using namespace metal;

using namespace raytracing;

constant unsigned int resourcesStride   [[function_constant(0)]];
constant bool useIntersectionFunctions  [[function_constant(1)]];
constant bool usePerPrimitiveData       [[function_constant(2)]];
//constant bool useResourcesBuffer = !usePerPrimitiveData
constant int max_bounce = 8;

constant unsigned int primes[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131,
    137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311,
    313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503,
    509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719,
    727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
    947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163,
    1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423,
    1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601, 1607, 1609, 1613, 1619,
    1621, 1627, 1637, 1657, 1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877,
    1879, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129,
    2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, 2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 2371, 2377,
    2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621, 2633, 2647, 2657,
    2659, 2663, 2671, 2677, 2683, 2687, 2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 2833, 2837, 2843, 2851, 2857, 2861,
    2879, 2887, 2897, 2903, 2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167,
    3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413,
    3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571, 3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, 3659, 3671,
    3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, 3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, 3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917, 3919, 3923,
    3929, 3931, 3943, 3947, 3967, 3989, 4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, 4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, 4153, 4157, 4159, 4177, 4201, 4211,
    4217, 4219, 4229, 4231, 4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, 4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481,
    4483, 4493, 4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751,
    4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, 4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 5009, 5011,
    5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279, 5281, 5297, 5303, 5309,
    5323, 5333, 5347, 5351, 5381, 5387, 5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, 5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563, 5569, 5573,
    5581, 5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, 5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, 5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849,
    5851, 5857, 5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, 5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133,
    6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, 6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, 6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 6373, 6379,
    6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 6679, 6689, 6691, 6701,
    6703, 6709, 6719, 6733, 6737, 6761, 6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, 6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 6947, 6949, 6959, 6961, 6967, 6971,
    6977, 6983, 6991, 6997, 7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253,
    7283, 7297, 7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, 7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561,
    7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, 7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 7841, 7853,
    7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919, 7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, 8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, 8117, 8123, 8147, 8161
};

#pragma mark - Generate Random Sample
// Returns the i'th element of the Halton sequence using the d'th prime number as a
// base. The Halton sequence is a low discrepency sequence: the values appear
// random, but are more evenly distributed than a purely random sequence. Each random
// value used to render the image uses a different independent dimension, `d`,
// and each sample (frame) uses a different index `i`. To decorrelate each pixel,
// you can apply a random offset to `i`.
float halton(unsigned int i, unsigned int d) {
    unsigned int b = primes[d % 1024];

    float f = 1.0f;
    float invB = 1.0f / b;

    float r = 0;

    while (i > 0) {
        f = f * invB;
        r = r + f * (i % b);
        i = i / b;
    }

    return r;
}

#pragma mark - Inline Utility Interpolate
// Interpolates the vertex attribute of an arbitrary type across the surface of a triangle
// given the barycentric coordinates and triangle index in an intersection structure.
template<typename T, typename IndexType>
inline T interpolateVertexAttribute(device T *attributes,
                                    IndexType i0,
                                    IndexType i1,
                                    IndexType i2,
                                    float2 uv) {
    // Look up value for each vertex.
    const T T0 = attributes[i0];
    const T T1 = attributes[i1];
    const T T2 = attributes[i2];

    // Compute the sum of the vertex attributes weighted by the barycentric coordinates.
    // The barycentric coordinates sum to one.
    return (1.0f - uv.x - uv.y) * T0 + uv.x * T1 + uv.y * T2;
}

template<typename T>
inline T interpolateVertexAttribute(thread T *attributes, float2 uv) {
    // Look up the value for each vertex.
    const T T0 = attributes[0];
    const T T1 = attributes[1];
    const T T2 = attributes[2];

    // Compute the sum of the vertex attributes weighted by the barycentric coordinates.
    // The barycentric coordinates sum to one.
    return (1.0f - uv.x - uv.y) * T0 + uv.x * T1 + uv.y * T2;
}

#pragma mark - Inline Utility Inverse

inline float4x4 inverse(float4x4 matrix) {
//    float det = determinant(matrix);
    float4x4 inv;
    
    float E_Matrix[4][4];
    float mik;
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            if(i == j)
                E_Matrix[i][j] = 1.f;
            else
                E_Matrix[i][j] = 0.f;
        }
    }
    float CalcuMatrix[4][8];
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            CalcuMatrix[i][j] = matrix[i][j];
        }
        for(int k = 4; k < 8; k++)
        {
            CalcuMatrix[i][k] = E_Matrix[i][k-4];
        }
    }

    for(int i = 1; i <= 4-1; i++)
    {
        for(int j = i+1; j <= 4; j++)
        {
            mik = CalcuMatrix[j-1][i-1]/CalcuMatrix[i-1][i-1];
            for(int k = i+1;k <= 8; k++)
            {
                CalcuMatrix[j-1][k-1] -= mik*CalcuMatrix[i-1][k-1];
            }
        }
    }
    for(int i=1;i<=4;i++)
    {
        float temp = CalcuMatrix[i-1][i-1];
        for(int j=1;j<=8;j++)
        {
            CalcuMatrix[i-1][j-1] = CalcuMatrix[i-1][j-1]/temp;
        }
    }
    for(int k=4-1;k>=1;k--)
    {
        for(int i=k;i>=1;i--)
        {
            mik = CalcuMatrix[i-1][k];
            for(int j=k+1;j<=8;j++)
            {
                CalcuMatrix[i-1][j-1] -= mik*CalcuMatrix[k][j-1];
            }
        }
    }
    float InverseMatrix[4][4];
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            InverseMatrix[i][j] = CalcuMatrix[i][j+4];
        }
    }

    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            inv[i][j] = InverseMatrix[i][j];
        }
    }

    return inv;
}

#pragma mark - Inline Utility Transform

__attribute__((always_inline))
float3 transformPoint(float3 p, float4x4 transform) {
    return (transform * float4(p.x, p.y, p.z, 1.0f)).xyz;
}

__attribute__((always_inline))
float3 transformDirection(float3 p, float4x4 transform) {
    return (transform * float4(p.x, p.y, p.z, 0.0f)).xyz;
}

#pragma mark - Inline Utility Sample
// Uses the inversion method to map two uniformly random numbers to a 3D
// unit hemisphere, where the probability of a given sample is proportional to the cosine
// of the angle between the sample direction and the "up" direction (0, 1, 0).
inline float3 sampleCosineWeightedHemisphere(float2 u) {
    float phi = 2.0f * M_PI_F * u.x;

    float cos_phi;
    float sin_phi = sincos(phi, cos_phi);

    float cos_theta = sqrt(u.y);
    float sin_theta = sqrt(1.0f - cos_theta * cos_theta);

    return float3(sin_theta * cos_phi, cos_theta, sin_theta * sin_phi);
}

// Sample uniform sphere
inline float3 sampleUniformSphere(float2 u) {
    float z = 1.0f - 2.0f * u.x;
    float r = sqrt(max(0.0f, 1.0f - z * z));
    float phi = 2.0f * M_PI_F * u.y;
    float x = r * cos(phi);
    float y = r * sin(phi);
    return float3(x, y, z);
}


// Sample cosine weighted hemisphere with power
inline float3 sampleCosineWeightedHemisphereWithPower(float2 u, float exponent) {
    float phi = 2.0f * M_PI_F * u.x;

    float cos_phi;
    float sin_phi = sincos(phi, cos_phi);

    float cos_theta = pow(u.y, 1/(exponent + 1));
    float sin_theta = sqrt(1.0f - cos_theta * cos_theta);

    return float3(sin_theta * cos_phi, cos_theta, sin_theta * sin_phi);
}

// Sample sphere cap
inline float3 sampleSphereCap(float2 u, float cos_theta_max) {
    float phi = 2.0f * M_PI_F * u.x;

    float cos_phi;
    float sin_phi = sincos(phi, cos_phi);

    float cos_theta = mix(cos_theta_max, 1.0f, u.y);
    float sin_theta = sqrt(1.0f - cos_theta * cos_theta);

    return float3(sin_theta * cos_phi, cos_theta, sin_theta * sin_phi);
}

// Maps two uniformly random numbers to the surface of a 2D area light
// source and returns the direction to this point, the amount of light that travels
// between the intersection point and the sample point on the light source, as well
// as the distance between these two points.

inline void sampleAreaLight(constant AreaLight & light,
                            float2 u,
                            float3 position,
                            thread float3 & lightDirection,
                            thread float3 & lightColor,
                            thread float & lightDistance)
{
    // Map to -1..1
    u = u * 2.0f - 1.0f;

    // Transform into the light's coordinate system.
    float3 samplePosition = light.position +
                            light.right * u.x +
                            light.up * u.y;

    // Compute the vector from sample point on  the light source to intersection point.
    lightDirection = samplePosition - position;

    lightDistance = length(lightDirection);

    float inverseLightDistance = 1.0f / max(lightDistance, 1e-3f);

    // Normalize the light direction.
    lightDirection *= inverseLightDistance;

    // Start with the light's color.
    lightColor = light.color;

    // Light falls off with the inverse square of the distance to the intersection point.
    lightColor *= (inverseLightDistance * inverseLightDistance);

    // Light also falls off with the cosine of the angle between the intersection point
    // and the light source.
    lightColor *= saturate(dot(-lightDirection, light.forward));
}

// Return the type for Scatter Record for different light.
struct EmitterRecord {
    bool accept;
    float3 out_direction;
    float3 emit;
    float pdf;
    float distance;
    bool test = false;
};

// Aligns a direction on the unit hemisphere such that the hemisphere's "up" direction
// (0, 1, 0) maps to the given surface normal direction.
inline float3 alignHemisphereWithNormal(float3 sample, float3 normal) {
    // Set the "up" vector to the normal
    float3 up = normal;

    // Find an arbitrary direction perpendicular to the normal, which becomes the
    // "right" vector.
    float3 right = normalize(cross(normal, float3(0.0072f, 1.0f, 0.0034f)));

    // Find a third vector perpendicular to the previous two, which becomes the
    // "forward" vector.
    float3 forward = cross(right, up);

    // Map the direction on the unit hemisphere to the coordinate system aligned
    // with the normal.
    return sample.x * right + sample.y * up + sample.z * forward;
}

// Return the type for a bounding box intersection function.
struct BoundingBoxIntersection {
    bool accept    [[accept_intersection]]; // Whether to accept or reject the intersection.
    float distance [[distance]];            // Distance from the ray origin to the intersection point.
};

// Sample a sphere light
inline EmitterRecord sampleSphereLight(Sphere sphere, float3 hit_point, float2 random_variable, float4x4 transform){
    EmitterRecord emitter_record;
    float3 origin = transformPoint(sphere.origin, transform);
    float distance = length(hit_point - origin);
    float radius = sphere.radius;
    if(distance < radius){
        emitter_record.accept = false;
        emitter_record.pdf = 0.f;
        return emitter_record;
    }
    float cos_max = sqrt(distance * distance - radius * radius) / distance;
    emitter_record.out_direction = alignHemisphereWithNormal(sampleSphereCap(random_variable, cos_max), normalize(origin - hit_point));
    float area = 2 * M_PI_F * (1-cos_max);
    emitter_record.pdf = abs(1.f / area);
    emitter_record.emit = sphere.material.color;
    emitter_record.distance = distance - radius;
    return emitter_record;
}

inline float spherePdf(Sphere sphere, float3 hit_point, float4x4 transform){
    float3 origin = transformPoint(sphere.origin, transform);
    float distance = length(hit_point - origin);
    float radius = sphere.radius;
    if(distance < radius){
        return 0.f;
    }
    float cos_max = sqrt(distance * distance - radius * radius) / distance;
    float area = 2 * M_PI_F * (1-cos_max);
    return abs(1.f / area);
}

// Sample a triangle light
inline EmitterRecord sampleTriangleLight(Triangle triangle, float3 hit_point, float2 random_variable, float4x4 transform){
    EmitterRecord emitter_record;
    float alpha = random_variable.x;
    float beta = random_variable.y;
    if (alpha + beta > 1){
        alpha = 1 - alpha;
        beta = 1 - beta;
    }
    float2 barycentric_coords = float2(alpha, beta);
    float3 objectSpaceSurfaceNormal = interpolateVertexAttribute(triangle.normals, barycentric_coords);
    float3 origin = interpolateVertexAttribute(triangle.positions, barycentric_coords);
    origin = transformPoint(origin, transform);
    float3 in_direction = hit_point - origin;
    
    emitter_record.out_direction = normalize(-in_direction);
    float3 worldSpaceSurfaceNormal = transformDirection(objectSpaceSurfaceNormal, transform);
    
    emitter_record.emit = triangle.material.color;
    
    float3 a = transformPoint(triangle.positions[0], transform) - transformPoint(triangle.positions[2], transform);
    float3 b = transformPoint(triangle.positions[1], transform) - transformPoint(triangle.positions[2], transform);
    float area = length(cross(a, b)) / 2.f;
    float cosine = abs(dot(emitter_record.out_direction, worldSpaceSurfaceNormal));
    emitter_record.distance = length(in_direction);
    emitter_record.pdf = length_squared(in_direction) / (cosine * area);
    return emitter_record;
}

inline float trianglePdf(Triangle triangle, float3 hit_point, float4x4 transform, float2 barycentric_coords){
    float3 objectSpaceSurfaceNormal = interpolateVertexAttribute(triangle.normals, barycentric_coords);
    float3 worldSpaceSurfaceNormal = transformDirection(objectSpaceSurfaceNormal, transform);
    float3 origin = interpolateVertexAttribute(triangle.positions, barycentric_coords);
    origin = transformPoint(origin, transform);
    float3 in_direction = hit_point - origin;
    float3 a = transformPoint(triangle.positions[0], transform) - transformPoint(triangle.positions[2], transform);
    float3 b = transformPoint(triangle.positions[1], transform) - transformPoint(triangle.positions[2], transform);
    float area = length(cross(a, b)) / 2.f;
    float cosine = abs(dot(in_direction, worldSpaceSurfaceNormal));
    return length_squared(in_direction) / (cosine * area);
}

# pragma mark - Materials' Scatter Methods

// Return the type for Scatter Record for different material.
struct ScatterRecord {
    bool accept;
    float3 out_direction;
    float3 attenuation;
    float pdf;
    bool test = false;
    // for glass
    bool is_refract = false;
};

// Material scatter function -> Metallic
inline ScatterRecord metallicScatter(float3 normal, float3 in_direction) {
    ScatterRecord scatter_record;
    
    scatter_record.out_direction = reflect(in_direction, normal);
    scatter_record.attenuation = float3(1.0f, 1.0f, 1.0f);
    scatter_record.accept = true;
    
    return scatter_record;
}

// Fresnel Dieletric
inline float fresnel_dielectric(float cos_theta_i, float eta_i, float eta_t)
{
    if (eta_t == eta_i)
        return 0.f;

    // Swap the indices of refraction if the interaction starts at the inside of the object
    bool entering = cos_theta_i > 0.0f;
    if (!entering)
    {
        float eta = 0.0f;
        eta = eta_i;
        eta_i = eta_t;
        eta_t = eta;
        cos_theta_i = -cos_theta_i;
    }

    // Using Sahl-Snell's law, calculate the squared sine of the angle between the normal and the transmitted ray
    float eta          = eta_i / eta_t;
    float sin_theta_t2 = eta * eta * (1 - cos_theta_i * cos_theta_i);

    // Total internal reflection!
    if (sin_theta_t2 > 1.0f)
        return 1.0f;

    float cos_theta_t = sqrt(1.0f - sin_theta_t2);

    float Rs = (eta_i * cos_theta_i - eta_t * cos_theta_t) / (eta_i * cos_theta_i + eta_t * cos_theta_t);
    float Rp = (eta_t * cos_theta_i - eta_i * cos_theta_t) / (eta_t * cos_theta_i + eta_i * cos_theta_t);

    return 0.5f * (Rs * Rs + Rp * Rp);
}


// Glossy scatter function -> Glossy
inline ScatterRecord glossyScatter(float3 normal, float3 in_direction, float random_variable) {
    ScatterRecord scatter_record;
    
    float ior = 1.5f;
    float cosine = dot(normal, in_direction);
    float reflected = fresnel_dielectric(cosine, ior, 1.0f);
    if (cosine > 0.0f){
        normal = -normal;
        ior = 1.0f / ior;
    }
    
    if(random_variable < reflected){
        scatter_record.out_direction = reflect(in_direction, normal);
        scatter_record.attenuation = float3(1.0f, 1.0f, 1.0f);
    } else {
        scatter_record.out_direction = refract(in_direction, normal, 1.0f/ior);
        scatter_record.attenuation = float3(1.0f, 1.0f, 1.0f);
        scatter_record.is_refract = true;
    }
    scatter_record.accept = true;
    
    return scatter_record;
}

// Phong Material
inline ScatterRecord phongScatter(float3 normal, float3 in_direction, float2 random_variable, Material material ){
    ScatterRecord scatter_record;
    
    float3 reflection_direction = reflect(in_direction, normal);
    float3 random_sample_direction = sampleCosineWeightedHemisphereWithPower(random_variable, material.exponent);
    scatter_record.out_direction = alignHemisphereWithNormal(random_sample_direction, reflection_direction);
    scatter_record.attenuation = material.color;
    scatter_record.accept = dot(in_direction, normal) < 0.f && dot(scatter_record.out_direction, normal) > 0.f;
    
    return scatter_record;
}

inline float phongPdf(float3 normal, float3 in_direction, float3 out_direction, Material material){
    if(dot(in_direction, normal) > 0.f || dot(out_direction, normal) < 0.f){
        return 0.f;
    }
    float3 mirror_direction = normalize(reflect(in_direction, normal));
    float cosine = saturate(dot(out_direction, mirror_direction));
    return ((material.exponent + 1.f)/(2.f*M_PI_F))*pow(cosine, material.exponent);
}

// Resources for a piece of triangle geometry.
struct TriangleResources {
//    device uint16_t *indices;
//    device float3 *vertexNormals;
//    device float3 *vertexColors;
//    device float3 *vertexPositions;
    device Triangle* triangles;
};

// Resources for a piece of sphere geometry.
struct SphereResources {
    device Sphere *spheres;
};

#pragma mark - Sphere Intersection
/*
 Custom sphere intersection function. The [[intersection]] keyword marks this as an intersection
 function. The [[bounding_box]] keyword means that this intersection function handles intersecting rays
 with bounding box primitives. To create sphere primitives, the sample creates bounding boxes that
 enclose the sphere primitives.

 The [[triangle_data]] and [[instancing]] keywords indicate that the intersector that calls this
 intersection function returns barycentric coordinates for triangle intersections and traverses
 an instance acceleration structure. These keywords must match between the intersection functions,
 intersection function table, intersector, and intersection result to ensure that Metal propagates
 data correctly between stages. Using fewer tags when possible may result in better performance,
 as Metal may need to store less data and pass less data between stages. For example, if you do not
 need barycentric coordinates, omitting [[triangle_data]] means Metal can avoid computing and storing
 them.

 The arguments to the intersection function contain information about the ray, primitive to be
 tested, and so on. The ray intersector provides this datas when it calls the intersection function.
 Metal provides other built-in arguments, but this sample doesn't use them.
 */
[[intersection(bounding_box, triangle_data, instancing)]]
BoundingBoxIntersection sphereIntersectionFunction(// Ray parameters passed to the ray intersector below
                                                   float3 origin                        [[origin]],
                                                   float3 direction                     [[direction]],
                                                   float minDistance                    [[min_distance]],
                                                   float maxDistance                    [[max_distance]],
                                                   // Information about the primitive.
                                                   unsigned int primitiveIndex          [[primitive_id]],
                                                   unsigned int geometryIndex           [[geometry_intersection_function_table_offset]],
                                                   // Custom resources bound to the intersection function table.
//                                                   device void *resources               [[buffer(0), function_constant(useResourcesBuffer)]]
//#if SUPPORTS_METAL_3
                                                   const device void* perPrimitiveData [[primitive_data]]
//#endif
                                                   )
{
    Sphere sphere;
    sphere = *(const device Sphere*)perPrimitiveData;
    
//#if SUPPORTS_METAL_3
//    // Look up the resources for this piece of sphere geometry.
//    if (usePerPrimitiveData) {
//        // Per-primitive data points to data from the specified buffer as was configured in the MTLAccelerationStructureBoundingBoxGeometryDescriptor.
//        sphere = *(const device Sphere*)perPrimitiveData;
//    } else
//#endif
//    {
//        device SphereResources& sphereResources = *(device SphereResources *)((device char *)resources + resourcesStride * geometryIndex);
//        // Get the actual sphere enclosed in this bounding box.
//        sphere = sphereResources.spheres[primitiveIndex];
//    }

    // Check for intersection between the ray and sphere mathematically.
    float3 oc = origin - sphere.origin;

    float a = dot(direction, direction);
    float b = 2 * dot(oc, direction);
    float c = dot(oc, oc) - sphere.radiusSquared;

    float disc = b * b - 4 * a * c;

    BoundingBoxIntersection ret;

    if (disc <= 0.0f) {
        // If the ray missed the sphere, return false.
        ret.accept = false;
    }
    else {
        // Otherwise, compute the intersection distance.
        ret.distance = (-b - sqrt(disc)) / (2 * a);

        // The intersection function must also check whether the intersection distance is
        // within the acceptable range. Intersection functions do not run in any particular order,
        // so the maximum distance may be different from the one passed into the ray intersector.
        
        // Modify to implement the transparent glass with refrection. Both root are needed to check.
        if (ret.distance > 1e-3f){
            ret.accept = ret.distance >= minDistance && ret.distance <= maxDistance;
        } else {
            ret.distance = (-b + sqrt(disc)) / (2 * a);
            ret.accept = ret.distance >= minDistance && ret.distance <= maxDistance;
        }
    }
    return ret;
}

#pragma mark - Generate Ray
#pragma region GenerateRay {

ray generate_ray(uint2 tid [[thread_position_in_grid]], constant Uniforms &uniforms [[buffer(0)]], texture2d<unsigned int> randomTex [[texture(0)]]) {
    // The ray to cast.
    ray ray;

    // Pixel coordinates for this thread.
    float2 pixel = (float2)tid;

    // Apply a random offset to the random number index to decorrelate pixels.
    unsigned int offset = randomTex.read(tid).x;

    // Add a random offset to the pixel coordinates for antialiasing.
    float2 r = float2(halton(offset + uniforms.frameIndex, 0),
                      halton(offset + uniforms.frameIndex, 1));

    pixel += r;

    // Map pixel coordinates to -1..1.
    float2 uv = (float2)pixel / float2(uniforms.width, uniforms.height);
    uv = uv * 2.0f - 1.0f;

    constant Camera & camera = uniforms.camera;

    // Rays start at the camera position.
    ray.origin = camera.position;

    // Map normalized pixel coordinates into camera's coordinate system.
    ray.direction = normalize(uv.x * camera.right +
                              uv.y * camera.up +
                              camera.forward);

    // Don't limit intersection distance.
    ray.max_distance = INFINITY;
    ray.min_distance = 1e-3;
    
    return ray;
}

#pragma endregion GenerateRay }

#pragma mark - Environment Map
#pragma region EnvironmentMap {

float3 get_env_color(float3 in_direction, texture2d<float> envmapTex [[texture(3)]]){
    float2 environmentMapUV = float2(atan2(in_direction.z, in_direction.x) / (2 * M_PI_F) + 0.5, 0.5 - asin(in_direction.y) / M_PI_F);
    constexpr sampler linearFilterSampler(coord::normalized, address::clamp_to_edge, filter::linear);
    float3 color = envmapTex.sample(linearFilterSampler, environmentMapUV).rgb;
    return color;
}

EmitterRecord uniform_sample_env(float2 random_variable, texture2d<float> envmapTex [[texture(3)]]){
    EmitterRecord emitter_record;
    
    float u = random_variable.x * 2 * M_PI_F;
    float v = acos(2 * random_variable.y - 1);
    
    emitter_record.out_direction = normalize(float3(-sin(v) * cos(u), cos(v), -sin(v) * sin(u)));
    emitter_record.emit = get_env_color(emitter_record.out_direction, envmapTex);
    emitter_record.distance = INFINITY;
    emitter_record.pdf = 1 / (4 * M_PI_F);
    
    return emitter_record;
}

//EmitterRecord uniform_sample_env(float2 random_variable, texture2d<float> envmapTex [[texture(3)]]){
//    EmitterRecord emitter_record;
//
//    float u = random_variable.x * 2 * M_PI_F;
//    float v = random_variable.y * M_PI_F;
//
//    constexpr sampler linearFilterSampler(coord::normalized, address::clamp_to_edge, filter::linear);
//    emitter_record.emit = envmapTex.sample(linearFilterSampler, random_variable).rgb * cos(v - M_PI_F / 2);
//    emitter_record.out_direction = normalize(float3(-sin(v) * cos(u), cos(v), -sin(v) * sin(u)));
//    emitter_record.distance = INFINITY;
//    emitter_record.pdf = 1 / (4 * M_PI_F);
//
//    return emitter_record;
//}

#pragma endregion EnvironmentMap }

#pragma mark - VDB Data
half4 get_vdb_data(float3 local_position, texture3d<half> vdbvolTex [[texture(4)]]){
    float3 uv = local_position + float3(0.5f);
    constexpr sampler linearFilterSampler(coord::normalized, address::clamp_to_edge, filter::linear);
    return vdbvolTex.sample(linearFilterSampler, uv);
}

#pragma mark - Ray Tracing Kernel - Mats

// Main ray tracing kernel.
kernel void raytracingKernelMats(
     uint2                                                  tid                       [[thread_position_in_grid]],
     constant Uniforms &                                    uniforms                  [[buffer(0)]],
     texture2d<unsigned int>                                randomTex                 [[texture(0)]],
     texture2d<float>                                       prevTex                   [[texture(1)]],
     texture2d<float, access::write>                        dstTex                    [[texture(2)]],
     texture2d<float>                                       envmapTex                 [[texture(3)]],
//     device void                                           *resources                 [[buffer(1), function_constant(useResourcesBuffer)]],
     constant MTLAccelerationStructureInstanceDescriptor   *instances                 [[buffer(2)]],
     constant AreaLight                                    *areaLights                [[buffer(3)]],
     instance_acceleration_structure                        accelerationStructure     [[buffer(4)]],
     intersection_function_table<triangle_data, instancing> intersectionFunctionTable [[buffer(5)]]
)
{
    // The sample aligns the thread count to the threadgroup size, which means the thread count
    // may be different than the bounds of the texture. Test to make sure this thread
    // is referencing a pixel within the bounds of the texture.
    if (tid.x > uniforms.width || tid.y > uniforms.height) return;
    
    // Apply a random offset to the random number index to decorrelate pixels.
    unsigned int offset = randomTex.read(tid).x;
    
    // Generate Ray.
    ray ray = generate_ray(tid, uniforms, randomTex);

    // Start with a fully white color. The kernel scales the light each time the
    // ray bounces off of a surface, based on how much of each light component
    // the surface absorbs.
    float3 color = float3(1.0f, 1.0f, 1.0f);

    float3 accumulatedColor = float3(0.0f, 0.0f, 0.0f);

    // Create an intersector to test for intersection between the ray and the geometry in the scene.
    intersector<triangle_data, instancing> i;

    // If the sample isn't using intersection functions, provide some hints to Metal for
    // better performance.
//    if (!useIntersectionFunctions) {
//        i.assume_geometry_type(geometry_type::triangle);
//        i.force_opacity(forced_opacity::opaque);
//    }

    typename intersector<triangle_data, instancing>::result_type intersection;
    
    ScatterRecord scatter_record;

    // Simulate up to three ray bounces. Each bounce propagates light backward along the
    // ray's path toward the camera.
    for (int bounce = 0; bounce < max_bounce + 1; bounce++) {
        // Get the closest intersection, not the first intersection. This is the default, but
        // the sample adjusts this property below when it casts shadow rays.
        i.accept_any_intersection(false);

        // Check for intersection between the ray and the acceleration structure. If the sample
        // isn't using intersection functions, it doesn't need to include one.
        if (useIntersectionFunctions)
            intersection = i.intersect(ray, accelerationStructure, RAY_MASK_PRIMARY, intersectionFunctionTable);
        else
            intersection = i.intersect(ray, accelerationStructure, bounce == 0 ? RAY_MASK_PRIMARY : RAY_MASK_SECONDARY);

        // Stop if the ray didn't hit anything and has bounced out of the scene.
        if (intersection.type == intersection_type::none) {
            accumulatedColor = get_env_color(ray.direction, envmapTex) * color;
            break;
        }
            
        unsigned int instanceIndex = intersection.instance_id;

        // Look up the mask for this instance, which indicates what type of geometry the ray hit.
        unsigned int mask = instances[instanceIndex].mask;

        // If the ray hit a light source, set the color to white, and stop immediately.
        if(mask & GEOMETRY_MASK_LIGHT){
            if(mask & GEOMETRY_MASK_SPHERE_LIGHT){
                Sphere area_light;
                area_light = *(const device Sphere*)intersection.primitive_data;
                float3 light_color = area_light.material.color;
                accumulatedColor = color * light_color;
                break;
            }
            if(mask & GEOMETRY_MASK_TRIANGLE_LIGHT){
                Triangle area_light;
                area_light = *(const device Triangle*)intersection.primitive_data;
                float3 light_color = area_light.material.color;
                accumulatedColor = color * light_color;
                break;
            }
        }


        // The ray hit something. Look up the transformation matrix for this instance.
        float4x4 objectToWorldSpaceTransform(1.0f);

        for (int column = 0; column < 4; column++)
            for (int row = 0; row < 3; row++)
                objectToWorldSpaceTransform[column][row] = instances[instanceIndex].transformationMatrix[column][row];

        // Compute the intersection point in world space.
        float3 worldSpaceIntersectionPoint = ray.origin + ray.direction * intersection.distance;

//        unsigned primitiveIndex = intersection.primitive_id;
//        unsigned int geometryIndex = instances[instanceIndex].accelerationStructureIndex;
        float2 barycentric_coords = intersection.triangle_barycentric_coord;

        float3 worldSpaceSurfaceNormal = 0.0f;
        float3 surfaceColor = 0.0f;
        Material material;

        if (mask & GEOMETRY_MASK_TRIANGLE) {
            Triangle triangle;
            
            float3 objectSpaceSurfaceNormal;
            
            triangle = *(const device Triangle*)intersection.primitive_data;
            
//#if SUPPORTS_METAL_3
//            if (usePerPrimitiveData) {
//                // Per-primitive data points to data from the specified buffer as was configured in the MTLAccelerationStructureTriangleGeometryDescriptor.
//                triangle = *(const device Triangle*)intersection.primitive_data;
//            } else
//#endif
//            {
//                // The ray hit a triangle. Look up the corresponding geometry's normal and UV buffers.
//                device TriangleResources & triangleResources = *(device TriangleResources *)((device char *)resources + resourcesStride * geometryIndex);
//
//                triangle.normals[0] =  triangleResources.vertexNormals[triangleResources.indices[primitiveIndex * 3 + 0]];
//                triangle.normals[1] =  triangleResources.vertexNormals[triangleResources.indices[primitiveIndex * 3 + 1]];
//                triangle.normals[2] =  triangleResources.vertexNormals[triangleResources.indices[primitiveIndex * 3 + 2]];
//
//                triangle.colors[0] =  triangleResources.vertexColors[triangleResources.indices[primitiveIndex * 3 + 0]];
//                triangle.colors[1] =  triangleResources.vertexColors[triangleResources.indices[primitiveIndex * 3 + 1]];
//                triangle.colors[2] =  triangleResources.vertexColors[triangleResources.indices[primitiveIndex * 3 + 2]];
//            }

            // Interpolate the vertex normal at the intersection point.
            objectSpaceSurfaceNormal = interpolateVertexAttribute(triangle.normals, barycentric_coords);

            // Interpolate the vertex color at the intersection point.
//            surfaceColor = interpolateVertexAttribute(triangle.colors, barycentric_coords);
            surfaceColor = triangle.material.color;
            material = triangle.material;

            // Transform the normal from object to world space.
            worldSpaceSurfaceNormal = normalize(transformDirection(objectSpaceSurfaceNormal, objectToWorldSpaceTransform));
        }
        else if (mask & GEOMETRY_MASK_SPHERE) {
            Sphere sphere;
            
            sphere = *(const device Sphere*)intersection.primitive_data;
            
//#if SUPPORTS_METAL_3
//            if (usePerPrimitiveData) {
//                // Per-primitive data points to data from the specified buffer as was configured in the MTLAccelerationStructureBoundingBoxGeometryDescriptor.
//                sphere = *(const device Sphere*)intersection.primitive_data;
//            } else
//#endif
//            {
//                // The ray hit a sphere. Look up the corresponding sphere buffer.
//                device SphereResources & sphereResources = *(device SphereResources *)((device char *)resources + resourcesStride * geometryIndex);
//                sphere = sphereResources.spheres[primitiveIndex];
//            }

            // Transform the sphere's origin from object space to world space.
            float3 worldSpaceOrigin = transformPoint(sphere.origin, objectToWorldSpaceTransform);

            // Compute the surface normal directly in world space.
            worldSpaceSurfaceNormal = normalize(worldSpaceIntersectionPoint - worldSpaceOrigin);

            // The sphere is a uniform color, so you don't need to interpolate the color across the surface.
            surfaceColor = sphere.color;
            surfaceColor = sphere.material.color;
            material = sphere.material;
        }

        // Choose a random point to sample on the light source.
        float2 r = float2(halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 1),
                          halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 2));
        
        ray.origin = worldSpaceIntersectionPoint;
        
        // Deal with scattering.
        if(material.is_metal){
            scatter_record = metallicScatter(worldSpaceSurfaceNormal, ray.direction);
            if(scatter_record.accept){
                ray.direction = scatter_record.out_direction;
                color *= material.color;
            }else{
                break;
            }
        } else if(material.is_glass){
            float random_variable = r.x;
            scatter_record = glossyScatter(worldSpaceSurfaceNormal, ray.direction, random_variable);
            if(scatter_record.accept){
                ray.direction = scatter_record.out_direction;
                color *= material.color;
            }else{
                break;
            }
        } else if(material.is_phong){
            float2 random_variable = r;
            scatter_record = phongScatter(worldSpaceSurfaceNormal, ray.direction, random_variable, material);
            if(scatter_record.accept){
                ray.direction = scatter_record.out_direction;
            }else{
                break;
            }
            color *= scatter_record.attenuation;
        } else if(material.is_contain_volume){
            
        } else{
            if(dot(worldSpaceSurfaceNormal, ray.direction) > 0){
//                worldSpaceSurfaceNormal = -worldSpaceSurfaceNormal;
            }
            float3 worldSpaceSampleDirection = sampleCosineWeightedHemisphere(r);
            worldSpaceSampleDirection = alignHemisphereWithNormal(worldSpaceSampleDirection, worldSpaceSurfaceNormal);
            ray.direction = worldSpaceSampleDirection;
            color *= material.color;
        }
        
        ray.direction = normalize(ray.direction);
    }

    // Average this frame's sample with all of the previous frames.
    if (uniforms.frameIndex > 0) {
        float3 prevColor = prevTex.read(tid).xyz;
        prevColor *= uniforms.frameIndex;

        accumulatedColor += prevColor;
        accumulatedColor /= (uniforms.frameIndex + 1);
    }

    dstTex.write(float4(accumulatedColor, 1.0f), tid);
}

#pragma mark - Rat Tracing Kernel - NEE

// Main ray tracing kernel.
kernel void raytracingKernelNEE(
     uint2                                                  tid                       [[thread_position_in_grid]],
     constant Uniforms &                                    uniforms                  [[buffer(0)]],
     texture2d<unsigned int>                                randomTex                 [[texture(0)]],
     texture2d<float>                                       prevTex                   [[texture(1)]],
     texture2d<float, access::write>                        dstTex                    [[texture(2)]],
     texture2d<float>                                       envmapTex                 [[texture(3)]],
     device void                                           *resources                 [[buffer(1)]],
     constant MTLAccelerationStructureInstanceDescriptor   *instances                 [[buffer(2)]],
     constant AreaLight                                    *areaLights                [[buffer(3)]],
     instance_acceleration_structure                        accelerationStructure     [[buffer(4)]],
     intersection_function_table<triangle_data, instancing> intersectionFunctionTable [[buffer(5)]],
     constant unsigned int                                 *lightIndexs               [[buffer(6)]],
     constant unsigned int                                 *lightCounts               [[buffer(7)]]
)
{
    // The sample aligns the thread count to the threadgroup size, which means the thread count
    // may be different than the bounds of the texture. Test to make sure this thread
    // is referencing a pixel within the bounds of the texture.
    if (tid.x > uniforms.width || tid.y > uniforms.height) return;
    
    // Apply a random offset to the random number index to decorrelate pixels.
    unsigned int offset = randomTex.read(tid).x;
    
    // Generate Ray.
    ray ray = generate_ray(tid, uniforms, randomTex);

    // Start with a fully white color. The kernel scales the light each time the
    // ray bounces off of a surface, based on how much of each light component
    // the surface absorbs.
    float3 color = float3(1.0f, 1.0f, 1.0f);

    float3 accumulatedColor = float3(0.0f, 0.0f, 0.0f);

    // Create an intersector to test for intersection between the ray and the geometry in the scene.
    intersector<triangle_data, instancing> i;

    typename intersector<triangle_data, instancing>::result_type intersection;
    
    ScatterRecord scatter_record;
    
    // Check if last is specular material
    bool no_light = uniforms.lightCount == 0;
    bool last_is_specular = false;

    // Simulate up to three ray bounces. Each bounce propagates light backward along the
    // ray's path toward the camera.
    for (int bounce = 0; bounce < max_bounce; bounce++) {
        // Get the closest intersection, not the first intersection. This is the default, but
        // the sample adjusts this property below when it casts shadow rays.
        i.accept_any_intersection(false);

        // Check for intersection between the ray and the acceleration structure. If the sample
        // isn't using intersection functions, it doesn't need to include one.
        intersection = i.intersect(ray, accelerationStructure, bounce == 0 || last_is_specular ? RAY_MASK_PRIMARY : RAY_MASK_SECONDARY, intersectionFunctionTable);

        // Stop if the ray didn't hit anything and has bounced out of the scene.
        if (intersection.type == intersection_type::none && (bounce == 0 || last_is_specular)) {
            accumulatedColor += get_env_color(ray.direction, envmapTex) * color;
            break;
        }
            
        unsigned int instanceIndex = intersection.instance_id;

        // Look up the mask for this instance, which indicates what type of geometry the ray hit.
        unsigned int mask = instances[instanceIndex].mask;

        // If the ray hit a light source, set the color to white, and stop immediately.
        if(mask & GEOMETRY_MASK_LIGHT){
            if(mask & GEOMETRY_MASK_SPHERE_LIGHT){
                Sphere area_light;
                area_light = *(const device Sphere*)intersection.primitive_data;
                float3 light_color = area_light.material.color;
                accumulatedColor += light_color;
                break;
            }
            if(mask & GEOMETRY_MASK_TRIANGLE_LIGHT){
                Triangle area_light;
                area_light = *(const device Triangle*)intersection.primitive_data;
                float3 light_color = area_light.material.color;
                accumulatedColor += light_color;
                break;
            }
        }

        // The ray hit something. Look up the transformation matrix for this instance.
        float4x4 objectToWorldSpaceTransform(1.0f);

        for (int column = 0; column < 4; column++)
            for (int row = 0; row < 3; row++)
                objectToWorldSpaceTransform[column][row] = instances[instanceIndex].transformationMatrix[column][row];

        // Compute the intersection point in world space.
        float3 worldSpaceIntersectionPoint = ray.origin + ray.direction * intersection.distance;

        float2 barycentric_coords = intersection.triangle_barycentric_coord;

        float3 worldSpaceSurfaceNormal = 0.0f;
        float3 surfaceColor = 0.0f;
        Material material;

        if (mask & GEOMETRY_MASK_TRIANGLE) {
            Triangle triangle;
            
            float3 objectSpaceSurfaceNormal;
            
            triangle = *(const device Triangle*)intersection.primitive_data;

            // Interpolate the vertex normal at the intersection point.
            objectSpaceSurfaceNormal = interpolateVertexAttribute(triangle.normals, barycentric_coords);

            surfaceColor = triangle.material.color;
            material = triangle.material;

            // Transform the normal from object to world space.
            worldSpaceSurfaceNormal = normalize(transformDirection(objectSpaceSurfaceNormal, objectToWorldSpaceTransform));
        }
        else if (mask & GEOMETRY_MASK_SPHERE) {
            Sphere sphere;
            
            sphere = *(const device Sphere*)intersection.primitive_data;

            // Transform the sphere's origin from object space to world space.
            float3 worldSpaceOrigin = transformPoint(sphere.origin, objectToWorldSpaceTransform);

            // Compute the surface normal directly in world space.
            worldSpaceSurfaceNormal = normalize(worldSpaceIntersectionPoint - worldSpaceOrigin);

            // The sphere is a uniform color, so you don't need to interpolate the color across the surface.
            surfaceColor = sphere.material.color;
            material = sphere.material;
        }
        
        // Choose a random direction to continue the path of the ray. This causes light to
        // bounce between surfaces. An app might evaluate a more complicated equation to
        // calculate the amount of light that reflects between intersection points.  However,
        // all the math in this kernel cancels out because this app assumes a simple diffuse
        // BRDF and samples the rays with a cosine distribution over the hemisphere (importance
        // sampling). This requires that the kernel only multiply the colors together. This
        // sampling strategy also reduces the amount of noise in the output image.
        float2 r = float2(halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 3),
                          halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 4));
        
        ray.origin = worldSpaceIntersectionPoint;
        float3 in_direction = ray.direction;
        
        // Deal with scattering.
        if(material.is_metal){
            scatter_record = metallicScatter(worldSpaceSurfaceNormal, ray.direction);
            if(scatter_record.accept){
                ray.direction = scatter_record.out_direction;
            }else{
                break;
            }
        } else if(material.is_glass){
            float random_variable = r.x;
            scatter_record = glossyScatter(worldSpaceSurfaceNormal, ray.direction, random_variable);
            
            if(scatter_record.accept){
                ray.direction = scatter_record.out_direction;
            }else{
                break;
            }
        } else if(material.is_phong){
            float2 random_variable = r;
            scatter_record = phongScatter(worldSpaceSurfaceNormal, ray.direction, random_variable, material);
            if(scatter_record.accept){
                ray.direction = scatter_record.out_direction;
            }else{
                break;
            }
        } else{
            if(dot(worldSpaceSurfaceNormal, ray.direction) > 0){
//                worldSpaceSurfaceNormal = -worldSpaceSurfaceNormal;
                break;
            }
            float3 worldSpaceSampleDirection = sampleCosineWeightedHemisphere(r);
            worldSpaceSampleDirection = alignHemisphereWithNormal(worldSpaceSampleDirection, worldSpaceSurfaceNormal);
            ray.direction = worldSpaceSampleDirection;
        }
        
        // Scale the ray color by the color of the surface to simulate the surface absorbing light.
        color *= surfaceColor;
        ray.direction = normalize(ray.direction);
        
        // Sample Environment map
        r = float2(halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 11),
                   halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 9));
        EmitterRecord envmap_record = uniform_sample_env(r, envmapTex);
        
        // Deal with Environment map
        if(material.is_metal){
            last_is_specular = true;
        }else if(material.is_glass){
            last_is_specular = true;
        }else{
            last_is_specular = false;
            // Sample light for diffuse material
            struct ray shadowRay;
            shadowRay.origin = worldSpaceIntersectionPoint + worldSpaceSurfaceNormal * 1e-3f;;
            shadowRay.direction = envmap_record.out_direction;
            shadowRay.max_distance = envmap_record.distance;
            i.accept_any_intersection(true);
            intersection = i.intersect(shadowRay, accelerationStructure, RAY_MASK_SHADOW, intersectionFunctionTable);
            if (intersection.type == intersection_type::none){
                float material_pdf;
                if(material.is_phong){
                    material_pdf = phongPdf(worldSpaceSurfaceNormal, in_direction, shadowRay.direction, material);
                }else{
                    material_pdf = saturate(dot(worldSpaceSurfaceNormal, shadowRay.direction)) / M_PI_F;
                }
                accumulatedColor += envmap_record.emit * color * material_pdf / envmap_record.pdf;
            }
        }
        
        // IF LIGHT_COUNT == 0, DO NOT NEED TO DO LIGHT SAMPLING
        if(no_light) continue;

        // Sample a random light instance
        EmitterRecord emitter_record;
        float light_instance_sample = halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 7);
        unsigned int light_instance_index = min((unsigned int)(light_instance_sample * uniforms.lightCount), uniforms.lightCount - 1);
        unsigned int light_instance = lightIndexs[light_instance_index];
        
        unsigned int geometries_index = instances[light_instance].accelerationStructureIndex;
        unsigned int light_mask = instances[light_instance_index].mask;
        
        float4x4 lightToWorldSpaceTransform(1.0f);

        for (int column = 0; column < 4; column++)
            for (int row = 0; row < 3; row++)
                lightToWorldSpaceTransform[column][row] = instances[instanceIndex].transformationMatrix[column][row];
        
        // Sample a random light geometry from one instance.
        float light_geometry_sample = halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 8);
        unsigned int light_geometry_index = min((unsigned int)(light_geometry_sample * lightCounts[light_instance_index]), lightCounts[light_instance_index] - 1);
        
        // Choose a random point to sample on the light source.
        r = float2(halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 1),
                          halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 2));
        
        if(light_mask & GEOMETRY_MASK_SPHERE_LIGHT){
            device SphereResources & sphereResources = *(device SphereResources *)((device char *)resources + resourcesStride * geometries_index);
            Sphere sphere = sphereResources.spheres[light_geometry_index];
            emitter_record = sampleSphereLight(sphere, worldSpaceIntersectionPoint, r, lightToWorldSpaceTransform);
        }
        if(light_mask & GEOMETRY_MASK_TRIANGLE_LIGHT){
            device TriangleResources & triangleResources = *(device TriangleResources *)((device char *)resources + resourcesStride * geometries_index);
            Triangle triangle = triangleResources.triangles[light_geometry_index];
            emitter_record = sampleTriangleLight(triangle, worldSpaceIntersectionPoint, r, lightToWorldSpaceTransform);
        }
        
        // Deal with light sampling
        if(material.is_metal){
            last_is_specular = true;
        }else if(material.is_glass){
            last_is_specular = true;
        }else{
            last_is_specular = false;
            // Sample light for diffuse material
            struct ray shadowRay;
            shadowRay.origin = worldSpaceIntersectionPoint + worldSpaceSurfaceNormal * 1e-3f;;
            shadowRay.direction = emitter_record.out_direction;
            shadowRay.max_distance = emitter_record.distance - 1e-3f;
            i.accept_any_intersection(true);
            intersection = i.intersect(shadowRay, accelerationStructure, RAY_MASK_SHADOW, intersectionFunctionTable);
            if (intersection.type == intersection_type::none){
                float material_pdf;
                if(material.is_phong){
                    material_pdf = phongPdf(worldSpaceSurfaceNormal, in_direction, shadowRay.direction, material);
                }else{
                    material_pdf = saturate(dot(worldSpaceSurfaceNormal, shadowRay.direction)) / M_PI_F;
                }
                accumulatedColor += emitter_record.emit * color * material_pdf / emitter_record.pdf * uniforms.totalLightCount;
            }
        }
    }

    // Average this frame's sample with all of the previous frames.
    if (uniforms.frameIndex > 0) {
        float3 prevColor = prevTex.read(tid).xyz;
        prevColor *= uniforms.frameIndex;

        accumulatedColor += prevColor;
        accumulatedColor /= (uniforms.frameIndex + 1);
    }

    dstTex.write(float4(accumulatedColor, 1.0f), tid);
}

#pragma mark - Rat Tracing Kernel - MIS

// Main ray tracing kernel.
kernel void raytracingKernelMIS(
     uint2                                                  tid                       [[thread_position_in_grid]],
     constant Uniforms &                                    uniforms                  [[buffer(0)]],
     texture2d<unsigned int>                                randomTex                 [[texture(0)]],
     texture2d<float>                                       prevTex                   [[texture(1)]],
     texture2d<float, access::write>                        dstTex                    [[texture(2)]],
     texture2d<float>                                       envmapTex                 [[texture(3)]],
     device void                                           *resources                 [[buffer(1)]],
     constant MTLAccelerationStructureInstanceDescriptor   *instances                 [[buffer(2)]],
     constant AreaLight                                    *areaLights                [[buffer(3)]],
     instance_acceleration_structure                        accelerationStructure     [[buffer(4)]],
     intersection_function_table<triangle_data, instancing> intersectionFunctionTable [[buffer(5)]],
     constant unsigned int                                 *lightIndexs               [[buffer(6)]],
     constant unsigned int                                 *lightCounts               [[buffer(7)]]
)
{
    
    // The sample aligns the thread count to the threadgroup size, which means the thread count
    // may be different than the bounds of the texture. Test to make sure this thread
    // is referencing a pixel within the bounds of the texture.
    if (tid.x > uniforms.width || tid.y > uniforms.height) return;
    
    // Apply a random offset to the random number index to decorrelate pixels.
    unsigned int offset = randomTex.read(tid).x;
    
    // Generate Ray.
    ray ray = generate_ray(tid, uniforms, randomTex);

    // Start with a fully white color. The kernel scales the light each time the
    // ray bounces off of a surface, based on how much of each light component
    // the surface absorbs.
    float3 color = float3(1.0f, 1.0f, 1.0f);

    float3 accumulatedColor = float3(0.0f, 0.0f, 0.0f);
    
    float Le_weight = 1.f;
    float nee_pdf = 1.f;
    float brdf_pdf = 1.f;
    float mis_power = 4.f;
    
    // Check if last is specular material
    bool last_is_specular = false;

    // Create an intersector to test for intersection between the ray and the geometry in the scene.
    intersector<triangle_data, instancing> i;

    typename intersector<triangle_data, instancing>::result_type intersection;
    
    ScatterRecord scatter_record;
    
    // IF LIGHT_COUNT == 0, DO NOT NEED TO DO LIGHT SAMPLING
    bool no_light = uniforms.lightCount == 0;

    // Simulate up to three ray bounces. Each bounce propagates light backward along the
    // ray's path toward the camera.
    for (int bounce = 0; bounce < max_bounce; bounce++) {
        i.accept_any_intersection(false);
        intersection = i.intersect(ray, accelerationStructure, RAY_MASK_PRIMARY, intersectionFunctionTable);

        // Stop if the ray didn't hit anything and has bounced out of the scene.
        if (intersection.type == intersection_type::none) {
            if(bounce > 0 && !last_is_specular){
                nee_pdf = 1 / (4 * M_PI_F);
                Le_weight = 1.f + pow(nee_pdf, mis_power) / pow(brdf_pdf, mis_power);
                Le_weight = 1.f / Le_weight;
            }else{
                Le_weight = 1.f;
            }
            accumulatedColor += Le_weight * get_env_color(ray.direction, envmapTex) * color;
            break;
        }
            
        unsigned int instanceIndex = intersection.instance_id;

        // Look up the mask for this instance, which indicates what type of geometry the ray hit.
        unsigned int mask = instances[instanceIndex].mask;

        // The ray hit something. Look up the transformation matrix for this instance.
        float4x4 objectToWorldSpaceTransform(1.0f);

        for (int column = 0; column < 4; column++)
            for (int row = 0; row < 3; row++)
                objectToWorldSpaceTransform[column][row] = instances[instanceIndex].transformationMatrix[column][row];
        
        // Compute the intersection point in world space.
        float3 worldSpaceIntersectionPoint = ray.origin + ray.direction * intersection.distance;
        
        // If the ray hit a light source, set the color to white, and stop immediately.
        if(mask & GEOMETRY_MASK_LIGHT){
            if(mask & GEOMETRY_MASK_SPHERE_LIGHT){
                Sphere area_light;
                area_light = *(const device Sphere*)intersection.primitive_data;
                float3 light_color = area_light.material.color;
                if(bounce > 0 && !last_is_specular){
                    nee_pdf = spherePdf(area_light, ray.origin, objectToWorldSpaceTransform);
//                    Le_weight = pow(brdf_pdf, mis_power) / (pow(nee_pdf, mis_power) + pow(brdf_pdf, mis_power));
                    Le_weight = 1.f + pow(nee_pdf, mis_power) / pow(brdf_pdf, mis_power);
                    Le_weight = 1.f / Le_weight;
                }else{
                    Le_weight = 1.f;
                }
                accumulatedColor += light_color * Le_weight * color;
            }
            if(mask & GEOMETRY_MASK_TRIANGLE_LIGHT){
                Triangle area_light;
                float2 barycentric_coords = intersection.triangle_barycentric_coord;
                area_light = *(const device Triangle*)intersection.primitive_data;
                float3 light_color = area_light.material.color;
                if(bounce > 0 && !last_is_specular){
                    nee_pdf = trianglePdf(area_light, ray.origin, objectToWorldSpaceTransform, barycentric_coords);
//                    Le_weight = pow(brdf_pdf, mis_power) / (pow(nee_pdf, mis_power) + pow(brdf_pdf, mis_power));
                    Le_weight = 1.f + pow(nee_pdf, mis_power) / pow(brdf_pdf, mis_power);
                    Le_weight = 1.f / Le_weight;
                }else{
                    Le_weight = 1.f;
                }
                accumulatedColor += light_color * Le_weight * color;
            }
            break;
        }

        float3 worldSpaceSurfaceNormal = 0.0f;
        float3 surfaceColor = 0.0f;
        Material material;

        if (mask & GEOMETRY_MASK_TRIANGLE) {
            Triangle triangle;
            
            float3 objectSpaceSurfaceNormal;
            
            triangle = *(const device Triangle*)intersection.primitive_data;
            float2 barycentric_coords = intersection.triangle_barycentric_coord;

            // Interpolate the vertex normal at the intersection point.
            objectSpaceSurfaceNormal = interpolateVertexAttribute(triangle.normals, barycentric_coords);

            surfaceColor = triangle.material.color;
            material = triangle.material;

            // Transform the normal from object to world space.
            worldSpaceSurfaceNormal = normalize(transformDirection(objectSpaceSurfaceNormal, objectToWorldSpaceTransform));
        }
        else if (mask & GEOMETRY_MASK_SPHERE) {
            Sphere sphere;
            
            sphere = *(const device Sphere*)intersection.primitive_data;

            // Transform the sphere's origin from object space to world space.
            float3 worldSpaceOrigin = transformPoint(sphere.origin, objectToWorldSpaceTransform);

            // Compute the surface normal directly in world space.
            worldSpaceSurfaceNormal = normalize(worldSpaceIntersectionPoint - worldSpaceOrigin);

            // The sphere is a uniform color, so you don't need to interpolate the color across the surface.
            surfaceColor = sphere.material.color;
            material = sphere.material;
        }
        
        // Choose a random direction to continue the path of the ray.
        float2 r = float2(halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 3),
                          halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 4));
        
        ray.origin = worldSpaceIntersectionPoint;
        float3 in_direction = ray.direction;
        
        // Deal with scattering.
        if(material.is_metal){
            scatter_record = metallicScatter(worldSpaceSurfaceNormal, ray.direction);
            if(scatter_record.accept){
                ray.direction = scatter_record.out_direction;
            }else{
                break;
            }
        } else if(material.is_glass){
            float random_variable = r.x;
            scatter_record = glossyScatter(worldSpaceSurfaceNormal, ray.direction, random_variable);
            if(scatter_record.accept){
                ray.direction = scatter_record.out_direction;
            }else{
                break;
            }
        } else if(material.is_phong){
            float2 random_variable = r;
            scatter_record = phongScatter(worldSpaceSurfaceNormal, ray.direction, random_variable, material);
            if(scatter_record.accept){
                ray.direction = scatter_record.out_direction;
            }else{
                break;
            }
        } else{
            if(dot(worldSpaceSurfaceNormal, ray.direction) > 0){
//                worldSpaceSurfaceNormal = -worldSpaceSurfaceNormal;
            }
            float3 worldSpaceSampleDirection = sampleCosineWeightedHemisphere(r);
            worldSpaceSampleDirection = alignHemisphereWithNormal(worldSpaceSampleDirection, worldSpaceSurfaceNormal);
            ray.direction = worldSpaceSampleDirection;
        }
        
        // Scale the ray color by the color of the surface to simulate the surface absorbing light.
        color *= surfaceColor;
        ray.direction = normalize(ray.direction);
        
        // Sample Environment map
        r = float2(halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 11),
                   halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 9));
        EmitterRecord envmap_record = uniform_sample_env(r, envmapTex);
        
        // Deal with Environment map
        if(material.is_metal){
            last_is_specular = true;
        }else if(material.is_glass){
            last_is_specular = true;
        }else{
            last_is_specular = false;
            // Sample light for diffuse material
            struct ray shadowRay;
            shadowRay.origin = worldSpaceIntersectionPoint + worldSpaceSurfaceNormal * 1e-3f;;
            shadowRay.direction = envmap_record.out_direction;
            shadowRay.max_distance = envmap_record.distance - 1e-3f;
            i.accept_any_intersection(true);
            intersection = i.intersect(shadowRay, accelerationStructure, RAY_MASK_SHADOW_ALL, intersectionFunctionTable);
            if (intersection.type == intersection_type::none){
                if(material.is_phong){
                    brdf_pdf = phongPdf(worldSpaceSurfaceNormal, in_direction, shadowRay.direction, material);
                }else{
                    brdf_pdf = saturate(dot(worldSpaceSurfaceNormal, shadowRay.direction)) / M_PI_F;
                }
                nee_pdf = envmap_record.pdf;
                //                Le_weight = pow(nee_pdf, mis_power) / (pow(nee_pdf, mis_power) + pow(brdf_pdf, mis_power));
                Le_weight = 1.f + pow(brdf_pdf, mis_power) / pow(nee_pdf, mis_power);
                Le_weight = 1.f / Le_weight;
                accumulatedColor += Le_weight * envmap_record.emit * color * brdf_pdf / nee_pdf;
            }
        }
        
        // IF LIGHT_COUNT == 0, DO NOT NEED TO DO LIGHT SAMPLING
        if(no_light) {
            // Calculate BRDF_PDF for next trace.
            if(material.is_metal){

            }else if(material.is_glass){

            }else if(material.is_phong){
                brdf_pdf = phongPdf(worldSpaceSurfaceNormal, in_direction, ray.direction, material);
            }else{
                brdf_pdf = max(0.f, dot(ray.direction, worldSpaceSurfaceNormal)) / M_PI_F;
            }
            continue;
        }

        // Sample a random light instance
        EmitterRecord emitter_record;
        float light_instance_sample = halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 7);
        unsigned int light_instance_index = min((unsigned int)(light_instance_sample * uniforms.lightCount), uniforms.lightCount - 1);
        unsigned int light_instance = lightIndexs[light_instance_index];
        
        unsigned int geometries_index = instances[light_instance].accelerationStructureIndex;
        unsigned int light_mask = instances[light_instance_index].mask;
        
        float4x4 lightToWorldSpaceTransform(1.0f);

        for (int column = 0; column < 4; column++)
            for (int row = 0; row < 3; row++)
                lightToWorldSpaceTransform[column][row] = instances[instanceIndex].transformationMatrix[column][row];
        
        // Sample a random light geometry from one instance.
        float light_geometry_sample = halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 8);
        unsigned int light_geometry_index = min((unsigned int)(light_geometry_sample * lightCounts[light_instance_index]), lightCounts[light_instance_index] - 1);
        
        // Choose a random point to sample on the light source.
        r = float2(halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 1),
                   halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 2));
        
        if(light_mask & GEOMETRY_MASK_SPHERE_LIGHT){
            device SphereResources & sphereResources = *(device SphereResources *)((device char *)resources + resourcesStride * geometries_index);
            Sphere sphere = sphereResources.spheres[light_geometry_index];
            emitter_record = sampleSphereLight(sphere, worldSpaceIntersectionPoint, r, lightToWorldSpaceTransform);
        }
        if(light_mask & GEOMETRY_MASK_TRIANGLE_LIGHT){
            device TriangleResources & triangleResources = *(device TriangleResources *)((device char *)resources + resourcesStride * geometries_index);
            Triangle triangle = triangleResources.triangles[light_geometry_index];
            emitter_record = sampleTriangleLight(triangle, worldSpaceIntersectionPoint, r, lightToWorldSpaceTransform);
        }
        
        // Deal with light sampling
        if(material.is_metal){
            last_is_specular = true;
        }else if(material.is_glass){
            last_is_specular = true;
        }else{
            last_is_specular = false;
            // Sample light for diffuse material
            struct ray shadowRay;
            shadowRay.origin = worldSpaceIntersectionPoint + worldSpaceSurfaceNormal * 1e-3f;;
            shadowRay.direction = emitter_record.out_direction;
            shadowRay.max_distance = emitter_record.distance - 1e-3f;
            i.accept_any_intersection(true);
            intersection = i.intersect(shadowRay, accelerationStructure, RAY_MASK_SHADOW, intersectionFunctionTable);
            if (intersection.type == intersection_type::none){
                if(material.is_phong){
                    brdf_pdf = phongPdf(worldSpaceSurfaceNormal, in_direction, shadowRay.direction, material);
                }else{
                    brdf_pdf = saturate(dot(worldSpaceSurfaceNormal, shadowRay.direction)) / M_PI_F;
                }
                nee_pdf = emitter_record.pdf;
//                Le_weight = pow(nee_pdf, mis_power) / (pow(nee_pdf, mis_power) + pow(brdf_pdf, mis_power));
                Le_weight = 1.f + pow(brdf_pdf, mis_power) / pow(nee_pdf, mis_power);
                Le_weight = 1.f / Le_weight;
                accumulatedColor += Le_weight * emitter_record.emit * color * brdf_pdf / nee_pdf * uniforms.totalLightCount;
            }
        }
                
        // Calculate BRDF_PDF for next trace.
        if(material.is_metal){

        }else if(material.is_glass){

        }else if(material.is_phong){
            brdf_pdf = phongPdf(worldSpaceSurfaceNormal, in_direction, ray.direction, material);
        }else{
            brdf_pdf = max(0.f, dot(ray.direction, worldSpaceSurfaceNormal)) / M_PI_F;
        }
    }

    // Average this frame's sample with all of the previous frames.
    if (uniforms.frameIndex > 0) {
        float3 prevColor = prevTex.read(tid).xyz;
        prevColor *= uniforms.frameIndex;

        accumulatedColor += prevColor;
        accumulatedColor /= (uniforms.frameIndex + 1);
    }

    dstTex.write(float4(accumulatedColor, 1.0f), tid);
}

#pragma mark - Ray Tracing Kernel - Vol

// Main ray tracing kernel.
kernel void raytracingKernelVOL(
     uint2                                                  tid                       [[thread_position_in_grid]],
     constant Uniforms &                                    uniforms                  [[buffer(0)]],
     texture2d<unsigned int>                                randomTex                 [[texture(0)]],
     texture2d<float>                                       prevTex                   [[texture(1)]],
     texture2d<float, access::write>                        dstTex                    [[texture(2)]],
     texture2d<float>                                       envmapTex                 [[texture(3)]],
     texture3d<half>                                        vdbvolTex                 [[texture(4)]],
//     device void                                           *resources                 [[buffer(1), function_constant(useResourcesBuffer)]],
     constant MTLAccelerationStructureInstanceDescriptor   *instances                 [[buffer(2)]],
     constant AreaLight                                    *areaLights                [[buffer(3)]],
     instance_acceleration_structure                        accelerationStructure     [[buffer(4)]],
     intersection_function_table<triangle_data, instancing> intersectionFunctionTable [[buffer(5)]],
     constant float                                        *maxDensity                [[buffer(8)]]
)
{
    // The sample aligns the thread count to the threadgroup size, which means the thread count
    // may be different than the bounds of the texture. Test to make sure this thread
    // is referencing a pixel within the bounds of the texture.
    if (tid.x > uniforms.width || tid.y > uniforms.height) return;
    
    // Apply a random offset to the random number index to decorrelate pixels.
    unsigned int offset = randomTex.read(tid).x;
    
    // Generate Ray.
    ray ray = generate_ray(tid, uniforms, randomTex);

    // Start with a fully white color. The kernel scales the light each time the
    // ray bounces off of a surface, based on how much of each light component
    // the surface absorbs.
    float3 color = float3(1.0f, 1.0f, 1.0f);

    float3 accumulatedColor = float3(0.0f, 0.0f, 0.0f);

    // Create an intersector to test for intersection between the ray and the geometry in the scene.
    intersector<triangle_data, instancing> i;

    // If the sample isn't using intersection functions, provide some hints to Metal for
    // better performance.
//    if (!useIntersectionFunctions) {
//        i.assume_geometry_type(geometry_type::triangle);
//        i.force_opacity(forced_opacity::opaque);
//    }

    typename intersector<triangle_data, instancing>::result_type intersection;
    
    ScatterRecord scatter_record;

    // Simulate up to three ray bounces. Each bounce propagates light backward along the
    // ray's path toward the camera.
    for (int bounce = 0; bounce < max_bounce + 1; bounce++) {
        // Get the closest intersection, not the first intersection. This is the default, but
        // the sample adjusts this property below when it casts shadow rays.
        i.accept_any_intersection(false);

        // Check for intersection between the ray and the acceleration structure. If the sample
        // isn't using intersection functions, it doesn't need to include one.
        intersection = i.intersect(ray, accelerationStructure, RAY_MASK_PRIMARY, intersectionFunctionTable);

        // Stop if the ray didn't hit anything and has bounced out of the scene.
        if (intersection.type == intersection_type::none) {
            accumulatedColor = get_env_color(ray.direction, envmapTex) * color;
            break;
        }
            
        unsigned int instanceIndex = intersection.instance_id;

        // Look up the mask for this instance, which indicates what type of geometry the ray hit.
        unsigned int mask = instances[instanceIndex].mask;

        // If the ray hit a light source, set the color to white, and stop immediately.
        if(mask & GEOMETRY_MASK_LIGHT){
            if(mask & GEOMETRY_MASK_SPHERE_LIGHT){
                Sphere area_light;
                area_light = *(const device Sphere*)intersection.primitive_data;
                float3 light_color = area_light.material.color;
                accumulatedColor = color * light_color;
                break;
            }
            if(mask & GEOMETRY_MASK_TRIANGLE_LIGHT){
                Triangle area_light;
                area_light = *(const device Triangle*)intersection.primitive_data;
                float3 light_color = area_light.material.color;
                accumulatedColor = color * light_color;
                break;
            }
        }


        // The ray hit something. Look up the transformation matrix for this instance.
        float4x4 objectToWorldSpaceTransform(1.0f);

        for (int column = 0; column < 4; column++)
            for (int row = 0; row < 3; row++)
                objectToWorldSpaceTransform[column][row] = instances[instanceIndex].transformationMatrix[column][row];

        // Compute the intersection point in world space.
        float3 worldSpaceIntersectionPoint = ray.origin + ray.direction * intersection.distance;

        float2 barycentric_coords = intersection.triangle_barycentric_coord;

        float3 worldSpaceSurfaceNormal = 0.0f;
        float3 surfaceColor = 0.0f;
        Material material;

        if (mask & GEOMETRY_MASK_TRIANGLE) {
            Triangle triangle;
            
            float3 objectSpaceSurfaceNormal;
            
            triangle = *(const device Triangle*)intersection.primitive_data;

            // Interpolate the vertex normal at the intersection point.
            objectSpaceSurfaceNormal = interpolateVertexAttribute(triangle.normals, barycentric_coords);

            // Interpolate the vertex color at the intersection point.
//            surfaceColor = interpolateVertexAttribute(triangle.colors, barycentric_coords);
            surfaceColor = triangle.material.color;
            material = triangle.material;

            // Transform the normal from object to world space.
            worldSpaceSurfaceNormal = normalize(transformDirection(objectSpaceSurfaceNormal, objectToWorldSpaceTransform));
        }
        else if (mask & GEOMETRY_MASK_SPHERE) {
            Sphere sphere;
            
            sphere = *(const device Sphere*)intersection.primitive_data;

            // Transform the sphere's origin from object space to world space.
            float3 worldSpaceOrigin = transformPoint(sphere.origin, objectToWorldSpaceTransform);

            // Compute the surface normal directly in world space.
            worldSpaceSurfaceNormal = normalize(worldSpaceIntersectionPoint - worldSpaceOrigin);

            // The sphere is a uniform color, so you don't need to interpolate the color across the surface.
            surfaceColor = sphere.color;
            surfaceColor = sphere.material.color;
            material = sphere.material;
        }

        // Choose a random point to sample on the light source.
        float2 r = float2(halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 1),
                          halton(offset + uniforms.frameIndex, 2 + bounce * 5 + 2));
        
        ray.origin = worldSpaceIntersectionPoint;
        bool in = false;
        
        // Deal with scattering.
        if(material.is_metal){
            scatter_record = metallicScatter(worldSpaceSurfaceNormal, ray.direction);
            if(scatter_record.accept){
                ray.direction = scatter_record.out_direction;
                color *= material.color;
            }else{
                break;
            }
        } else if(material.is_glass){
            float random_variable = r.x;
            scatter_record = glossyScatter(worldSpaceSurfaceNormal, ray.direction, random_variable);
            if(scatter_record.accept){
                if(scatter_record.is_refract) in=true;
                ray.direction = scatter_record.out_direction;
                color *= material.color;
            }else{
                break;
            }
        } else if(material.is_phong){
            float2 random_variable = r;
            scatter_record = phongScatter(worldSpaceSurfaceNormal, ray.direction, random_variable, material);
            if(scatter_record.accept){
                ray.direction = scatter_record.out_direction;
            }else{
                break;
            }
            color *= scatter_record.attenuation;
        } else if(material.is_contain_volume){
            in = true;
        } else{
            if(dot(worldSpaceSurfaceNormal, ray.direction) > 0){
//                worldSpaceSurfaceNormal = -worldSpaceSurfaceNormal;
            }
            float3 worldSpaceSampleDirection = sampleCosineWeightedHemisphere(r);
            worldSpaceSampleDirection = alignHemisphereWithNormal(worldSpaceSampleDirection, worldSpaceSurfaceNormal);
            ray.direction = worldSpaceSampleDirection;
            color *= material.color;
        }
        
        ray.direction = normalize(ray.direction);
        
//        // For Homogeneous Medium (with or without glass)
        if(in && material.is_contain_volume && *maxDensity == 0.f){
            bool out = false;
//            ray.origin = ray.origin + 1e-3f;
            for (int volume_bounce = 0; volume_bounce < max_bounce + 1; volume_bounce++){
                r = float2(halton(offset + uniforms.frameIndex, 2 + bounce * 5 + volume_bounce * 3 + 1),
                                                      halton(offset + uniforms.frameIndex, 2 + bounce * 5 + volume_bounce * 3 + 2));
        
                float t = -log(1.f - halton(offset + uniforms.frameIndex, 2 + bounce * 5 + volume_bounce * 3 + 3)) / material.density;
                float3 sigma_t = float3(material.density, material.density, material.density);
                float3 sigma_s = material.albedo * material.density;
                float3 sigma_a = sigma_t - sigma_s;
                ray.max_distance = t;
                i.accept_any_intersection(true);
                intersection = i.intersect(ray, accelerationStructure, GEOMETRY_MASK_VOLUME_CONTAINER, intersectionFunctionTable);
                if(intersection.type == intersection_type::none){

                    ray.origin = ray.origin + ray.direction * t;
                    accumulatedColor += sigma_a/sigma_t * material.emission * color;
                    color *= sigma_s/sigma_t;
                    ray.direction = sampleUniformSphere(r); // CHANGE PHASE FUNCTION HERE
                }else{
                    ray.origin = ray.origin + ray.direction * (intersection.distance);
                    if(material.is_glass){
                        float random_variable = r.x;
                        if (mask & GEOMETRY_MASK_TRIANGLE) {
                            Triangle triangle;
                            float3 objectSpaceSurfaceNormal;
                            triangle = *(const device Triangle*)intersection.primitive_data;
                            objectSpaceSurfaceNormal = interpolateVertexAttribute(triangle.normals, barycentric_coords);
                            worldSpaceSurfaceNormal = normalize(transformDirection(objectSpaceSurfaceNormal, objectToWorldSpaceTransform));
                        }
                        else if (mask & GEOMETRY_MASK_SPHERE) {
                            Sphere sphere;
                            sphere = *(const device Sphere*)intersection.primitive_data;
                            float3 worldSpaceOrigin = transformPoint(sphere.origin, objectToWorldSpaceTransform);
                            worldSpaceSurfaceNormal = normalize(worldSpaceIntersectionPoint - worldSpaceOrigin);
                        }
                        scatter_record = glossyScatter(worldSpaceSurfaceNormal, ray.direction, random_variable);
                        out = scatter_record.is_refract;
                        ray.direction = scatter_record.out_direction;
                        color *= material.color;
                        if(out) {
                            ray.max_distance = INFINITY;
                            break;
                        }
                    }else{
                        //                    color *= exp(-material.density * intersection.distance);
                        ray.max_distance = INFINITY;
                        out = true;
                        break;
                    }
                }
            }
            if(!out) break;
        }
        
        if(in && material.is_contain_volume && *maxDensity != 0.f){
            bool out = false;
            float4x4 worldToObjectSpaceTransform = inverse(objectToWorldSpaceTransform);
            for (int volume_bounce = 0; volume_bounce < max_bounce*32 + 1; volume_bounce++){
                r = float2(halton(offset + uniforms.frameIndex, 2 + bounce * 5 + volume_bounce * 3 + 1),
                           halton(offset + uniforms.frameIndex, 2 + bounce * 5 + volume_bounce * 3 + 2));
                
                float t = -log(1.f - halton(offset + uniforms.frameIndex, 2 + bounce * 5 + volume_bounce * 3 + 3)) / *maxDensity;
                float3 localSpaceIntersectionPoint = transformPoint(ray.origin, worldToObjectSpaceTransform);
                half4 vdb_data = get_vdb_data(localSpaceIntersectionPoint, vdbvolTex);
                float density = (float)vdb_data.a;
                float3 albedo = (float3)vdb_data.rgb;
                
                float3 sigma_t = float3(density);
                float3 sigma_s = albedo * sigma_t;
//                float3 sigma_a = sigma_t - sigma_s;
                float prob_null = (*maxDensity - (float)vdb_data.a) / *maxDensity;
                float prob_s = sigma_s.x / *maxDensity;
//                float prob_a = sigma_a.x / *maxDensity;
                
                ray.max_distance = t;
                i.accept_any_intersection(false);
                intersection = i.intersect(ray, accelerationStructure, GEOMETRY_MASK_VOLUME_CONTAINER, intersectionFunctionTable);
                
                if(intersection.type == intersection_type::none){
                    float random_variable = halton(offset + uniforms.frameIndex, 2 + bounce * 5 + volume_bounce * 3 + 4);
                    ray.origin = ray.origin + ray.direction * t;
                    if(random_variable <= prob_null){

                    }else if( prob_null < random_variable <= prob_null + prob_s){
                        color *= albedo;
                        ray.direction = sampleUniformSphere(r);
                    }else {
//                        accumulatedColor += sigma_a/sigma_t * material.emission * color;
//                        color *= sigma_a/sigma_t * exp(-sigma_a * t);
                        color *= albedo;
                    }
                }else{
                    ray.origin = ray.origin + ray.direction * (intersection.distance);
                    if(material.is_glass){
                        float random_variable = r.x;
                        if (mask & GEOMETRY_MASK_TRIANGLE) {
                            Triangle triangle;
                            float3 objectSpaceSurfaceNormal;
                            triangle = *(const device Triangle*)intersection.primitive_data;
                            objectSpaceSurfaceNormal = interpolateVertexAttribute(triangle.normals, barycentric_coords);
                            worldSpaceSurfaceNormal = normalize(transformDirection(objectSpaceSurfaceNormal, objectToWorldSpaceTransform));
                        }
                        else if (mask & GEOMETRY_MASK_SPHERE) {
                            Sphere sphere;
                            sphere = *(const device Sphere*)intersection.primitive_data;
                            float3 worldSpaceOrigin = transformPoint(sphere.origin, objectToWorldSpaceTransform);
                            worldSpaceSurfaceNormal = normalize(worldSpaceIntersectionPoint - worldSpaceOrigin);
                        }
                        scatter_record = glossyScatter(worldSpaceSurfaceNormal, ray.direction, random_variable);
                        out = scatter_record.is_refract;
                        ray.direction = scatter_record.out_direction;
                        color *= material.color;
                        if(out) {
                            ray.max_distance = INFINITY;
                            break;
                        }
                    }else{
                        //                    color *= exp(-material.density * intersection.distance);
                        ray.max_distance = INFINITY;
                        out = true;
                        break;
                    }
                }
            }
            if(!out) break;
        }
    }

    // Average this frame's sample with all of the previous frames.
    if (uniforms.frameIndex > 0) {
        float3 prevColor = prevTex.read(tid).xyz;
        prevColor *= uniforms.frameIndex;

        accumulatedColor += prevColor;
        accumulatedColor /= (uniforms.frameIndex + 1);
    }

    dstTex.write(float4(accumulatedColor, 1.0f), tid);
}


// Screen filling quad in normalized device coordinates.
constant float2 quadVertices[] = {
    float2(-1, -1),
    float2(-1,  1),
    float2( 1,  1),
    float2(-1, -1),
    float2( 1,  1),
    float2( 1, -1)
};

struct CopyVertexOut {
    float4 position [[position]];
    float2 uv;
};

// Simple vertex shader that passes through NDC quad positions.
vertex CopyVertexOut copyVertex(unsigned short vid [[vertex_id]]) {
    float2 position = quadVertices[vid];

    CopyVertexOut out;

    out.position = float4(position, 0, 1);
    out.uv = position * 0.5f + 0.5f;

    return out;
}

// Simple fragment shader that copies a texture and applies a simple tonemapping function.
fragment float4 copyFragment(CopyVertexOut in [[stage_in]],
                             texture2d<float> tex)
{
    constexpr sampler sam(min_filter::nearest, mag_filter::nearest, mip_filter::none);

    float3 color = tex.sample(sam, in.uv).xyz;

    // Apply a simple tonemapping function to reduce the dynamic range of the
    // input image into a range which the screen can display.
    color = color / (1.0f + color);

    return float4(color, 1.0f);
}
