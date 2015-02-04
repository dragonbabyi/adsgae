/**
 * Real-time Realistic Ocean Lighting using Seamless Transitions from Geometry to BRDF
 * Copyright (c) 2009 INRIA
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holders nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/**
 * Authors: Eric Bruneton & Jonathan Dupuy
 */

////////////////////////////////////////////////////////////////////
precision highp float;

layout(std140) uniform Matrices {
    mat4 view;
    mat4 proj;
    vec4 sunpos;
};

uniform sampler2DArray fftWavesSampler;	// ocean surface

uniform vec4 GRID_SIZES;
uniform float hdrExposure;
uniform vec3 seaColor; // sea bottom color
//uniform float normals;
uniform float choppy;
uniform vec4 choppy_factor;
uniform float jacobian_scale;
uniform vec2 gridSize;

////////////////////////////////////////////////////////////////////
#ifdef _VERTEX_

#define LAYER_HEIGHT		0.0
const float M_PI = 3.141592657;

in vec4 position;   //vertex position
 
out vec2 u; 	// horizontal coordinates in world space used to compute P(u)
out vec3 P; 	// wave point P(u) in world space

mat4 projMatrix = transpose(proj);
mat4 modelView = transpose(view);

mat4 screenToCamera = inverse(projMatrix);
mat4 cameraToWorld = inverse(modelView);
vec3 worldCamera = vec3(cameraToWorld[3][0], cameraToWorld[3][1], cameraToWorld[3][2]); 

vec2 oceanPos(vec4 vertex) { 
    vec3 cameraDir = normalize((screenToCamera * vertex).xyz);
    vec3 worldDir = (cameraToWorld * vec4(cameraDir, 0.0)).xyz;
    float t = -worldCamera.z / worldDir.z;
    return worldCamera.xy + t * worldDir.xy;
}

void main() {
 
    u = oceanPos(position);
    vec2 ux = oceanPos(position + vec4(gridSize.x, 0.0, 0.0, 0.0));
    vec2 uy = oceanPos(position + vec4(0.0, gridSize.y, 0.0, 0.0));
    vec2 dux = ux - u;
    vec2 duy = uy - u;

    // sum altitudes (use grad to get correct mipmap level)
    vec3 dP = vec3(0.0);
    dP.z += textureGrad(fftWavesSampler, vec3(u / GRID_SIZES.x, LAYER_HEIGHT), dux / GRID_SIZES.x, duy / GRID_SIZES.x).x;
    dP.z += textureGrad(fftWavesSampler, vec3(u / GRID_SIZES.y, LAYER_HEIGHT), dux / GRID_SIZES.y, duy / GRID_SIZES.y).y;
    dP.z += textureGrad(fftWavesSampler, vec3(u / GRID_SIZES.z, LAYER_HEIGHT), dux / GRID_SIZES.z, duy / GRID_SIZES.z).z;
    dP.z += textureGrad(fftWavesSampler, vec3(u / GRID_SIZES.w, LAYER_HEIGHT), dux / GRID_SIZES.w, duy / GRID_SIZES.w).w;

    // choppy
    if (choppy > 0.0) {
        dP.xy += choppy_factor.x*textureGrad(fftWavesSampler, vec3(u / GRID_SIZES.x, 3.0), dux / GRID_SIZES.x, duy / GRID_SIZES.x).xy;
        dP.xy += choppy_factor.y*textureGrad(fftWavesSampler, vec3(u / GRID_SIZES.y, 3.0), dux / GRID_SIZES.y, duy / GRID_SIZES.y).zw;
        dP.xy += choppy_factor.z*textureGrad(fftWavesSampler, vec3(u / GRID_SIZES.z, 4.0), dux / GRID_SIZES.z, duy / GRID_SIZES.z).xy;
        dP.xy += choppy_factor.w*textureGrad(fftWavesSampler, vec3(u / GRID_SIZES.w, 4.0), dux / GRID_SIZES.w, duy / GRID_SIZES.w).zw;
    }
    P = vec3(u + dP.xy, dP.z);   // varying

    // float d = (projMatrix * modelView * vec4(P, 1.0)).x;   // will be negative in the opposite direction
    // float angle = -R(d);  //R = 0.4773;  //0.4776959488
    float angle = -0.4773;
    // rotate -R along axis y
    float ca = cos(angle * M_PI / 180.0);
    float sa = sin(angle * M_PI / 180.0);
    vec4 pos;
    pos.x = ca * P.x + sa * P.z;
    pos.y = P.y;
    pos.z = -sa * P.x + ca * P.z;
    pos.w = 1.0;           
	// Final position
    gl_Position = projMatrix * modelView * pos;
    // gl_Position = projMatrix * modelView * vec4(P, 1.0);
	//gl_Position = projMatrix * modelView * vec4(u, dP.z, 1.0);

}

#endif

////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////


#ifdef _FRAGMENT_


#define SEA_CONTRIB 
#define SUN_CONTRIB 
#define SKY_CONTRIB 
#define HARDWARE_ANISTROPIC_FILTERING
#define FOAM_CONTRIB

#define LAYER_JACOBIAN_XX 	5.0
#define LAYER_JACOBIAN_YY	6.0
#define LAYER_JACOBIAN_XY	7.0


#define TRANSMITTANCE_NON_LINEAR

const float SCALE = 1000.0;

const float Rg = 6360.0 * SCALE;
const float Rt = 6420.0 * SCALE;
const float RL = 6421.0 * SCALE;

// Rayleigh
const float HR = 8.0 * SCALE;
const vec3 betaR = vec3(5.8e-3, 1.35e-2, 3.31e-2) / SCALE;

// Mie
// DEFAULT
const float HM = 1.2 * SCALE;
const vec3 betaMSca = vec3(4e-3) / SCALE;
const vec3 betaMEx = betaMSca / 0.9;
const float mieG = 0.8;
// CLEAR SKY
/*const float HM = 1.2 * SCALE;
const vec3 betaMSca = vec3(20e-3) / SCALE;
const vec3 betaMEx = betaMSca / 0.9;
const float mieG = 0.76;*/
// PARTLY CLOUDY
/*const float HM = 3.0 * SCALE;
const vec3 betaMSca = vec3(3e-3) / SCALE;
const vec3 betaMEx = betaMSca / 0.9;
const float mieG = 0.65;*/
 
uniform sampler3D slopeVarianceSampler; 
uniform sampler2D skySampler;
 
in vec2 u; 	// horizontal coordinates in world space used to compute P(u)
in vec3 P; 	// wave point P(u) in world space

out vec4 FragColor;

const float M_PI = 3.141592657;
const vec3 earthPos = vec3(0.0, 0.0, 6360010.0); 
const float SUN_INTENSITY = 100.0;
  
mat4 projMatrix = transpose(proj);
mat4 modelView = transpose(view);

mat4 cameraToWorld = inverse(modelView);
vec3 worldCamera = vec3(cameraToWorld[3][0], cameraToWorld[3][1], cameraToWorld[3][2]); 
// vec3 worldSunDir = vec3(sunpos.x, sunpos.y, sunpos.z);

vec3 worldSunDir = normalize(vec3(sunpos.x, sunpos.y, sunpos.z));

////////
// optical depth for ray (r,mu) of length d, using analytic formula
// (mu=cos(view zenith angle)), intersections with ground ignored
// H=height scale of exponential density function
float opticalDepth(float H, float r, float mu, float d) {
    float a = sqrt((0.5/H)*r);
    vec2 a01 = a*vec2(mu, mu + d / r);
    vec2 a01s = sign(a01);
    vec2 a01sq = a01*a01;
    float x = a01s.y > a01s.x ? exp(a01sq.x) : 0.0;
    vec2 y = a01s / (2.3193*abs(a01) + sqrt(1.52*a01sq + 4.0)) * vec2(1.0, exp(-d/H*(d/(2.0*r)+mu)));
    return sqrt((6.2831*H)*r) * exp((Rg-r)/H) * (x + dot(y, vec2(1.0, -1.0)));
}

// transmittance(=transparency) of atmosphere for ray (r,mu) of length d
// (mu=cos(view zenith angle)), intersections with ground ignored
// uses analytic formula instead of transmittance texture
vec3 analyticTransmittance(float r, float mu, float d) {
    return exp(- betaR * opticalDepth(HR, r, mu, d) - betaMEx * opticalDepth(HM, r, mu, d));
}

vec3 transmittanceWithShadow(float r, float mu) {
    return mu < -sqrt(1.0 - (Rg / r) * (Rg / r)) ? vec3(0.0) : analyticTransmittance(r, mu, Rg);
}

vec3 sunRadiance(float r, float muS) {
     return transmittanceWithShadow(r, muS) * SUN_INTENSITY;
//    return transmittanceWithShadow(r, muS) * SUN_INTENSITY * cos(sunpos.w);
}

void sunRadianceAndSkyIrradiance(vec3 worldP, vec3 worldS, out vec3 sunL)
{
    vec3 worldV = normalize(worldP); // vertical vector
    float r = length(worldP);
    float muS = dot(worldV, worldS);
	sunL = sunRadiance(r, muS);
}

vec3 hdr(vec3 L) {
    L = L * hdrExposure;
    L.r = L.r < 1.413 ? pow(L.r * 0.38317, 1.0 / 2.2) : 1.0 - exp(-L.r);
    L.g = L.g < 1.413 ? pow(L.g * 0.38317, 1.0 / 2.2) : 1.0 - exp(-L.g);
    L.b = L.b < 1.413 ? pow(L.b * 0.38317, 1.0 / 2.2) : 1.0 - exp(-L.b);
    return L;
}

// ---------------------------------------------------------------------------
// REFLECTED SUN RADIANCE
// ---------------------------------------------------------------------------

// assumes x>0
float erfc(float x) {
	return 2.0 * exp(-x * x) / (2.319 * x + sqrt(4.0 + 1.52 * x * x));
}

float erf(float x) {
	float a  = 0.140012;
	float x2 = x*x;
	float ax2 = a*x2;
	return sign(x) * sqrt( 1.0 - exp(-x2*(4.0/M_PI + ax2)/(1.0 + ax2)) );
}

float Lambda(float cosTheta, float sigmaSq) {
	float v = cosTheta / sqrt((1.0 - cosTheta * cosTheta) * (2.0 * sigmaSq));
    return max(0.0, (exp(-v * v) - v * sqrt(M_PI) * erfc(v)) / (2.0 * v * sqrt(M_PI)));
	//return (exp(-v * v)) / (2.0 * v * sqrt(M_PI)); // approximate, faster formula
}

// L, V, N, Tx, Ty in world space
float reflectedSunRadiance(vec3 L, vec3 V, vec3 N, vec3 Tx, vec3 Ty, vec2 sigmaSq) {
    vec3 H = normalize(L + V);
    float zetax = dot(H, Tx) / dot(H, N);
    float zetay = dot(H, Ty) / dot(H, N);

    float zL = dot(L, N); // cos of source zenith angle
    float zV = dot(V, N); // cos of receiver zenith angle
    float zH = dot(H, N); // cos of facet normal zenith angle
    float zH2 = zH * zH;

    float p = exp(-0.5 * (zetax * zetax / sigmaSq.x + zetay * zetay / sigmaSq.y))
                / (2.0 * M_PI * sqrt(sigmaSq.x * sigmaSq.y));

    float tanV = atan(dot(V, Ty), dot(V, Tx));
    float cosV2 = 1.0 / (1.0 + tanV * tanV);
    float sigmaV2 = sigmaSq.x * cosV2 + sigmaSq.y * (1.0 - cosV2);

    float tanL = atan(dot(L, Ty), dot(L, Tx));
    float cosL2 = 1.0 / (1.0 + tanL * tanL);
    float sigmaL2 = sigmaSq.x * cosL2 + sigmaSq.y * (1.0 - cosL2);

    float fresnel = 0.02 + 0.98 * pow(1.0 - dot(V, H), 5.0);

    zL = max(zL, 0.01);
    zV = max(zV, 0.01);

    return fresnel * p  / ((1.0 + Lambda(zL, sigmaL2) + Lambda(zV, sigmaV2)) * zV * zH2 * zH2 * 4.0);
}

// ---------------------------------------------------------------------------
// REFLECTED SKY RADIANCE
// ---------------------------------------------------------------------------

// V, N, Tx, Ty in world space
vec2 U(vec2 zeta, vec3 V, vec3 N, vec3 Tx, vec3 Ty) {
    vec3 f = normalize(vec3(-zeta, 1.0)); // tangent space
    vec3 F = f.x * Tx + f.y * Ty + f.z * N; // world space
    vec3 R = 2.0 * dot(F, V) * F - V;
    return R.xy / (1.0 + R.z);
}

float meanFresnel(float cosThetaV, float sigmaV) {
	return pow(1.0 - cosThetaV, 5.0 * exp(-2.69 * sigmaV)) / (1.0 + 22.7 * pow(sigmaV, 1.5));
}

// V, N in world space
float meanFresnel(vec3 V, vec3 N, vec2 sigmaSq) {
    vec2 v = V.xy; // view direction in wind space
    vec2 t = v * v / (1.0 - V.z * V.z); // cos^2 and sin^2 of view direction
    float sigmaV2 = dot(t, sigmaSq); // slope variance in view direction
    return meanFresnel(dot(V, N), sqrt(sigmaV2));
}

// CIE XYZ to sRGB
vec3 XYZ2RGB(vec3 xyz) {
  vec3 rgb;
  float X = xyz[0];
  float Y = xyz[1];
  float Z = xyz[2];
  rgb[0] = 3.2410*X - 1.5374*Y - 0.4986*Z;
  rgb[1] = -0.9692*X + 1.8760*Y + 0.0416*Z;
  rgb[2] = 0.0556*X - 0.2040*Y + 1.5070*Z; 
  return rgb;
}

// manual anisotropic filter
vec4 myTexture2DGrad(sampler2D tex, vec2 u, vec2 s, vec2 t)
{
    const float TEX_SIZE = 1024.0; // 'tex' size in pixels
    const int N = 1; // use (2*N+1)^2 samples
    vec4 r = vec4(0.0);
    float l = max(0.0, log2(max(length(s), length(t)) * TEX_SIZE) - 0.0);
    for (int i = -N; i <= N; ++i) {
        for (int j = -N; j <= N; ++j) {
            r += textureLod(tex, u + (s * float(i) + t * float(j)) / float(N), l);
        }
    }
    return r / pow(2.0 * float(N) + 1.0, 2.0);
}

// V, N, Tx, Ty in world space;
vec3 meanSkyRadiance(vec3 V, vec3 N, vec3 Tx, vec3 Ty, vec2 sigmaSq) {
    vec4 result = vec4(0.0);

    const float eps = 0.001;
    vec2 u0 = U(vec2(0.0), V, N, Tx, Ty);
    vec2 dux = 2.0 * (U(vec2(eps, 0.0), V, N, Tx, Ty) - u0) / eps * sqrt(sigmaSq.x);
    vec2 duy = 2.0 * (U(vec2(0.0, eps), V, N, Tx, Ty) - u0) / eps * sqrt(sigmaSq.y);

//    result = textureGrad(skySampler, u0 * (0.5 / 1.1) + 0.5, dux * (0.5 / 1.1), duy * (0.5 / 1.1));
//    result = textureGrad(skySampler, u0 * (0.5 / 1.1) + 0.5, dux, duy);

//    result = myTexture2DGrad(skySampler, u0 * (0.5 / 1.1) + 0.5, dux * (0.5 / 1.1), duy * (0.5 / 1.1));
    result = texture(skySampler, u0 * (0.5 / 1.1) + 0.5);

    result.rgb = abs(result.rgb / result.w);

    return result.rgb;
   
}


void main() {

	vec3 V = normalize(worldCamera - P);

	vec2 slopes = texture(fftWavesSampler, vec3(u / GRID_SIZES.x, 1.0)).xy;
	slopes += texture(fftWavesSampler, vec3(u / GRID_SIZES.y, 1.0)).zw;
	slopes += texture(fftWavesSampler, vec3(u / GRID_SIZES.z, 2.0)).xy;
	slopes += texture(fftWavesSampler, vec3(u / GRID_SIZES.w, 2.0)).zw;

	if(choppy > 0.0)
	{
		float Jxx, Jxy, Jyy;
		vec4 lambda = choppy_factor;
		// Jxx1..4 : partial Jxx
		float Jxx1 = texture(fftWavesSampler, vec3(u / GRID_SIZES.x, LAYER_JACOBIAN_XX)).r;
		float Jxx2 = texture(fftWavesSampler, vec3(u / GRID_SIZES.y, LAYER_JACOBIAN_XX)).g;
		float Jxx3 = texture(fftWavesSampler, vec3(u / GRID_SIZES.z, LAYER_JACOBIAN_XX)).b;
		float Jxx4 = texture(fftWavesSampler, vec3(u / GRID_SIZES.w, LAYER_JACOBIAN_XX)).a;
		Jxx = dot((lambda), vec4(Jxx1,Jxx2,Jxx3,Jxx4));

		// Jyy1..4 : partial Jyy
		float Jyy1 = texture(fftWavesSampler, vec3(u / GRID_SIZES.x, LAYER_JACOBIAN_YY)).r;
		float Jyy2 = texture(fftWavesSampler, vec3(u / GRID_SIZES.y, LAYER_JACOBIAN_YY)).g;
		float Jyy3 = texture(fftWavesSampler, vec3(u / GRID_SIZES.z, LAYER_JACOBIAN_YY)).b;
		float Jyy4 = texture(fftWavesSampler, vec3(u / GRID_SIZES.w, LAYER_JACOBIAN_YY)).a;
		Jyy = dot((lambda), vec4(Jyy1,Jyy2,Jyy3,Jyy4));

		slopes /= (1.0 + vec2(Jxx, Jyy));
	}

	vec3 N = normalize(vec3(-slopes.x, -slopes.y, 1.0));
	if (dot(V, N) < 0.0) {
		N = reflect(N, V); // reflects backfacing normals
	}

	float Jxx = dFdx(u.x);
	float Jxy = dFdy(u.x);
	float Jyx = dFdx(u.y);
	float Jyy = dFdy(u.y);
	float A = Jxx * Jxx + Jyx * Jyx;
	float B = Jxx * Jxy + Jyx * Jyy;
	float C = Jxy * Jxy + Jyy * Jyy;
	const float SCALE = 10.0;
	float ua = pow(A / SCALE, 0.25);
	float ub = 0.5 + 0.5 * B / sqrt(A * C);
	float uc = pow(C / SCALE, 0.25);
	vec2 sigmaSq = texture(slopeVarianceSampler, vec3(ua, ub, uc)).xw;

	sigmaSq = max(sigmaSq, 2e-5);

	vec3 Ty = normalize(vec3(0.0, N.z, -N.y));
	vec3 Tx = cross(Ty, N);

	vec3 Rf = vec3(0.0);
	vec3 Rs = vec3(0.0);
	vec3 Ru = vec3(0.0);

#if defined(SEA_CONTRIB) || defined(SKY_CONTRIB)
	float fresnel = 0.02 + 0.98 * meanFresnel(V, N, sigmaSq);
#endif

	vec3 Lsun;
	sunRadianceAndSkyIrradiance(worldCamera + earthPos, worldSunDir, Lsun);

	FragColor = vec4(0.0);

// // #ifdef SUN_CONTRIB
	Rs += reflectedSunRadiance(worldSunDir, V, N, Tx, Ty, sigmaSq) * Lsun / 10.0;
	FragColor.rgb = Rs;
// // #endif

// #ifdef SKY_CONTRIB
	vec3 Lsky = meanSkyRadiance(V, N, Tx, Ty, sigmaSq) / 100.0;
	Rs += fresnel * Lsky;
	FragColor.rgb = Rs;
// #endif

// #ifdef SEA_CONTRIB

// Esky is the total irradiance
	float omega = 2 * M_PI;
    float sunSr = 6.87e-5;        // http://en.wikipedia.org/wiki/Solid_angle
	vec3 Esky = Lsun * sunSr + Lsky * omega;
	vec3 Lsea = seaColor * Esky / M_PI;
	Ru += (1.0 - fresnel) * Lsea;
	FragColor.rgb += Ru;
// #endif

	FragColor.rgb = hdr(FragColor.rgb); 

}

#endif