
#ifdef _VERTEX_

layout(std140) uniform Matrices {
    mat4 view;
    mat4 proj;
    vec4 sunpos;
};

in vec3 vPosition;

out vec4 p;
// out float solarElevation;
out vec4 sun;


mat4 projMatrix = transpose(proj);
mat4 modelView = transpose(view);

void main()  {
     
    p = vec4(vPosition, 1.0);
    // solarElevation = sunpos.w;
    sun = sunpos;
    gl_Position = projMatrix * modelView * vec4(vPosition*1000.0, 1.0);   // radius
}


#endif

#ifdef _FRAGMENT_ 
precision highp float;

// uniform sampler1D DatasetTexture; 
// uniform float albedo;
// uniform float turbidity;
// //uniform float solarElevation;

uniform vec3 HosekRadiances;
uniform vec3 XHosekABC;
uniform vec3 XHosekDEF;
uniform vec3 XHosekGHI;
uniform vec3 YHosekABC;
uniform vec3 YHosekDEF;
uniform vec3 YHosekGHI;
uniform vec3 ZHosekABC;
uniform vec3 ZHosekDEF;
uniform vec3 ZHosekGHI;
 
in vec4 sun;
in vec4 p;
out vec4 FragColor;

#ifndef NIL
#define NIL                         0
#endif

#ifndef MATH_PI 
#define MATH_PI                     3.141592653589793
#endif

#ifndef MATH_DEG_TO_RAD
#define MATH_DEG_TO_RAD             ( MATH_PI / 180.0 )
#endif

#ifndef MATH_RAD_TO_DEG
#define MATH_RAD_TO_DEG             ( 180.0 / MATH_PI )
#endif

#ifndef DEGREES
#define DEGREES                     * MATH_DEG_TO_RAD
#endif

#ifndef TERRESTRIAL_SOLAR_RADIUS
#define TERRESTRIAL_SOLAR_RADIUS    ( ( 0.51 DEGREES ) / 2.0 )
#endif
 

float HosekWilkie(vec3 ABC, vec3 DEF, vec3 GHI, float theta, float gamma)
{   // F(theta, r)
    float expM = exp(DEF.y * gamma);
    float rayM = cos(gamma)*cos(gamma);
    float mieM = (1.0 + cos(gamma)*cos(gamma)) / pow((1.0 + GHI.z*GHI.z - 2.0*GHI.z*cos(gamma)), 1.5);
    float zenith = sqrt(cos(theta));

    return (1.0 + ABC.x * exp(ABC.y / (cos(theta) + 0.01))) 
           * (ABC.z + DEF.x * expM + DEF.z * rayM + GHI.x * mieM + GHI.y * zenith);

}

vec3 XYZtoxyY(vec3 XYZ)
{
    vec3 xyY;
    xyY.x = XYZ.x / (XYZ.x + XYZ.y + XYZ.z);
    xyY.y = XYZ.y / (XYZ.x + XYZ.y + XYZ.z);
    xyY.z = XYZ.y;

    return xyY;
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

// CIE XYZ to CIE RGB  
// vec3 xyz2rgb(vec3 xyz) {
//   vec3 rgb;
//   float X = xyz[0];
//   float Y = xyz[1];
//   float Z = xyz[2];
//   rgb[0] = 0.41847*X - 0.15866*Y - 0.082835*Z;
//   rgb[1] = -0.091169*X + 0.25243*Y + 0.015708*Z;
//   rgb[2] = 0.0009209*X - 0.0025498*Y + 0.1786*Z;
//   return rgb;
// }

vec4 xyYtoRGB(vec3 xyY)
{
    vec3 XYZ;
    vec4 RGB;

    // Convert xyY to XYZ
    XYZ.x = xyY.x * (xyY.z / xyY.y);
    XYZ.y = xyY.z;
    XYZ.z = (1.0 - xyY.x - xyY.y) * (xyY.z / xyY.y);
     
    // Convert XYZ to RGB
    RGB.xyz = XYZ2RGB(XYZ);

    // Deal with negative values
    float w = min(0.0, RGB.x);
    w = min(w, RGB.y);
    w = min(w, RGB.z);
    RGB += -w;

    // Scale down if necessary, preserving color
    float f = max(RGB.y, RGB.z);
    f = max(RGB.x, f); // f is now largest color component
    f = max(1.0, f); // if f is less than 1.0, set to 1.0 (to do nothing)
    RGB.xyz = RGB.xyz / f;

    // Apply gamma correction
    float oneOverGamma = 1.0/2.2;
    RGB = pow(RGB, vec4(oneOverGamma, oneOverGamma, oneOverGamma, oneOverGamma));

    RGB.w = 1.0;

    RGB = clamp(RGB, 0.0, 1.0);

    return RGB;
}


vec3 toneMap(vec3 xyY, float mC, float mR, float Ldmax)
{
    // This is based on Durand's operator, which is based on Ferwerda
    // which is based on Ward...

    // Convert Kcd/m2 to cd/m2
    float Y = xyY.z * 1000.0;

    // deal with negative luminances (nonsensical)
    // max() doesn't work here for some reason...
    if (Y < 0.0)
    {
        Y = 0.0;
    }

    vec3 XYZ;
    float R;
    // Convert xyY to XYZ
    XYZ.x = xyY.x * (xyY.z / xyY.y);
    XYZ.y = xyY.z;
    XYZ.z = (1.0 - xyY.x - xyY.y) * (xyY.z / xyY.y);

    const vec3 scotopic = vec3(-0.702, 1.039, 0.433);
    R = dot(XYZ, scotopic);

    // Apply perceptual blue-shift in scotopic conditions (Durand)
    // const vec3 blue = vec3(0.3, 0.3, 1.0);
    // xyY = mix(xyY, blue, k);   // blue-shift for viewing night scenes

    // Apply the photopic and scotopic scales that were precomputed in
    // the app per Ferwerda and Durand's operator
    float Ldp = Y * mC;
    float Lds = R * mR;

    Y = Ldp + 0.5 * Lds;
    // Normalize to display range
    Y = Y / Ldmax;

    xyY.z = Y;

    return xyY;
}

//wiki: Atmosphere_model. The U.S. Standard Atmosphere.
float P(float altitude) {    
    return 101.0;
}

float T(float altitude) {
    return 15.0;
}


void main() 
{
    vec3 XYZ;
    //
    // if(p.z < 0) {
    //     FragColor = vec4(1.0, 0.0, 0.0, 1.0);
    //     return;
    // }
    
    float solarElevation = MATH_PI/2.0 - sun.w;
    vec4 debug;
    // ArHosekSkyModelState skymodel_state = arhosek_xyz_skymodelstate_alloc_init(solarElevation, debug);
   
    float nr = sqrt(p.x*p.x + p.y*p.y)/sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
    float th = asin(nr); 
    float ph = atan(p.y, p.x);

    // angle between viewpoint to the sun and lookat direction, computed from law of cosines  
    float cangle = cos(solarElevation) * sin(th) * cos(ph) + sin(solarElevation) * cos(th);
    if(cangle > 1.0)  { cangle = 1.0; } 
    if(cangle < -1.0) { cangle = -1.0;}   
    float gamma = acos(cangle);
    float theta = th;

    // Update: moved the state conputation to CPU code
    XYZ[0] = HosekWilkie(XHosekABC, XHosekDEF, XHosekGHI, theta, gamma) * HosekRadiances.x;
    XYZ[1] = HosekWilkie(YHosekABC, YHosekDEF, YHosekGHI, theta, gamma) * HosekRadiances.y;
    XYZ[2] = HosekWilkie(ZHosekABC, ZHosekDEF, ZHosekGHI, theta, gamma) * HosekRadiances.z;
 
/*
    if(abs(sin(solarElevation+th)) > 1.000)   {
    //because the float point error, this case happens and lead to the black spot close to sun!!!!!
    //change to 1.0001 then won't happen 
        FragColor = vec4(0.0, 0.0, abs(sin(solarElevation+th)) - 1.0, 1.0);
    }
*/

    // XYZ to xyY
    vec3 xyYsHW = XYZtoxyY(XYZ);
    // float isothermalEffect = 2.5 - theta; // add isothermal effects (sky darkens with altitude)
    // xyYsHW.z *= isothermalEffect;
    // tone mapping
    float mC = 0.75;
    float mR = 0.19;
    // float k = 0.75;   // blue-shift
    float Ldmax = theta*80.0;   //maximum image luminance
    vec3 xyYs = toneMap(xyYsHW, mC, mR, Ldmax);

    vec4 skyColor = xyYtoRGB(xyYs);

    FragColor = skyColor; 
    

    //************** add sun color  *****************

    float sun_radius = 695800;   //(km)
    
    vec3 viewpoint = vec3(0.0, 0.0, 3.5);  // scene->camera.xyz
    // apparent altitude
    float r = sqrt((p.x - viewpoint.x)*(p.x - viewpoint.x)+(p.y - viewpoint.y)*(p.y - viewpoint.y));
    float ha = atan((p.z - viewpoint.z)/r) * 180.0 / MATH_PI;  // in degree
  
    //refraction index formula  
    float P = 101;    // pressure P(h)
    float T = 15;     // temperature 
    float prefix = (P/101.0) * 283.0/(273.0 + T); 
    
    float R = prefix/tan((ha + 7.31/(ha + 4.4)) * MATH_PI / 180.0);   //in arcminutes
    float h = ha - R/60.0;   // true altitude in degree

    // tracing the "true ray" and see whether there is a intersection with the sphere
    // if has intersection, mark as sun color
    vec3 viewdir = vec3(p.x, p.y, r*tan(h * MATH_PI / 180.0) + viewpoint.z);
    vec3 dc = viewpoint - sun.xyz;

    // solve p=r.start-center + t*r.direction; p*p -radius^2=0  
    float a = dot(viewdir, viewdir);
    float b = 2*dot(viewdir, dc);
    float c = dot(dc, dc) - sun_radius * sun_radius;

    float discriminant = b*b - 4*a*c;
    if (discriminant > 0)   // intersect with real sun
    {
        FragColor = vec4(1.0, 1.0, 1.0, 1.0);
    }

}
 
#endif