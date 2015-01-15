// sky dome pass 2
// tone mapping

#ifdef _VERTEX_

layout(std140) uniform Matrices {
    mat4 view;
    mat4 proj;
    vec4 sunpos;
};

in vec3 vPosition;
// in vec4 position;  

out vec4 p; 
out vec4 sun;

mat4 projMatrix = transpose(proj);
mat4 modelView = transpose(view);

void main() {
    p = vec4(vPosition*10.0, 1.0);
    sun = sunpos;
    gl_Position = projMatrix * modelView * vec4(vPosition*1000.0, 1.0); 
}


#endif

#ifdef _FRAGMENT_ 

#ifndef MATH_PI 
#define MATH_PI                     3.141592653589793
#endif
precision highp float;

uniform sampler2D skySample2D; 
uniform sampler1D htrueSample1D;
  
in vec4 p;
in vec4 sun; 
// in vec2 texcd;
out vec4 FragColor;


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

    RGB = clamp(RGB, 0.0, 1.0);
    RGB.w = 1.0;

    return RGB;
}

// interactive calibration
vec3 toneMap5(vec3 XYZ, float xmean) 
{
    float a = 0.0;   //aperture
    float c = 50.0;
     
    vec3 skyColor = XYZ;
  
    float s = pow(2, a+1) * xmean /(1.0 + c);  
    float e = s * c;
    float n = 0.5;
    float u = (pow(1, n) - pow(0.02, n))/(pow(e, n) - pow(s, n));
    float v = pow(1, n) - u*pow(e, n);

    for (int i = 0; i < 3; ++i)
    {
        clamp(skyColor[i], s, e); 
        skyColor[i] = pow( u * pow(skyColor[i], n) + v, 1.0/n);
    }

    return skyColor;
}


// vec3 hdr(vec3 L) {
//     // L = L * hdrExposure;
//     L.r = L.r < 1.413 ? pow(L.r * 0.38317, 1.0 / 2.2) : 1.0 - exp(-L.r);
//     L.g = L.g < 1.413 ? pow(L.g * 0.38317, 1.0 / 2.2) : 1.0 - exp(-L.g);
//     L.b = L.b < 1.413 ? pow(L.b * 0.38317, 1.0 / 2.2) : 1.0 - exp(-L.b);
//     return L;
// }


// const float SCALE = 1000.0;
// const float Rg = 6360.0; // (m)
//  // Rayleigh
// const float HR = 8.0 * SCALE;
// const vec3 betaR = vec3(5.8e-3, 1.35e-2, 3.31e-2) / SCALE;

// // Mie
// // DEFAULT
// const float HM = 1.2 * SCALE;
// const vec3 betaMSca = vec3(4e-3) / SCALE;
// const vec3 betaMEx = betaMSca / 0.9;

// ////////
// // optical depth for ray (r,mu) of length d, using analytic formula
// // (mu=cos(view zenith angle)), intersections with ground ignored
// // H=height scale of exponential density function
// float opticalDepth(float H, float r, float mu, float d) {
//     float a = sqrt((0.5/H)*r);
//     vec2 a01 = a*vec2(mu, mu + d / r);
//     vec2 a01s = sign(a01);
//     vec2 a01sq = a01*a01;
//     float x = a01s.y > a01s.x ? exp(a01sq.x) : 0.0;
//     vec2 y = a01s / (2.3193*abs(a01) + sqrt(1.52*a01sq + 4.0)) * vec2(1.0, exp(-d/H*(d/(2.0*r)+mu)));
//     return sqrt((6.2831*H)*r) * exp((Rg-r)/H) * (x + dot(y, vec2(1.0, -1.0)));
// }

// // transmittance(=transparency) of atmosphere for ray (r,mu) of length d
// // (mu=cos(view zenith angle)), intersections with ground ignored
// // uses analytic formula instead of transmittance texture
// vec3 analyticTransmittance(float r, float mu, float d) {
//     return exp(- betaR * opticalDepth(HR, r, mu, d) - betaMEx * opticalDepth(HM, r, mu, d));
// }
 

void main()
{

    // float sinSkyTheta = sqrt(p.x*p.x + p.y*p.y)/sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
    // float SkyTheta = asin(sinSkyTheta);    // sky theta  
    float tanSkyTheta = sqrt(p.x*p.x + p.y*p.y)/p.z;
    float SkyTheta = atan(tanSkyTheta);
    float nr = abs(SkyTheta) / MATH_PI ;   // map th(radians) 0 ~ PI/2, scale 1/2
    float ph = atan(p.y, p.x);   //  -PI ~ PI
    vec2 texcoord = vec2(cos(ph)*nr + 0.5, sin(ph)*nr + 0.5);
 
	vec4 XYZW = texture(skySample2D, texcoord.st);
    vec3 XYZ = XYZW.xyz / XYZW.w; 
    //************* tone mapping *************
 
    vec3 rgb = XYZ2RGB(XYZ);
    // vec3 skyc = mix(pow(0.38317*rgb,vec3(1.0/2.2)), 1.-exp(-rgb), step(1.413,rgb));
    vec3 avg = textureLod (skySample2D, vec2(0.5, 0.5), 10.0).xyz;
    float hdrExposure = 5.0 - 2.75 * cos(sun.w);
    // float hdrExposure = 5.0;
    vec3 Lrgb = XYZ2RGB(avg);
    float avglum = 0.2126 * Lrgb.x + 0.7152 * Lrgb.y + 0.0722 * Lrgb.z;
    avglum *= hdrExposure;

    vec3 skyc = toneMap5(rgb, avglum); 

    FragColor = vec4(skyc, 1.0);

	//************** add sun  *****************

    float sun_radius = 695800;   //(km)
    
    vec3 viewpoint = vec3(0.0, 0.0, 0.035);  // scene->camera.xyz
    // apparent altitude
    float d = sqrt((p.x - viewpoint.x)*(p.x - viewpoint.x)+(p.y - viewpoint.y)*(p.y - viewpoint.y));
    float ha = atan((p.z - viewpoint.z)/d) * 180.0 / MATH_PI;  // in degree
   
    // //refraction index formula  
    float P = 101;    // pressure P(h)
    float T = 25;     // temperature 
    float prefix = (P/101.0) * 283.0/(273.0 + T); 
    
    float R, h, index;
    // R = prefix/tan((ha + 7.31/(ha + 4.4)) * MATH_PI / 180.0);   //in arcminutes
    // h = ha - R/60.0;   // true altitude in degree

    //////////////////////
    // if (ha > 0.2)
    // {
        R = prefix/tan((ha + 7.31/(ha + 4.4)) * MATH_PI / 180.0);   //in arcminutes
        h = ha - R/60.0;   // true altitude in degree
    // }
    // else {// using the precomputed true altitude
    //     // int index = int((ha*60.0 + 4.0)/0.25);
    //     // h = texelFetch(htrueSample1D, index, 0).r;
    //     index = (ha * 60.0 + 8.0) * 8.0;
    //     if (index < 0.0)
    //     {
    //         index = 0.0;
    //     }
    //     // float h1 = texelFetch(htrueSample1D, int(floor(index)), 0).r;
    //     // float h2 = texelFetch(htrueSample1D, int(ceil(index)), 0).r;
    //     // // interpolation
    //     // h = (index - floor(index)) * h1/60.0 + (ceil(index) - index) * h2/60.0;  //in degree
    //     // h = 0.0;

    //    h = texture(htrueSample1D, index/256.0).r/60.0;
    //}

    //////////////////////
    float distance = 149600000.0; // (km)
    const float sunSolidAngle = 0.0000687;
    float halfangle = acos(1.0 - (sunSolidAngle / MATH_PI / 2.0));  
    float suncenterTheta = sun.w + halfangle;   // 
 
    float Sun_ha = (MATH_PI/2.0 - suncenterTheta) * 180.0 / MATH_PI;     // apparent altitude in degree
    float Sun_R = prefix/tan((Sun_ha + 7.31/(Sun_ha + 4.4)) * MATH_PI / 180.0); //in arcminutes
    
    // tracing the "true ray" and see whether there is a intersection with the sphere
    // if has intersection, mark as sun color
    // vec3 viewdir = vec3(p.x, p.y, d*tan(h * MATH_PI / 180.0) + viewpoint.z);
    vec3 viewdir = normalize(vec3(p.x, p.y, d*tan(h * MATH_PI / 180.0) + viewpoint.z));  
    float Sun_theta = suncenterTheta + (Sun_R/60.0) * MATH_PI / 180.0;  // sun center true theta in radians
    // sun_phi = 0.0
    vec3 realSun = vec3(sin(Sun_theta)*distance, 0.0, cos(Sun_theta)*distance);
    vec3 dc = viewpoint - realSun;

    // vec3 dc = viewpoint - sun.xyz;
    // vec3 dc = - sun.xyz;

    // solve p=r.start-center + t*r.direction; p*p -radius^2=0  
    float a = dot(viewdir, viewdir);
    float b = 2*dot(viewdir, dc);
    float c = dot(dc, dc) - sun_radius * sun_radius;

    float discriminant = b*b - 4*a*c;
    if (discriminant > 0)   // intersect with real sun
    {
        //worldCamera + earthPos, worldSunDir
     // //  vec3 worldV = vec3(0.0, 0.0, 1.0); // vertical vector
        // vec3 worldSunDir = vec3(sin(suncenterTheta), 0.0, cos(suncenterTheta));
        // vec3 worldSunDir = normalize(vec3(sin(suncenterTheta)*distance, 0.0, cos(suncenterTheta)*distance - viewpoint.z));
        vec3 worldSunDir = normalize(vec3(sin(Sun_theta)*distance, 0.0, cos(Sun_theta)*distance - viewpoint.z));
        // vec3 worldSunDir = normalize(vec3(sun.x, sun.y, sun.z));

        // float r = 6360045.0;   // (m)
        // float muS = dot(worldV, worldSunDir); 
        // vec3 sunColor = analyticTransmittance(r, muS, Rg) ;   // radiance

        vec3 sunColor = vec3(253, 14, 19) / 255.0;
        // // limb darkening 
        float u = 0.3;  
        float d = acos(dot(viewdir, worldSunDir)) / halfangle;   // distance from the center of the sun 
        // avoid middle black
        if ( d > 0.1 )
        {
            float limbdarkening = 1.0 - u *(1.0 - sqrt(1.0 - d*d));  // /(sun_radius*sun_radius)
            sunColor *= limbdarkening;
        }
        // float limbdarkening = 1.0 - u *(1.0 - sqrt(1.0 - d*d));  // /(sun_radius*sun_radius)
        // sunColor *= limbdarkening;
         
        // FragColor = vec4(limbdarkening, 0.0, 0.0, 1.0);
        FragColor = vec4(sunColor, 1.0);
        // FragColor = vec4( texelFetch(htrueSample1D, 0.3, 0).r, 0.0, 0.0, 1.0);

    }



}




#endif
