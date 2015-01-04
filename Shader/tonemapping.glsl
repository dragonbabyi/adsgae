// sky dome pass 2
// tone mapping

#ifdef _VERTEX_

layout(std140) uniform Matrices {
    mat4 view;
    mat4 proj;
    vec4 sunpos;
};

// in vec3 vPosition;
in vec4 position;  

// out vec4 p; 
out vec4 sun;
//test
out vec2 texcd;

mat4 projMatrix = transpose(proj);
mat4 modelView = transpose(view);

void main() {
 //    p = vec4(vPosition, 1.0);
     sun = sunpos;
 //    gl_Position = projMatrix * modelView * vec4(vPosition*1000.0, 1.0); 

    texcd = position.zw;
    gl_Position = vec4(position.xy, 0.0, 1.0);
}


#endif

#ifdef _FRAGMENT_ 

#ifndef MATH_PI 
#define MATH_PI                     3.141592653589793
#endif

precision highp float;

uniform sampler2D skySample2D; 
  
// in vec4 p;
in vec4 sun; 
in vec2 texcd;
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

    RGB.w = 1.0;

    RGB = clamp(RGB, 0.0, 1.0);

    return RGB;
}

//Durand's operator
vec3 toneMap(vec3 xyY, float mC, float mR, float k, float Ldmax)
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

    Y = Ldp + k * Lds;
    // Normalize to display range
    Y = Y / Ldmax;

    xyY.z = Y;

    return xyY;
}

vec3 Uncharted2Tonemap(vec3 x)
{
    float A = 0.15;
    float B = 0.50;
    float C = 0.10;
    float D = 0.20;
    float E = 0.02;
    float F = 0.30;
    float W = 11.2;
    return ((x*(A*x+C*B)+D*E)/(x*(A*x+B)+D*F))-E/F;
}

// Uncharted 2 operator by Jahn Hable
vec3 toneMap2(vec3 rbg)
{
    float W = 11.2;

    rbg *= 16;   //exposure_adjustment
    float ExposureBias = 2.0f;
    vec3 curr = Uncharted2Tonemap(ExposureBias * rbg);

    vec3 whiteScale = 1.0f/Uncharted2Tonemap(vec3(W, W, W));
    vec3 RGB = curr*whiteScale;
   // Apply gamma correction
    float oneOverGamma = 1.0/2.2;
    vec3 res = pow(RGB, vec3(oneOverGamma, oneOverGamma, oneOverGamma));

    return res;
}

//Jim Hejl and Richard Burgess-Dawson
vec3 toneMap3(vec3 color, float Ldmax)
{
    // color *= exposure_adjustment;
    
    float x = max(0,color.x-0.004);
    float y = max(0,color.y-0.004);
    float z = max(0,color.z-0.004);

    vec3 retColor;
    retColor.x = (x*(6.2*x+.5))/(x*(6.2*x+1.7)+0.06);
    retColor.y = (y*(6.2*y+.5))/(y*(6.2*y+1.7)+0.06);
    retColor.z = (z*(6.2*z+.5))/(z*(6.2*z+1.7)+0.06);
    return retColor;
}

//Reinhard, stark, shirley, ferwerda
float toneMap4(vec3 XYZ, float avglog)
{
    float lum = dot(XYZ, vec3(0.2126, 0.7152, 0.0722));

    int number_of_pixel = 256;
    float image_key = exp(number_of_pixel * avglog) / number_of_pixel;     //eq1
     
    float middle_grey = 0.18;
    float scaledL = lum * middle_grey / image_key;    //ea2
    float Y = scaledL / (1+scaledL);        //eq3
     
    return Y;
}


// interactive calibration
vec3 toneMap5(vec3 color, float avglum) 
{

    return color;
}

void main()
{
    // get the true value from texture
    // p, sun --> texcoord
/*
    float nr = sqrt(p.x*p.x + p.y*p.y)/sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
    float th = asin(nr); 
    float r = th / MATH_PI / 2.0;
    float ph = atan(p.y, p.x);
    vec2 texcoord = vec2(cos(ph)*r, sin(ph)*r);

	vec3 XYZ = texture(skySample2D, texcoord.st).xyz;
*/
	vec3 XYZ = texture(skySample2D, texcd.st).xyz;

    

    //************* tone mapping *************

	// XYZ to xyY
    vec3 xyYsHW = XYZtoxyY(XYZ);
    float isothermalEffect = MATH_PI/2 - sun.w; // add isothermal effects (sky darkens with altitude)
    xyYsHW.z *= isothermalEffect;

    float mC = 0.05;    //  m = Ld/Lw // display/real
    float mR = 4.0;
    float k = 0.2;   // blue-shift

	vec2 LmaxTexcoord = vec2(1.0 - sin(sun.w)/2.0, 0.5); 
    vec3 Lrgb = texture(skySample2D, LmaxTexcoord.st).xyz;
     Lrgb = XYZ2RGB(Lrgb);
    float Ldmax = 0.2126 * Lrgb[0] + 0.7152 * Lrgb[1] + 0.0722 * Lrgb[2];
    
    vec3 xyYs = toneMap(xyYsHW, mC, mR, k, Ldmax*100.0);
    vec4 skyColor = xyYtoRGB(xyYs);
    FragColor = skyColor;
    
    //************** 
    // vec4 skyColor = xyYtoRGB(xyYsHW);
    // FragColor = skyColor; 
    //************** 
    // vec3 avg = texture(skySample2D, vec2(0.5, 0.5), 8.0).xyz;
    // float avglum = 0.2126 * avg.x + 0.7152 * avg.y + 0.0722 * avg.z;

    // vec3 rgb = XYZ2RGB(XYZ);
    // vec3 skyc = toneMap4(rgb, avglum);   // 50 * sun.w
    // FragColor = vec4(skyc, 1.0); 

    // FragColor = vec4(XYZ, 1.0);
    //************** 
    // float avglog = texture(skySample2D, vec2(0.5, 0.5), 8.0).w;
    // float scaledL = toneMap4(XYZ, avglog);    
    // vec3 xyYs = XYZtoxyY(XYZ);
    // xyYs.z = scaledL; 
    // vec4 skyColor = xyYtoRGB(xyYs);
    // FragColor = skyColor; 
   
	//************** add sun color  *****************
/*
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
*/

}




#endif
