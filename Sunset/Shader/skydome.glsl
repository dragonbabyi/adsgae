// skydome pass 1
// render to Spherical sky texture

#ifdef _VERTEX_

layout(std140) uniform Matrices {
    mat4 view;
    mat4 proj;
    vec4 sunpos;
};

in vec4 position;

out vec4 sun;
out vec2 vTexcoord;
 
void main()  {
    
    sun = sunpos;
	vTexcoord = position.zw;
    gl_Position = vec4(position.xy, 0.0, 1.0); 
}

#endif

#ifdef _FRAGMENT_ 
precision highp float;

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
in vec2 vTexcoord;
 
layout(location = 0) out vec4 FragColor;

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
 

float HosekWilkie(vec3 ABC, vec3 DEF, vec3 GHI, float costheta, float gamma, float cosGamma)
{   // F(theta, r)
    float expM = exp(DEF.y * gamma);
    float rayM = cosGamma * cosGamma;
    float mieM = (1.0 + cosGamma*cosGamma) / pow((1.0 + GHI.z*GHI.z - 2.0*GHI.z*cosGamma), 1.5);
    float zenith = sqrt(costheta);

    return (1.0 + ABC.x * exp(ABC.y / (costheta + 0.01))) 
           * (ABC.z + DEF.x * expM + DEF.z * rayM + GHI.x * mieM + GHI.y * zenith);

}


void main() 
{
    vec3 XYZ;
    
    float solarElevation = MATH_PI/2.0 - sun.w;
    
    float rx = vTexcoord.x - 0.5;
    float ry = vTexcoord.y - 0.5;
    float nr = sqrt(rx*rx + ry*ry)/0.5;
    float th = nr * 0.5 * MATH_PI;    // linear store the color by degree
    float theta = th; //asin(nr);     
    float ph = atan(ry, rx);

    if(nr > 1.0) {
    	// discard;
        FragColor = vec4(0.0, 0.0, 0.0, 0.0);
    }else{
        // angle between viewpoint to the sun and lookat direction, computed from law of cosines  
        float cosGamma = cos(solarElevation) * sin(theta) * cos(ph) + sin(solarElevation) * cos(theta);
        if(cosGamma > 1.0)  { cosGamma = 1.0; } 
        if(cosGamma < -1.0) { cosGamma = -1.0;}   
        float gamma = acos(cosGamma);
        float costheta = cos(theta);

        XYZ[0] = HosekWilkie(XHosekABC, XHosekDEF, XHosekGHI, costheta, gamma, cosGamma) * HosekRadiances.x;
        XYZ[1] = HosekWilkie(YHosekABC, YHosekDEF, YHosekGHI, costheta, gamma, cosGamma) * HosekRadiances.y;
        XYZ[2] = HosekWilkie(ZHosekABC, ZHosekDEF, ZHosekGHI, costheta, gamma, cosGamma) * HosekRadiances.z;
         
        FragColor = vec4(XYZ, 1.0);  
    }

}
 
#endif
