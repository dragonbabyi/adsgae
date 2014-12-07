//
//  Skydome.h
//  Sunset
//


#ifndef __Sunset__Skydome__
#define __Sunset__Skydome__

#include <iostream>
#include <vector>
#include "vec4.h"
#include "mat4.h"
#include "Program.h" 

// using core modern OpenGL
#include <GL/glew.h>
#include <GLFW/glfw3.h>


//unsigned int skyTexSize = 256;
// struct vec3 {
//     float x,y,z;
//     vec3(float X=0, float Y=0, float Z=0):x(X),y(Y),z(Z) {};
// };

class Skydome {
private:
    int radius;
    
    float solarElevation;
    float turbidity;
    float albedo;   //albedo ->= 0, because open ocean albedo is near 0.
    
    float HosekRadiances[3];
    float HosekConfig[3][9];
    
    // GL vertex array object IDs
    enum {VARRAY, NUM_VARRAYS};
    unsigned int varrayIDs[NUM_VARRAYS];
    
    // GL buffer object IDs
    enum {POSITION_BUFFER, INDEX_BUFFER, NUM_BUFFERS};  //ADD NORMAL
    unsigned int bufferIDs[NUM_BUFFERS];
    
    // vertex attribute IDs
    unsigned int positionAttrib;
//    unsigned int normalAttrib;
    
//    GLfloat *data;
//    GLuint DatasetTexture;
    
    Program* skydomeShader[1];
    
public:
    Skydome();
    ~Skydome();
    
    void loadProgram();
    void generateTexture();
    void updateShader();
    void HosekSkyModel_Configuration();
    void draw( GLFWwindow *win, float sunTheta);
    
};



#endif /* defined(__Sunset__Skydome__) */
