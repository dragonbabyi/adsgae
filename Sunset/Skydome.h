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
 

class Skydome {
private:
    int radius;
    
    float solarElevation;
    float turbidity;
    float albedo;    
    
    float HosekRadiances[3];
    float HosekConfig[3][9];
    
    // GL vertex array object IDs
    enum {VARRAY, QUAD_VERTEX, NUM_VARRAYS};
    unsigned int varrayIDs[NUM_VARRAYS];
    
    // GL buffer object IDs
    enum {POSITION_BUFFER, INDEX_BUFFER, QUAD_VERTEX_BUFFER, QUAD_INDEX_BUFFER, NUM_BUFFERS};  //ADD NORMAL
    unsigned int bufferIDs[NUM_BUFFERS];
    
    GLuint depthrenderbuffer;
    GLuint skyframebuffer;
    
    // vertex attribute IDs
    unsigned int positionAttrib, posQuadAttrib;
    
    Program* skydomeShader[2];   //1 scattering 2 tone mapping
    
public:
    GLuint htex;  // htrue--happ
    GLuint skytexture;
    unsigned int skyTexSize;
    bool init;
    
  public:
    Skydome();
    ~Skydome();
    
    void loadProgram();
    
    void updateShader();
    void HosekSkyModel_Configuration(float solarAngle);
    void draw( GLFWwindow *win, float sunTheta);
    
};



#endif /* defined(__Sunset__Skydome__) */
