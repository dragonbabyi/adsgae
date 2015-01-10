//
//  Ocean.h
//  Sunset
//
//  Created by Yuping on 5/4/14.
//  Copyright (c) 2014 Yuping. All rights reserved.
//

#ifndef __Sunset__Ocean__
#define __Sunset__Ocean__

#include <iostream>
#include "vec4.h"
#include "mat4.h"
#include "Program.h"

// using core modern OpenGL
#include <GL/glew.h>
#include <GLFW/glfw3.h>

enum {
    // textures
    TEXTURE_SKY,
    TEXTURE_SPECTRUM12,
    TEXTURE_SPECTRUM34,
    TEXTURE_SLOPE_VARIANCE,
    TEXTURE_FFT_PING,
    TEXTURE_FFT_PONG,
    TEXTURE_BUTTERFLY,
    TEXTURE_GAUSSZ,
    TEXTURE_COUNT,
    
    // GL vertex array object IDs
    VARRAY_MESH = 0,
    VARRAY_QUAD,
    VARRAY_COUNT,
    
    // buffers
    BUFFER_GRID_INDEX = 0,
    BUFFER_GRID_VERTEX,
    BUFFER_QUAD_INDEX,
    BUFFER_QUAD_VERTEX,
    BUFFER_COUNT,
    
    // renderbuffers
    RENDERBUFFER_DEPTH = 0,
    RENDERBUFFER_COUNT,
    
    // framebuffers
    FRAMEBUFFER_FFT0 = 0,
    FRAMEBUFFER_FFT1,
    FRAMEBUFFER_SKY,
    FRAMEBUFFER_VARIANCES,
    FRAMEBUFFER_GAUSS,
    FRAMEBUFFER_COUNT,
    
    // programs
    PROGRAM_RENDER = 0,
//    PROGRAM_SHOW_SPECTRUM,
    PROGRAM_INIT,
    PROGRAM_VARIANCES,
    PROGRAM_FFTX,
    PROGRAM_FFTY,
//    PROGRAM_WHITECAP_PRECOMPUTE,
    PROGRAM_COUNT
};

class Ocean
{
private:
    
    GLuint      renderbuffers[RENDERBUFFER_COUNT];
    GLuint      framebuffers[FRAMEBUFFER_COUNT];
    GLuint      textures[TEXTURE_COUNT];
    GLuint      buffers[BUFFER_COUNT];
    GLuint      varrays[VARRAY_COUNT];
    
    Program*    programs[PROGRAM_COUNT];
    
    struct WinDow {
        int width, height;
    }window;
    
    //double updateTime;
    // vertex attribute IDs
    unsigned int positionAttrib;  //mesh
    unsigned int posQuadAttrib[PROGRAM_COUNT];
    
    int vboSize;
    int vboVertices;
    float *spectrum12;
    float *spectrum34;
    
public:
    
    Ocean(GLFWwindow *win);
    ~Ocean();
    void generateMesh();
    void loadPrograms();
    void generateWavesSpectrum();
    void simulateFFTWaves(float t);
    void computeSlopeVarianceTex();
    void drawQuad(int programindex);
    
    void updateShader();
    void draw(GLFWwindow *win);
    
    //debug
    void drawToScreen(GLFWwindow *win);
};



#endif /* defined(__Sunset__Ocean__) */
