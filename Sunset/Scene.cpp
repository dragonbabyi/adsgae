//
//  Scene.cpp
//  Sunset
//
//  Created by Yuping on 5/3/14.
//  Copyright (c) 2014 Yuping. All rights reserved.
//

#include "AppContext.h"
#include "Scene.h"


// using core modern OpenGL
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <math.h>

#ifndef F_PI
#define F_PI 3.1415926f
#endif

Scene::Scene(GLFWwindow *win)
{
    glGenBuffers(NUM_BUFFERS, bufferIDs);
    glBindBuffer(GL_UNIFORM_BUFFER, bufferIDs[MATRIX_BUFFER]);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(UniformMatrix), 0, GL_STREAM_DRAW);
    glBindBufferBase(GL_UNIFORM_BUFFER, AppContext::MATRIX_UNIFORMS, bufferIDs[MATRIX_BUFFER]);
    
    camera.z 	    = 5.0f;  //3.5f;
    camera.velx		= 0.0f;
    camera.vely		= 0.0f;
    camera.velz		= 0.0f;
    camera.x		= 0.0f;
    camera.y		= 0.0f;
    camera.theta 	= -2.0f;   //-1.0f;
    camera.phi 		= 2.f;
    camera.fovy 	= 45.0f;//atan(14/len_mm)  // 22mm lens ~ 32.0f  200mm ~~ 4.0f
    camera.vel		= 2.0f;
    
    sunTheta = 89.0 * M_PI / 180.0;
    sunPhi = 0.0;
    sunThetaVel = 0;
    
    // update view
    viewport(win);
    view();
    sunmove();
}

void Scene::viewport(GLFWwindow *win)
{
    // get window dimensions
    glfwGetFramebufferSize(win, &width, &height);
    
    // this viewport makes a 1 to 1 mapping of physical pixels to GL
    // "logical" pixels
    glViewport(0, 0, width, height);
    
    // adjust 3D projection into this window
    float ch = camera.z;
    
    uMatrix.proj = mat4f::perspectiveProjection(camera.fovy,
                                                      float(width)/float(height),
                                                      0.1 * ch,
                                                      300000.0 * ch);
}


void Scene::view()
{
    // update view matrix
    uMatrix.view = mat4f(
                       0.0, -1.0, 0.0, -camera.x,
                       0.0, 0.0, 1.0, -camera.z,
                       -1.0, 0.0, 0.0, -camera.y,
                       0.0, 0.0, 0.0, 1.0
                       );
    
    //inverse
//    right up forward position
//    1	0	0	-1	0
//    2	-1	0	0	0
//    3	0	1	0	3.5
//    4	0	0	0	1
   
    uMatrix.view = mat4f::rotatey(camera.phi) * uMatrix.view;
//	uMatrix.view = mat4f::rotatex(-camera.theta) * uMatrix.view;   //rotatex(angle in degree)
    uMatrix.view = mat4f::rotatex(camera.theta) * uMatrix.view;
}

void Scene::sunmove()
{
    //update sun position
//    float sunRadius = 695800.0;
    float sunDistance = 149600000.0;   //(km)
    uMatrix.sunpos = vec4f(sin(sunTheta) * cos(sunPhi)*sunDistance , sin(sunTheta) * sin(sunPhi)*sunDistance , cos(sunTheta)*sunDistance , sunTheta);
 
}

//
// call before drawing each frame to update per-frame scene state
//
void Scene::update()
{
    view();
    sunmove();
    
    // update uniform block
    glBindBuffer(GL_UNIFORM_BUFFER, bufferIDs[MATRIX_BUFFER]);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(UniformMatrix), &uMatrix);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);

}

