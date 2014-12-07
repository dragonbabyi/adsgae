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
    
    camera.z 	    = 3.5f;
    camera.velx		= 0.0f;
    camera.vely		= 0.0f;
    camera.velz		= 0.0f;
    camera.x		= 0.0f;
    camera.y		= 0.0f;
    camera.theta 	= 27.0f; //27.0f  in degree    +: looking up
    camera.phi 		= 2.f;  //-625.0f
    camera.fovy 	= 30.0f;
    camera.vel		= 2.0f;
    
    sunTheta = 80.0 * M_PI / 180.0;
    sunPhi = 0.0;
    
//    if(vboParams.x != width ||
//       vboParams.y != height ||
//       vboParams.z != gridSize ||
//       vboParams.w != camera.theta)
//    {
//        generateMesh();
//        
//        vboParams.x = width;
//        vboParams.y = height;
//        vboParams.z = gridSize;
//        vboParams.w = camera.theta;
//    }
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
   
    uMatrix.view = mat4f::rotatey(camera.phi) * uMatrix.view;
	uMatrix.view = mat4f::rotatex(-camera.theta) * uMatrix.view;   //rotatex(angle in degree)
    
}

void Scene::sunmove()
{
    //update sun position
//    float sunRadius = 695800.0;
    float sunDistance = 149600000;   //(km)
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
