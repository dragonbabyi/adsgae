//
//  Scene.h
//  Sunset
//
//  Created by Yuping on 5/3/14.
//  Copyright (c) 2014 Yuping. All rights reserved.
//

#ifndef __Sunset__Scene__
#define __Sunset__Scene__

#include <iostream>
#include "vec4.h"
#include "mat4.h"


// using core modern OpenGL
#include <GL/glew.h>
#include <GLFW/glfw3.h>


class Scene {
public:
    struct UniformMatrix {
        mat4f view;
        mat4f proj;
        vec4f sunpos;
    } uMatrix;
    
    // GL uniform buffer IDs
    enum {MATRIX_BUFFER, NUM_BUFFERS};
    unsigned int bufferIDs[NUM_BUFFERS];
 
    int width, height;         // current window dimensions
    
    struct Cam {
        float velx, vely, velz, x, y, z, velt,velp, theta, phi, fovy, velf, vel;
    }camera;
    
//    float solarElevation;
    float sunPhi, sunTheta;
    float sunThetaVel;   // velocity
    
    vec4f vboParams;
    
public:
    Scene(GLFWwindow *win);
    
    // set up new window viewport and projection
    void viewport(GLFWwindow *win);
    
    // set view using orbitAngle
    void view();
    
    // update sun position
    void sunmove();

    // update shader uniform state each frame
    void update();

};

#endif /* defined(__Sunset__Scene__) */
