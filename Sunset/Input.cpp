//
//  Input.cpp
//  Sunset
//
//  Created by Yuping on 5/4/14.
//  Copyright (c) 2014 Yuping. All rights reserved.
//

#include "Input.h"
#include "Scene.h"
#include "Skydome.h"
#include "ocean.h"

//// using core modern OpenGL
//#include <GL/glew.h>
//#include <GLFW/glfw3.h>

#include <math.h>
#include <algorithm>

#ifndef F_PI
#define F_PI 3.1415926f
#endif

//
// called when a mouse button is pressed.
// Remember where we were, and what mouse button it was.
//
void Input::mousePress(GLFWwindow *win, int b, int action)
{
    if (action == GLFW_PRESS) {
        // hide cursor, record button
        glfwSetInputMode(win, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
        button = b;
    }
    else {
        // display cursor, update button state
        glfwSetInputMode(win, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
        button = -1;       // no button
    }
}

//
// called when the mouse moves
// use difference between oldX,oldY and x,y to define a rotation
//
void Input::mouseMove(GLFWwindow *win, Scene *scene, double x, double y)
{
    // only update view after at least one old position is stored
    if (button == GLFW_MOUSE_BUTTON_LEFT && button == oldButton) {
        // record differences & update last position
        float dx = float(x - oldX);
        float dy = float(y - oldY);
          
        scene->camera.theta += dy/40.0;
        scene->camera.phi += dx/40.0;
        
        if(scene->camera.theta  < -87.f)
            scene->camera.theta = -87.f;
        if(scene->camera.theta  > 89.f)
            scene->camera.theta = 89.f;
        
    }
    
    // update prior mouse state
    oldButton = button;
    oldX = x;
    oldY = y;
}

//
// called when any key is pressed
//
void Input::keyPress(GLFWwindow *win, int key, Scene *scene, Skydome *skydome, Ocean *ocean)
{
    switch (key) {
        case 'B':
            scene->sunThetaVel = 1.0;
            break;
        case 'N':
            scene->sunThetaVel = -1.0;
            break;
        case 'C':
            scene->sunThetaVel = 0.01;
            break;
        case 'V':
            scene->sunThetaVel = -0.01;
            break;
            
        //camera.fovy
        case 'F':
            scene->camera.velf = 1.0;
            break;
        case 'G':
            scene->camera.velf = -1.0;
            break;
        case 'D':
            scene->camera.fovy = 5.0;
            break;
        case 'H':
            scene->camera.fovy = 45.0;
            break;
        case 'O':
            scene->camera.phi = 2.0;   // in angle
            break;
        case 'P':
            scene->camera.phi = -182.0;   // opposite side
            break;
     
        //reload the shaders
        case 'R':   //int = 82
            skydome->updateShader();
            //ocean->updateShader();
            break;
        // Escape: exit
        case GLFW_KEY_ESCAPE:
            glfwSetWindowShouldClose(win, true);
            break;
            
    }
    
}

//
// called when any key is released
//
void Input::keyRelease(GLFWwindow *win, int key, Scene *scene)
{
    switch (key) {
        case 'B':
            scene->sunThetaVel = 0.0;
            break;
        case 'N':
            scene->sunThetaVel = 0.0;
            break;
        case 'C':
            scene->sunThetaVel = 0.0;
            break;
        case 'V':
            scene->sunThetaVel = 0.0;
            break;
            
        //camera.fovy
        case 'F':
            scene->camera.velf = 0.0;
            break;
        case 'G':
            scene->camera.velf = 0.0;
            break;
            
    }
}

//
// update view if necessary based on a/d keys
//
void Input::keyUpdate(Scene *scene)
{
    if (scene->camera.velf != 0 || scene->sunThetaVel != 0) {

        double now = glfwGetTime();
        double delta = (now - updateTime);
        
        scene->camera.fovy += scene->camera.velf * delta;

        if (scene->sunThetaVel != 0 &&scene->sunThetaVel <= 1.1 && scene->sunThetaVel >= -1.1 && scene->sunTheta > 0.0 && scene->sunTheta < 1.65) {
//            printf("sunTheta %f  sunThetaVel %f  delta %f\n", scene->sunTheta, scene->sunThetaVel, delta);
            scene->sunTheta += scene->sunThetaVel * delta/10.0;
        }
        
        // remember time for next update
        updateTime = now;
        
    }
}
