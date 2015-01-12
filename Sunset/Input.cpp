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

// using core modern OpenGL
#include <GL/glew.h>
#include <GLFW/glfw3.h>

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
        
//      scene->sunPhi += dx / 400.0;
//		scene->sunTheta += dy / 400.0;
        // set limit: sunPhi [0, 2*M_PI], sunTheta [-M_PI/2, M_PI/2]  // * M_PI / 180.0
        // don't allow move across the zenith, allow move below the horizon
        if (scene->sunTheta + dy/400.0 > 0.01 && scene->sunTheta + dy/400.0 < 1.6) {
            scene->sunTheta += dy / 400.0;
        }
        printf("%f\n", scene->sunTheta);
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

        //camera position
        case 'Z':
            scene->camera.vely = fmax(-1.0f, scene->camera.vely - 1.0f);
            break;
        case 'S':
            scene->camera.vely = fmin(1.0f, scene->camera.vely + 1.0f);
            break;
        case 'X':
            scene->camera.velx = fmax(-1.0f, scene->camera.velx - 1.0f);
            break;
        case 'D':
            scene->camera.velx = fmin(1.0f, scene->camera.velx + 1.0f);
            break;
        case 'A':
            scene->camera.velz = fmax(-1.0f, scene->camera.velz - 1.0f);
            break;
        case 'W':
            scene->camera.velz = fmin(1.0f, scene->camera.velz + 1.0f);
            break;
            
        case 'B':
            scene->sunThetaVel = 1;
            // scene->sunTheta += 0.01;
            break;
        case 'N':
            scene->sunThetaVel = -1;
            // scene->sunTheta -= 0.01;
            break;
            
        //lookat direction
        case GLFW_KEY_UP:
            //scene->camera.theta = fmin(scene->camera.theta + 1.0f, 90.0f - 0.001f);
            scene->camera.velt = fmin(1.0f, scene->camera.velt + 1.0f);
            break;
        case GLFW_KEY_DOWN:
            //scene->camera.theta = fmax(scene->camera.theta - 1.0f, -45.0f);
            scene->camera.velt = fmax(-1.0f, scene->camera.velt - 1.0f);
            break;
        case GLFW_KEY_LEFT:
            //scene->camera.phi -= 5.0;
            scene->camera.velp = fmax(-5.0f, scene->camera.velp - 5.0f);
            break;
        case GLFW_KEY_RIGHT:
            //scene->camera.phi += 5.0;
            scene->camera.velp = fmin(5.0f, scene->camera.velp + 5.0f);
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
        case 'Z':
            scene->camera.vely = fmin(1.0f, scene->camera.vely + 1.0f);
            break;
        case 'S':
            scene->camera.vely = fmax(-1.0f, scene->camera.vely - 1.0f);
            break;
        case 'X':
            scene->camera.velx = fmin(1.0f, scene->camera.velx + 1.0f);
            break;
        case 'D':
            scene->camera.velx = fmax(-1.0f, scene->camera.velx - 1.0f);
            break;
        case 'A':
            scene->camera.velz = fmin(1.0f, scene->camera.velz + 1.0f);
            break;
        case 'W':
            scene->camera.velz = fmax(-1.0f, scene->camera.velz - 1.0f);
            break;
            
        case 'B':
            scene->sunThetaVel = 0;
            break;
        case 'N':
            scene->sunThetaVel = 0;
            break;
            
        //lookat direction
        case GLFW_KEY_UP:
            //scene->camera.theta = fmin(scene->camera.theta + 1.0f, 90.0f - 0.001f);
            scene->camera.velt = fmax(-1.0f, scene->camera.velt - 1.0f);
            break;
        case GLFW_KEY_DOWN:
            //scene->camera.theta = fmax(scene->camera.theta - 1.0f, -45.0f);
            scene->camera.velt = fmin(1.0f, scene->camera.velt + 1.0f);
            break;
        case GLFW_KEY_LEFT:
            //scene->camera.phi -= 5.0;
            scene->camera.velp = fmin(5.0f, scene->camera.velp + 5.0f);
            break;
        case GLFW_KEY_RIGHT:
            //scene->camera.phi += 5.0;
            scene->camera.velp = fmax(-5.0f, scene->camera.velp - 5.0f);
            break;

    }
}

//
// update view if necessary based on a/d keys
//
void Input::keyUpdate(Scene *scene)
{
    if (scene->camera.velx != 0 || scene->camera.vely != 0 || scene->camera.velz != 0 || scene->camera.velt != 0|| scene->camera.velp != 0 || scene->sunThetaVel != 0) {
        double now = glfwGetTime();
        double delta = (now - updateTime);
        
        // update pan based on time elapsed since last update
        // ensures uniform rate of change
        
        scene->camera.x += scene->camera.vel * delta * scene->camera.velx * fmax(scene->camera.z*0.5f,1.0f);
        scene->camera.y += scene->camera.vel * delta * scene->camera.vely * fmax(scene->camera.z*0.5f,1.0f);
        scene->camera.z += scene->camera.vel * delta * scene->camera.velz * fmax(scene->camera.z*0.5f,1.0f);
        
        scene->camera.theta += scene->camera.vel * delta * scene->camera.velt * 5.0;
        scene->camera.phi += scene->camera.vel * delta * scene->camera.velp * 5.0;
        
        //        printf("%f\n", scene->camera.theta);
        if(scene->camera.theta  < -27.f)
            scene->camera.theta = -27.f;
        if(scene->camera.theta  > 180.f)
            scene->camera.theta = 180.f;
        
        if(scene->camera.z  < 0.2f)
            scene->camera.z = 0.2f;
        
        if (scene->sunThetaVel != 0 && scene->sunTheta > 0.0 && scene->sunTheta < 1.65) {
            scene->sunTheta += scene->sunThetaVel * delta/1000.0;
//            printf("%f\n", scene->sunTheta);
        }
        
        // remember time for next update
        updateTime = now;
        
    }
}
