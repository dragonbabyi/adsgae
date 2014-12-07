//
//  SunsetDemo.cpp
//  Sunset
//
//  Created by Yuping on 4/20/14.
//  Copyright (c) 2014 Yuping. All rights reserved.
//


#include "AppContext.h"
#include "Scene.h"
#include "Input.h"
#include "Skydome.h"
#include "Ocean.h"
#include "Program.h"

// using core modern OpenGL
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <stdio.h>

///////
// Clean up any context data
AppContext::~AppContext()
{
    // if any are NULL, deleting a NULL pointer is OK
    delete scene;
    delete input;
    delete ocean;
    delete sky;
}

///////
// GLFW callbacks must use extern "C"
extern "C" {
    
    //
    // called for GLFW error
    //
    void winError(int error, const char *description)
    {
        fprintf(stderr, "GLFW error: %s\n", description);
    }
    
    //
    // called whenever the window size changes
    //
    void reshape(GLFWwindow *win, int width, int height)
    {
        AppContext *appctx = (AppContext*)glfwGetWindowUserPointer(win);
        
        appctx->scene->viewport(win);
        appctx->input->redraw = true;
    }
    
    //
    // called when mouse button is pressed
    //
    void mousePress(GLFWwindow *win, int button, int action, int mods)
    {
        AppContext *appctx = (AppContext*)glfwGetWindowUserPointer(win);
        
        appctx->input->mousePress(win, button, action);
    }
    
    //
    // called when mouse is moved
    //
    void mouseMove(GLFWwindow *win, double x, double y)
    {
        AppContext *appctx = (AppContext*)glfwGetWindowUserPointer(win);
        
        appctx->input->mouseMove(win, appctx->scene, x,y);
    }
    
    //
    // called on any keypress
    //
    void keyPress(GLFWwindow *win, int key, int scancode, int action, int mods)
    {
        AppContext *appctx = (AppContext*)glfwGetWindowUserPointer(win);
        
        if (action == GLFW_PRESS)
            appctx->input->keyPress(win, key, appctx->scene, appctx->sky, appctx->ocean);
        else if (action == GLFW_RELEASE)
            appctx->input->keyRelease(win, key, appctx->scene);
    }
}

// initialize GLFW - windows and interaction
GLFWwindow *initGLFW(AppContext *appctx)
{
    // set error callback before init
    glfwSetErrorCallback(winError);
    if (! glfwInit())
        return 0;
    
    // ask for a window with dimensions 843 x 480 (HD 480p)
    // using core OpenGL
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);
    GLFWwindow *win = glfwCreateWindow(820, 600, "OpenGL Demo", 0, 0);
    if (! win) {
        glfwTerminate();
        return 0;
    }
    glfwMakeContextCurrent(win);
    
    // use GLEW on windows to access modern OpenGL functions
    glewExperimental = true;
    glewInit();

    printf("1# %d​\n", glGetError());

	// store context pointer to access application data
    glfwSetWindowUserPointer(win, appctx);
    
    // set callback functions to be called by GLFW
    glfwSetFramebufferSizeCallback(win, reshape);
    glfwSetKeyCallback(win, keyPress);
    glfwSetMouseButtonCallback(win, mousePress);
    glfwSetCursorPosCallback(win, mouseMove);
    
    // set OpenGL state
    glEnable(GL_DEPTH_TEST);      // tell OpenGL to handle overlapping surfaces
    return win;
    
}

int main(int argc, char *argv[])
{
    // collected data about application for use in callbacks
    AppContext appctx;
  
    GLFWwindow *win = initGLFW(&appctx);
    if (! win) return 1;
   
    // initialize context (after GLFW)
    appctx.scene = new Scene(win);
    appctx.input = new Input;
    appctx.ocean = new Ocean(win);
    appctx.sky = new Skydome();
    
    //all program should be initialized
    while (!glfwWindowShouldClose(win)) {
        
        // check for continuous key updates to view
        appctx.input->keyUpdate(appctx.scene);
        appctx.scene->update();        //scene should update every frame
        
        glClearColor(1.0,0.8,0.5,1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glDisable(GL_DEPTH_TEST);
        
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        
        //////////////////////////////////////
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
        
        appctx.ocean->draw(win);
        appctx.sky->draw(win, appctx.scene->sunTheta);
        
        // show what we drew
        glfwSwapBuffers(win);
 
        glfwPollEvents();

    }

    
    glfwDestroyWindow(win);
    glfwTerminate();
    
    return 0;
    
}