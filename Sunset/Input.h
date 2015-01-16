//
//  Input.h
//  Sunset
//
//  Created by Yuping on 5/4/14.
//  Copyright (c) 2014 Yuping. All rights reserved.
//

#ifndef __Sunset__Input__
#define __Sunset__Input__


// using core modern OpenGL
#include <GL/glew.h>
#include <GLFW/glfw3.h>


struct GLFWwindow;
class Scene;
class Ocean;
class Skydome;
class Sun;

class Input {
    // private data
private:
    int button, oldButton;      // which mouse button was pressed?
    double oldX, oldY;          // location of mouse at last event
    double updateTime;          // time (in seconds) of last update
    
    // public data
public:
    bool redraw;                // true if we need to redraw
    bool redrawsky;         //only when sun changes
    
    // public methods
public:
    // initialize
    Input() : button(-1), oldButton(-1), oldX(0), oldY(0), redraw(true), redrawsky(false) {
        updateTime = glfwGetTime();
    }
    
    // handle mouse press / release
    void mousePress(GLFWwindow *win, int button, int action);
    
    // handle mouse motion
    void mouseMove(GLFWwindow *win, Scene *scene, double x, double y);
    
    // handle key press
    void keyPress(GLFWwindow *win, int key, Scene *,  Skydome *, Ocean *);
    
    // handle key release
    void keyRelease(GLFWwindow *win, int key, Scene *);
    
    // update view (if necessary) based on key input
    void keyUpdate(Scene *scene);
};

#endif /* defined(__Sunset__Input__) */
