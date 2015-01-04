//
//  AppContext.h
//  Sunset
//
//  Created by Yuping on 5/3/14.
//  Copyright (c) 2014 Yuping. All rights reserved.
//

#ifndef Sunset_AppContext_h
#define Sunset_AppContext_h

struct AppContext {
    class Scene *scene;
    class Input *input;
    class Ocean *ocean;
    class Skydome *sky;
    
    // uniform matrix block indices
    enum { MATRIX_UNIFORMS, TEXTURE_SKY };
    
    AppContext() : scene(0), input(0), ocean(0), sky(0) {}
    
    // clean up any context data
    ~AppContext();
};

#endif
