//
//  Skydome.cpp
//  Sunset
//
//  Created by Yuping on 3/24/14.
//  Copyright (c) 2014 Yuping. All rights reserved.
//

#include "AppContext.h"
#include "Skydome.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <sstream>
#include <cmath>
#include <new>
#include <exception>

//data file
#include "ArHosekSkyModelData_CIEXYZ.h"

// using core modern OpenGL
#include <GL/glew.h>
#include <GLFW/glfw3.h>


void ArHosekSkyModel_CookConfiguration(double* dataset, float *HosekConfigChannel, float turbidity, float albedo, float theta);
void ArHosekSkyModel_CookRadianceConfiguration(double* dataset, float &HosekRadiancesChannel, float turbidity, float albedo, float theta);

Skydome::Skydome() {
//    Initialize program
 
    skydomeShader[0] = NULL;
    skydomeShader[1] = NULL;
    
    skyframebuffer = 0;
    glGenFramebuffers(1, &skyframebuffer);
    
    //add sky texture
    skyTexSize = 1024;
    glGenTextures(1, &skytexture);
    glActiveTexture(GL_TEXTURE0 + skytexture);
	glBindTexture(GL_TEXTURE_2D, skytexture);
    
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, skyTexSize, skyTexSize, 0, GL_RGBA, GL_FLOAT, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glGenerateMipmapEXT(GL_TEXTURE_2D);
    
    //depth
    glGenRenderbuffers(1, &depthrenderbuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, depthrenderbuffer);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, 1024, 768);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthrenderbuffer);
    
    glGenBuffers(NUM_BUFFERS, bufferIDs);
    glGenVertexArrays(NUM_VARRAYS, varrayIDs);
    
    //draw to texture
    vec4f *vertexdata = new vec4f[4];
    vertexdata[0] = vec4f(-1.0, -1.0, 0.0, 0.0);
    vertexdata[1] = vec4f(+1.0, -1.0, 1.0, 0.0);
    vertexdata[2] = vec4f(-1.0, +1.0, 0.0, 1.0);
    vertexdata[3] = vec4f(+1.0, +1.0, 1.0, 1.0);
    glBindBuffer(GL_ARRAY_BUFFER, bufferIDs[QUAD_VERTEX_BUFFER]);
    glBufferData(GL_ARRAY_BUFFER, 4 * sizeof(vec4f), vertexdata, GL_STATIC_DRAW);
    
    GLuint indicesdata[6] = {0, 1, 2, 2, 1, 3};
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufferIDs[QUAD_INDEX_BUFFER]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 6 * sizeof(GLuint), indicesdata, GL_STATIC_DRAW);
    
    //default
    solarElevation = M_PI/2.0;
    turbidity = 8.0;
    albedo = 0.03;
    init = false;
    
    //////////////////////  create meshes (hemishpere)  //////////////////////////
    
    std::vector<GLfloat> vert;
    std::vector<GLint> indices;
    
    float radius = 1.0;    // have to be the same for the sun shader part
    int iFactor = 180;
    int h = 90;
    vert.resize(iFactor*h*3);
    indices.resize(iFactor*h*6);
    
    float phi, theta;
    std::vector<GLfloat>::iterator v = vert.begin();
    for (int i = 0; i < h ; i++) {
        theta = (M_PI*i)/h;
        
        for (int j = 0; j<iFactor; j++) {
            phi = (2*M_PI*j)/iFactor;
            
            *v++ = (float)(sin(theta)*cos(phi)*radius);
            *v++ = (float)(sin(theta)*sin(phi)*radius);
            *v++ = (float)(cos(theta)*radius);
            
        }
    }
    
    // build index array linking sets of three vertices into triangles
    std::vector<GLint>::iterator k = indices.begin();
    for (int i = 0; i<h; i++) {
        for (int j = 0; j<iFactor; j++) {
            *k++ = i*iFactor + j;
            *k++ = (i+1)*iFactor + j;
            *k++ = i*iFactor + (j + 1);
            *k++ = (i+1)*iFactor + j;
            *k++ = (i+1)*iFactor + (j + 1);
            *k++ = i*iFactor + (j + 1);
            
        }
    }
    /////////
    
    glBindBuffer(GL_ARRAY_BUFFER, bufferIDs[POSITION_BUFFER]);
    glBufferData(GL_ARRAY_BUFFER, iFactor*h*3*sizeof(float), &vert.front(), GL_STATIC_DRAW);
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufferIDs[INDEX_BUFFER]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, iFactor*h*6*sizeof(int), &indices.front(), GL_STATIC_DRAW);
    
    loadProgram();
//    updateShader();
}


Skydome::~Skydome() {
//    glDeleteProgram(skydomeShader);
    glDeleteBuffers(NUM_BUFFERS, bufferIDs);
    glDeleteTextures(1, &skytexture);
}

void Skydome::loadProgram() {
    char* files[1];
//    files[0] = "Shader/skydome1.glsl";
    files[0] = "Shader/skydome.glsl";
//    if (skydomeShader[0] != NULL) {
//        delete skydomeShader[0];
//        skydomeShader[0] = NULL;
//    }
    skydomeShader[0] = new Program(1, files);
    
    // second shader
//    files[0] = "Shader/tonemapping.glsl";
    files[0] = "Shader/sun.glsl";
    
    skydomeShader[1] = new Program(1, files);

}


void Skydome::updateShader()
{
    glUseProgram(skydomeShader[0]->program);
    glBindVertexArray(varrayIDs[QUAD_VERTEX]);
    posQuadAttrib = glGetAttribLocation(skydomeShader[0]->program, "position");
    glEnableVertexAttribArray(posQuadAttrib);
    glBindBuffer(GL_ARRAY_BUFFER, bufferIDs[QUAD_VERTEX_BUFFER]);
    glVertexAttribPointer(posQuadAttrib, 4, GL_FLOAT, GL_FALSE, 0, 0);   // 4!!
    
    glUniformBlockBinding(skydomeShader[0]->program,
                          glGetUniformBlockIndex(skydomeShader[0]->program,"Matrices"),
                          AppContext::MATRIX_UNIFORMS);
    
    glUniform3fv(glGetUniformLocation(skydomeShader[0]->program, "HosekRadiances"), 1, HosekRadiances);
    glUniform3f(glGetUniformLocation(skydomeShader[0]->program, "XHosekABC"), HosekConfig[0][0], HosekConfig[0][1], HosekConfig[0][2]);
    glUniform3f(glGetUniformLocation(skydomeShader[0]->program, "XHosekDEF"), HosekConfig[0][3], HosekConfig[0][4], HosekConfig[0][5]);
    glUniform3f(glGetUniformLocation(skydomeShader[0]->program, "XHosekGHI"), HosekConfig[0][6], HosekConfig[0][7], HosekConfig[0][8]);
    glUniform3f(glGetUniformLocation(skydomeShader[0]->program, "YHosekABC"), HosekConfig[1][0], HosekConfig[1][1], HosekConfig[1][2]);
    glUniform3f(glGetUniformLocation(skydomeShader[0]->program, "YHosekDEF"), HosekConfig[1][3], HosekConfig[1][4], HosekConfig[1][5]);
    glUniform3f(glGetUniformLocation(skydomeShader[0]->program, "YHosekGHI"), HosekConfig[1][6], HosekConfig[1][7], HosekConfig[1][8]);
    glUniform3f(glGetUniformLocation(skydomeShader[0]->program, "ZHosekABC"), HosekConfig[2][0], HosekConfig[2][1], HosekConfig[2][2]);
    glUniform3f(glGetUniformLocation(skydomeShader[0]->program, "ZHosekDEF"), HosekConfig[2][3], HosekConfig[2][4], HosekConfig[2][5]);
    glUniform3f(glGetUniformLocation(skydomeShader[0]->program, "ZHosekGHI"), HosekConfig[2][6], HosekConfig[2][7], HosekConfig[2][8]);
    

    glUseProgram(skydomeShader[1]->program);
    // enable vertex arrays
    glBindVertexArray(varrayIDs[VARRAY]);
    positionAttrib = glGetAttribLocation(skydomeShader[1]->program, "vPosition");
    glEnableVertexAttribArray(positionAttrib);
    glBindBuffer(GL_ARRAY_BUFFER, bufferIDs[POSITION_BUFFER]);
    glVertexAttribPointer(positionAttrib, 3, GL_FLOAT, GL_FALSE, 0, 0);
    
    glUniformBlockBinding(skydomeShader[1]->program,
                          glGetUniformBlockIndex(skydomeShader[1]->program,"Matrices"),
                          AppContext::MATRIX_UNIFORMS);
    
    //skytexture
    glBindTexture(GL_TEXTURE_2D, skytexture);
    glUniform1i(glGetUniformLocation(skydomeShader[1]->program, "skySample2D"), skytexture);

    // turn off everything we enabled
//    glBindTexture(GL_TEXTURE_2D, 0);  //0
    glBindVertexArray(0);
    glUseProgram(0);

}

void Skydome::draw( GLFWwindow *win, float theta) {
    
    // update HosekSkyModel_Configuration parameters
    if ( !init ) {
        solarElevation = M_PI/2.0 - theta;
        HosekSkyModel_Configuration( solarElevation );
        init = true;
    }
    else if ( abs(theta + solarElevation - M_PI/2.0) < 1e-8) {   //recompute the configuration if the sunTheta changes
        solarElevation = M_PI/2.0 - theta;
        HosekSkyModel_Configuration( solarElevation );
    }
    
    updateShader();
    
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    // render to skytexture
    glEnable(GL_DEPTH_TEST);
    
    glBindTexture(GL_TEXTURE_2D, skytexture);
    
    glBindFramebuffer(GL_FRAMEBUFFER, skyframebuffer);
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, skytexture, 0);
    
    GLenum DrawBuffers[1] = {GL_COLOR_ATTACHMENT0};
    glDrawBuffers(1, DrawBuffers);
    
    //check
    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    {
        printf("%d\n", glCheckFramebufferStatus(GL_FRAMEBUFFER));
    }
    
    glViewport(0, 0, skyTexSize, skyTexSize);
    
    glUseProgram(skydomeShader[0]->program);
    
    glBindVertexArray(varrayIDs[QUAD_VERTEX]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufferIDs[QUAD_INDEX_BUFFER]);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
    
//    glBindTexture(GL_TEXTURE_2D, 0);
    glUseProgram(0);
    
    //debug
//    glGetError();
//    float pixeldata[16] = {0};
//    glReadPixels( 500, 500, 2, 2, GL_RGBA, GL_FLOAT, pixeldata);
    
    ///////////////////////////////////////////////////////////
    // final render
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glDrawBuffer(GL_BACK);
    
//    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
//    glEnable(GL_DEPTH_TEST);
//    glEnable(GL_CULL_FACE);
//    glCullFace(GL_BACK);
    
// test tone mapping    /////debugging
/*    glViewport(0, 0, skyTexSize, skyTexSize);
    
    glUseProgram(skydomeShader[1]->program);
    
    glBindVertexArray(varrayIDs[QUAD_VERTEX]);
    positionAttrib = glGetAttribLocation(skydomeShader[1]->program, "position");
    glEnableVertexAttribArray(positionAttrib);
    glBindBuffer(GL_ARRAY_BUFFER, bufferIDs[QUAD_VERTEX_BUFFER]);
    glVertexAttribPointer(positionAttrib, 4, GL_FLOAT, GL_FALSE, 0, 0);   // 4!!
    
    glUniformBlockBinding(skydomeShader[1]->program,
                          glGetUniformBlockIndex(skydomeShader[1]->program,"Matrices"),
                          AppContext::MATRIX_UNIFORMS);

    //skytexture
    glActiveTexture(GL_TEXTURE0 + skytexture);
    glGenerateMipmap(GL_TEXTURE_2D);
    
    glBindTexture(GL_TEXTURE_2D, skytexture);
    glUniform1i(glGetUniformLocation(skydomeShader[1]->program, "skySample2D"), skytexture);
    
    glBindVertexArray(varrayIDs[QUAD_VERTEX]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufferIDs[QUAD_INDEX_BUFFER]);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
*/
    
//    // get window dimensions
    int width, height;
    glfwGetFramebufferSize(win, &width, &height);
    glViewport(0, 0, width, height);

    glUseProgram(skydomeShader[1]->program);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); //GL_LINE
  
    //skytexture
    glActiveTexture(GL_TEXTURE0 + skytexture);
    glGenerateMipmap(GL_TEXTURE_2D);
    
//    glBindTexture(GL_TEXTURE_2D, skytexture);
    glUniform1i(glGetUniformLocation(skydomeShader[1]->program, "skySample2D"), skytexture);
    
    glBindVertexArray(varrayIDs[VARRAY]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufferIDs[INDEX_BUFFER]);
    glDrawElements(GL_TRIANGLES, 97200, GL_UNSIGNED_INT, 0);

//    glDisable(GL_CULL_FACE);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
//    glBindTexture(GL_TEXTURE_2D, 0);
    glBindVertexArray(0);
    glUseProgram(0);

}


void Skydome::HosekSkyModel_Configuration(float solarAngle) {
    for( unsigned int channel = 0; channel < 3; ++channel )
    {
        ArHosekSkyModel_CookConfiguration(datasetsXYZ[channel], HosekConfig[channel], turbidity, albedo, solarAngle);
        ArHosekSkyModel_CookRadianceConfiguration(datasetsXYZRad[channel], HosekRadiances[channel], turbidity, albedo, solarAngle );
    }
}


void ArHosekSkyModel_CookConfiguration(double* dataset, float *config,float turbidity,float albedo,float solar_elevation )
{
    double  * elev_matrix;
    
    int     int_turbidity = (int)turbidity;
    double  turbidity_rem = turbidity - (double)int_turbidity;
    
    solar_elevation = pow(solar_elevation / (M_PI / 2.0), (1.0 / 3.0));
    
    // alb 0 low turb
    
    elev_matrix = dataset + ( 9 * 6 * (int_turbidity-1) );
    
    
    for( unsigned int i = 0; i < 9; ++i )
    {
        //(1-t).^3* A1 + 3*(1-t).^2.*t * A2 + 3*(1-t) .* t .^ 2 * A3 + t.^3 * A4;
        // config[i] =
        // (1.0-albedo) * (1.0 - turbidity_rem)
        // *
        config[i] = (1.0-albedo) * (1.0 - turbidity_rem)
        *( pow(1.0-solar_elevation, 5.0) * elev_matrix[i]  +
                     5.0  * pow(1.0-solar_elevation, 4.0) * solar_elevation * elev_matrix[i+9] +
                     10.0*pow(1.0-solar_elevation, 3.0)*pow(solar_elevation, 2.0) * elev_matrix[i+18] +
                     10.0*pow(1.0-solar_elevation, 2.0)*pow(solar_elevation, 3.0) * elev_matrix[i+27] +
                     5.0*(1.0-solar_elevation)*pow(solar_elevation, 4.0) * elev_matrix[i+36] +
                     pow(solar_elevation, 5.0)  * elev_matrix[i+45]);
        
        // printf("%f  %f  %f  \n", solar_elevation, elev_matrix[i], config[i]);
    }
    
    // alb 1 low turb
    elev_matrix = dataset + (9*6*10 + 9*6*(int_turbidity-1));
    for(unsigned int i = 0; i < 9; ++i)
    {
        //(1-t).^3* A1 + 3*(1-t).^2.*t * A2 + 3*(1-t) .* t .^ 2 * A3 + t.^3 * A4;
        config[i] +=
        (albedo) * (1.0 - turbidity_rem)
        * ( pow(1.0-solar_elevation, 5.0) * elev_matrix[i]  +
           5.0  * pow(1.0-solar_elevation, 4.0) * solar_elevation * elev_matrix[i+9] +
           10.0*pow(1.0-solar_elevation, 3.0)*pow(solar_elevation, 2.0) * elev_matrix[i+18] +
           10.0*pow(1.0-solar_elevation, 2.0)*pow(solar_elevation, 3.0) * elev_matrix[i+27] +
           5.0*(1.0-solar_elevation)*pow(solar_elevation, 4.0) * elev_matrix[i+36] +
           pow(solar_elevation, 5.0)  * elev_matrix[i+45]);
    }
    
    if(int_turbidity == 10)
        return;
    
    // alb 0 high turb
    elev_matrix = dataset + (9*6*(int_turbidity));
    for(unsigned int i = 0; i < 9; ++i)
    {
        //(1-t).^3* A1 + 3*(1-t).^2.*t * A2 + 3*(1-t) .* t .^ 2 * A3 + t.^3 * A4;
        config[i] +=
        (1.0-albedo) * (turbidity_rem)
        * ( pow(1.0-solar_elevation, 5.0) * elev_matrix[i]  +
           5.0  * pow(1.0-solar_elevation, 4.0) * solar_elevation * elev_matrix[i+9] +
           10.0*pow(1.0-solar_elevation, 3.0)*pow(solar_elevation, 2.0) * elev_matrix[i+18] +
           10.0*pow(1.0-solar_elevation, 2.0)*pow(solar_elevation, 3.0) * elev_matrix[i+27] +
           5.0*(1.0-solar_elevation)*pow(solar_elevation, 4.0) * elev_matrix[i+36] +
           pow(solar_elevation, 5.0)  * elev_matrix[i+45]);
    }
    
    // alb 1 high turb
    elev_matrix = dataset + (9*6*10 + 9*6*(int_turbidity));
    for(unsigned int i = 0; i < 9; ++i)
    {
        //(1-t).^3* A1 + 3*(1-t).^2.*t * A2 + 3*(1-t) .* t .^ 2 * A3 + t.^3 * A4;
        config[i] +=
        (albedo) * (turbidity_rem)
        * ( pow(1.0-solar_elevation, 5.0) * elev_matrix[i]  +
           5.0  * pow(1.0-solar_elevation, 4.0) * solar_elevation * elev_matrix[i+9] +
           10.0*pow(1.0-solar_elevation, 3.0)*pow(solar_elevation, 2.0) * elev_matrix[i+18] +
           10.0*pow(1.0-solar_elevation, 2.0)*pow(solar_elevation, 3.0) * elev_matrix[i+27] +
           5.0*(1.0-solar_elevation)*pow(solar_elevation, 4.0) * elev_matrix[i+36] +
           pow(solar_elevation, 5.0)  * elev_matrix[i+45]);
    }

}

void ArHosekSkyModel_CookRadianceConfiguration(double* dataset, float &HosekRadiancesChannel,float turbidity,float albedo,float solar_elevation)
{
    double* elev_matrix;
    
    int int_turbidity = (int)turbidity;
    double turbidity_rem = turbidity - (double)int_turbidity;
    double res;
    solar_elevation = pow(solar_elevation / (M_PI / 2.0), (1.0 / 3.0));
    
    // alb 0 low turb
    elev_matrix = dataset + (6*(int_turbidity-1));
    //(1-t).^3* A1 + 3*(1-t).^2.*t * A2 + 3*(1-t) .* t .^ 2 * A3 + t.^3 * A4;
    res = (1.0-albedo) * (1.0 - turbidity_rem) *
    ( pow(1.0-solar_elevation, 5.0) * elev_matrix[0] +
     5.0*pow(1.0-solar_elevation, 4.0)*solar_elevation * elev_matrix[1] +
     10.0*pow(1.0-solar_elevation, 3.0)*pow(solar_elevation, 2.0) * elev_matrix[2] +
     10.0*pow(1.0-solar_elevation, 2.0)*pow(solar_elevation, 3.0) * elev_matrix[3] +
     5.0*(1.0-solar_elevation)*pow(solar_elevation, 4.0) * elev_matrix[4] +
     pow(solar_elevation, 5.0) * elev_matrix[5]);
    
    // alb 1 low turb
    elev_matrix = dataset + (6*10 + 6*(int_turbidity-1));
    //(1-t).^3* A1 + 3*(1-t).^2.*t * A2 + 3*(1-t) .* t .^ 2 * A3 + t.^3 * A4;
    res += (albedo) * (1.0 - turbidity_rem) *
    ( pow(1.0-solar_elevation, 5.0) * elev_matrix[0] +
     5.0*pow(1.0-solar_elevation, 4.0)*solar_elevation * elev_matrix[1] +
     10.0*pow(1.0-solar_elevation, 3.0)*pow(solar_elevation, 2.0) * elev_matrix[2] +
     10.0*pow(1.0-solar_elevation, 2.0)*pow(solar_elevation, 3.0) * elev_matrix[3] +
     5.0*(1.0-solar_elevation)*pow(solar_elevation, 4.0) * elev_matrix[4] +
     pow(solar_elevation, 5.0) * elev_matrix[5]);
    if(int_turbidity == 10) {
        HosekRadiancesChannel = res;
        return ;
    }
    
    // alb 0 high turb
    elev_matrix = dataset + (6*(int_turbidity));
    //(1-t).^3* A1 + 3*(1-t).^2.*t * A2 + 3*(1-t) .* t .^ 2 * A3 + t.^3 * A4;
    res += (1.0-albedo) * (turbidity_rem) *
    ( pow(1.0-solar_elevation, 5.0) * elev_matrix[0] +
     5.0*pow(1.0-solar_elevation, 4.0)*solar_elevation * elev_matrix[1] +
     10.0*pow(1.0-solar_elevation, 3.0)*pow(solar_elevation, 2.0) * elev_matrix[2] +
     10.0*pow(1.0-solar_elevation, 2.0)*pow(solar_elevation, 3.0) * elev_matrix[3] +
     5.0*(1.0-solar_elevation)*pow(solar_elevation, 4.0) * elev_matrix[4] +
     pow(solar_elevation, 5.0) * elev_matrix[5]);
    
    // alb 1 high turb
    elev_matrix = dataset + (6*10 + 6*(int_turbidity));
    //(1-t).^3* A1 + 3*(1-t).^2.*t * A2 + 3*(1-t) .* t .^ 2 * A3 + t.^3 * A4;
    res += (albedo) * (turbidity_rem) *
    ( pow(1.0-solar_elevation, 5.0) * elev_matrix[0] +
     5.0*pow(1.0-solar_elevation, 4.0)*solar_elevation * elev_matrix[1] +
     10.0*pow(1.0-solar_elevation, 3.0)*pow(solar_elevation, 2.0) * elev_matrix[2] +
     10.0*pow(1.0-solar_elevation, 2.0)*pow(solar_elevation, 3.0) * elev_matrix[3] +
     5.0*(1.0-solar_elevation)*pow(solar_elevation, 4.0) * elev_matrix[4] +
     pow(solar_elevation, 5.0) * elev_matrix[5]);
    
    HosekRadiancesChannel = res;
}


