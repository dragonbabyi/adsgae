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
//    data = NULL;
    skydomeShader[0] = NULL;
    
    glGenBuffers(NUM_BUFFERS, bufferIDs);
    glGenVertexArrays(NUM_VARRAYS, varrayIDs);    
 
    //default
    solarElevation = M_PI/2.0;
    turbidity = 6.7;
    albedo = 0.03;
    HosekSkyModel_Configuration();
    
    //////////////////////  draw  hemishpere  //////////////////////////
    
    std::vector<GLfloat> vert;
    std::vector<GLint> indices;
    // create meshes
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
    
    generateTexture();
    
    updateShader();
}


Skydome::~Skydome() {
    
//    glDeleteProgram(skydomeShader);
    glDeleteBuffers(NUM_BUFFERS, bufferIDs);
//    glDeleteTextures(1, &DatasetTexture);
}

void Skydome::loadProgram() {
    char* files[1];
    files[0] = "Shader/skydome.glsl";
//    if (skydomeShader[0] != NULL) {
//        delete skydomeShader[0];
//        skydomeShader[0] = NULL;
//    }
    skydomeShader[0] = new Program(1, files);

}

//Todo:  render to texture
void Skydome::generateTexture() {
 
 
}

void Skydome::updateShader()
{
    loadProgram();

    glUseProgram(skydomeShader[0]->program);
    
    // enable vertex arrays
    glBindVertexArray(varrayIDs[VARRAY]);
    positionAttrib = glGetAttribLocation(skydomeShader[0]->program, "vPosition");
    glEnableVertexAttribArray(positionAttrib);
    glBindBuffer(GL_ARRAY_BUFFER, bufferIDs[POSITION_BUFFER]);
    glVertexAttribPointer(positionAttrib, 3, GL_FLOAT, GL_FALSE, 0, 0);
 
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
    
    // turn off everything we enabled
//    glBindTexture(GL_TEXTURE_1D, 0);  //0
    glBindVertexArray(0);
    glUseProgram(0);

}

void Skydome::draw( GLFWwindow *win, float theta) {
    //recompute the configuration if the sunTheta changes
    if( abs(theta + solarElevation - M_PI/2.0) < 1e-8) {
        solarElevation = M_PI/2.0 - theta;
        HosekSkyModel_Configuration();
    }
    
    //
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glDrawBuffer(GL_BACK);
    
    // get window dimensions
    int width, height;
    glfwGetFramebufferSize(win, &width, &height);
    glViewport(0, 0, width, height);
	   
    glEnable(GL_DEPTH_TEST);
//	glEnable(GL_CULL_FACE);
//	glCullFace(GL_BACK);
    
    glUseProgram(skydomeShader[0]->program);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    
//    ///DatasetTexture
//    glActiveTexture(GL_TEXTURE0 + DatasetTexture);
//    glBindTexture(GL_TEXTURE_1D, DatasetTexture);
    
    /////////// debug ///////////////////
//    GLfloat Pixels[64*3];
//    glGetTexImage(GL_TEXTURE_1D, 0, GL_RGB, GL_FLOAT, &Pixels);
    
    glBindVertexArray(varrayIDs[VARRAY]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bufferIDs[INDEX_BUFFER]);
//    glDrawElements(GL_TRIANGLES, iFactor*h*6, GL_UNSIGNED_INT, 0);
    glDrawElements(GL_TRIANGLES, 97200, GL_UNSIGNED_INT, 0);
    
//    glDisable(GL_CULL_FACE);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
//    glBindTexture(GL_TEXTURE_1D, 0);
    glBindVertexArray(0);
    glUseProgram(0);

}


void Skydome::HosekSkyModel_Configuration() {
    for( unsigned int channel = 0; channel < 3; ++channel )
    {
        ArHosekSkyModel_CookConfiguration(datasetsXYZ[channel], HosekConfig[channel], turbidity, albedo, solarElevation);
        ArHosekSkyModel_CookRadianceConfiguration(datasetsXYZRad[channel], HosekRadiances[channel], turbidity, albedo, solarElevation );
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


