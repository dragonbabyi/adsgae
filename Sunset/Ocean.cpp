//
//  Ocean.cpp
//  Sunset
//
//  Created by Yuping on 5/4/14.
//  Copyright (c) 2014 Yuping. All rights reserved.
//

#include "Ocean.h"
#include "AppContext.h"

// using core modern OpenGL
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <math.h>

// constant
unsigned int skyTexSize = 1024;
////////////////////////
float gridSize = 4.0f;
float hdrExposure = 1.05;
bool grid = false; //false
bool animate = true;
bool seaContrib = true;
bool sunContrib = true;
bool skyContrib = true;
bool foamContrib = false;
bool manualFilter = false;
bool show_spectrum = false;
float show_spectrum_zoom = 1.0;
bool show_spectrum_linear = false;
bool normals = false;
bool choppy = false; //true
float choppy_factor0 = 2.3f;	// Control Choppiness
float choppy_factor1 = 2.1f;	// Control Choppiness
float choppy_factor2 = 1.3f;	// Control Choppiness
float choppy_factor3 = 0.9f;	// Control Choppiness

// WAVES SPECTRUM
const int N_SLOPE_VARIANCE = 4; // size of the 3d texture containing precomputed filtered slope variances
float GRID1_SIZE = 893.0; // size in meters (i.e. in spatial domain) of the first grid
float GRID2_SIZE = 101.0; // size in meters (i.e. in spatial domain) of the second 
float GRID3_SIZE = 21.0; //51 // size in meters (i.e. in spatial domain) of the third grid
float GRID4_SIZE = 11.0; // size in meters (i.e. in spatial domain) of the fourth grid
float WIND = 12.0; //12.0 wind speed in meters per second (at 10m above surface)
float OMEGA = 2.0f; // sea state (inverse wave age)
bool propagate = true; // wave propagation?
float A = 1.0; // wave amplitude factor
const float cm = 0.23; // Eq 59
const float km = 370.0; // Eq 59
float speed = 0.6f;

// FFT WAVES
const int PASSES = 8; // number of passes needed for the FFT 6 -> 64, 7 -> 128, 8 -> 256, etc
const int FFT_SIZE = 1 << PASSES; // size of the textures storing the waves in frequency and spatial domains

// Foam
float jacobian_scale = 0.2f;


////////////////////////////////////////////////////////////////////////
float *computeButterflyLookupTexture();

Ocean::Ocean(GLFWwindow *win)
{
    vboSize = 0;
    vboVertices = 0;
    spectrum12 = NULL;
    spectrum34 = NULL;
    
    //Initialize programs
	for(GLuint i = 0; i < PROGRAM_COUNT; ++i)
		programs[i] = NULL;
    
    // Gen GL Objects
	glGenFramebuffers(FRAMEBUFFER_COUNT, framebuffers);
	glGenTextures(TEXTURE_COUNT, textures);
	glGenRenderbuffers(RENDERBUFFER_COUNT, renderbuffers);
	glGenBuffers(BUFFER_COUNT, buffers);
    // add vertex array
    glGenVertexArrays(VARRAY_COUNT, varrays);
    // get window dimensions
    glfwGetFramebufferSize(win, &window.width, &window.height);
    
    glBindRenderbuffer(GL_RENDERBUFFER, renderbuffers[RENDERBUFFER_DEPTH]);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, window.width, window.height);
	glBindRenderbuffer(GL_RENDERBUFFER, 0);
    
    float maxAnisotropy = 1.0f;
    
    ///debug texture///////////////
    glActiveTexture(GL_TEXTURE0 + TEXTURE_DEBUG);
    glBindTexture(GL_TEXTURE_2D, textures[TEXTURE_DEBUG]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, 820, 600, 0, GL_RGBA, GL_FLOAT, NULL);
    
    //////////////////
    
	glActiveTexture(GL_TEXTURE0 + TEXTURE_SPECTRUM12);
	glBindTexture(GL_TEXTURE_2D, textures[TEXTURE_SPECTRUM12]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, FFT_SIZE, FFT_SIZE, 0, GL_RGBA, GL_FLOAT, NULL);
    
 	glActiveTexture(GL_TEXTURE0 + TEXTURE_SPECTRUM34);
	glBindTexture(GL_TEXTURE_2D, textures[TEXTURE_SPECTRUM34]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, FFT_SIZE, FFT_SIZE, 0, GL_RGBA, GL_FLOAT, NULL);
    
  	glActiveTexture(GL_TEXTURE0 + TEXTURE_SLOPE_VARIANCE);
	glBindTexture(GL_TEXTURE_3D, textures[TEXTURE_SLOPE_VARIANCE]);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA16F, N_SLOPE_VARIANCE, N_SLOPE_VARIANCE, N_SLOPE_VARIANCE, 0, GL_RGBA, GL_FLOAT, NULL);
    
	glActiveTexture(GL_TEXTURE0 + TEXTURE_FFT_PING);
	glBindTexture(GL_TEXTURE_2D_ARRAY, textures[TEXTURE_FFT_PING]);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAX_ANISOTROPY_EXT, maxAnisotropy);
    glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, GL_RGBA32F, FFT_SIZE, FFT_SIZE, 10, 0, GL_RGBA, GL_FLOAT, NULL); // 8 = 1 for y + 2 for slope + 2 for D + 3 for Jacobians (Jxx, Jyy, Jxy)
    glGenerateMipmap(GL_TEXTURE_2D_ARRAY);
    
	glActiveTexture(GL_TEXTURE0 + TEXTURE_FFT_PONG);
	glBindTexture(GL_TEXTURE_2D_ARRAY, textures[TEXTURE_FFT_PONG]);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAX_ANISOTROPY_EXT, maxAnisotropy);
    glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, GL_RGBA32F, FFT_SIZE, FFT_SIZE, 10, 0, GL_RGBA, GL_FLOAT, NULL);
    glGenerateMipmap(GL_TEXTURE_2D_ARRAY);
    
	glActiveTexture(GL_TEXTURE0 + TEXTURE_BUTTERFLY);
	glBindTexture(GL_TEXTURE_2D, textures[TEXTURE_BUTTERFLY]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    //generate data and initial texture
    float *data = computeButterflyLookupTexture();
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, FFT_SIZE, PASSES, 0, GL_RGBA, GL_FLOAT, data);
	delete[] data;

    generateWavesSpectrum();   //create TEXTURE_SPECTRUM12 TEXTURE_SPECTRUM34
    
	// back to default framebuffer
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glDrawBuffer(GL_BACK);
    
	// Grid
	generateMesh();
    
    // Programs
    loadPrograms();
    
	// Slope
	computeSlopeVarianceTex();  //textures[TEXTURE_SLOPE_VARIANCE]

}

Ocean::~Ocean()
{
    for (int i=0; i<PROGRAM_COUNT; i++) {
        glDeleteProgram(programs[i]->program);
    }
    
    glDeleteTextures(TEXTURE_COUNT, textures);
    glDeleteBuffers(BUFFER_COUNT, buffers);
  
	glDeleteFramebuffers(FRAMEBUFFER_COUNT, framebuffers);
	glDeleteRenderbuffers(RENDERBUFFER_COUNT, renderbuffers);

}

void Ocean::updateShader()
{
    loadPrograms();
    computeSlopeVarianceTex();
    
}

void Ocean::draw( GLFWwindow *win, unsigned int skytex)
{
    static double now = glfwGetTime();
    static double update = glfwGetTime();
    
    update = glfwGetTime();
    float delta = update - now;
    
    static double t = 0.0;
    if(animate)
		t += delta*speed;
    
    ///////////////////////////////////////////////////////
    // solve fft
    simulateFFTWaves(t);
	glActiveTexture(GL_TEXTURE0 + TEXTURE_FFT_PING);
    glGenerateMipmap(GL_TEXTURE_2D_ARRAY);

    //////////////////////////////////////////
	// Final Rendering
    
    /////// debug /////
    glBindTexture(GL_TEXTURE_2D, textures[TEXTURE_DEBUG]);
    glBindFramebuffer(GL_FRAMEBUFFER, framebuffers[debugframebuffer]);
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, TEXTURE_DEBUG, 0);
    
    GLenum DrawBuffers[1] = {GL_COLOR_ATTACHMENT0};
    glDrawBuffers(1, DrawBuffers);
    
    //check
    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
    {
        printf("%d\n", glCheckFramebufferStatus(GL_FRAMEBUFFER));
    }
    
    
//    glBindFramebuffer(GL_FRAMEBUFFER, 0);
//	glDrawBuffer(GL_BACK);
    glViewport(0, 0, window.width, window.height);
    
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

    //update shader
	glUseProgram(programs[PROGRAM_RENDER]->program);
    glActiveTexture(GL_TEXTURE0 + skytex);
    glGenerateMipmap(GL_TEXTURE_2D);
    
    glBindTexture(GL_TEXTURE_2D, skytex);
    glUniform1i(glGetUniformLocation(programs[PROGRAM_RENDER]->program, "skySampler"), skytex);
    
    glActiveTexture(GL_TEXTURE0 + TEXTURE_SLOPE_VARIANCE);   ///// ??
    glBindTexture(GL_TEXTURE_2D, TEXTURE_SLOPE_VARIANCE);
    glUniform1i(glGetUniformLocation(programs[PROGRAM_RENDER]->program, "slopeVarianceSampler"), TEXTURE_SLOPE_VARIANCE);
    
	glBindTexture(GL_TEXTURE_2D_ARRAY, textures[TEXTURE_FFT_PING]);
	glUniform1i(glGetUniformLocation(programs[PROGRAM_RENDER]->program, "fftWavesSampler"), TEXTURE_FFT_PING);
    
    glUniformBlockBinding(programs[PROGRAM_RENDER]->program,
                          glGetUniformBlockIndex(programs[PROGRAM_RENDER]->program,"Matrices"),
                          AppContext::MATRIX_UNIFORMS);
    
    // bind vertices position in shader
    glBindVertexArray(varrays[VARRAY_MESH]);
    positionAttrib = glGetAttribLocation(programs[PROGRAM_RENDER]->program, "position");
    glEnableVertexAttribArray(positionAttrib);
    glBindBuffer(GL_ARRAY_BUFFER, buffers[BUFFER_GRID_VERTEX]);
    glVertexAttribPointer(positionAttrib, 4, GL_FLOAT, GL_FALSE, 0, 0);

    if (grid)
	{
		glPolygonMode(GL_FRONT, GL_LINE);
		glPolygonMode(GL_BACK, GL_LINE);
	}
	else
	{
		glPolygonMode(GL_FRONT, GL_FILL);
		glPolygonMode(GL_BACK, GL_FILL);
	}
    
    //draw mesh
    /////////////////////////////////////////////////////////////////////////////////////////
    glUseProgram(programs[PROGRAM_RENDER]->program);
    glBindVertexArray(varrays[VARRAY_MESH]);
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffers[BUFFER_GRID_INDEX]);
    glDrawElements(GL_TRIANGLES, vboSize, GL_UNSIGNED_INT, 0);

    now = update;
    
    
    /////// debug /////
    glBindTexture(GL_TEXTURE_2D, textures[TEXTURE_DEBUG]);
//    GLfloat *pixels = (GLfloat*) malloc(123,000 * sizeof(GLfloat) * 4);
    GLfloat *pixels = new GLfloat[600000];
    glReadPixels(0, 0, 820, 600, GL_RGBA, GL_FLOAT, pixels);
    
    for (int i=12300; i >0; i--) {
        printf("%f  %f  %f  %f \n", pixels[4*i], pixels[4*i+1], pixels[4*i+2], pixels[4*i+3]);
    }
////////////////////
   
    glDisable(GL_CULL_FACE);
 
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    glUseProgram(0);

}
///////////////////////////////////////////

float sqr(float x)
{
    return x * x;
}

float omega(float k)
{
    return sqrt(9.81 * k * (1.0 + sqr(k / km))); // Eq 24
}

// 1/kx and 1/ky in meters
float spectrum(float kx, float ky, bool omnispectrum = false)
{
    float U10 = WIND;
    float Omega = OMEGA;
    
    // phase speed
    float k = sqrt(kx * kx + ky * ky);
    float c = omega(k) / k;
    
    // spectral peak
    float kp = 9.81 * sqr(Omega / U10); // after Eq 3
    float cp = omega(kp) / kp;
    
    // friction velocity
    float z0 = 3.7e-5 * sqr(U10) / 9.81 * pow(U10 / cp, 0.9f); // Eq 66
    float u_star = 0.41 * U10 / log(10.0 / z0); // Eq 60
    
    float Lpm = exp(- 5.0 / 4.0 * sqr(kp / k)); // after Eq 3
    float gamma = Omega < 1.0 ? 1.7 : 1.7 + 6.0 * log(Omega); // after Eq 3 // log10 or log??
    float sigma = 0.08 * (1.0 + 4.0 / pow(Omega, 3.0f)); // after Eq 3
    float Gamma = exp(-1.0 / (2.0 * sqr(sigma)) * sqr(sqrt(k / kp) - 1.0));
    float Jp = pow(gamma, Gamma); // Eq 3
    float Fp = Lpm * Jp * exp(- Omega / sqrt(10.0) * (sqrt(k / kp) - 1.0)); // Eq 32
    float alphap = 0.006 * sqrt(Omega); // Eq 34
    float Bl = 0.5 * alphap * cp / c * Fp; // Eq 31
    
    float alpham = 0.01 * (u_star < cm ? 1.0 + log(u_star / cm) : 1.0 + 3.0 * log(u_star / cm)); // Eq 44
    float Fm = exp(-0.25 * sqr(k / km - 1.0)); // Eq 41
    float Bh = 0.5 * alpham * cm / c * Fm; // Eq 40
    
    Bh *= Lpm;
    
    if (omnispectrum)
    {
        return A * (Bl + Bh) / (k * sqr(k)); // Eq 30
    }
    
    float a0 = log(2.0) / 4.0;
    float ap = 4.0;
    float am = 0.13 * u_star / cm; // Eq 59
    float Delta = tanh(a0 + ap * pow(c / cp, 2.5f) + am * pow(cm / c, 2.5f)); // Eq 57
    
    float phi = atan2(ky, kx);
    
    if (propagate)
    {
        if (kx < 0.0)
        {
            return 0.0;
        }
        else
        {
            Bl *= 2.0;
            Bh *= 2.0;
        }
    }
    
	// remove waves perpendicular to wind dir
	float tweak = sqrt(fmax(kx/sqrt(kx*kx+ky*ky),0.0f));
	tweak = 1.0f;
    return A * (Bl + Bh) * (1.0 + Delta * cos(2.0 * phi)) / (2.0 * M_PI * sqr(sqr(k))) * tweak; // Eq 67
}

void Ocean::drawQuad(int programindex)
{
    vec4f *vertexdata = new vec4f[4];
    vertexdata[0] = vec4f(-1.0, -1.0, 0.0, 0.0);
    vertexdata[1] = vec4f(+1.0, -1.0, 1.0, 0.0);
    vertexdata[2] = vec4f(-1.0, +1.0, 0.0, 1.0);
    vertexdata[3] = vec4f(+1.0, +1.0, 1.0, 1.0);
    glBindBuffer(GL_ARRAY_BUFFER, buffers[BUFFER_QUAD_VERTEX]);
    glBufferData(GL_ARRAY_BUFFER, 4 * sizeof(vec4f), vertexdata, GL_STATIC_DRAW);
    
    GLuint indices[6] = {0, 1, 2, 2, 1, 3};
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffers[BUFFER_QUAD_INDEX]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 6 * sizeof(GLuint), indices, GL_STATIC_DRAW);
    
    ////////////////////////////////////////////////////////////////////////////////
    glBindVertexArray(varrays[VARRAY_QUAD]);
    posQuadAttrib[programindex] = glGetAttribLocation(programs[programindex]->program, "position");
    
    glBindBuffer(GL_ARRAY_BUFFER, buffers[BUFFER_QUAD_VERTEX]);
    glVertexAttribPointer(posQuadAttrib[programindex], 4, GL_FLOAT, GL_FALSE, 0, 0);
    
    glEnableVertexAttribArray(posQuadAttrib[programindex]);

    ////////////////////////////////////////////////////////////////////////////////
    glUseProgram(programs[programindex]->program);
    glBindVertexArray(varrays[VARRAY_QUAD]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffers[BUFFER_QUAD_INDEX]);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
    
    delete[] vertexdata;

}

// ----------------------------------------------------------------------------
// PROGRAM RELOAD
// ----------------------------------------------------------------------------

void Ocean::loadPrograms()
{
	char* files[2];
	char options[512];
    //	files[0] = "Shader/atmosphere.glsl";
    //	files[1] = "Shader/ocean.glsl";
    files[0] = "Shader/ocean.glsl";
    sprintf(options, "#define %sSEA_CONTRIB\n#define %sSUN_CONTRIB\n#define %sSKY_CONTRIB\n#define %sHARDWARE_ANISTROPIC_FILTERING\n#define %sFOAM_CONTRIB\n",
	        seaContrib ? "" : "NO_", sunContrib ? "" : "NO_", skyContrib ? "" : "NO_", manualFilter ? "NO_" : "", foamContrib ? "" : "NO_");
	
	if (programs[PROGRAM_RENDER] != NULL)
	{
		delete programs[PROGRAM_RENDER];
		programs[PROGRAM_RENDER] = NULL;
	}
    programs[PROGRAM_RENDER] = new Program(1, files, options);
	glUseProgram(programs[PROGRAM_RENDER]->program);
    
    glUniform1f(glGetUniformLocation(programs[PROGRAM_RENDER]->program, "hdrExposure"), hdrExposure);
	glUniform1f(glGetUniformLocation(programs[PROGRAM_RENDER]->program, "jacobian_scale"), jacobian_scale);
    glUniform3f(glGetUniformLocation(programs[PROGRAM_RENDER]->program, "seaColor"), 10.0/255.0, 40.0/255.0, 120.0/255.0);
	glUniform4f(glGetUniformLocation(programs[PROGRAM_RENDER]->program, "GRID_SIZES"), GRID1_SIZE, GRID2_SIZE, GRID3_SIZE, GRID4_SIZE);
	glUniform2f(glGetUniformLocation(programs[PROGRAM_RENDER]->program, "gridSize"), gridSize/float(window.width), gridSize/float(window.height));
	glUniform1f(glGetUniformLocation(programs[PROGRAM_RENDER]->program, "spectrum"), show_spectrum);
	glUniform1f(glGetUniformLocation(programs[PROGRAM_RENDER]->program, "normals"), normals);
	glUniform1f(glGetUniformLocation(programs[PROGRAM_RENDER]->program, "choppy"), choppy);
	glUniform4f(glGetUniformLocation(programs[PROGRAM_RENDER]->program, "choppy_factor"),choppy_factor0,choppy_factor1,choppy_factor2,choppy_factor3);
    
	files[0] = "Shader/init.glsl";
	if (programs[PROGRAM_INIT] != NULL)
	{
		delete programs[PROGRAM_INIT];
		programs[PROGRAM_INIT] = NULL;
	}
	programs[PROGRAM_INIT] = new Program(1, files);
	glUseProgram(programs[PROGRAM_INIT]->program);
	glUniform1i(glGetUniformLocation(programs[PROGRAM_INIT]->program, "spectrum_1_2_Sampler"), TEXTURE_SPECTRUM12);
	glUniform1i(glGetUniformLocation(programs[PROGRAM_INIT]->program, "spectrum_3_4_Sampler"), TEXTURE_SPECTRUM34);
    glUniform1f(glGetUniformLocation(programs[PROGRAM_INIT]->program, "FFT_SIZE"),FFT_SIZE);
	glUniform4f(glGetUniformLocation(programs[PROGRAM_INIT]->program, "INVERSE_GRID_SIZES"),
		        2.0 * M_PI * FFT_SIZE / GRID1_SIZE,
		        2.0 * M_PI * FFT_SIZE / GRID2_SIZE,
		        2.0 * M_PI * FFT_SIZE / GRID3_SIZE,
		        2.0 * M_PI * FFT_SIZE / GRID4_SIZE);
    
	files[0] = "Shader/variances.glsl";
	if (programs[PROGRAM_VARIANCES] != NULL)
	{
		delete programs[PROGRAM_VARIANCES];
		programs[PROGRAM_VARIANCES] = NULL;
	}
	programs[PROGRAM_VARIANCES] = new Program(1, files);
	glUseProgram(programs[PROGRAM_VARIANCES]->program);
	glUniform1f(glGetUniformLocation(programs[PROGRAM_VARIANCES]->program, "N_SLOPE_VARIANCE"), N_SLOPE_VARIANCE);
	glUniform1i(glGetUniformLocation(programs[PROGRAM_VARIANCES]->program, "spectrum_1_2_Sampler"), TEXTURE_SPECTRUM12);
	glUniform1i(glGetUniformLocation(programs[PROGRAM_VARIANCES]->program, "spectrum_3_4_Sampler"), TEXTURE_SPECTRUM34);
	glUniform1i(glGetUniformLocation(programs[PROGRAM_VARIANCES]->program, "FFT_SIZE"), FFT_SIZE);
    
    files[0] = "Shader/fftx.glsl";
	if (programs[PROGRAM_FFTX] != NULL)
	{
		delete programs[PROGRAM_FFTX];
		programs[PROGRAM_FFTX] = NULL;
	}
	programs[PROGRAM_FFTX] = new Program(1, files);
    glUseProgram(programs[PROGRAM_FFTX]->program);
//	glUniform1i(glGetUniformLocation(programs[PROGRAM_FFTX]->program, "butterflySampler"), TEXTURE_BUTTERFLY);
    glUniform1i(glGetUniformLocation(programs[PROGRAM_FFTX]->program, "nLayers"), choppy ? 8 : 3);
	glUniform1i(glGetUniformLocation(programs[PROGRAM_FFTX]->program, "sLayer"), 0);

    
	files[0] = "Shader/ffty.glsl";
	if (programs[PROGRAM_FFTY] != NULL)
	{
		delete programs[PROGRAM_FFTY];
		programs[PROGRAM_FFTY] = NULL;
	}
	programs[PROGRAM_FFTY] = new Program(1, files);
    glUseProgram(programs[PROGRAM_FFTY]->program);
//	glUniform1i(glGetUniformLocation(programs[PROGRAM_FFTY]->program, "butterflySampler"), TEXTURE_BUTTERFLY);
    glUniform1i(glGetUniformLocation(programs[PROGRAM_FFTX]->program, "nLayers"), choppy ? 8 : 3);
	glUniform1i(glGetUniformLocation(programs[PROGRAM_FFTX]->program, "sLayer"), 0);

    
	// Back to default pipeline
	glUseProgram(0);
}

// ----------------------------------------------------------------------------
// MESH GENERATION
// ----------------------------------------------------------------------------

float frandom(long *seed);

void Ocean::generateMesh()
{
    if (vboSize != 0)
    {
        glDeleteBuffers(1, &buffers[BUFFER_GRID_VERTEX]);
        glDeleteBuffers(1, &buffers[BUFFER_GRID_INDEX]);
    }
    glGenBuffers(1, &buffers[BUFFER_GRID_VERTEX]);
    glBindBuffer(GL_ARRAY_BUFFER, buffers[BUFFER_GRID_VERTEX]);
    
    //float horizon = tan(camera.theta / 180.0 * M_PI);
//    float horizon = tan(27.0 / 180.0 * M_PI);
    float horizon = tan(-1.0 / 180.0 * M_PI);
    float s = fmin(1.1f, 0.5f + horizon * 0.5f);
    
    float vmargin = 0.1;
    float hmargin = 0.1;
    
    vec4f *data = new vec4f[int(ceil(window.height * (s + vmargin) / gridSize) + 5) * int(ceil(window.width * (1.0 + 2.0 * hmargin) / gridSize) + 5)];
//    vec4f *data = new vec4f[int(ceil(window.height * 1.2 / gridSize) + 5) * int(ceil(window.width * 1.2 / gridSize) + 5)];
//debug1
//    printf("%d \n", int(ceil(window.height * (s + vmargin) / gridSize) + 5) * int(ceil(window.width * (1.0 + 2.0 * hmargin) / gridSize) + 5));
    
    int n = 0;
    int nx = 0;
    
    for (float j = window.height * s - 0.1/* - gridSize*/; j > -window.height * vmargin - gridSize; j -= gridSize)
//    for (float j = window.height * s - 0.1; j > 0; j -= gridSize)
    {
        nx = 0;
        for (float i = -window.width * hmargin; i < window.width * (1.0 + hmargin) + gridSize; i += gridSize)
        {
            data[n++] = vec4f(-1.0 + 2.0 * i / window.width, -1.0 + 2.0 * j / window.height, 0.0, 1.0);
            nx++;
        }
//        if (n>3) {
//            break;
//        }
        
    }
	vboVertices = n;
    glBufferData(GL_ARRAY_BUFFER, n * sizeof(vec4f), data, GL_STATIC_DRAW);
    delete[] data;
    
    glGenBuffers(1, &buffers[BUFFER_GRID_INDEX]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffers[BUFFER_GRID_INDEX]);
    
    vboSize = 0;
    GLuint *indices = new GLuint[6 * int(ceil(window.height * (s + vmargin) / gridSize) + 4) * int(ceil(window.width * (1.0 + 2.0 * hmargin) / gridSize) + 4)];
    
    int nj = 0;
    for (float j = window.height * s - 0.1; j > -window.height * vmargin; j -= gridSize)
    {
        int ni = 0;
        for (float i = -window.width * hmargin; i < window.width * (1.0 + hmargin); i += gridSize)
        {
            indices[vboSize++] = ni + (nj + 1) * nx;
            indices[vboSize++] = (ni + 1) + (nj + 1) * nx;
            indices[vboSize++] = (ni + 1) + nj * nx;
            
            indices[vboSize++] = ni + nj * nx;
            indices[vboSize++] = ni + (nj + 1) * nx;
            indices[vboSize++] = (ni + 1) + nj * nx;
            ni++;
        }
        nj++;
    }
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, vboSize * 4, indices, GL_STATIC_DRAW);
    delete[] indices;
    
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

// ----------------------------------------------------------------------------
// WAVES SPECTRUM GENERATION
// ----------------------------------------------------------------------------

long lrandom(long *seed)
{
    *seed = (*seed * 1103515245 + 12345) & 0x7FFFFFFF;
    return *seed;
}

float frandom(long *seed)
{
    long r = lrandom(seed) >> (31 - 24);
    return r / (float)(1 << 24);
}

inline float grandom(float mean, float stdDeviation, long *seed)
{
    float x1, x2, w, y1;
    static float y2;
    static int use_last = 0;
    
    if (use_last)
    {
        y1 = y2;
        use_last = 0;
    }
    else
    {
        do
        {
            x1 = 2.0f * frandom(seed) - 1.0f;
            x2 = 2.0f * frandom(seed) - 1.0f;
            w  = x1 * x1 + x2 * x2;
        }
        while (w >= 1.0f);
        w  = sqrt((-2.0f * log(w)) / w);
        y1 = x1 * w;
        y2 = x2 * w;
        use_last = 1;
    }
    return mean + y1 * stdDeviation;
}

void getSpectrumSample(int i, int j, float lengthScale, float kMin, float *result)
{
    static long seed = 1234;
    float dk = 2.0 * M_PI / lengthScale;
    float kx = i * dk;
    float ky = j * dk;
    if (abs(kx) < kMin && abs(ky) < kMin)
    {
        result[0] = 0.0;
        result[1] = 0.0;
    }
    else
    {
        float S = spectrum(kx, ky);
        float h = sqrt(S / 2.0) * dk;
        float phi = frandom(&seed) * 2.0 * M_PI;
        result[0] = h * cos(phi);
        result[1] = h * sin(phi);
    }
}

// generates the waves spectrum
void Ocean::generateWavesSpectrum()
{
    if (spectrum12 != NULL)
    {
        delete[] spectrum12;
        delete[] spectrum34;
    }
    spectrum12 = new float[FFT_SIZE * FFT_SIZE * 4];
    spectrum34 = new float[FFT_SIZE * FFT_SIZE * 4];
    
    for (int y = 0; y < FFT_SIZE; ++y)
    {
        for (int x = 0; x < FFT_SIZE; ++x)
        {
            int offset = 4 * (x + y * FFT_SIZE);
            int i = x >= FFT_SIZE / 2 ? x - FFT_SIZE : x;
            int j = y >= FFT_SIZE / 2 ? y - FFT_SIZE : y;
            getSpectrumSample(i, j, GRID1_SIZE, M_PI / GRID1_SIZE, spectrum12 + offset);
            getSpectrumSample(i, j, GRID2_SIZE, M_PI * FFT_SIZE / GRID1_SIZE, spectrum12 + offset + 2);
            getSpectrumSample(i, j, GRID3_SIZE, M_PI * FFT_SIZE / GRID2_SIZE, spectrum34 + offset);
            getSpectrumSample(i, j, GRID4_SIZE, M_PI * FFT_SIZE / GRID3_SIZE, spectrum34 + offset + 2);
        }
    }
    
    glActiveTexture(GL_TEXTURE0 + TEXTURE_SPECTRUM12);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, FFT_SIZE, FFT_SIZE, 0, GL_RGBA, GL_FLOAT, spectrum12);
    glActiveTexture(GL_TEXTURE0 + TEXTURE_SPECTRUM34);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, FFT_SIZE, FFT_SIZE, 0, GL_RGBA, GL_FLOAT, spectrum34);
    
}

float getSlopeVariance(float kx, float ky, float *spectrumSample)
{
    float kSquare = kx * kx + ky * ky;
    float real = spectrumSample[0];
    float img = spectrumSample[1];
    float hSquare = real * real + img * img;
    return kSquare * hSquare * 2.0;
}

// precomputes filtered slope variances in a 3d texture, based on the wave spectrum
void Ocean::computeSlopeVarianceTex()
{
    // slope variance due to all waves, by integrating over the full spectrum
    float theoreticSlopeVariance = 0.0;
    float k = 5e-3;
    while (k < 1e3)
    {
        float nextK = k * 1.001;
        theoreticSlopeVariance += k * k * spectrum(k, 0, true) * (nextK - k);
        k = nextK;
    }
    
    // slope variance due to waves, by integrating over the spectrum part
    // that is covered by the four nested grids. This can give a smaller result
    // than the theoretic total slope variance, because the higher frequencies
    // may not be covered by the four nested grid. Hence the difference between
    // the two is added as a "delta" slope variance in the "variances" shader,
    // to be sure not to lose the variance due to missing wave frequencies in
    // the four nested grids
    float totalSlopeVariance = 0.0;
    for (int y = 0; y < FFT_SIZE; ++y)
    {
        for (int x = 0; x < FFT_SIZE; ++x)
        {
            int offset = 4 * (x + y * FFT_SIZE);
            float i = 2.0 * M_PI * (x >= FFT_SIZE / 2 ? x - FFT_SIZE : x);
            float j = 2.0 * M_PI * (y >= FFT_SIZE / 2 ? y - FFT_SIZE : y);
            totalSlopeVariance += getSlopeVariance(i / GRID1_SIZE, j / GRID1_SIZE, spectrum12 + offset);
            totalSlopeVariance += getSlopeVariance(i / GRID2_SIZE, j / GRID2_SIZE, spectrum12 + offset + 2);
            totalSlopeVariance += getSlopeVariance(i / GRID3_SIZE, j / GRID3_SIZE, spectrum34 + offset);
            totalSlopeVariance += getSlopeVariance(i / GRID4_SIZE, j / GRID4_SIZE, spectrum34 + offset + 2);
        }
    }
    
    glBindFramebuffer(GL_FRAMEBUFFER, framebuffers[FRAMEBUFFER_VARIANCES]);
//    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    glViewport(0, 0, N_SLOPE_VARIANCE, N_SLOPE_VARIANCE);
    
    glUseProgram(programs[PROGRAM_VARIANCES]->program);
    glUniform1i(glGetUniformLocation(programs[PROGRAM_VARIANCES]->program, "spectrum_1_2_Sampler"), TEXTURE_SPECTRUM12);
    glUniform1i(glGetUniformLocation(programs[PROGRAM_VARIANCES]->program, "spectrum_3_4_Sampler"), TEXTURE_SPECTRUM34);
    glUniform4f(glGetUniformLocation(programs[PROGRAM_VARIANCES]->program, "GRID_SIZES"), GRID1_SIZE, GRID2_SIZE, GRID3_SIZE, GRID4_SIZE);
    glUniform1f(glGetUniformLocation(programs[PROGRAM_VARIANCES]->program, "slopeVarianceDelta"), theoreticSlopeVariance - totalSlopeVariance);
    
    for (int layer = 0; layer < N_SLOPE_VARIANCE; ++layer)
    {
        glFramebufferTexture3D(GL_FRAMEBUFFER,
		                          GL_COLOR_ATTACHMENT0,
		                          GL_TEXTURE_3D,
		                          textures[TEXTURE_SLOPE_VARIANCE],
		                          0,
		                          layer);
//        GLenum check = glCheckFramebufferStatus(GL_FRAMEBUFFER);
//        std::cout << "11#" << check << "\n";
        
		glUniform1f(glGetUniformLocation(programs[PROGRAM_VARIANCES]->program, "c"), layer);
        drawQuad(PROGRAM_VARIANCES);
        
	}
    
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}


// ----------------------------------------------------------------------------
// WAVES GENERATION AND ANIMATION (using FFT on GPU)
// ----------------------------------------------------------------------------

int bitReverse(int i, int N)
{
	int j = i;
	int M = N;
	int Sum = 0;
	int W = 1;
	M = M / 2;
	while (M != 0)
	{
		j = (i & M) > M - 1;
        
        
		Sum += j * W;
        
        
        
		W *= 2;
		M = M / 2;
	}
	return Sum;
}

void computeWeight(int N, int k, float &Wr, float &Wi)
{
    Wr = cosl(2.0 * M_PI * k / float(N));
    Wi = sinl(2.0 * M_PI * k / float(N));
}

float *computeButterflyLookupTexture()
{
	float *data = new float[FFT_SIZE * PASSES * 4];
    
	for (int i = 0; i < PASSES; i++)
	{
		int nBlocks  = (int) powf(2.0, float(PASSES - 1 - i));
		int nHInputs = (int) powf(2.0, float(i));
		for (int j = 0; j < nBlocks; j++)
		{
		    for (int k = 0; k < nHInputs; k++)
		    {
		        int i1, i2, j1, j2;
		        if (i == 0)
		        {
		            i1 = j * nHInputs * 2 + k;
		            i2 = j * nHInputs * 2 + nHInputs + k;
		            j1 = bitReverse(i1, FFT_SIZE);
		            j2 = bitReverse(i2, FFT_SIZE);
		        }
		        else
		        {
		            i1 = j * nHInputs * 2 + k;
		            i2 = j * nHInputs * 2 + nHInputs + k;
		            j1 = i1;
		            j2 = i2;
		        }
                
		        float wr, wi;
		        computeWeight(FFT_SIZE, k * nBlocks, wr, wi);
                
		        int offset1 = 4 * (i1 + i * FFT_SIZE);
		        data[offset1 + 0] = (j1 + 0.5) / FFT_SIZE;
		        data[offset1 + 1] = (j2 + 0.5) / FFT_SIZE;
		        data[offset1 + 2] = wr;
		        data[offset1 + 3] = wi;
                
		        int offset2 = 4 * (i2 + i * FFT_SIZE);
		        data[offset2 + 0] = (j1 + 0.5) / FFT_SIZE;
		        data[offset2 + 1] = (j2 + 0.5) / FFT_SIZE;
		        data[offset2 + 2] = -wr;
		        data[offset2 + 3] = -wi;
		    }
		}
	}
    
	return data;
}

void Ocean::simulateFFTWaves(float t)
{
	// init
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffers[FRAMEBUFFER_FFT0]);
    
	for (int i = 0; i < 8; ++i)
	{
		glFramebufferTextureLayer(   GL_FRAMEBUFFER,
                                     GL_COLOR_ATTACHMENT0 + i,
                                     textures[TEXTURE_FFT_PING],
                                     0,
                                     i);
	}
    
	GLenum drawBuffers[8] =
	{
		GL_COLOR_ATTACHMENT0,
		GL_COLOR_ATTACHMENT1,
		GL_COLOR_ATTACHMENT2,
		GL_COLOR_ATTACHMENT3,
		GL_COLOR_ATTACHMENT4,
		GL_COLOR_ATTACHMENT5,
		GL_COLOR_ATTACHMENT6,
		GL_COLOR_ATTACHMENT7
	};
    
	glDrawBuffers(choppy ? 8 : 3, drawBuffers);

 	glViewport(0, 0, FFT_SIZE, FFT_SIZE);
	glUseProgram(programs[PROGRAM_INIT]->program);
    
    glBindTexture(GL_TEXTURE_2D, textures[TEXTURE_SPECTRUM12]);
	glUniform1i(glGetUniformLocation(programs[PROGRAM_INIT]->program, "spectrum_1_2_Sampler"), TEXTURE_SPECTRUM12);
    
    glBindTexture(GL_TEXTURE_2D, textures[TEXTURE_SPECTRUM34]);
	glUniform1i(glGetUniformLocation(programs[PROGRAM_INIT]->program, "spectrum_3_4_Sampler"), TEXTURE_SPECTRUM34);
	glUniform1f(glGetUniformLocation(programs[PROGRAM_INIT]->program, "t"), t);

	drawQuad(PROGRAM_INIT);
    
    // fft passes
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffers[FRAMEBUFFER_FFT1]);
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, textures[TEXTURE_FFT_PING], 0);
    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, textures[TEXTURE_FFT_PONG], 0);
    
	glUseProgram(programs[PROGRAM_FFTX]->program);

	glUniform1i(glGetUniformLocation(programs[PROGRAM_FFTX]->program, "butterflySampler"), TEXTURE_BUTTERFLY);
		for (int i = 0; i < PASSES; ++i)
	{
		glUniform1f(glGetUniformLocation(programs[PROGRAM_FFTX]->program, "pass"), float(i + 0.5) / PASSES);
		if (i%2 == 0)
		{
		    glUniform1i(glGetUniformLocation(programs[PROGRAM_FFTX]->program, "imgSampler"), TEXTURE_FFT_PING);
		    glDrawBuffer(GL_COLOR_ATTACHMENT1);
		}
		else
		{
		    glUniform1i(glGetUniformLocation(programs[PROGRAM_FFTX]->program, "imgSampler"), TEXTURE_FFT_PONG);
		    glDrawBuffer(GL_COLOR_ATTACHMENT0);
		}
 
		drawQuad(PROGRAM_FFTX);
	}
    
	glUseProgram(programs[PROGRAM_FFTY]->program);
	glUniform1i(glGetUniformLocation(programs[PROGRAM_FFTY]->program, "butterflySampler"), TEXTURE_BUTTERFLY);
    
	for (int i = PASSES; i < 2 * PASSES; ++i)
	{
		glUniform1f(glGetUniformLocation(programs[PROGRAM_FFTY]->program, "pass"), float(i - PASSES + 0.5) / PASSES);
		if (i%2 == 0)
		{
		    glUniform1i(glGetUniformLocation(programs[PROGRAM_FFTY]->program, "imgSampler"), TEXTURE_FFT_PING);
		    glDrawBuffer(GL_COLOR_ATTACHMENT1);
		}
		else
		{
		    glUniform1i(glGetUniformLocation(programs[PROGRAM_FFTY]->program, "imgSampler"), TEXTURE_FFT_PONG);
		    glDrawBuffer(GL_COLOR_ATTACHMENT0);
		}
		drawQuad(PROGRAM_FFTY);
	}

	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

