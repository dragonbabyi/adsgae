#ifdef _VERTEX_
in vec4 position;
out vec2 uvIn;

void main() {
    uvIn = position.zw;
    gl_Position = vec4(position.xy, 0.0, 1.0);
}

#endif

#ifdef _GEOMETRY_

layout(triangles) in;
layout(triangle_strip, max_vertices = 24) out;

uniform int nLayers;
uniform int sLayer;

in vec2 uvIn[];
out vec2 uv;

void main() {
    for (int i = sLayer; i < sLayer + nLayers; ++i) {
        gl_Layer = i;
        gl_PrimitiveID = i;
        gl_Position = gl_in[0].gl_Position;
        uv = uvIn[0];
        EmitVertex();
        gl_Position = gl_in[1].gl_Position;
        uv = uvIn[1];
        EmitVertex();
        gl_Position = gl_in[2].gl_Position;
        uv = uvIn[2];
        EmitVertex();
        EndPrimitive();
    }
}

#endif

#ifdef _FRAGMENT_ 

uniform sampler2D butterflySampler;
uniform sampler2DArray imgSampler; // 2 complex inputs (= 4 values) per layer

uniform float pass;

in vec2 uv;
out vec4 FragColor;

// performs two FFTs on two inputs packed in a single texture
// returns two results packed in a single vec4
vec4 fft2(int layer, vec2 i, vec2 w) {
    vec4 input1 = textureLod(imgSampler, vec3(uv.x, i.x, layer), 0.0);
    vec4 input2 = textureLod(imgSampler, vec3(uv.x, i.y, layer), 0.0);
    float res1x = w.x * input2.x - w.y * input2.y;
    float res1y = w.y * input2.x + w.x * input2.y;
    float res2x = w.x * input2.z - w.y * input2.w;
    float res2y = w.y * input2.z + w.x * input2.w;
    return input1 + vec4(res1x, res1y, res2x, res2y);
}

void main() {
    vec4 data = textureLod(butterflySampler, vec2(uv.y, pass), 0.0);
    vec2 i = data.xy;
    vec2 w = data.zw;

    FragColor = fft2(gl_PrimitiveID, i, w);
}

#endif
