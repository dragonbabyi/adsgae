 
#ifdef _VERTEX_
in vec4 position;
out vec2 uv;	// texcoords

void main() {
    uv = position.zw;
    gl_Position = vec4(position.xy, 0.0, 1.0);
}

#endif

#ifdef _FRAGMENT_

#define LAYER_JACOBIAN_XX 	5.0
#define LAYER_JACOBIAN_YY	6.0
#define LAYER_JACOBIAN_XY	7.0

uniform sampler2DArray fftWavesSampler;
uniform vec4 choppy;

in vec2 uv;
out vec4 FragColor[4];

void main() {

	// fftWavesSampler has 4 height values
	// heights.r : Lx, Lz = x0
	// heights.b : Lx, Lz = x1
	// heights.b : Lx, Lz = x2
	// heights.a : Lx, Lz = x3
	// with x0 < x1 < x2 < x3
	vec4 heights = texture(fftWavesSampler, vec3(uv, 0.0));

	FragColor[0] = vec4(heights.x, heights.x * heights.x, heights.y, heights.y * heights.y);
	FragColor[1] = vec4(heights.z, heights.z * heights.z, heights.w, heights.w * heights.w);

	// store Jacobian coeff value and variance
	vec4 Jxx = choppy*texture(fftWavesSampler, vec3(uv, LAYER_JACOBIAN_XX));
	vec4 Jyy = choppy*texture(fftWavesSampler, vec3(uv, LAYER_JACOBIAN_YY));
	vec4 Jxy = choppy*choppy*texture(fftWavesSampler, vec3(uv, LAYER_JACOBIAN_XY));

	// Store partial jacobians
	vec4 res = 0.25 + Jxx + Jyy + choppy*Jxx*Jyy - Jxy*Jxy;
	vec4 res2 = res*res;
	FragColor[2] = vec4(res.x, res2.x, res.y, res2.y);
	FragColor[3] = vec4(res.z, res2.z, res.w, res2.w);

}

#endif
