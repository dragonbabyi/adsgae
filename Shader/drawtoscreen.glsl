
////////////////////////////////////////////////////////////////////
#ifdef _VERTEX_

in vec3 position;   //vertex position
in vec2 Texcoord;

out vec2 vTexcoord;

void main() {
	vTexcoord = Texcoord;
	gl_Position = vec4(position, 1.0);
}

#endif

////////////////////////////////////////////////////////////////////

#ifdef _FRAGMENT_

uniform sampler2D texturesample;

in vec2 vTexcoord;
out vec4 FragColor;

void main()
{
	FragColor = texture(texturesample, vTexcoord.st);
	FragColor = vec4(vTexcoord, 0.0, 1.0);
	//debug
    //FragColor = vec4(1.0, 0.0, 0.0, 1.0);
}


#endif