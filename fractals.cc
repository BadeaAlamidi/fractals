#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592654

#define MAX_LEVEL 5
#define SIGMA 2
//#define SIGMA 5

// stolen from http://c-faq.com/lib/gaussian.html solution 2
double gaussrand()
{
	static double U, V;
	static int phase = 0;
	double Z;

	if(phase == 0) {
		U = (rand() + 1.) / (RAND_MAX + 2.);
		V = rand() / (RAND_MAX + 1.);
		Z = sqrt(-2 * log(U)) * sin(2 * PI * V);
	} else
		Z = sqrt(-2 * log(U)) * cos(2 * PI * V);

	phase = 1 - phase;

	return Z;
}

float f3(float delta, float x0, float x1, float x2) {
	return (x0 + x1 + x2) / 3 + delta * gaussrand();
}

float f4(float delta, float x0, float x1, float x2, float x3) {
	return (x0 + x1 + x2 + x3) / 4 + delta * gaussrand();
}

int main(int argc, char** argv){
	if (argc < 2){
		std::cerr << "Usage: ./fractals seed" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	int seed = atoi(argv[1]);

	srand(seed);
    //std::cout << gaussrand() << std::endl;
	int i, Stage;
	float delta; // "standard deviation for current level"
	int x, y, y0, D, d; // "indexing variables"

	// unused variables in psuedo code
	(void) i; (void) y0;

	// BEGIN:
	const int N = pow(2,MAX_LEVEL);
	float X[33][33]; // N + 1...
	float Colors[33][33][3]; // rgb values pertaining to each vertex in array X
	float Normals[32][32][3]; // per-vertex normal for each vertex in terrain that is surrounded by other vertices from above and to the right
	// "set the initial random corners"
	delta = SIGMA;
	X[0][0] = delta * gaussrand();
	X[0][N] = delta * gaussrand();
	X[N][0] = delta * gaussrand();
	X[N][N] = delta * gaussrand();
	// float hausdorff = 3 - (log(N)/log(Stage));
	float hausdorff = 0.84;

	D = N;
	d = N/2;

	for (Stage = 1; Stage < MAX_LEVEL + 1; ++Stage){
		// "going from grid type I to type II"
		delta = delta * pow(0.5, 0.5 * hausdorff);
		
		// "interpolate and offset mid points"
			// it is ambiguous whether this is supposed to run on the 1st iteration, hence the +1 in the stop condition
		//std::cout << "1\n";
		for (x = d; x < N - d + 1; x+=D){ 
			for (y = d; y < N - d + 1; y+=D){ 
				X[x][y] =f4(delta, X[x+d][y+d],X[x+d][y-d],X[x-d][y+d], X[x-d][y-d]);
			}
		}
		// "displace existing points"
		//std::cout << "2\n";
		for (x = 0; x < N+1; x+=D){
			for (y = 0; y < N+1; y+=D){
				X[x][y] = X[x][y] + delta * gaussrand();
			}
		}
		// "going from grid type II to type I"
		delta = delta * pow(0.5, 0.5 * hausdorff);
		// "interpolate and offset mid points at boundary"
		//std::cout << "3\n";
		for (x = d; x < N - d + 1; x+=D){
			X[x][0] = f3 (delta , X[x+d][0],X[x-d][0],X[x][d]); 
			X[x][N] = f3 (delta , X[x+d][N],X[x-d][N],X[x][N-d]); 
			X[0][x] = f3 (delta , X[0][x+d],X[0][x-d],X[d][x]);
			X[N][x] = f3 (delta , X[N][x+d],X[N][x-d],X[N-d][x]);
		}
		// "interpolate and offset mid points in interior"
		//std::cout << "4\n";
		for (x = d; x < N - d + 1; x+=D){ 
			for (y = D; y < N - d + 1; y+=D){ 
				X[x][y] =f4(delta, X[x][y+d], X[x][y-d],X[x+d][y],X[x-d][y]);
			}
		}
		//std::cout << "5\n";
		for (x = D; x < N - d + 1; x+=D){ 
			for (y = d; y < N - d + 1; y+=D){ 
				X[x][y] =f4(delta, X[x][y+d], X[x][y-d],X[x+d][y],X[x-d][y]);
			}
		}

		// "displace existing points"
		//std::cout << "6\n";
		for (x =0; x < N+1; x+=D){
			for (y =0; y < N+1; y+=D){
				X[x][y] = X[x][y] + delta *gaussrand();
			}
		}
		//std::cout << "7\n";
		for (x = d; x < N - d + 1; x+=D){
			for (y = d; y < N - d + 1; y +=D){
				X[x][y] = X[x][y] + delta * gaussrand();
			}
		}
		// " prepare for next level"
		D = D/2;
		d = d/2;
		//std::cout << Stage << std::endl;
	}
	//std::cout << std::fixed;
	//std::cout << std::setprecision(3);
	// for (int i = 0; i < 33; ++i)
	// {
	// 	for (int j = 0; j < 33; ++j)
	// 		//std::cout << X[i][j] << ' ';
	// 	//std::cout << std::endl;
	// }
	
	// colorization for standard deviation of 2:
	for (size_t i=0; i < 33; ++i)
		for (size_t j=0; j < 33; ++j){
			if (X[i][j] <= -4 ){
				Colors[i][j][0] = 0; //R
				Colors[i][j][1] = 0; //G
				Colors[i][j][2] = 1; //B
			}
			else if(X[i][j] <= -2){
				Colors[i][j][0] = 0.3; //R
				Colors[i][j][1] = 0.3; //G
				Colors[i][j][2] = 1; //B
			}
			else if (X[i][j] <= -1){
				Colors[i][j][0] = 0.6; //R
				Colors[i][j][1] = 0.6; //G
				Colors[i][j][2] = 1; //B
			}
			else if (X[i][j] <= 0){
				Colors[i][j][0] = 0.75; //R
				Colors[i][j][1] = 0.75; //G
				Colors[i][j][2] = 1; //B
			}
			else if (X[i][j] < 1){ //sand
				Colors[i][j][0] = 0.96; //R
				Colors[i][j][1] = 0.93; //G
				Colors[i][j][2] = 0.81; //B
			}
			else{ // green to white linear interpolation:
				float rb_color;
				if (X[i][j]>6 ) rb_color = 6; //limit coloration to cover the positive 0.1% of binomial distribution (SIGMA * 3)
				else rb_color = X[i][j];
				rb_color = (X[i][j] - 1) / 5; // turning Z value to a fraction to be used in linear interpolation
				Colors[i][j][0] = rb_color; // R
				Colors[i][j][1] = 1; // G
				Colors[i][j][2] = rb_color; // B
			}
		}

	// normal calculation:
	for (int i = 0; i < 32; ++i)
		for (int j = 0; j < 32; ++j){
				//I:0 J:1 K1
				//I:1 J:0 K2
				//K = 0 - 1
				//J = +K1 - 0
				//I = +K2 - 0
				float K1 = X[i][j] - X[i][j+1];
				float K2 = X[i][j] - X[i+1][j];
				
				Normals[i][j][0]=K2;
				Normals[i][j][1]=K1;
				Normals[i][j][2]=-1;
		}

	unsigned int vertex_count = 0;
	std::cout << "Display \"display\" \"Screen\" \"rgbdouble\"" << std::endl;
	std::cout << "Format 1280 960" << std::endl;
	std::cout << "CameraEye   0.0 0.0 14.0" << std::endl;
	std::cout << "CameraAt    16.0 16.0 3.0" << std::endl;
	std::cout << "CameraUp   0 0 1" << std::endl;
	std::cout << "Background 0 0.3 0.7" << std::endl;
	std::cout << "WorldBegin" << std::endl;
	//std::cout << "PointLight 16 16 20 1.0 1.0 1.0  100" << std::endl;
	std::cout << "FarLight 16 16 20  1.0  1.0  1.0  0.1" << std::endl;
	std::cout << "Ka 0.3" << std::endl;
	std::cout << "Kd 0.7" << std::endl;
	std::cout << "Color 1 0 0 " << std::endl;
	std::cout << "PolySet \"PCN\"" << std::endl;
	std::cout << pow(N,2) * 4 << std::endl << pow(N,2)*2 << std::endl;

	for (size_t i = 0; i < 32; ++i) // size is hardcoded as reminder that H is also hardcoded in main logic.
							 // if you want to change this, be sure that H and the size of the array N are also changed and not hardcoded
	{
		for (size_t j = 0; j < 32; ++j){
			std::cout << i << ' ' << j << ' ' << X[i][j] 	     << ' ' << Colors[i][j][0] << ' ' << Colors[i][j][1] << ' ' << Colors[i][j][2] 		 	   << ' ' << Normals[i][j][0] <<  ' ' << Normals[i][j][1] <<  ' ' << Normals[i][j][2] << std::endl;
			std::cout << i << ' ' << j+1 << ' ' << X[i][j+1] 	 << ' ' << Colors[i][j+1][0] << ' ' << Colors[i][j+1][1] << ' ' << Colors[i][j+1][2] 	   << ' ' << Normals[i][j+1][0] <<  ' ' << Normals[i][j+1][1] <<  ' ' << Normals[i][j+1][2] << std::endl;
			std::cout << i+1 << ' ' << j+1 << ' ' << X[i+1][j+1] << ' ' << Colors[i+1][j+1][0] << ' ' << Colors[i+1][j+1][1] << ' ' << Colors[i+1][j+1][2] << ' ' << Normals[i+1][j+1][0] <<  ' ' << Normals[i+1][j+1][1] <<  ' ' << Normals[i+1][j+1][2] << std::endl;
			std::cout << i+1 << ' ' << j << ' ' << X[i+1][j] 	 << ' ' << Colors[i+1][j][0] << ' ' << Colors[i+1][j][1] << ' ' << Colors[i+1][j][2] 	   << ' ' << Normals[i+1][j][0] <<  ' ' << Normals[i+1][j][1] <<  ' ' << Normals[i+1][j][2] << std::endl;
			vertex_count +=4;
		}
	}
	// face description for rd_view scene
//	for (size_t i = 0; i < vertex_count; ++i){
//		std::cout << i << ' ';
//		if ((i+1)%4 == 0)
//		std::cout << -1 << std::endl;
//	}
	for (size_t i = 0; i < vertex_count; i+=4)
	{
		std::cout << i << ' ' << i+1 << ' ' << i+3 << ' ' << -1 << std::endl;
		std::cout << i+1 << ' ' << i+2 << ' ' << i+3 << ' ' << -1 << std::endl;
	}
	std::cout << "WorldEnd" << std::endl;
}