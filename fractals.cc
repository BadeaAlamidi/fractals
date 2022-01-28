#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define PI 3.141592654

#define MAX_LEVEL 5
#define SIGMA 2

// used to dereference one-dimensional array as a 2d array. the width is equal to how many elements exist in a row
inline int index2D(int column, int row, int width){
	return column * width + row;
}

// used to dereference one-dimensional arrays as a 3d array. the width is equal to how many elements exist in a row
inline int index3D(int column, int row, int plane, int column_width, int row_width){
//	return (column*column_width) + (row*row_width) + plane;
	return (row_width * column_width * plane) + (column*column_width + row);

}
// [1][2][2] 33 + 5 + 2 = 40

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
	int max_level = 5; // determines the level of recursion of main logic @ "BEGIN"
	bool normals = false; // used to conditionally print normal vectors for each vertex
	bool seed_provided = false; // used to ensure that the seed argument was provided in invocation
	
	//else if (argc >= 3){
	//	max_level = atoi(argv[2]);;;;;;;;;;;;
	//}
	int seed, opt;
	int deviation = SIGMA;
	while ((opt = getopt(argc, argv, "ns:l:d:")) != -1){
		switch(opt) {
			case 'n' : 
				normals = true;
			break;
			case 'l' : 
				max_level = atoi(optarg);
				std::cout << "#max level = "<< max_level << std::endl;
				if (max_level < 2 || max_level > 7){
					max_level = 5;
					std::cout << "unrecommended value for max_level. values must be in the range 2-7" << std::endl;
				}
			break;
			case 's' : 
				seed_provided = true;
				seed = atoi(optarg);
				std::cout << "#seed level = "<< seed << std::endl;
			break;
			case 'd' : 
				deviation = atoi(optarg);
				if (!deviation) deviation = SIGMA; // check if atoi returns 0
			break;
			default : 
				std::cerr << "Usage: ./fractals -s (seed. required) -l(max_level, default: 5) -n (normal flag, optional) -d(optional standard_deviation value. default: 2)\n
							  values for max_level must range from 2-7\n
					  		  value for standard deviation are recommended to be in the range of 2-5" << std::endl;
				exit(EXIT_FAILURE);
			break;
		}
	}
	if (!seed_provided){
		std::cerr << "Usage: ./fractals -s (seed. required) -l(max_level, default: 5) -n (normal flag, optional) -d(optional standard_deviation value. default: 2)\n
					  values for max_level must range from 2-7\n
					  value for standard deviation are recommended to be in the range of 2-5" << std::endl;
		exit(EXIT_FAILURE);
	}
	// int seed = atoi(argv[1]);

	srand(seed);
    //std::cout << gaussrand() << std::endl;
	int i, Stage;
	float delta; // "standard deviation for current level"
	size_t x, y, y0, D, d; // "indexing variables"

	// unused variables in psuedo code
	(void) i; (void) y0;

	// BEGIN:
	const size_t N = pow(2,max_level);
	// float X[33][33]; // N + 1...
	float* X = new float[(N+1) * (N+1)];
	//float Colors[33][33][3]; // rgb values pertaining to each vertex in array X
	float * Colors = new float[(N+1) * (N+1) * 3];
	//float Normals[33][33][3]; // per-vertex normal for each vertex in terrain that is surrounded by other vertices from above and to the right
	float * Normals = new float[(N+1) * (N+1) * 3];

	// "set the initial random corners"
	delta = deviation;
	X[index2D(0,0,N+1)] = delta * gaussrand();
	X[index2D(0,N,N+1)] = delta * gaussrand();
	X[index2D(N,0,N+1)] = delta * gaussrand();
	X[index2D(N,N,N+1)] = delta * gaussrand();
	
	float hausdorff = 3 - (log(N)/log(max_level));
	//float hausdorff = 0.84;

	D = N;
	d = N/2;

	for (Stage = 1; Stage < max_level + 1; ++Stage){
		// "going from grid type I to type II"
		delta = delta * pow(0.5, 0.5 * hausdorff);
		
		// "interpolate and offset mid points"
			// it is ambiguous whether this is supposed to run on the 1st iteration, hence the +1 in the stop condition
		//std::cout << "1\n";
		for (x = d; x < N - d + 1; x+=D){ 
			for (y = d; y < N - d + 1; y+=D){ 
				X[index2D(x,y, N+1)] =f4(delta, X[index2D(x+d,y+d,N+1)],X[index2D(x+d,y-d,N+1)],X[index2D(x-d,y+d,N+1)], X[index2D(x-d,y-d,N+1)]);
			}
		}
		// "displace existing points"
		//std::cout << "2\n";
		for (x = 0; x < N+1; x+=D){
			for (y = 0; y < N+1; y+=D){
				X[index2D(x,y,N+1)] = X[index2D(x,y,N+1)] + delta * gaussrand();
			}
		}
		// "going from grid type II to type I"
		delta = delta * pow(0.5, 0.5 * hausdorff);
		// "interpolate and offset mid points at boundary"
		//std::cout << "3\n";
		for (x = d; x < N - d + 1; x+=D){
			X[index2D(x,0,N+1)] = f3 (delta , X[index2D(x+d,0,N+1)],X[index2D(x-d,0,N+1)],X[index2D(x,d  , N+1)]); 
			X[index2D(x,N,N+1)] = f3 (delta , X[index2D(x+d,N,N+1)],X[index2D(x-d,N,N+1)],X[index2D(x,N-d, N+1)]); 
			X[index2D(0,x,N+1)] = f3 (delta , X[index2D(0,x+d,N+1)],X[index2D(0,x-d,N+1)],X[index2D(d,x,N+1)]);
			X[index2D(N,x,N+1)] = f3 (delta , X[index2D(N,x+d,N+1)],X[index2D(N,x-d,N+1)],X[index2D(N-d,x,N+1)]);
		}
		// "interpolate and offset mid points in interior"
		//std::cout << "4\n";
		for (x = d; x < N - d + 1; x+=D){ 
			for (y = D; y < N - d + 1; y+=D){ 
				X[index2D(x,y,N+1)] =f4(delta, X[index2D(x,y+d,N+1)], X[index2D(x,y-d,N+1)],X[index2D(x+d,y,N+1)],X[index2D(x-d,y,N+1)]);
			}
		}
		//std::cout << "5\n";
		for (x = D; x < N - d + 1; x+=D){ 
			for (y = d; y < N - d + 1; y+=D){ 
				X[index2D(x,y,N+1)] =f4(delta, X[index2D(x,y+d,N+1)], X[index2D(x,y-d,N+1)],X[index2D(x+d,y,N+1)],X[index2D(x-d,y,N+1)]);
			}
		}

		// "displace existing points"
		//std::cout << "6\n";
		for (x =0; x < N+1; x+=D){
			for (y =0; y < N+1; y+=D){
				X[index2D(x,y,N+1)] = X[index2D(x,y,N+1)] + delta *gaussrand();
			}
		}
		//std::cout << "7\n";
		for (x = d; x < N - d + 1; x+=D){
			for (y = d; y < N - d + 1; y +=D){
				X[index2D(x,y,N+1)] = X[index2D(x,y,N+1)] + delta * gaussrand();
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
	for (size_t i=0; i < N+1; ++i)
		for (size_t j=0; j < N+1; ++j){
			if (X[index2D(i, j, N + 1)] <= -4 ){
				Colors[index3D(i,j,0,N+1,N + 1)] = 0; //R
				Colors[index3D(i,j,1, N+1,N + 1)] = 0; //G
				Colors[index3D(i,j,2, N+1,N + 1)] = 1; //B
			}
			else if(X[index2D(i,j,N+1)] <= -2){
				Colors[index3D(i,j,0,N+1,N + 1)] = 0.3; //R
				Colors[index3D(i,j,1,N+1,N + 1)] = 0.3; //G
				Colors[index3D(i,j,2,N+1,N + 1)] = 1; //B
			}
			else if (X[index2D(i,j, N+1)] <= -1){
				Colors[index3D(i,j,0,N+1,N + 1)] = 0.6; //R
				Colors[index3D(i,j,1,N+1,N + 1)] = 0.6; //G
				Colors[index3D(i,j,2,N+1,N + 1)] = 1; //B
			}
			else if (X[index2D(i, j, N + 1)] <= 0){
				Colors[index3D(i,j,0,N+1,N + 1)] = 0.75; //R
				Colors[index3D(i,j,1,N+1,N + 1)] = 0.75; //G
				Colors[index3D(i,j,2,N+1,N + 1)] = 1; //B
			}
			else if (X[index2D(i, j, N + 1)] < 1){ //sand
				Colors[index3D(i,j,0,N+1,N + 1)] = 0.96; //R
				Colors[index3D(i,j,1,N+1,N + 1)] = 0.93; //G
				Colors[index3D(i,j,2,N+1,N + 1)] = 0.81; //B
			}
			else{ // green to white linear interpolation:
				float rb_color;
				if (X[index2D(i, j, N + 1)]>6 ) rb_color = 6; //limit coloration to cover the positive 0.1% of binomial distribution (SIGMA * 3)
				else rb_color = X[index2D(i, j, N + 1)];
				rb_color = (X[index2D(i, j, N + 1)] - 1) / 5; // turning Z value to a fraction to be used in linear interpolation
				Colors[index3D(i,j,0,N+1,N + 1)] = rb_color; // R
				Colors[index3D(i,j,1,N+1,N + 1)] = 1; // G
				Colors[index3D(i,j,2,N+1,N + 1)] = rb_color; // B
			}
		}

	// normal calculation:
//	for (int i = 0; i < 32; ++i)
//		for (int j = 0; j < 32; ++j){
//				//I:0 J:1 K1
//				//I:1 J:0 K2
//				//K = 0 - 1
//				//J = +K1 - 0
//				//I = +K2 - 0
//				float K1 = X[i][j+1] - X[i][j];
//				float K2 = X[i+1][j] - X[i][j];
//				
//				Normals[i][j][0]=K2;
//				Normals[i][j][1]=K1;
//				Normals[i][j][2]=-1;
//		}
	// "normal" calculation (uses the same logic from marching cubes)
	if (normals){
		for (size_t i = 1; i < N; ++i){
			for (size_t j = 1; j < N; ++j){
				Normals[index3D(i,j,0,N+1,N + 1)] = (X[index2D(i+1,j,N+1)] - X[index2D(i-1,j,N+1)]) / 2;
				Normals[index3D(i,j,1,N+1,N + 1)] = (X[index2D(i,j+1,N+1)] - X[index2D(i,j-1,N+1)]) / 2;
				Normals[index3D(i,j,2,N+1,N + 1)] = -1;
			}
			// sides
			Normals[index3D(i,0,0,N+1,N + 1)]  = (X[index2D(i+1,0,N+1)] - X[index2D(i-1,0,N+1)]) / 2;
			Normals[index3D(i,0,1,N+1,N + 1)]  = X[index2D(i,1,N+1)] - X[index2D(i,0,N+1)];
			Normals[index3D(i,0,2,N+1,N + 1)]  = -1;

			Normals[index3D(i,N,0,N+1,N + 1)] = (X[index2D(i+1,N,N+1)] - X[index2D(i-1,N,N+1)]) / 2;
			Normals[index3D(i,N,1,N+1,N + 1)] = X[index2D(i,N,N+1)] - X[index2D(i,31,N+1)];
			Normals[index3D(i,N,2,N+1,N + 1)] = -1;

			Normals[index3D(0,i,0,N+1,N + 1)]  = X[index2D(1,i,N+1)] - X[index2D(0,i,N+1)];
			Normals[index3D(0,i,1,N+1,N + 1)]  = (X[index2D(0,i+1,N+1)] - X[index2D(0,i-1,N+1)]) / 2;
			Normals[index3D(0,i,2,N+1,N + 1)]  = -1;

			Normals[index3D(N,i,0,N+1,N + 1)] =  X[index2D(N,i,N+1)] -   X[index2D(N,i,N+1)];
			Normals[index3D(N,i,1,N+1,N + 1)] = (X[index2D(N,i+1,N+1)] - X[index2D(N,i-1,N+1)]) / 2;
			Normals[index3D(N,i,2,N+1,N + 1)] = -1;
		}
	}

	unsigned int vertex_count = 0;
	std::cout << "Display \"display\" \"Screen\" \"rgbdouble\"" << std::endl;
	std::cout << "Format 1280 960" << std::endl;
	std::cout << "CameraEye   0.0 0.0 14.0" << std::endl;
	// std::cout << "CameraAt    16.0 16.0 3.0" << std::endl;
	std::cout << "CameraAt    " << N/2 << ' ' << N/2 << ' ' << 3.0 << std::endl;
	std::cout << "CameraUp   0 0 1" << std::endl;
	std::cout << "Background 0 0.3 0.7" << std::endl;
	std::cout << "WorldBegin" << std::endl;
	std::cout << "FarLight " << N/2 << ' ' << N/2 << " 20 1.0  1.0  1.0  1" << std::endl;
	//std::cout << "FarLight 16 16 20  1.0  1.0  1.0  1" << std::endl;
	std::cout << "Ka 0.3" << std::endl;
	std::cout << "Kd 0.7" << std::endl;
	std::cout << "Color 1 0 0 " << std::endl;
	std::cout << "PolySet \"PC" << (normals ? "N\"" : "\"") << std::endl;
	std::cout << pow(N,2) * 4 << std::endl << pow(N,2)*2 << std::endl;

	if (normals){
		for (size_t i = 0; i < N; ++i){
			for (size_t j = 0; j < N; ++j){
				std::cout << i << ' ' << j << ' ' << X[index2D(i,j,N+1)] 	     << ' ' << Colors[index3D(i,j,0,N+1,N + 1)    ]<< ' ' << Colors[index3D(i,j,1,N+1,N + 1)    ]<< ' ' << Colors[index3D(i,j,2,N+1,N + 1)    ]<< ' ' << Normals[index3D(i,j,0,N+1,N + 1)    ]<< ' ' << Normals[index3D(i,j,1,N+1,N + 1)    ]<< ' ' << Normals[index3D(i,j,2,N+1,N + 1)    ]<< std::endl;
				std::cout << i << ' ' << j+1 << ' ' << X[index2D(i,j+1,N+1)] 	 << ' ' << Colors[index3D(i,j+1,0,N+1,N + 1)  ]<< ' ' << Colors[index3D(i,j+1,1,N+1,N + 1)  ]<< ' ' << Colors[index3D(i,j+1,2,N+1,N + 1)  ]<< ' ' << Normals[index3D(i,j+1,0,N+1,N + 1)  ]<< ' ' << Normals[index3D(i,j+1,1,N+1,N + 1)  ]<< ' ' << Normals[index3D(i,j+1,2,N+1,N + 1)  ]<< std::endl;
				std::cout << i+1 << ' ' << j+1 << ' ' << X[index2D(i+1,j+1,N+1)] << ' ' << Colors[index3D(i+1,j+1,0,N+1,N + 1)]<< ' ' << Colors[index3D(i+1,j+1,1,N+1,N + 1)]<< ' ' << Colors[index3D(i+1,j+1,2,N+1,N + 1)]<< ' ' << Normals[index3D(i+1,j+1,0,N+1,N + 1)]<< ' ' << Normals[index3D(i+1,j+1,1,N+1,N + 1)]<< ' ' << Normals[index3D(i+1,j+1,2,N+1,N + 1)]<< std::endl;
				std::cout << i+1 << ' ' << j << ' ' << X[index2D(i+1,j,N+1)] 	 << ' ' << Colors[index3D(i+1,j,0,N+1,N + 1)  ]<< ' ' << Colors[index3D(i+1,j,1,N+1,N + 1)  ]<< ' ' << Colors[index3D(i+1,j,2,N+1,N + 1)  ]<< ' ' << Normals[index3D(i+1,j,0,N+1,N + 1)  ]<< ' ' << Normals[index3D(i+1,j,1,N+1,N + 1)  ]<< ' ' << Normals[index3D(i+1,j,2,N+1,N + 1)  ]<< std::endl;
				vertex_count +=4;
			}
		}
	}
	else {
		for (size_t i = 0; i < N; ++i){
			for (size_t j = 0; j < N; ++j){
				std::cout << i << ' ' << j << ' ' << X[index2D(i,j,N+1)] 	     << ' ' << Colors[index3D(i,j,0,N+1,N + 1)    ]<< ' ' << Colors[index3D(i,j,1,N+1,N + 1)    ]<< ' ' << Colors[index3D(i,j,2,N+1,N + 1)    ]<< std::endl;
				std::cout << i << ' ' << j+1 << ' ' << X[index2D(i,j+1,N+1)] 	 << ' ' << Colors[index3D(i,j+1,0,N+1,N + 1)  ]<< ' ' << Colors[index3D(i,j+1,1,N+1,N + 1)  ]<< ' ' << Colors[index3D(i,j+1,2,N+1,N + 1)  ]<< std::endl;
				std::cout << i+1 << ' ' << j+1 << ' ' << X[index2D(i+1,j+1,N+1)] << ' ' << Colors[index3D(i+1,j+1,0,N+1,N + 1)]<< ' ' << Colors[index3D(i+1,j+1,1,N+1,N + 1)]<< ' ' << Colors[index3D(i+1,j+1,2,N+1,N + 1)]<< std::endl;
				std::cout << i+1 << ' ' << j << ' ' << X[index2D(i+1,j,N+1)] 	 << ' ' << Colors[index3D(i+1,j,0,N+1,N + 1)  ]<< ' ' << Colors[index3D(i+1,j,1,N+1,N + 1)  ]<< ' ' << Colors[index3D(i+1,j,2,N+1,N + 1)  ]<< std::endl;
				vertex_count +=4;
			}
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
	delete[] X;
	delete[] Normals;
	delete[] Colors;
}