#include <iostream>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592654

#define MAX_LEVEL 5
#define SIGMA 2

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

int main(){
	srand(10);
    std::cout << gaussrand() << std::endl;
	int i, Stage;
	float delta; // "standard deviation for current level"
	int x, y, y0, D, d; // "indexing variables"

	// BEGIN:
	const int N = pow(2,MAX_LEVEL);
	float X[33][33]; // N + 1...
	// "set the initial random corners"
	delta = SIGMA;
	X[0][0] = delta * gaussrand();
	X[0][N] = delta * gaussrand();
	X[N][0] = delta * gaussrand();
	X[N][N] = delta * gaussrand();

	D = N;
	d = N/2;

	for (Stage = 1; Stage < MAX_LEVEL + 1; ++Stage){
		// "going from grid type I to type II"
		float hausdorff = 3 - (log(N)/log(Stage));
		delta = delta * pow(0.5, 0.5 * hausdorff);
		
		// "interpolate and offset mid points"
			// it is ambiguous whether this is supposed to run on the 1st iteration, hence the +1 in the stop condition
		for (x = d; x < N - d + 1; x+=D){ 
			for (y = d; y < N - d + 1; y+=D){ 
				X[x][y] =f4(delta, X[x+d][y+d],X[x+d][y-d],X[x-d][y+d], X[x-d][y-d]);
			}
		}
		// "displace existing points"
		for (x = 0; x < N+1; x+=D){
			for (y = 0; y < N+1; y+=D){
				X[x][y] = X[x][y] + delta * gaussrand();
			}
		}
		// "going from grid type II to type I"
		delta = delta * pow(0.5, 0.5 * hausdorff);
		// "interpolate and offset mid points at boundary"
		for (x = d; x < N - d + 1; x+=D){
			X[x][0] = f3 (delta , X[x+d][0],X[x-d][0],X[x][d]); 
			X[x][N] = f3 (delta , X[x+d][N],X[x-d][N],X[x][N-d]); 
			X[0][x] = f3 (delta , X[0][x+d],X[0][x-d],X[d][x]);
			X[N][x] = f3 (delta , X[N][x+d],X[N][x-d],X[N-d][x]);
		}
		// "interpolate and offset mid points in interior"
		for (x = d; x < N - d + 1; x+=D){ 
			for (y = D; y < N - d + 1; y+=D){ 
				X[x][y] =f4(delta, X[x][y+d], X[x][y-d],X[x+d][y],X[x-d][y]);
			}
		}

		for (x = D; x < N - d + 1; x+=D){ 
			for (y = d; y < N - d + 1; y+=D){ 
				X[x][y] =f4(delta, X[x][y+d], X[x][y-d],X[x+d][y],X[x-d][y]);
			}
		}

		// "displace existing points"
		for (x =0; x < N+1; x+=D){
			for (y =0; y < N+1; x+=D){
				X[x][y] = X[x][y] + delta *gaussrand();
			}
		}
		for (x = d; x < N - d + 1; x+=D){
			for (y = d; y < N - d + 1; y +=D){
				X[x][y] = X[x][y] + delta * gaussrand();
			}
		}
		// " prepare for next level"
		D = D/2;
		d = d/2;
	}
}