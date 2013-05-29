#ifndef _DYNAMICS_
#define _DYNAMICS_


#include <vector>
#include <map>
#include <queue>
#include <list>
#include <time.h>
#include <float.h>
#include <fstream>
#include <windows.h>
#include <algorithm>
#include <ppl.h>
#include <omp.h>
#include "Dynamics.h"
#include "callisto.h"
#include "matrix.h"
#include "utils.h"
#include "Primitive.h"


class Dynamics{

public:
	double tau;
	Dynamics();
	double vmax;
	double phimax;
	double r;

	Dynamics(const double& dt){
		r = 0.35;
		tau = dt;
		vmax = 0.5;
		phimax = 3.14159 / 3;
	}

	Matrix<3> dynamics_zero_noise(const Matrix<3>& X, const Matrix<2>& u){
		
		int seg = 5;
		double t = tau * 1.0 / seg;
		Matrix<3> x_old = X;
		Matrix<3> x_new = x_old;

		for(int i = 1; i <= seg; i++){
			x_new[0] = x_old[0] + t * u[0] * cos(x_old[2]);
			x_new[1] = x_old[1] + t * u[0] * sin(x_old[2]);
			x_new[2] = x_old[2] + t * u[0] * tan(u[1]) / (r * tan(phimax)); 

			if(x_new[2] > 2*3.14159)
				x_new[2] = x_new[2] - 2 * 3.1415926;
			if(x_new[2] < 0)
				x_new[2] = x_new[2] + 2 * 3.1415926;

			x_old = x_new;

		}
		return x_new;
	}

	Matrix<3> dynamics_noise(const Matrix<3>& X, const Matrix<2>& u, const Matrix<2>& noise){
		int seg = 5;
		double t = tau * 1.0 / seg;
		Matrix<3> x_old = X;
		Matrix<3> x_new = x_old;

		Matrix<2> corruptu = u + noise;

		for(int i = 1; i <= seg; i++){
			x_new[0] = x_old[0] + t * corruptu[0] * cos(x_old[2]);
			x_new[1] = x_old[1] + t * corruptu[0] * sin(x_old[2]);
			x_new[2] = x_old[2] + t * corruptu[0] * tan(corruptu[1]) / (r * tan(phimax)); 

			if(x_new[2] > 2*3.14159)
				x_new[2] = x_new[2] - 2 * 3.1415926;
			if(x_new[2] < 0)
				x_new[2] = x_new[2] + 2 * 3.1415926;

			x_old = x_new;

		}
		return x_new;
	}

};



#endif