#ifndef _LQGMP_
#define _LQGMP_

#define _CRT_RAND_S
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
#include "Dynamics.h"
#include "RRT.h"


class LQGMP{

public:
	struct PathNodeLQG{
		Matrix<3> T;
		Matrix<2> u; 

		Matrix<3,3> A;
		Matrix<3,2> B;
		Matrix<3,2> V;
		Matrix<2,3> L;
		Matrix<3,2> K;

		Matrix<6,6> F;
		Matrix<6,4> G;
		Matrix<6,6> Sigma;

		Matrix<6> y;
		Matrix<6,6> R;

		PathNodeLQG(){
			T.reset();
			u.reset();
			A.reset();
			B.reset();
			V.reset();
			L.reset();
			K.reset();
			F.reset();
			G.reset();
			Sigma.reset();
			y.reset();
			R.reset();
		}
	};

	double dt;
	Matrix<3,3> P0; //initial covariance.
	Matrix<2,2> N; //sense noise
	Matrix<2,2> M; //process noise;
	Matrix<3,3> C;
	Matrix<2,2> D;
	Matrix<2,3> H;
	
	Primitive Prim;

	std::vector<PathNodeLQG> m_rrtpathLQG;
	std::vector<PathNodeLQG> pathlqg;
	std::vector<RRT::PathNode> m_rrtpath;
	std::vector<PathNodeLQG> m_solutionLQG;

	Matrix<2> goal;
	double goal_radius;

	double r;
	double phimax;
	double vmax;
	int cal_obstacles;

	LQGMP(){};
	LQGMP(const std::vector<RRT::PathNode>& rawpath, const double& timestep, const Matrix<3,3>& initialCov, Primitive& pm, const int& cal_obs){
		//m_rrtpath.clear();
		m_rrtpath = rawpath;
		dt = timestep;
		P0 = initialCov;
		
		Prim = pm;

		int size = (int)m_rrtpath.size();
		pathlqg.clear();

		M = identity<2>() * 0.05*0.05;
		N = identity<2>() * 0.08*0.08;

		H.reset();
		H(0,0) = 1; H(0,1) = 0; H(0,2) = 0; 
		H(1,0) = 0; H(1,1) = 1; H(1,2) = 0;
		
		C = identity<3>();
		D = identity<2>();

		Dynamics dy(dt);
		r = dy.r; vmax = dy.vmax; phimax = dy.phimax;

		cal_obstacles = cal_obs;
	}

	void createABVLK();
	void setGoal(const Matrix<2>& gl, const double& gr);
	double computeConfidence(const Matrix<2>& pos, const Matrix<2,2>& R, const int& cal_obstacles, const int& cal_environment, const int& cal_point);
	double computeProbability(const Matrix<3,3>& initialCov, const int& cal_obstacles, const int& cal_environment, const int& cal_point);
	void draw_prior_distribution(const int& cal_ellipse);
	bool LQGSimulate(const Matrix<3,3>& initialCov,  Primitive& PrimCollision, const int& cal_lqg);

	//truncation
	void query(const Matrix<2>& pos, const Matrix<2,2>& R, std::vector<std::pair<Matrix<2,1>, double>>& cvx, const int& cal_obstacles, const int& cal_environment, const int& cal_point, std::vector<std::pair<Matrix<2,1>, double>>& cvxprob);
	double computeStatics(const Matrix<2>& pos, const Matrix<2,2>& R, const std::vector<std::pair<Matrix<2>, double>>& cvxprob);
	double boolprobsuccess(const int& cal_obstacles, const int& cal_environment, const int& cal_point, const Matrix<3,3>& initialCov);
	void lqgmpTruncation( Matrix<2*3>& X,  Matrix<2*3, 2*3>& S, std::vector<std::pair<Matrix<2,1>, double>>& cvx);
	double computeLQGMPTruncate(const int& cal_obstacles, const int& cal_environment, const int& cal_point, const Matrix<3,3>& initialCov);
	void draw_truncate_distribution(const int& cal_ellipse_trunc);
};

#endif