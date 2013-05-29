#ifndef _RRT_
#define _RRT_
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


class RRT{

public:

	struct TreeNode{
		Matrix<3> T; //x, y, theta, v
		Matrix<2> u; //v, phi

		int bp;
		bool marked;
		int attempts;

		double depth;
		double clearance;

		std::vector<int> children;

		TreeNode()
		{
			T.reset();
			u.reset();
			bp = -1;
			marked = false;
			depth = 0.0;
			clearance = 9999999;
			attempts = 0;
		}
	};

	struct PathNode{
		Matrix<3> T;
		Matrix<2> u;
	};

	Matrix<3> start;
	Matrix<2> goal;
	double dt;
	double rthreshold;

	double factor;
	double planbias;
	double maxstep;

	int maxchildren;
	int maxattempts;
	int maxiter;
	int maxNNiter;
	int maxtries;
	double plantime;
	double K;
	std::vector<PathNode> rrtpath;
	std::vector<TreeNode> rrttree;
	std::vector<int> rrtpaths;
	std::vector<std::vector<PathNode>> pathSet;

	double vmax;
	double vmin;
	double phimax;
	double phimin;

	int cal_obstacles;

	RRT(){};
	RRT(const Matrix<3>& Start, const Matrix<2>& gl, const double& timestep, const double& ptime, const double& gr, const int& cal_obs){
		start = Start;
		goal = gl;
		dt = timestep;
		plantime = ptime;
		rthreshold = gr;
		pathSet.clear();

		Dynamics dy(dt);
		vmax = dy.vmax;
		vmin = 0;
		phimax = dy.phimax;
		phimin = -phimax;
		cal_obstacles = cal_obs;

	}

	void setPlannerDefaults();
	void initTree(std::vector<TreeNode>& tree, const Matrix<3>& pose);
	double dist(const Matrix<3>& p1, const Matrix<3>& p2);
	int nearestNeighbor(const std::vector<TreeNode>& tree, const Matrix<3>& point);
	bool rrtStep(std::vector<TreeNode>& tree, std::vector<int>& paths, bool& found,  const clock_t& sc);
	bool buildRRT(std::vector<TreeNode>& tree, std::vector<int>& paths, const clock_t& sc);
	bool executeRRT(const clock_t& sc);
	bool Plan_K_Seconds();
	void showPath(const int& cal_plan);
};


#endif