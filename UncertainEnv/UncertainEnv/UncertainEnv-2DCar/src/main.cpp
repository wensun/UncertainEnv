///////////////////
//This code is used to implement the idea of real-time replanning. 
//namely, replanning ahead. 
//assume the time for replanning is also dt
///////////////////////////////////

#define _CRT_RAND_S

#define DIM 2
#define X_DIM 4  //X = (x,y,theta, v)  theta: orientation, v speed.
#define U_DIM 2  //a, phi, a: acceleration, phi: steering wheel angel
#define Z_DIM 3  //it can observe position x and position y.  
#define INFTY 9e9


//define constratins:


const static double goal_radius = 0.3;
const static double plan_goal_radius = 0.6;
const static double dt = 0.5;  // asume the time for replan is also dt. 
const static double car_l = 1.0;
const static double MaxPlanTime = 1.0;


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
#include "LQGMP.h"

using namespace Concurrency;
//using namespace std;

static Matrix<2> goal;
static Matrix<3> start;
static Matrix<3, 3> P0; //covariance

/***********************Environment settings***************************/
int cal_environment;
int cal_line;
int cal_start;
int cal_obstacles;
int cal_goal;
int cal_rrt;
int cal_paths;
int cal_ellipse;
int cal_ellipse_trunc;
int cal_point;
int cal_cienvironment, cal_ciobstacles, cal_cipoint;
int cal_plan;
int cal_execute;
int cal_robot;
static std::vector<Matrix<2,1>> obstaclesSet; //ax + by = b
/******************End of Environment settings************************/


double Random()
{
	double a;
	a = rand()%1000;
	return a / 1000;
}

/**********************init environment***************************/
void initEnvironment() 
{
	goal[0] = 4.5;
	goal[1] = 3;

	start[0] = 0.5; start[1] = 0.5; start[2] = 1.0/2.0 * M_PI;
	P0 = identity<3>() * 0.001;
	// callisto params
	CAL_SetViewParams(0, 0, 0, 14, 0, 0, 0, 0, 1, 0); 

	CAL_CreateGroup(&cal_environment, 0, true, "Environment");
	CAL_CreateGroup(&cal_obstacles, cal_environment, true, "Obstacles");
	CAL_SetGroupColor(cal_obstacles, 1, 0, 0, 0.5);

	CAL_CreateBox(cal_obstacles, 0.5, 6, 0.2, -0.25, 3, 0); 	
	CAL_CreateBox(cal_obstacles, 0.5, 6, 0.2, 5.25, 3, 0);
	CAL_CreateBox(cal_obstacles, 6, 0.5, 0.2, 2.5, 6.25, 0);
	CAL_CreateBox(cal_obstacles, 6, 0.5, 0.2, 2.5, -0.25, 0);
	CAL_CreateBox(cal_obstacles, 2, 4, 0.2, 2, 3, 0);
	CAL_CreateBox(cal_obstacles, 1, 2, 0.2, 4.5, 1,0);
	CAL_CreateGroup(&cal_rrt, 0, false, "RRT");
	CAL_SetGroupColor(cal_rrt, 1, 0, 1);
	
	CAL_CreateGroup(&cal_line, 0, false, "Line");
	CAL_SetGroupColor(cal_line, 0, 0, 0);
	
	CAL_CreateGroup(&cal_paths, 0, false, "Paths");
	CAL_SetGroupColor(cal_paths, 0, 0, 1);

	CAL_CreateGroup(&cal_point, 0, true, "Point");
	CAL_CreateSphere(cal_point, 0.001, 0, 0, 0);

	CAL_CreateGroup(&cal_robot, 0, true, "Robot");
	CAL_SetGroupColor(cal_robot, 1, 0, 0);

	CAL_CreateGroup(&cal_ellipse, 0, false, "Ellipse");
	CAL_SetGroupColor(cal_ellipse, 0, 1, 0, 0.5);

	CAL_CreateGroup(&cal_ellipse_trunc, 0, false, "Ellipse_trunc");
	CAL_SetGroupColor(cal_ellipse_trunc, 1, 0, 0);

	CAL_CreateGroup(&cal_plan, 0, false, "newplan");
	CAL_SetGroupColor(cal_plan, 0, 1, 0);

	CAL_CreateGroup(&cal_execute, 0, false, "execute");
	CAL_SetGroupColor(cal_execute, 0, 1, 0);
	
	CAL_CreateGroup(&cal_goal, 0, false, "Goal region");
	CAL_SetGroupColor(cal_goal, 0, 1, 1, 0.5);
	CAL_CreateCylinder(cal_goal, (float) goal_radius, 0.05f, 0, 0, 0);
	CAL_SetGroupPosition(cal_goal, (float) goal[0], (float) goal[1], -0.025f);
	CAL_SetGroupOrientation(cal_goal, (float) (M_PI*0.5), 0, 0);

	CAL_CreateGroup(&cal_start, 0, false, "start");
	CAL_SetGroupColor(cal_start, 0, 1, 0);
	CAL_CreateSphere(cal_start, 0.04, start[0], start[1], 0);
	
}


void uncertainEnv()
{
	Primitive prim;
	for(int i = 0; i < 10; i++){
		prim.CreateVertex();
		float line[6] = {prim.randomvertex[0][0], prim.randomvertex[0][1], 0, prim.randomvertex[1][0], prim.randomvertex[1][1], 0};
		int np[1] = {2};
		CAL_CreatePolyline(cal_line, 1, np, line);

		float line1[6] = {prim.randomvertex[1][0], prim.randomvertex[1][1], 0, prim.randomvertex[2][0], prim.randomvertex[2][1], 0};
		CAL_CreatePolyline(cal_line, 1, np, line1);
		float line2[6] = {prim.randomvertex[2][0], prim.randomvertex[2][1], 0, prim.randomvertex[3][0], prim.randomvertex[3][1], 0};
		CAL_CreatePolyline(cal_line, 1, np, line2);
		float line3[6] = {prim.randomvertex[3][0], prim.randomvertex[3][1], 0, prim.randomvertex[0][0], prim.randomvertex[0][1], 0};
		CAL_CreatePolyline(cal_line, 1, np, line3);

		float line4[6] = {prim.randomvertex[4][0], prim.randomvertex[4][1], 0, prim.randomvertex[5][0], prim.randomvertex[5][1], 0};
		CAL_CreatePolyline(cal_line, 1, np, line4);
		float line5[6] = {prim.randomvertex[5][0], prim.randomvertex[5][1], 0, prim.randomvertex[6][0], prim.randomvertex[6][1], 0};
		CAL_CreatePolyline(cal_line, 1, np, line5);
		float line6[6] = {prim.randomvertex[6][0], prim.randomvertex[6][1], 0, prim.randomvertex[7][0], prim.randomvertex[7][1], 0};
		CAL_CreatePolyline(cal_line, 1, np, line6);
		float line7[6] = {prim.randomvertex[7][0], prim.randomvertex[7][1], 0, prim.randomvertex[4][0], prim.randomvertex[4][1], 0};
		CAL_CreatePolyline(cal_line, 1, np, line7);
		float line8[6] = {prim.randomvertex[8][0], prim.randomvertex[8][1], 0, prim.randomvertex[9][0], prim.randomvertex[9][1], 0};
		CAL_CreatePolyline(cal_line, 1, np, line8);
		float line9[6] = {prim.randomvertex[8][0], prim.randomvertex[8][1], 0, prim.randomvertex[10][0], prim.randomvertex[10][1], 0};
		CAL_CreatePolyline(cal_line, 1, np, line9);

	}

	prim.CreateSegments(); //create segments for background environment.
	//for(int i = 0; i < prim.seg.size(); i++){
	//	std::cout<<"("<<prim.seg[i].S1(0,0)<<" "<<prim.seg[i].S1(1,1)<<") "<<"("<<prim.seg[i].S2(0,0)<<" "<<prim.seg[i].S2(1,1)<<") "<<std::endl;
	//}

}



/*********************End of environment init**************************/

void showPath(const int& cal_plan, const std::vector<RRT::PathNode>& rrtpath)
{
	Matrix<3> finalnode = zeros<3,1>();
	for(int t = 0; t < (int)rrtpath.size() - 1; t++){
		double tau = dt*1.0/5;
		Dynamics dyn(tau);
		Matrix<3> oldtmp = rrtpath[t].T;
		Matrix<3> newtmp = oldtmp;
		Matrix<2> u = rrtpath[t].u;
		for(int k = 0; k < 5; k++){
			newtmp = dyn.dynamics_zero_noise(oldtmp, u);
			finalnode = newtmp;
			float line[6] = {oldtmp[0], oldtmp[1], 0, newtmp[0], newtmp[1], 0};
			int np[1] = {2};
			CAL_CreatePolyline(cal_plan, 1, np, line);
			oldtmp = newtmp;
		}
	}
}



int main()
{

	srand(1000);
	
	CAL_Initialisation (true, true, true);
	initEnvironment();
	uncertainEnv();

	double plantime = 2;
	double dt = 0.5;
	Matrix<3,3> P0 = identity<3>() * 0.1*0.1;


	RRT rrt(start, goal, dt, plantime, goal_radius, cal_obstacles);
	rrt.setPlannerDefaults();
	rrt.Plan_K_Seconds();
	showPath(cal_rrt, rrt.pathSet[0]);
	//rrt.showPath(cal_rrt);

	LQGMP lqgmp(rrt.pathSet[0], dt, P0);
	double prob = exp(lqgmp.computeProbability(P0, cal_obstacles, cal_environment, cal_point));
	lqgmp.draw_prior_distribution(cal_ellipse);

	double ps = lqgmp.boolprobsuccess(cal_obstacles, cal_environment, cal_point, P0);
	std::cout<<ps<<std::endl;

	double pstrunc = lqgmp.computeLQGMPTruncate(cal_obstacles, cal_environment, cal_point, P0);
	std::cout<<pstrunc<<std::endl;
	lqgmp.draw_truncate_distribution(cal_ellipse_trunc);

	//for(int i = 0; i < 10; i++){
	//	Primitive collision;
	//	collision.CreateVertex();
	//	collision.CreateRandomSegs();
	//	LQGMP lqgmp(rrt.pathSet[0], dt, P0);
	//	lqgmp.LQGSimulate(P0, collision, cal_execute);
	//}


	int num;
	std::cin>>num;

	CAL_End();
	return 0;
}