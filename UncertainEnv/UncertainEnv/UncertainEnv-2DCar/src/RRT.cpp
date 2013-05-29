#include "RRT.h"


void RRT::setPlannerDefaults()
{
	factor = 2.0;
	K = 0.35 * 3.1415926;
	//w_max = k / ( v + 1).

	planbias = 0.1;
	maxiter = 1000;
	maxNNiter = 1000;

	maxtries = 1000;

	maxstep = 1;
	maxchildren = 15;
	maxattempts = 10;
}

void RRT::initTree(std::vector<TreeNode>& tree, const Matrix<3>& pose)
{
	TreeNode n;
	n.T = pose;

	n.bp = -1;

	n.depth = 0.0;
	n.marked = false;
	n.attempts = 0;

	tree.push_back(n);
}

double RRT::dist(const Matrix<3>& p1, const Matrix<3>& p2)
{
	double dis = 0;
	dis = sqrt(tr(~(p1 - p2) * (p1 - p2)));

	return dis;
}

int RRT::nearestNeighbor(const std::vector<TreeNode>& tree, const Matrix<3>& point)
{
	int closest = -1;
	double mindist = 10000;

	for(int i = 0; i < tree.size(); i++){
		Matrix<3> tmpnode = tree[i].T;
		double tmpdis = dist(tmpnode, point);
		if(tmpdis < mindist && tree[i].attempts < maxattempts && tree[i].marked == false){
			closest = i;
			mindist = tmpdis;
		}
	}

	return closest;
}


bool RRT::rrtStep(std::vector<TreeNode>& tree, std::vector<int>& paths, bool& found, const clock_t& sc){
	
	found = false;
	if((double)((clock() - sc)/CLOCKS_PER_SEC) > plantime)
		return false;

	int node = -1;
	int tries = 0;

	Matrix<3> point;
	point.reset();
	double rgoal = rthreshold;

	tries = 0;
	do{
		if(tries > 25)
			rgoal = 0.5;

		if(random() < planbias){
			point[0] = goal[0] + 0.1*(2.0*random() - 1.0);
			point[1] = goal[1] + 0.1*(2.0*random() - 1.0);
			point[2] = 0 + random() * (2 * 3.1416 - 0); 
		}
		else{
			point[0] = 0 + random() * 5;
			point[1] = 0 + random() * 6;
			point[2] = 0 + random() * (2 * 3.1416 - 0);
		}

		int col = -1;
		CAL_CheckPointCollision(cal_obstacles, point[0], point[1], 0.0, false, &col);

		if(col == 0){
			node = nearestNeighbor(tree, point);
			if(node != -1){
				tree[node].attempts ++;
			}
		}

	} while ((node == -1) && (++tries < maxNNiter));
	
	if(tries == maxtries){
			return false;  //return false when cannot sample
	}

	if(node == -1)
		return false;

	Matrix<3> x_new, x_old;
	x_new.reset(); x_old.reset();
	Matrix<2> rancontrol;
	rancontrol.reset();
	x_old = tree[node].T;

	rancontrol[0] = vmin + random() * (vmax - vmin);
	rancontrol[1] = phimin + random() * (phimax - phimin);

	bool valid = false;
	int col = -1;
	double tau = dt * 1.0 / 5.0;
	for(int seg = 0; seg < 5; seg++){
		Dynamics dyn(tau);
		x_new = dyn.dynamics_zero_noise(x_old, rancontrol);
		
		CAL_CheckLineCollision(cal_obstacles, x_old[0], x_old[1], 0.0, x_new[0], x_new[1], 0.0, false, &col);
		if(col != 0)
			break;
		if(x_new[2] > 2* 3.1415926)
			x_new[2] = x_new[2] - 2 * 3.1415926;
		if(x_new[2] < 0)
			x_new[2] = x_new[2] + 2 * 3.1415926;

		x_old = x_new;
	}
	
	if(col == 0){
		TreeNode newnode;
		newnode.T = x_new;
		newnode.bp = node;
		newnode.u = rancontrol;
		
		newnode.marked = false;
		newnode.attempts = 0;

		int newid = (int) tree.size();
		tree[node].children.push_back(newid);
		tree.push_back(newnode);
		
		double dg = tr(~(x_new.subMatrix<2,1>(0,0) - goal) * (x_new.subMatrix<2,1>(0,0) - goal));
		if(dg < rthreshold * rthreshold)
		{
			TreeNode& tnode = tree[newid];
			tnode.marked = true;
			paths.push_back(newid);
			found = true;
			return true;
		}
	}
	return true;
}


bool RRT::buildRRT(std::vector<TreeNode>& tree, std::vector<int>& paths, const clock_t& sc)
{
	bool stepRet = true;
	// Initialize tree
	initTree(tree, start);

	bool found = false;
	for (int i = 0; stepRet ; ++i) {
		stepRet = rrtStep(tree, paths, found, sc);
		if(found == true)
			break;
	}
	if (stepRet && !paths.empty()) 
	{
		//std::cout << "\nRRT build: " << (double) (clock() - startClk) / CLOCKS_PER_SEC << std::endl;
		int numpaths = (int)paths.size();
		//std::cout << "Num paths found: " << numpaths << std::endl;

	} 
	else {
		std::cout << "Unable to find solution, reusing previous solution" << std::endl;
	}

	return (stepRet && !paths.empty());
}


bool RRT::executeRRT(const clock_t& sc)
{
	rrttree.clear();
	rrtpath.clear();
	rrtpaths.clear();

	return buildRRT(rrttree, rrtpaths, sc);
}

bool RRT::Plan_K_Seconds()
{
	clock_t startCLK = clock();
	
	while(((double)(clock() - startCLK) / CLOCKS_PER_SEC) < plantime){
	
		bool s = executeRRT(startCLK);
		if(s == true){
			
			std::vector<PathNode> onePath;
			onePath.clear();
			int tmp = rrtpaths[0];
			PathNode tmpnode;
			tmpnode.T = rrttree[tmp].T;
			tmpnode.u = rrttree[tmp].u;
			onePath.push_back(tmpnode);

			int parent = -1;
			while(parent != 0){
				parent = rrttree[tmp].bp;
				tmpnode.T = rrttree[parent].T;
				tmpnode.u = rrttree[parent].u;
				onePath.push_back(tmpnode);
				tmp = parent;
			}
			std::reverse(onePath.begin(), onePath.end());

			for(int i = 0; i < (int) onePath.size() - 1; i++){
				onePath[i].u = onePath[i+1].u;
			}
			onePath[(int)onePath.size() - 1].u[0] = 0;
			onePath[(int)onePath.size() - 1].u[1] = 0;

			pathSet.push_back(onePath);
			onePath.clear();
		}
	}

	if(pathSet.size() == 0){
		std::cout<<"No Path find in "<<plantime<<" seconds"<<std::endl;
		return false;
	}
	//std::cout<<pathSet.size()<<std::endl;
	return true;
}


void RRT::showPath(const int& cal_plan)
{

	Matrix<3> finalnode = zeros<3,1>();
	for(int i = 0; i < (int)pathSet.size(); i++){
		for(int t = 0; t < (int)pathSet[i].size() - 1; t++){
			double tau = dt*1.0/5;
			Dynamics dyn(tau);
			Matrix<3> oldtmp = pathSet[i][t].T;
			Matrix<3> newtmp = oldtmp;
			Matrix<2> u = pathSet[i][t].u;
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
}