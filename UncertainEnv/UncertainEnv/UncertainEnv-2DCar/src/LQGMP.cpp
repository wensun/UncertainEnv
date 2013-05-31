#include "LQGMP.h"

void LQGMP::createABVLK()
{
	int ns = (int) m_rrtpath.size();	
	pathlqg.clear();
	for(int i = 0; i < ns; i++){
		PathNodeLQG tmpnode;
		tmpnode.T = m_rrtpath[i].T;
		tmpnode.u = m_rrtpath[i].u;
		pathlqg.push_back(tmpnode);
	}

	for(int i = 0; i < ns; i++){
		RRT::PathNode tmp = m_rrtpath[i];
		Matrix<2> u = tmp.u;
		Matrix<3> X = tmp.T;

		Matrix<3,3> A;
		A.reset();
		Matrix<3,2> B;
		B.reset();
		Matrix<3,2> V;
		V.reset();

		double car_l = r * tan(phimax);

		A(0,0) = 1; A(0,1) = 0; A(0,2) = -dt*u[0]*sin(X[2]);
		A(1,0) = 0; A(1,1) = 1; A(1,2) = dt*u[0]*cos(X[2]);
		A(2,0) = 0; A(2,1) = 0; A(2,2) = 1; 
		
		B(0,0) = dt*cos(X[2]); B(0,1) = 0;
		B(1,0) = dt*sin(X[2]); B(1,1) = 0;
		B(2,0) = dt*tan(u[1])/car_l; B(2,1) = dt*u[0]*(1 + tan(u[1])*tan(u[1])) / car_l;

		V = B;
		pathlqg[i].A = A;
		pathlqg[i].B = B;
		pathlqg[i].V = V;
		
	}

	Matrix<3,3> S = C;
	int length = (int)pathlqg.size() - 1;
	for(int k = length - 1; k != -1; k--){
		pathlqg[k].L = -!(~pathlqg[k].B*S*pathlqg[k].B + D)*~pathlqg[k].B*S*pathlqg[k].A;
		S = C + ~pathlqg[k].A*S*pathlqg[k].A + ~pathlqg[k].A*S*pathlqg[k].B*pathlqg[k].L;
	}

	Matrix<3,3> P = P0;
	for(int k = 1; k <= length; k++){
		P = pathlqg[k-1].A*P*~pathlqg[k-1].A + pathlqg[k-1].V * M *~pathlqg[k-1].V;
		
		pathlqg[k].K = P * ~H * !(H*P*~H + N);
		P = (identity<3>() - pathlqg[k].K * H) * P;
	}

	pathlqg[0].Sigma = zeros<6,6>();
	pathlqg[0].Sigma.insert(0,0, P0);

	for(int k = 1; k <= length; k++){
		Matrix<6,6> F;
		F.insert(0,0, pathlqg[k-1].A);
		F.insert(0,3, pathlqg[k-1].B * pathlqg[k-1].L);
		F.insert(3,0, pathlqg[k].K * H * pathlqg[k-1].A);
		F.insert(3,3, pathlqg[k-1].A + pathlqg[k-1].B * pathlqg[k-1].L - pathlqg[k].K * H * pathlqg[k-1].A);
		pathlqg[k-1].F = F;

		Matrix<6, 4> G;
		G.insert(0,0, pathlqg[k-1].V);
		G.insert(0,2, zeros<3,2>());
		G.insert(3,0, pathlqg[k].K * H * pathlqg[k-1].V);
		G.insert(3,2, pathlqg[k].K);
		pathlqg[k-1].G = G;
		
		Matrix<4,4> Q;
		Q.reset();
		Q.insert(0,0, M);
		Q.insert(2,2, N);
		pathlqg[k].Sigma = F * pathlqg[k-1].Sigma * ~F + G * Q * ~G;
	}
}

void LQGMP::setGoal(const Matrix<2>& gl, const double& gr)
{
	goal = gl;
	goal_radius = gr;
}

double LQGMP::computeConfidence(const Matrix<2>& pos, const Matrix<2,2>& R, const int& cal_obstacles, const int& cal_environment, const int& cal_point){
	
	Matrix<2,2> EVec, EVal;
	jacobi(R, EVec, EVal);
	Matrix<3,3> Temp = identity<3>();
	Temp.insert(0,0, ~EVec);
	Matrix<4,1> q = quatFromRot(Temp);

	CAL_SetGroupQuaternion(cal_obstacles, (float) q(0,0), (float) q(1,0), (float) q(2,0), (float) q(3,0));
	CAL_SetGroupScaling(cal_environment, 1/(float)sqrt(EVal(0,0)), 1/(float)sqrt(EVal(1,1)), 1);

	Matrix<2,2> invScale = zeros<2,2>();
	invScale(0,0) = 1/(float)sqrt(EVal(0,0));
	invScale(1,1) = 1/(float)sqrt(EVal(1,1));
	Matrix<2> transPos =  invScale * ~EVec * pos;

	CAL_SetGroupPosition(cal_point, (float) transPos[0], (float) transPos[1], 0);

	int num_pairs;
	CAL_GetClosestPairs(cal_point, cal_environment, &num_pairs);
	SCALResult* results = new SCALResult[num_pairs];
	CAL_GetResults(results);
	double distance = results[0].distance;
	delete [] results;

	CAL_SetGroupQuaternion(cal_obstacles, 0,0,0,1);
	CAL_SetGroupScaling(cal_environment, 1,1,1);

	return distance;
}

double LQGMP::computeProbability(const Matrix<3,3>& initialCov, const int& cal_obstacles, const int& cal_environment, const int& cal_point)
{
	P0 = initialCov;
	createABVLK();
	double log_prob = 0;
	for(int i = 0; i < (int)pathlqg.size(); i++)
	{
		Matrix<2,1> pos = pathlqg[i].T.subMatrix<2,1>(0,0);
		Matrix<2,2> R = pathlqg[i].Sigma.subMatrix<2,2>(0,0);
		double conf = computeConfidence(pos, R, cal_obstacles, cal_environment, cal_point);
		double p = log(incompletegamma(1, 0.5*conf*conf));
		log_prob += p;
	}

	return log_prob;
}

void LQGMP::draw_prior_distribution(const int& cal_ellipse){

	//createABVLK();
	for(int i = 0; i < (int)pathlqg.size(); i+=3){
		drawEllipse2d(pathlqg[i].T.subMatrix<2,1>(0,0), pathlqg[i].Sigma.subMatrix<2,2>(0,0), cal_ellipse, false);
	}
}


bool LQGMP::LQGSimulate(const Matrix<3,3>& initialCov, Primitive& PrimCollision, const int& cal_lqg){
	Dynamics dyn(dt);
	P0 = initialCov;
	createABVLK();

	Matrix<3> x = pathlqg[0].T;
	Matrix<3> x_true_old = x + sampleGaussian(zeros<3,1>(), P0);
	Matrix<3> x_true_new = x_true_old;
	Matrix<3> x_est_d = x - pathlqg[0].T;

	for(int i = 1; i < (int)pathlqg.size(); i++){
		Matrix<2> u = pathlqg[i-1].u + pathlqg[i-1].L*x_est_d;
		x_true_new = dyn.dynamics_noise(x_true_old, u, sampleGaussian(zeros<2,1>(), M));

		bool flag = false;
		flag = PrimCollision.checkline(x_true_new.subMatrix<2,1>(0,0), x_true_old.subMatrix<2,1>(0,0));
		if(flag == false){
			//std::cout<<"collision happen"<<std::endl;
			return false;
		}
		float line[6] = {x_true_new[0], x_true_new[1], 0.0, x_true_old[0], x_true_old[1], 0.0};
		int np[1] = {2};
		CAL_CreatePolyline(cal_lqg, 1, np, line);

		x_est_d = pathlqg[i-1].A * x_est_d + pathlqg[i-1].B * pathlqg[i-1].L * x_est_d;
		Matrix<2> z = H * x_true_new - H * pathlqg[i].T + sampleGaussian(zeros<2,1>(), N);
		x_est_d = x_est_d + pathlqg[i].K * (z - H * x_est_d);
		x_true_old = x_true_new;
	}
	return true;
}


//truncation
void LQGMP::query(const Matrix<2>& pos, const Matrix<2,2>& R, std::vector<std::pair<Matrix<2,1>, double>>& cvx, const int& cal_obstacles, const int& cal_environment, const int& cal_point, std::vector<std::pair<Matrix<2,1>, double>>& cvxprob)
{
	Matrix<2> transPos = pos;
	CAL_SetGroupPosition(cal_point, (float) transPos[0], (float) transPos[1], 0);
	int col  = -1;
	CAL_CheckPointCollision(cal_environment, transPos[0], transPos[1], 0, false, &col);
	if(col != 0){
		std::cout<<"collision happen in query"<<std::endl;
		return;
	}

	int num_pairs = 0;
	if(CAL_GetClosestPairs(cal_point, cal_obstacles, &num_pairs) != 0){
		cvx.clear();
		cvxprob.clear();
		return;
	}

	SCALResult* results = new SCALResult[num_pairs];
	CAL_GetResults(results);

	//store all the contacts.
	std::vector<Matrix<2>> contacts;
	contacts.resize(num_pairs);
	for(int i = 0; i < num_pairs; i++){
		Matrix<2> p = zeros<2,1>();
		p[0] = (double)results[i].vector1[0];
		p[1] = (double)results[i].vector1[1];
		contacts[i] = p;

		Matrix<2> po = zeros<2,1>();
		po[0] = (double)results[i].vector0[0];
		po[1] = (double)results[i].vector0[1];
		Matrix<2> n = (p - po);
		n = n / sqrt(tr(~n*n));
		double b = tr(~n * p);
		cvxprob.push_back(std::make_pair(n, b));
	}

	//compute local convex region.
	for(int i = 0; i < num_pairs; i++){
		Matrix<2> p = zeros<2,1>();
		p[0] = (double)results[i].vector1[0];
		p[1] = (double)results[i].vector1[1];

		int j = 0;
		for(j = 0; j < (int)contacts.size(); j++){
			if(sqrt(tr(~(contacts[j] - p)*(contacts[j] - p))) <= 0.005)
				break;
		}
		if(j == (int)contacts.size())
			;
		else if( j < (int)contacts.size()){
			Matrix<2> po = zeros<2,1>();
			po[0] = (double)results[i].vector0[0];
			po[1] = (double)results[i].vector0[1];

			Matrix<2> n = (p - po);
			n = n/sqrt(tr(~n*n));
			double b = tr(~n*p);
			cvx.push_back(std::make_pair(n, b));

			for(int k = i+1; k < (int)contacts.size(); k++){
				Matrix<2> pp = contacts[k];
				if(tr(~n * pp) > b){
					contacts.erase(contacts.begin() + k);
					k--;
				}
			}

		}
	}
	delete []results;
}

//compute the probablity of success for a step. there is no truncation.
double LQGMP::computeStatics(const Matrix<2>& pos, const Matrix<2,2>& R, const std::vector<std::pair<Matrix<2>, double>>& cvxprob)
{
	double ps = 0.0;
	int cvxlen = (int)cvxprob.size();
	for(int i = 0; i < cvxlen; i++){
		Matrix<2> a = cvxprob[i].first;
		double b = cvxprob[i].second;

		int s = 0;
		for(s = 0; s < (int)Prim.seg.size(); s++){
			if(sqrt(tr(~(a-Prim.seg[s].a)*(a-Prim.seg[s].a))) + abs(b - Prim.seg[s].b) < 0.0005)
				break; //find the corresponding line segment. 
		}
		if(s < (int)Prim.seg.size()){ //means find the line segment
			Primitive::segments ssegment = Prim.seg[s];
			//start compute the approximated probablity of success.
			Matrix<2,2> A = ssegment.RotationM;
			Matrix<2,1> x = pos;
			Matrix<2,1> p1 = ssegment.p1; Matrix<2,1> p2 = ssegment.p2;
			Matrix<2,1> D = A*p1 - A*p2; Matrix<2,1> E = (~A)*x - 2*(~A)*p1 + A*p2;
			Matrix<2,1> F = -(~A)*x + (~A)*p1;
			double mean = tr(~p1*~A*x) - tr(~p2*~A*x) - tr(~p1*~A*p1) + tr(~p2*~A*p1);
			double var = tr(~D*R*D) + tr(~E*ssegment.S1*E) + tr(~F*ssegment.S2*F);
			double alpha = (0 - mean) / sqrt(var);
			double p = 1 - cdf(alpha);
			ps += p;
		}
		else if(s == (int)Prim.seg.size()){ //did not find line segment, then it should be a vertex.
			//search for vertex;
			int v = 0;
			for(v = 0; v < (int)Prim.points.size(); v++){
				Matrix<2,1> contact = Prim.points[v].first;
				if(abs(tr(~a*contact)-b) < 0.0005)
					break; //find the vertex as the contact point.
			}
			if(v < (int)Prim.points.size()){
				Matrix<2,1> c = Prim.points[v].first;
				Matrix<2,2> ccov = Prim.points[v].second;
				Matrix<2,1> x = pos;
				double mean = 2*tr(~c*x) - tr(~x*x) - tr(~c*c);
				Matrix<2> D = 2*c - 2*x; Matrix<2> E = 2*x - 2*c;
				double var = tr(~D*R*D) + tr(~E*ccov*E);
				double alpha = (0 - mean) / sqrt(var);
				double p = 1 - cdf(alpha);
				ps += p;
			}
			else
				std::cout<<"did not find a vertex nor a line segments"<<std::endl;
		}
	}
	return 1 - ps;
}


//return probablity of success measured by using boole inequality.
double LQGMP::boolprobsuccess(const int& cal_obstacles, const int& cal_environment, const int& cal_point, const Matrix<3,3>& initialCov)
{
	P0 = initialCov;
	createABVLK();
	double probs = 1.0;
	for(int i = 0; i < (int)pathlqg.size(); i++)
	{
		Matrix<2,1> pos = pathlqg[i].T.subMatrix<2,1>(0,0);
		Matrix<2,2> R = pathlqg[i].Sigma.subMatrix<2,2>(0,0);
		std::vector<std::pair<Matrix<2,1>, double>> cvx;
		std::vector<std::pair<Matrix<2,1>, double>> cvxprob;
		cvx.clear(); cvxprob.clear();
		query(pos, R, cvx, cal_obstacles, cal_environment, cal_point, cvxprob);

		//double stepp = computeStatics(pos, R, cvxprob);
		double stepp = computeStatics(pos, R, cvx);
		probs *= stepp;
	}
	return probs;
}


//truncate the gaussian distribution. Current ignoring truncating the vertex in the environments
void LQGMP::lqgmpTruncation(Matrix<2*3>& X, Matrix<2*3, 2*3>& S, std::vector<std::pair<Matrix<2,1>, double>>& cvx)
{
	Matrix<6> xDelta; 
	Matrix<6,6> sDelta;
	double yMean = 0, yVar = 0, yNewMean = 0, yNewVar = 0;
	xDelta.reset();
	sDelta.reset();


	int ncvx = (int)cvx.size();
	for(int i = 0; i < ncvx; i++){
		Matrix<2> a = cvx[i].first;
		double b = cvx[i].second;

		int s = 0;
		for(s = 0; s < (int)Prim.seg.size(); s++){
			if(sqrt(tr(~(a-Prim.seg[s].a)*(a-Prim.seg[s].a))) + abs(b - Prim.seg[s].b) < 0.0001)
				break; //find the corresponding line segment. 
		}
		if(s < (int)Prim.seg.size()){ //means find the line segment
			//std::cout<<"aaaa"<<std::endl;
			Primitive::segments ssegment = Prim.seg[s];
			//start compute the approximated probablity of success.
			Matrix<2,2> A = ssegment.RotationM;
			Matrix<2,1> x = X.subMatrix<2,1>(0,0);
			Matrix<2,1> p1 = ssegment.p1; Matrix<2,1> p2 = ssegment.p2;
			Matrix<2,1> D = A*p1 - A*p2; Matrix<2,1> E = (~A)*x - 2*(~A)*p1 + A*p2;
			Matrix<2,1> F = -(~A)*x + (~A)*p1;
			double M = tr(~p1*~A*x) - tr(~p2*~A*x) - tr(~p1*~A*p1) + tr(~p2*~A*p1) - tr(~D*x)- tr(~E*p1) - tr(~F*p2);
			Matrix<10,1> aa = zeros<10,1>();
			aa.insert<2,1>(0,0,D); aa.insert<2,1>(6,0, E); aa.insert<2,1>(8,0, F);
			double bb = -M;
			Matrix<10,1> augmentX = zeros<10,1>(); augmentX.insert<6,1>(0,0,X); augmentX.insert<2,1>(6,0, p1); augmentX.insert<2,1>(8,0, p2);
			Matrix<10,10> augmentCov = zeros<10,10>(); augmentCov.insert<6,6>(0,0,S); augmentCov.insert<2,2>(6,6, ssegment.S1); augmentCov.insert<2,2>(8,8, ssegment.S2);
			yMean = tr(~aa*augmentX);
			yVar = tr(~aa*augmentCov*aa);
			truncate(bb, yMean, yVar, yNewMean, yNewVar);
			Matrix<10,1> L = augmentCov*aa / yVar;
			Matrix<10,1> xDeviation = L*(yMean - yNewMean);
			Matrix<10,10> SDeviation = (yVar - yNewVar)*(L*~L);
			xDelta += xDeviation.subMatrix<6,1>(0,0);
			sDelta += SDeviation.subMatrix<6,6>(0,0);

			int index1 = ssegment.index1; int index2 = ssegment.index2;
			Prim.tmppoints[index1].first -= xDeviation.subMatrix<2,1>(6,0);
			Prim.tmppoints[index2].first -= xDeviation.subMatrix<2,1>(8,0);
			Prim.tmppoints[index1].second -= SDeviation.subMatrix<2,2>(6,6);
			Prim.tmppoints[index2].second -= SDeviation.subMatrix<2,2>(8,8);

		}
		else if(s == (int)Prim.seg.size()){
			//search for vertex;
			int v = 0;
			for(v = 0; v < (int)Prim.points.size(); v++){
				Matrix<2,1> contact = Prim.points[v].first;
				if(abs(tr(~a*contact)-b) < 0.0001)
					break; //find the vertex as the contact point.
			}
			if(v < (int)Prim.points.size()){
				Matrix<2,1> c = Prim.points[v].first;
				Matrix<2,2> ccov = Prim.points[v].second;
				Matrix<2,1> x = X.subMatrix<2,1>(0,0);
				Matrix<2> D = 2*c - 2*x; Matrix<2> E = 2*x - 2*c;
				double M = 2*tr(~c*x) - tr(~x*x) - tr(~c*c) - tr(~D*x) - tr(~E*c);
				Matrix<8,1> aa = zeros<8,1>(); aa.insert<2,1>(0,0,D); aa.insert<2,1>(6,0, E);
				double bb = -M;
				Matrix<8,1> augmentX = zeros<8,1>(); augmentX.insert<6,1>(0,0, X); augmentX.insert<2,1>(6,0, c);
				Matrix<8,8> augmentCov = zeros<8,8>(); augmentCov.insert<6,6>(0,0, S); augmentCov.insert<2,2>(6,6, ccov);
				yMean = tr(~aa*augmentX);
				yVar = tr(~aa*augmentCov*aa);
				truncate(bb, yMean, yVar, yNewMean, yNewVar);
				Matrix<8,1> L = augmentCov*aa / yVar;
				Matrix<8,1> xDeviation = L*(yMean - yNewMean);
				Matrix<8,8> SDeviation = (yVar - yNewVar)*(L*~L);
				xDelta += xDeviation.subMatrix<6,1>(0,0);
				sDelta += SDeviation.subMatrix<6,6>(0,0);

				Prim.tmppoints[v].first -= xDeviation.subMatrix<2,1>(6,0);
				Prim.tmppoints[v].second -= SDeviation.subMatrix<2,2>(6,6);
			}
		}
	}
	X -= xDelta;
	S -= sDelta;
	Prim.UpdateEnvironment(cal_obstacles);
}


double LQGMP::computeLQGMPTruncate(const int& cal_obstacles, const int& cal_environment, const int& cal_point, const Matrix<3,3>& initialCov){

	double ps = 1.0;
	P0 = initialCov;
	createABVLK();

	Matrix<4,4> Q;
	Q.reset();
	Q.insert(0,0, M);
	Q.insert(2,2, N);

	pathlqg[0].y = zeros<6,1>();
	pathlqg[0].R = zeros<6,6>(); pathlqg[0].R.insert<3,3>(0,0, P0);
	Matrix<6> X = zeros<6,1>(); X.insert<3,1>(0,0, pathlqg[0].T); X.insert<3,1>(3,0, pathlqg[0].T);
	X += pathlqg[0].y;
	for(int i =1; i < (int)pathlqg.size(); i++){
		Matrix<2> pos = X.subMatrix<2,1>(0,0);
		Matrix<2,2> cov = pathlqg[i-1].R.subMatrix<2,2>(0,0);
		
		std::vector<std::pair<Matrix<2,1>, double>> cvx;
		std::vector<std::pair<Matrix<2,1>, double>> cvxprob;
		cvx.clear(); cvxprob.clear();
		query(pos, cov, cvx, cal_obstacles, cal_environment, cal_point, cvxprob);
		//double ps_i = computeStatics(pos, cov, cvxprob); 
		double ps_i = computeStatics(pos, cov, cvx);
		ps *= ps_i;

		lqgmpTruncation(X, pathlqg[i-1].R, cvx);
		Matrix<6> pathx = zeros<6,1>(); pathx.insert<3,1>(0,0, pathlqg[i-1].T); pathx.insert<3,1>(3,0, pathlqg[i-1].T);
		pathlqg[i-1].y = X - pathx;

		pathlqg[i].y = pathlqg[i-1].F * pathlqg[i-1].y;
		pathlqg[i].R = pathlqg[i-1].F * pathlqg[i-1].R * ~pathlqg[i-1].F + pathlqg[i-1].G * Q * ~pathlqg[i-1].G;

		Matrix<6> nextpathx  = zeros<6,1>(); nextpathx.insert<3,1>(0,0, pathlqg[i].T); nextpathx.insert<3,1>(3,0, pathlqg[i].T);
		X = pathlqg[i].y + nextpathx;
	}

	return ps;
}

void LQGMP::draw_truncate_distribution(const int& cal_ellipse_trunc)
{
	for(int i = 0; i < (int)pathlqg.size(); i+=3){
		Matrix<2> pathx = pathlqg[i].T.subMatrix<2,1>(0,0);
		Matrix<2> pos = pathx + pathlqg[i].y.subMatrix<2,1>(0,0);
		Matrix<2,2> Cov = pathlqg[i].R.subMatrix<2,2>(0,0);
		drawEllipse2d(pos, Cov, cal_ellipse_trunc, false);
	}

	for(int i = 0; i < (int)Prim.points.size(); i++){
		std::cout<<i<<": ("<<Prim.points[i].first[0]<<" "<<Prim.points[i].first[1]<<")"<<std::endl;	
	}
}