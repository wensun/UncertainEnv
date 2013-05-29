#ifndef _PRIMITIVE_
#define _PRIMITIVE_

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

class Primitive{

public:
	int num;
	
	struct segments{
		Matrix<2> a;
		double b;
		Matrix<2> p1;
		Matrix<2,2> S1;
		Matrix<2> p2;
		Matrix<2,2> S2;
		Matrix<2,2> RotationM;

		segments(){
			a = zeros<2,1>();
			b = 0;
			p1 = zeros<2,1>();
			p2 = zeros<2,1>();
			S1 =zeros<2,2>();
			S2 = zeros<2,2>();
			RotationM = zeros<2,2>();
		}
	};

	std::vector<std::pair<Matrix<2>, Matrix<2,2>>> points;   //pair: mean and Covariance.
	std::vector<segments> seg;
	std::vector<Matrix<2>> randomvertex;
	std::vector<segments> randomseg;

	Primitive(){
		num = 11;
		points.resize(num);
		seg.clear();
		randomseg.clear();
		randomvertex.clear();

		Matrix<2> p1 = zeros<2,1>(); p1(0,0) = 0; p1(1,0) = 6;
		Matrix<2,2> S1 = 0.03 * identity<2>();
		points[0] = std::make_pair(p1, S1);

		Matrix<2> p2 = zeros<2,1>(); p2(0,0) = 5; p2(1,0) = 6;
		Matrix<2,2> S2 = 0.03 * identity<2>();
		points[1] = std::make_pair(p2, S2);

		Matrix<2> p3 = zeros<2,1>(); p3(0,0) = 5; p3(1,0) = 0;
		Matrix<2,2> S3 = 0.01 * identity<2>();
		points[2] = std::make_pair(p3, S3);

		Matrix<2> p4 = zeros<2,1>(); p4(0,0) = 0; p4(1,0) = 0;
		Matrix<2,2> S4 = 0.01 * identity<2>();
		points[3] = std::make_pair(p4, S4);

		Matrix<2> p5 = zeros<2,1>(); p5(0,0) = 1; p5(1,0) = 5;
		Matrix<2,2> S5 = 0.005 * identity<2>();
		points[4] = std::make_pair(p5, S5);
		
		Matrix<2> p6 = zeros<2,1>(); p6(0,0) = 3; p6(1,0) = 5;
		Matrix<2,2> S6 = 0.005 * identity<2>();
		points[5] = std::make_pair(p6, S6);

		Matrix<2> p7 = zeros<2,1>(); p7(0,0) = 3; p7(1,0) = 1;
		Matrix<2,2> S7 = 0.01 * identity<2>();
		points[6] = std::make_pair(p7, S7);

		Matrix<2> p8 = zeros<2,1>(); p8(0,0) = 1; p8(1,0) = 1;
		Matrix<2,2> S8 = 0.01 * identity<2>();
		points[7] = std::make_pair(p8, S8);

		Matrix<2> p9 = zeros<2,1>(); p9(0,0) = 4; p9(1,0) = 2;
		Matrix<2,2> S9 = 0.02 * identity<2>();
		points[8] = std::make_pair(p9, S9);

		Matrix<2> p10 = zeros<2,1>(); p10(0,0) = 5; p10(1,0) = 2;
		Matrix<2,2> S10 = 0.001 * identity<2>();
		S10(0,0) = 0; //no uncertainty in the x direction.
		points[9] = std::make_pair(p10, S10);

		Matrix<2> p11 = zeros<2,1>(); p11(0,0) = 4; p11(1,0) = 0;
		Matrix<2,2> S11 = 0.001 * identity<2>();
		S11(1,1) = 0;
		points[10] = std::make_pair(p11, S11);

	}
	
	void CreateVertex()
	{
		double trunc = 0.5;
		randomvertex.resize((int)points.size());
		for(int i = 0; i < (int)points.size(); i++){
			//implement truncated gaussians. 
			while(true){
				randomvertex[i] = sampleGaussian(points[i].first, points[i].second);
				if(tr(~(randomvertex[i] - points[i].first) * (randomvertex[i] - points[i].first)) < trunc)
					break;
			}
		}

	}

	void CreateRandomSegs()
	{
		Matrix<2,2> rotation = zeros<2,2>();
		rotation(0, 1) = 1;
		rotation(1,0) = -1;

		std::vector<Matrix<2>> vertex;
		vertex.clear();
		for(int i = 0; i < (int)randomvertex.size(); i++)
			vertex.push_back(randomvertex[i]);

		double len = 0;

		Matrix<2> a = rotation * (vertex[0] - vertex[1]);
		double b = tr(~a * vertex[0]);
		len = sqrt(tr(~a*a));
		segments line;
		line.a = a/len; line.b = b/len; line.p1 = vertex[0]; line.p2 = vertex[1];
		randomseg.push_back(line);

		a = rotation * (vertex[1] - vertex[2]);
		b = tr(~a * vertex[1]);
		len = sqrt(tr(~a*a));
		line.a = a/len; line.b = b/len; line.p1 = vertex[1]; line.p2 = vertex[2];
		randomseg.push_back(line);

		a = rotation * (vertex[2] - vertex[3]);
		b = tr(~a * vertex[2]);
		len = sqrt(tr(~a*a));
		line.a = a/len; line.b = b/len; line.p1 = vertex[2]; line.p2 = vertex[3];
		randomseg.push_back(line);
		
		a = rotation * (vertex[3] - vertex[0]);
		b = tr(~a * vertex[3]);
		len = sqrt(tr(~a*a));
		line.a = a/len; line.b = b/len; line.p1 = vertex[3]; line.p2 = vertex[0];
		randomseg.push_back(line);


		rotation(0,1) = -1;
		rotation(1,0) = 1;

		a = rotation * (vertex[4] - vertex[5]);
		b = tr(~a * vertex[4]);
		len = sqrt(tr(~a*a));
		line.a = a/len; line.b = b/len; line.p1 = vertex[4]; line.p2 = vertex[5];
		randomseg.push_back(line);

		a = rotation * (vertex[5] - vertex[6]);
		b = tr(~a * vertex[5]);
		len = sqrt(tr(~a*a));
		line.a = a/len; line.b = b/len; line.p1 = vertex[5]; line.p2 = vertex[6];
		randomseg.push_back(line);

		a = rotation * (vertex[6] - vertex[7]);
		b = tr(~a * vertex[6]);
		len = sqrt(tr(~a*a));
		line.a = a/len; line.b = b/len; line.p1 = vertex[6]; line.p2 = vertex[7];
		randomseg.push_back(line);

		a = rotation * (vertex[7] - vertex[4]);
		b = tr(~a * vertex[7]);
		len = sqrt(tr(~a*a));
		line.a = a/len; line.b = b/len; line.p1 = vertex[7]; line.p2 = vertex[4];
		randomseg.push_back(line);

		a = rotation * (vertex[8] - vertex[9]);
		b = tr(~a * vertex[8]);
		len = sqrt(tr(~a*a));
		line.a = a/len; line.b = b/len; line.p1 = vertex[8]; line.p2 = vertex[9];
		randomseg.push_back(line);

		a = rotation * (vertex[10] - vertex[8]);
		b = tr(~a * vertex[8]);
		len = sqrt(tr(~a*a));
		line.a = a/len; line.b = b/len; line.p1 = vertex[10]; line.p2 = vertex[8];
		randomseg.push_back(line);
	}

	void CreateSegments()
	{
		//CreateVertex();
		Matrix<2,2> rotation = zeros<2,2>();
		rotation(0, 1) = 1;
		rotation(1,0) = -1;

		std::vector<Matrix<2>> vertex;
		vertex.clear();
		for(int i = 0; i < (int)points.size(); i++)
			vertex.push_back(points[i].first);

		double len = 0;

		Matrix<2> a = rotation * (vertex[0] - vertex[1]);
		double b = tr(~a * vertex[0]);
		len = sqrt(tr(~a*a));
		segments line;
		line.a = a/len; line.b = b/len; line.p1 = vertex[0]; line.p2 = vertex[1];
		line.S1 = points[0].second; line.S2 = points[1].second;
		line.RotationM = rotation;
		seg.push_back(line);

		a = rotation * (vertex[1] - vertex[2]);
		b = tr(~a * vertex[1]);
		len = sqrt(tr(~a*a));
		line.a = a/len; line.b = b/len; line.p1 = vertex[1]; line.p2 = vertex[2];
		line.S1 = points[1].second; line.S2 = points[2].second;
		line.RotationM = rotation;
		seg.push_back(line);

		a = rotation * (vertex[2] - vertex[3]);
		b = tr(~a * vertex[2]);
		len = sqrt(tr(~a*a));
		line.a = a/len; line.b = b/len; line.p1 = vertex[2]; line.p2 = vertex[3];
		line.S1 = points[2].second; line.S2 = points[3].second;
		line.RotationM = rotation;
		seg.push_back(line);
		
		a = rotation * (vertex[3] - vertex[0]);
		b = tr(~a * vertex[3]);
		len = sqrt(tr(~a*a));
		line.a = a/len; line.b = b/len; line.p1 = vertex[3]; line.p2 = vertex[0];
		line.S1 = points[3].second; line.S2 = points[0].second;
		line.RotationM = rotation;
		seg.push_back(line);


		rotation(0,1) = -1;
		rotation(1,0) = 1;

		a = rotation * (vertex[4] - vertex[5]);
		b = tr(~a * vertex[4]);
		len = sqrt(tr(~a*a));
		line.a = a/len; line.b = b/len; line.p1 = vertex[4]; line.p2 = vertex[5];
		line.S1 = points[4].second; line.S2 = points[5].second;
		line.RotationM = rotation;
		seg.push_back(line);

		a = rotation * (vertex[5] - vertex[6]);
		b = tr(~a * vertex[5]);
		len = sqrt(tr(~a*a));
		line.a = a/len; line.b = b/len; line.p1 = vertex[5]; line.p2 = vertex[6];
		line.S1 = points[5].second; line.S2 = points[6].second;
		line.RotationM = rotation;
		seg.push_back(line);

		a = rotation * (vertex[6] - vertex[7]);
		b = tr(~a * vertex[6]);
		len = sqrt(tr(~a*a));
		line.a = a/len; line.b = b/len; line.p1 = vertex[6]; line.p2 = vertex[7];
		line.S1 = points[6].second; line.S2 = points[7].second;
		line.RotationM = rotation;
		seg.push_back(line);

		a = rotation * (vertex[7] - vertex[4]);
		b = tr(~a * vertex[7]);
		len = sqrt(tr(~a*a));
		line.a = a/len; line.b = b/len; line.p1 = vertex[7]; line.p2 = vertex[4];
		line.S1 = points[7].second; line.S2 = points[4].second;
		line.RotationM = rotation;
		seg.push_back(line);

		a = rotation * (vertex[8] - vertex[9]);
		b = tr(~a * vertex[8]);
		len = sqrt(tr(~a*a));
		line.a = a/len; line.b = b/len; line.p1 = vertex[8]; line.p2 = vertex[9];
		line.S1 = points[8].second; line.S2 = points[9].second;
		line.RotationM = rotation;
		seg.push_back(line);

		a = rotation * (vertex[10] - vertex[8]);
		b = tr(~a * vertex[8]);
		len = sqrt(tr(~a*a));
		line.a = a/len; line.b = b/len; line.p1 = vertex[10]; line.p2 = vertex[8];
		line.S1 = points[10].second; line.S2 = points[8].second;
		line.RotationM = rotation;
		seg.push_back(line);

	}

	bool checkCollision(const Matrix<2>& pos){
		Matrix<2> x = pos;
		for(int i = 0; i < 4; i++){
			Matrix<2> a = randomseg[i].a;
			double b = randomseg[i].b;
			if(tr(~a*x) > b)
				return false;
		}

		bool flag = true;
		for(int i = 4; i < 8; i++){
			Matrix<2> a = randomseg[i].a;
			double b = randomseg[i].b;
			if(tr(~a*x) <= b)
				flag = false;
		}
		if(flag == true)
			return false;

		flag = true;
		for(int i = 8; i < 10; i++){
			Matrix<2> a = randomseg[i].a;
			double b = randomseg[i].b;
			if(tr(~a*x) <= b)
				flag = false;
		}
		if(flag == true)
			return false;

		return true;
	}

	bool checkline(const Matrix<2>& p1, const Matrix<2>& p2){
		for(int i = 0; i <= 10; i++){
			double alpha = i*1.0 / 10;
			Matrix<2> tmpp = alpha * p1 + (1 - alpha)*p2;
			bool tmpflag = checkCollision(tmpp);
			if(tmpflag == false)
				return false;
		}
		return true;
	}

};


#endif