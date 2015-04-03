// This file is part of Diamond1D
// Copyright (C) 2015 Qiqi Wang, qiqi@mit.edu AND Maitham Alhubail, hubailmm@mit.edu
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef H_DIAMONDISCRETIZATION1D
#define H_DIAMONDISCRETIZATION1D

#include "DiamondTop.h"
#include "DiamondBottom.h"

class DiamondDiscretization1D
{
private:
	DiamondTop    *dt1;
	DiamondTop    *dt2;
	DiamondBottom *db1;
	DiamondBottom *db2;	
	int executeFnc;
	int substeps;
	int levels;
	int foundationLength;
	int dataPointSize;
	double *foundation;
	double *staging;
	double *localGhost;
	double *remoteGhost;
	double *memory;
	PointStruct *topPoints;
	PointStruct *bottomPoints;
	int myrank;
	int Lrank;
	int Rrank;
	int bufferMode;
	MPI_Request *reqs,*recvReq;
	bool needReqFree;
	

public:
	DiamondDiscretization1D(int foundationLength,int dataPointSize,int substeps,int Lrank,int myrank,int Rrank,int bufferMode);
	void setInitialValue(double *initData);
	void writeOutput(char *filename);
	void getValues(double *values);	
	void printOutput();
	void registerStepFunction(void fnc(PointStruct *));
	void calculate(int cycles);
	void calculateInlined(int cycles);

	~DiamondDiscretization1D();
};

#endif 
