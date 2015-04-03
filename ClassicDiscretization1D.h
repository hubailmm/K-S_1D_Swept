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

#ifndef H_CLASSICDISCRETIZATION1D
#define H_CLASSICDISCRETIZATION1D

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string.h>
#include "mpi.h"
#include "DiamondGlobals.h"
#include "Substeps.h"

class ClassicDiscretization1D
{
private:		
	int substeps;
	int foundationLength;
	int dataPointSize;
	int Lrank;
	int Rrank;
	double *foundation;
	double *staging;
	PointStruct *points;
	MPI_Request *reqs;
	void initPointStructs();
	void initGhostData();

public:
	ClassicDiscretization1D(int foundationLength,int dataPointSize,int substeps,int Lrank,int Rrank);
	void setInitialValue(double *initData);	
	void writeOutput(char *filename);	
	void calculate(int substepsToPerform);
	void printOutput();
	~ClassicDiscretization1D();

};

#endif

