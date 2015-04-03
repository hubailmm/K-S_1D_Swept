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

#ifndef H_DIAMONDHALF
#define H_DIAMONDHALF

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string.h>
#include "mpi.h"
#include "DiamondGlobals.h"
#include "Substeps.h"

using namespace std;

typedef void (*fp)(double**,double**,double**,double*);
class DiamondHalf
{

protected:
	int foundationLength;
	int dataPointSize;
	int diamondType;
	int levels;
	int totalDataPoints;
	int Lrank;
	int Rrank;
	int type;
	int substeps;
	double *foundation;
	double *staging;
	double *localGhost;
	PointStruct *points;
	int bufferMode;
	vector< void (*)(PointStruct *)> stepFunctions;
	//fp fncs[10];

public:
	DiamondHalf(int foundationLength,int dataPointSize,int substeps,int type,int Lrank,int Rrank,double *foundation,double *staging,double *localGhost,int bufferMode);	
	~DiamondHalf();
};

#endif
