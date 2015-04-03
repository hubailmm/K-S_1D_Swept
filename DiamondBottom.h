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

#ifndef H_DIAMONDBOTTOM
#define H_DIAMONDBOTTOM

#include "DiamondHalf.h"
#include <ctime>
class DiamondBottom : public DiamondHalf
{
private:
	double *rightGhost;
	double *leftGhost;
	void prepareInputPtrs(int i,int j,PointStruct *point);
	void getGhostValues(int level);
	void adjustGhostValues(int level);

public:
	DiamondBottom(int foundationLength,int dataPointSize,int substeps,int type,int Lrank,int Rrank,double *foundation,double *staging,double *localGhost,double *remoteGhost,PointStruct *points,int bufferMode);	
	double *getFoundationPtr();
	void getFoundationValues(double *foundationValues);	
	void initPointStructs();
	int calculate(int executeFnc);
	~DiamondBottom();
};

#endif

