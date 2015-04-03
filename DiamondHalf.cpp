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

#include "DiamondHalf.h"


DiamondHalf::DiamondHalf(int foundationLength,int dataPointSize,int substeps,int type,int Lrank,int Rrank,double *foundation,double *staging,double *localGhost,int bufferMode)
{
	this->foundationLength = foundationLength;
	this->dataPointSize    = dataPointSize;
	this->levels           = foundationLength/2;
	this->totalDataPoints  = 0;	
	this->Lrank            = Lrank;
	this->Rrank            = Rrank;
	this->substeps         = substeps;
	this->type             = type;
	this->localGhost       = localGhost;
	this->foundation       = foundation;
	this->staging          = staging;

	for(int i=foundationLength;i>=2;i=i-2)this->totalDataPoints += i;
	this->bufferMode = bufferMode;
	
}

DiamondHalf::~DiamondHalf()
{
	this->stepFunctions.clear();
}
