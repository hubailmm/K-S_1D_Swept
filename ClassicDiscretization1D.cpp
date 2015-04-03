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

#include "ClassicDiscretization1D.h"

extern int myrank;
ClassicDiscretization1D::ClassicDiscretization1D(int foundationLength,int dataPointSize,int substeps,int Lrank,int Rrank)
{
	foundation    = new double[(foundationLength+2) * dataPointSize * substeps];
	staging       = new double[(foundationLength+2) * dataPointSize * substeps];

	points     = (PointStruct*) malloc(6 * sizeof(PointStruct));
	reqs = (MPI_Request*) malloc(2 * sizeof(MPI_Request));
	this->substeps = substeps;
	this->dataPointSize = dataPointSize;
	this->foundationLength = foundationLength;
	this->Lrank = Lrank;
	this->Rrank = Rrank;	
}

void ClassicDiscretization1D::setInitialValue(double *initData)
{
	int index = dataPointSize*substeps;
	for(int i=0;i<foundationLength;i++)
	{
		foundation[index] = initData[i];
		
		for(int s=0;s<substeps;s++)
		{
			index++;			
			foundation[index] = 0;
		}
	}
	this->initPointStructs();
	this->initGhostData();
}

void ClassicDiscretization1D::printOutput()
{
	int index = 0;
	for(int i=0;i<foundationLength;i++)
	{
		printf("u[%d] = %f\n",i+(myrank*foundationLength),foundation[index+(dataPointSize*substeps)]);
		//index+=3;
		if(i != foundationLength - 1)
		for(int s=0;s<substeps;s++)
		{
			index++;			
		}
	}	
}
void ClassicDiscretization1D::initPointStructs()
{
	points[0].output  = &staging[1*dataPointSize*substeps];
	points[0].ulInput = &foundation[0];
	points[0].uInput = &foundation[1*dataPointSize*substeps];
	points[0].urInput = &foundation[2*dataPointSize*substeps];

	points[1].output  = &staging[foundationLength*dataPointSize*substeps];
	points[1].ulInput = &foundation[(foundationLength-1)*dataPointSize*substeps];
	points[1].uInput = &foundation[foundationLength*dataPointSize*substeps];
	points[1].urInput = &foundation[(foundationLength+1)*dataPointSize*substeps];

	points[2].output  = &staging[2*dataPointSize*substeps];
	points[2].ulInput = &foundation[1*dataPointSize*substeps];
	points[2].uInput = &foundation[2*dataPointSize*substeps];
	points[2].urInput = &foundation[3*dataPointSize*substeps];

	points[3].output  = &foundation[1*dataPointSize*substeps];
	points[3].ulInput = &staging[0];
	points[3].uInput = &staging[1*dataPointSize*substeps];
	points[3].urInput = &staging[2*dataPointSize*substeps];

	points[4].output  = &foundation[foundationLength*dataPointSize*substeps];
	points[4].ulInput = &staging[(foundationLength-1)*dataPointSize*substeps];
	points[4].uInput = &staging[foundationLength*dataPointSize*substeps];
	points[4].urInput = &staging[(foundationLength+1)*dataPointSize*substeps];

	points[5].output  = &foundation[2*dataPointSize*substeps];
	points[5].ulInput = &staging[1*dataPointSize*substeps];
	points[5].uInput = &staging[2*dataPointSize*substeps];
	points[5].urInput = &staging[3*dataPointSize*substeps];

}
void ClassicDiscretization1D::initGhostData()
{
	MPI_Isend(&foundation[1*dataPointSize*substeps],dataPointSize*substeps,MPI_DOUBLE,Lrank,1,MPI_COMM_WORLD,&reqs[0]);
	MPI_Isend(&foundation[foundationLength*dataPointSize*substeps],dataPointSize*substeps,MPI_DOUBLE,Rrank,2,MPI_COMM_WORLD,&reqs[1]);
	MPI_Recv(&foundation[0],dataPointSize*substeps,MPI_DOUBLE,Lrank,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Recv(&foundation[(foundationLength+1)*dataPointSize*substeps],dataPointSize*substeps,MPI_DOUBLE,Rrank,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
}
void ClassicDiscretization1D::calculate(int substepsToPerform)
{
	int fncIndex = 0;
	for(int i=1;i<=substepsToPerform;i++)
	{
		int indexShift = 0;
		if(i % 2 == 0)
		{
			indexShift = 3;
		}
		MPI_Request_free(&reqs[0]);
		MPI_Request_free(&reqs[1]);
		
		PointStruct point;
		point = points[indexShift+0];
		memcpy(point.output+(1*dataPointSize),point.uInput,(substeps-1) * dataPointSize * sizeof(double));
		SubSteps::executeStepFnc(fncIndex,&point);
		MPI_Isend(point.output,dataPointSize*substeps,MPI_DOUBLE,Lrank,1,MPI_COMM_WORLD,&reqs[0]);
		point = points[indexShift+1];
		memcpy(point.output+(1*dataPointSize),point.uInput,(substeps-1) * dataPointSize * sizeof(double));
		SubSteps::executeStepFnc(fncIndex,&point);
		MPI_Isend(point.output,dataPointSize*substeps,MPI_DOUBLE,Rrank,2,MPI_COMM_WORLD,&reqs[1]);
		
		point = points[indexShift+2];
		for(int j=1;j<foundationLength-1;j++)
		{
			memcpy(point.output+(1*dataPointSize),point.uInput,(substeps-1) * dataPointSize * sizeof(double));
			SubSteps::executeStepFnc(fncIndex,&point);
			point.output  += dataPointSize*substeps;
			point.ulInput += dataPointSize*substeps;
			point.uInput  += dataPointSize*substeps;
			point.urInput += dataPointSize*substeps;
		}
		if(i % 2 == 0)
		{
			MPI_Recv(&foundation[0],dataPointSize*substeps,MPI_DOUBLE,Lrank,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			MPI_Recv(&foundation[(foundationLength+1)*dataPointSize*substeps],dataPointSize*substeps,MPI_DOUBLE,Rrank,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		else
		{
			MPI_Recv(&staging[0],dataPointSize*substeps,MPI_DOUBLE,Lrank,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			MPI_Recv(&staging[(foundationLength+1)*dataPointSize*substeps],dataPointSize*substeps,MPI_DOUBLE,Rrank,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		
		fncIndex++;if(fncIndex == this->substeps) fncIndex = 0;
	}
}

ClassicDiscretization1D::~ClassicDiscretization1D()
{
	free(reqs);
	free(points);
	delete[] foundation;
    delete[] staging;
}

