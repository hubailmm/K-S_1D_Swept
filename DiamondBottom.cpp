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

#include "DiamondBottom.h"

DiamondBottom::DiamondBottom(int foundationLength,int dataPointSize,int substeps,int type,int Lrank,int Rrank,double *foundation,double *staging,double *localGhost,double *remoteGhost,PointStruct *points,int bufferMode) 
	      : DiamondHalf(foundationLength,dataPointSize,substeps,type,Lrank,Rrank,foundation,staging,localGhost,bufferMode)
{
	if(type == 1)
	{
		this->rightGhost = localGhost;
		this->leftGhost  = remoteGhost;
	}
	else if(type == 2)
	{
		this->rightGhost = remoteGhost;
		this->leftGhost  = localGhost;
	}

	this->points = points;

}

inline void DiamondBottom::prepareInputPtrs(int i,int j,PointStruct *point)
{
	double *sourceData;
	double *destData;
		
	if(i%2 != 0)
	{
		sourceData = foundation;
		destData   = staging;
	}
	else
	{
		sourceData = staging;
		destData   = foundation;
	}

	int dim =0;
	int index,level = i-1;
	int elementIndex = ((levels-i)*dataPointSize*substeps)+((j)*dataPointSize*substeps);
	if(i == levels) 
	{
		if(j == 0)
		{
			point->ulInput = &this->leftGhost[(foundationLength * dataPointSize * substeps) - (2 * dataPointSize * substeps)];
			point->uInput  = &sourceData[elementIndex];
			point->urInput = &sourceData[elementIndex+(dataPointSize*substeps)];
		}
		else if(j == levels*2-1)
		{
			point->ulInput = &sourceData[elementIndex-(dataPointSize*substeps)];
			point->uInput  = &sourceData[elementIndex];
			point->urInput = &this->rightGhost[(foundationLength * dataPointSize * substeps) - (dataPointSize * substeps)];
		}
		else
		{
			point->ulInput = &sourceData[elementIndex-(dataPointSize*substeps)];
			point->uInput  = &sourceData[elementIndex];
			point->urInput = &sourceData[elementIndex+(dataPointSize*substeps)];
		}
	}
	else 
	{
		point->ulInput = &sourceData[elementIndex-(dataPointSize*substeps)];
		point->uInput  = &sourceData[elementIndex];
		point->urInput = &sourceData[elementIndex+(dataPointSize*substeps)];
		
	}
	point->output = &destData[elementIndex];

}

inline void DiamondBottom::getGhostValues(int level)
{
	if(type == 1)
	{
		int index = (level-1)*dataPointSize*2*substeps;
		MPI_Recv(&this->leftGhost[index],2*dataPointSize*substeps,MPI_DOUBLE,Lrank,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		
	}
	else if(type == 2)
	{
		int index = (level-1)*dataPointSize*2*substeps;
		MPI_Recv(&this->rightGhost[index],2*dataPointSize*substeps,MPI_DOUBLE,Rrank,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		
	}	
}

void DiamondBottom::adjustGhostValues(int level)
{
	int index = (level-1)*dataPointSize*2*substeps;
	double *destData;
	double *leftLocation;
	double *rightLocation;
		
	if(level%2 != 0)
	{
		destData   = foundation;
	}
	else
	{
		destData   = staging;
	}
	leftLocation  = destData;
	rightLocation = destData + (this->levels * dataPointSize * substeps);
	if(level != levels)
	{
		leftLocation  += (levels-level-1) * dataPointSize * substeps;		
		rightLocation += (level-1) * dataPointSize * substeps;
		memcpy(leftLocation ,&leftGhost[index] ,2 * dataPointSize * substeps * sizeof(double));
		memcpy(rightLocation,&rightGhost[index],2 * dataPointSize * substeps * sizeof(double));
	}
	else
	{
		rightLocation += (levels-1) * dataPointSize * substeps;
		memcpy(leftLocation ,&leftGhost[(foundationLength * dataPointSize * substeps) - (dataPointSize * substeps)] ,1 * dataPointSize * substeps * sizeof(double));
		memcpy(rightLocation,&rightGhost[(foundationLength * dataPointSize * substeps) - (2 * dataPointSize * substeps)],1 * dataPointSize * substeps * sizeof(double));
	}

}

void DiamondBottom::initPointStructs()
{
	int pointIndex = 0;

	for(int i=1;i<this->levels;i++)
	{
		this->prepareInputPtrs(i,0,&points[pointIndex]);				
		pointIndex++;
	}
	this->prepareInputPtrs(levels,1,&points[pointIndex]);				
	pointIndex++;
	this->prepareInputPtrs(levels,0,&points[pointIndex]);				
	pointIndex++;	
	this->prepareInputPtrs(levels,(levels*2)-1,&points[pointIndex]);				
	pointIndex++;
}

int DiamondBottom::calculate(int executeFnc)
{
	int fncIndex = executeFnc;
	int totalExecuted =0;
	if(this->bufferMode == 1)
	{
		double start = MPI_Wtime();
		std::clock_t startTime = std::clock();
		if(type == 1)
			MPI_Recv(leftGhost,2*dataPointSize*substeps*levels,MPI_DOUBLE,Lrank,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		else if(type == 2)
			MPI_Recv(rightGhost,2*dataPointSize*substeps*levels,MPI_DOUBLE,Rrank,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		
		std::clock_t endTime = std::clock();
		double totalTime = (endTime - startTime) / (double)CLOCKS_PER_SEC;
		double ellapsed = MPI_Wtime() - start;		
	}
	else if(this->bufferMode == 2)
	{
		double start = MPI_Wtime();
		std::clock_t startTime = std::clock();
		if(type == 1)
			MPI_Recv(leftGhost,2*dataPointSize*substeps*levels/2,MPI_DOUBLE,Lrank,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		else if(type == 2)
			MPI_Recv(rightGhost,2*dataPointSize*substeps*levels/2,MPI_DOUBLE,Rrank,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		
		std::clock_t endTime = std::clock();
		double totalTime = (endTime - startTime) / (double)CLOCKS_PER_SEC;
		double ellapsed = MPI_Wtime() - start;
	}
	int pointIndex = 0;
	PointStruct point;
	
	for(int i=1;i<this->levels;i++)
	{
		point  = points[pointIndex];
		pointIndex++;
		
		if(this->bufferMode != 1 && this->bufferMode != 2)
		{
			double start = MPI_Wtime();
			this->getGhostValues(i);
			double ellapsed = MPI_Wtime() - start;
	                
		}
		else if(this->bufferMode == 2)
		{
			if(i==levels/2)
			{
				double start = MPI_Wtime();
				std::clock_t startTime = std::clock();
				if(type == 1)
					MPI_Recv(leftGhost+(2*dataPointSize*substeps*levels/2),2*dataPointSize*substeps*levels/2,MPI_DOUBLE,Lrank,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);					
				else if(type == 2)
					MPI_Recv(rightGhost+(2*dataPointSize*substeps*levels/2),2*dataPointSize*substeps*levels/2,MPI_DOUBLE,Rrank,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		
				std::clock_t endTime = std::clock();
				double totalTime = (endTime - startTime) / (double)CLOCKS_PER_SEC;
				double ellapsed = MPI_Wtime() - start;
			}
		}

		adjustGhostValues(i);
		for(int j=0;j<i*2;j++)
		{
			memcpy(point.output+(1*dataPointSize),point.uInput,(substeps-1) * dataPointSize * sizeof(double));
			SubSteps::executeStepFnc(fncIndex,&point);
			point.output  += dataPointSize*substeps;
			point.ulInput += dataPointSize*substeps;
			point.uInput  += dataPointSize*substeps;
			point.urInput += dataPointSize*substeps;
			
		}		
		fncIndex++;
		if(fncIndex==substeps)fncIndex=0;		
	}
	if(this->bufferMode != 1 && this->bufferMode != 2)
	{
			double start = MPI_Wtime();
			this->getGhostValues(levels);
			double ellapsed = MPI_Wtime() - start;
	                
	}
	adjustGhostValues(levels);
	point = points[pointIndex];
	pointIndex++;
	for(int j=1;j<levels*2-1;j++)
	{
		memcpy(point.output+(1*dataPointSize),point.uInput,(substeps-1) * dataPointSize * sizeof(double));
		SubSteps::executeStepFnc(fncIndex,&point);
		point.output  += dataPointSize*substeps;
		point.ulInput += dataPointSize*substeps;
		point.uInput  += dataPointSize*substeps;
		point.urInput += dataPointSize*substeps;
	}
	point = points[pointIndex];
	if(this->type == 2)
	{
		//A fix for DiamondBottom Type 2:
		point.ulInput = &this->leftGhost[(foundationLength * dataPointSize * substeps) - (2 * dataPointSize * substeps)];
		//Fix End
	}
	pointIndex++;
	memcpy(point.output+(1*dataPointSize),point.uInput,(substeps-1) * dataPointSize * sizeof(double));
	SubSteps::executeStepFnc(fncIndex,&point);
	point = points[pointIndex];
	if(this->type == 2)
	{
		//A fix for DiamondBottom Type 2:
		point.urInput = &this->rightGhost[(foundationLength * dataPointSize * substeps) - (1 * dataPointSize * substeps)];
		//Fix End
	}
	memcpy(point.output+(1*dataPointSize),point.uInput,(substeps-1) * dataPointSize * sizeof(double));
	SubSteps::executeStepFnc(fncIndex,&point);
	fncIndex++;
	if(fncIndex==substeps)fncIndex=0;		
	return fncIndex;
}

DiamondBottom::~DiamondBottom()
{
}
