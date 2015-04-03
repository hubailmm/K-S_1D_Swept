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

#include "DiamondTop.h"

DiamondTop::DiamondTop(int foundationLength,int dataPointSize,int substeps,int type,int Lrank,int Rrank,double *foundation,double *staging,double *localGhost,double *remoteGhost,PointStruct *points,int bufferMode) 
	   :DiamondHalf(foundationLength,dataPointSize,substeps,type,Lrank,Rrank,foundation,staging,localGhost,bufferMode)
{
	if(this->bufferMode == 1)
		reqs = (MPI_Request*) malloc(1 * sizeof(MPI_Request));
	else if(this->bufferMode == 2)
		reqs = (MPI_Request*) malloc(2 * sizeof(MPI_Request));		
	else
		reqs = (MPI_Request*) malloc(levels * sizeof(MPI_Request));	

	needReqFree = false;
	
	this->remoteGhost = remoteGhost;
	this->points      = points;
}

inline void DiamondTop::prepareInputPtrs(int i,int j,PointStruct *point)
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
		
	int elementIndex = ((j)*dataPointSize*(substeps)) + ((i)*dataPointSize*(substeps));
	
	point->ulInput = &sourceData[elementIndex-(1*dataPointSize*(substeps))];
	point->uInput  = &sourceData[elementIndex];
	point->urInput = &sourceData[elementIndex+(1*dataPointSize*(substeps))];
	point->output  = &destData[elementIndex];

}

void DiamondTop::initGhostData()
{
	int Rindex = foundationLength*dataPointSize*(substeps)-(2*dataPointSize*(substeps));
	if(needReqFree)
	{
		freeRequests();		
	}

	if(type == 1)
	{
		if(this->bufferMode == 1 || this->bufferMode == 2)
		{
			for(int i=0;i<2*dataPointSize*substeps;i++)
			{
				remoteGhost[i] = foundation[Rindex+i];
			}
		}
		else
		{
			MPI_Isend(&foundation[Rindex],2*dataPointSize*substeps,MPI_DOUBLE,Rrank,2,MPI_COMM_WORLD,&reqs[0]);
		}		
		for(int i=0;i<2*dataPointSize*substeps;i++)
		{
			localGhost[i] = foundation[i];
		}
	}
	else if(type == 2)
	{
		for(int i=0;i<2*dataPointSize*substeps;i++)
		{
			localGhost[i] = foundation[Rindex+i];
		}
		if(this->bufferMode == 1 || this->bufferMode == 2)
		{
			for(int i=0;i<2*dataPointSize*substeps;i++)
			{
				remoteGhost[i] = foundation[i];
			}
		}
		else
		{
			MPI_Isend(&foundation[0],2*dataPointSize*substeps,MPI_DOUBLE,Lrank,1,MPI_COMM_WORLD,&reqs[0]);
		}
	}
	//printf("Init Ghost Sent type=%d!\n",type);	
}

void DiamondTop::initPointStructs()
{
	int pointIndex = 0;
	for(int i=1;i<levels;i++)	
	{
		int j = ((levels-i)*2)-1;
		prepareInputPtrs(i,j,&points[pointIndex]);				
		pointIndex++;
	}
}
int DiamondTop::calculate(int executeFnc)
{
	int fncIndex = executeFnc;
	int pointIndex = 0;
	initGhostData();
	PointStruct point;
	double *output;
	double **ulInput;
	double **uInput;
	double **urInput;
	//Loop for each Level of this DiamondTop
	for(int i=1;i<levels;i++)	
	{
		point  = points[pointIndex];		
		pointIndex++;

		for(int j=((levels-i)*2)-1;j>=0;j--)
		{
			memcpy(point.output+(1*dataPointSize),point.uInput,(substeps-1) * dataPointSize * sizeof(double));
			SubSteps::executeStepFnc(fncIndex,&point);

			if(j==foundationLength-(i*2)-2)
			{
				if(type == 1)
				{
					int index = i*2*dataPointSize*substeps;
					if(this->bufferMode == 1)
					{
						for(int c=0;c<2*dataPointSize*substeps;c++)
						{
							remoteGhost[index] = point.output[c];							
							index++;
						}
						if(this->bufferMode == 2 && i+1 == levels/2)
						{
							MPI_Isend(remoteGhost,2*dataPointSize*substeps*levels/2,MPI_DOUBLE,Rrank,2,MPI_COMM_WORLD,&reqs[0]);
						}					
					}
					else if(this->bufferMode == 2)
					{
						for(int c=0;c<2*dataPointSize*substeps;c++)
						{
							remoteGhost[index] = point.output[c];							
							index++;
						}
						if(i+1 == levels/2)
						{
							MPI_Isend(remoteGhost,2*dataPointSize*substeps*levels/2,MPI_DOUBLE,Rrank,2,MPI_COMM_WORLD,&reqs[0]);
						}					
					}
					else
					{
						MPI_Isend(point.output,2*dataPointSize*substeps,MPI_DOUBLE,Rrank,2,MPI_COMM_WORLD,&reqs[i]);
					}					
					
					if(i==levels-1)
					{
						int index = i*2*dataPointSize*substeps;
						for(int c=0;c<2*dataPointSize*substeps;c++)
						{
							localGhost[index] = point.output[c];
							index++;
						}						
					}						
				}
				else if(type == 2)
				{
					int index = i*2*dataPointSize*substeps;
					for(int c=0;c<2*dataPointSize*substeps;c++)
					{
						localGhost[index] = point.output[c];
						index++;
					}
					if(i==levels-1)
					{
						int index = i*2*dataPointSize*substeps;
						if(this->bufferMode == 1)
						{
							for(int c=0;c<2*dataPointSize*substeps;c++)
							{
								remoteGhost[index] = point.output[c];;
								index++;
							}
							if(this->bufferMode == 2 && i+1 == levels/2)
							{
								MPI_Isend(remoteGhost,2*dataPointSize*substeps*levels/2,MPI_DOUBLE,Lrank,1,MPI_COMM_WORLD,&reqs[0]);
							}							
						}
						else if(this->bufferMode == 2)
						{
							for(int c=0;c<2*dataPointSize*substeps;c++)
							{
								remoteGhost[index] = point.output[c];;
								index++;
							}
							if(i+1 == levels/2)
							{
								MPI_Isend(remoteGhost,2*dataPointSize*substeps*levels/2,MPI_DOUBLE,Lrank,1,MPI_COMM_WORLD,&reqs[0]);
							}							
						}

						else
						{
							MPI_Isend(point.output,2*dataPointSize*substeps,MPI_DOUBLE,Lrank,1,MPI_COMM_WORLD,&reqs[i]);
						}						
					}
				}
			}
			else if(j == 0)
			{
				//point.output[0] = 111;
				if(type == 1)
				{
					int index = i*2*dataPointSize*substeps;
					for(int c=0;c<2*dataPointSize*substeps;c++)
					{
						localGhost[index] = point.output[c];
						index++;
					}			
				}
				else if(type == 2)
				{
					if(this->bufferMode == 1)
					{
						int index = i*2*dataPointSize*substeps;
						for(int c=0;c<2*dataPointSize*substeps;c++)
						{
							remoteGhost[index] = point.output[c];
							index++;
						}						
					}
					else if(this->bufferMode == 2)
					{
						int index = i*2*dataPointSize*substeps;
						for(int c=0;c<2*dataPointSize*substeps;c++)
						{
							remoteGhost[index] = point.output[c];
							index++;
						}
						if(i+1 == levels/2)
						{
							MPI_Isend(remoteGhost,2*dataPointSize*substeps*levels/2,MPI_DOUBLE,Lrank,1,MPI_COMM_WORLD,&reqs[0]);
						}	
					}
					else
					{
						MPI_Isend(point.output,2*dataPointSize*substeps,MPI_DOUBLE,Lrank,1,MPI_COMM_WORLD,&reqs[i]);	
					}
				}
			}

			point.output  -= dataPointSize*substeps;
			point.ulInput -= dataPointSize*substeps;
			point.uInput  -= dataPointSize*substeps;
			point.urInput -= dataPointSize*substeps;

		}
		fncIndex++;
		if(fncIndex==substeps)fncIndex=0;
	}	
	needReqFree = true;
	if(this->bufferMode == 1)
	{
		if(this->type == 1)
			MPI_Isend(remoteGhost,2*dataPointSize*substeps*levels,MPI_DOUBLE,Rrank,2,MPI_COMM_WORLD,&reqs[0]);	
		else
			MPI_Isend(remoteGhost,2*dataPointSize*substeps*levels,MPI_DOUBLE,Lrank,1,MPI_COMM_WORLD,&reqs[0]);
	}
	else if(this->bufferMode == 2)
	{
		if(this->type == 1)
			MPI_Isend(remoteGhost+2*dataPointSize*substeps*levels/2,2*dataPointSize*substeps*levels/2,MPI_DOUBLE,Rrank,2,MPI_COMM_WORLD,&reqs[1]);	
		else
			MPI_Isend(remoteGhost+2*dataPointSize*substeps*levels/2,2*dataPointSize*substeps*levels/2,MPI_DOUBLE,Lrank,1,MPI_COMM_WORLD,&reqs[1]);
	}
	//printf("TOP type %d returning %d total execuions: %d\n",type,fncIndex,totalExecuted);
	return fncIndex;
}


void DiamondTop::freeRequests()
{
	if(this->bufferMode == 1)	
		MPI_Request_free(&reqs[0]);
	else if(this->bufferMode == 2)
	{
		MPI_Request_free(&reqs[0]);
		MPI_Request_free(&reqs[1]);
	}
	else
		for(int i=0;i<levels;i++)
			MPI_Request_free(&reqs[i]);
	//printf("All MPI_Requests Were Cleared!!\n");
}

DiamondTop::~DiamondTop()
{	
	free(reqs);	
}
