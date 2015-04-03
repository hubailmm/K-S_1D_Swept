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

#include "DiamondDiscretization1D.h"

DiamondDiscretization1D::DiamondDiscretization1D(int foundationLength,int dataPointSize,int substeps,int Lrank,int myrank,int Rrank,int bufferMode)
{
	localGhost    = new double[foundationLength * dataPointSize * substeps];
	remoteGhost   = new double[foundationLength * dataPointSize * substeps];
	foundation    = new double[foundationLength * dataPointSize * substeps];
	staging       = new double[foundationLength * dataPointSize * substeps];

	topPoints     = (PointStruct*) malloc(foundationLength/2 * sizeof(PointStruct));
	bottomPoints  = (PointStruct*) malloc(((foundationLength/2) + 2) * sizeof(PointStruct));

	dt1 = new DiamondTop(foundationLength,dataPointSize,substeps,1,myrank,Rrank,foundation,staging,localGhost,remoteGhost,topPoints,bufferMode);
	dt1->initPointStructs();
	dt2 = new DiamondTop(foundationLength,dataPointSize,substeps,2,Lrank,myrank,foundation,staging,localGhost,remoteGhost,topPoints,bufferMode);

	db1 = new DiamondBottom(foundationLength,dataPointSize,substeps,1,Lrank,myrank,foundation,staging,localGhost,remoteGhost,bottomPoints,bufferMode);
	db1->initPointStructs();
	db2 = new DiamondBottom(foundationLength,dataPointSize,substeps,2,myrank,Rrank,foundation,staging,localGhost,remoteGhost,bottomPoints,bufferMode);

	this->substeps = substeps;
	this->dataPointSize = dataPointSize;
	this->foundationLength = foundationLength;
	this->myrank = myrank;
	this->Lrank  = Lrank;
	this->Rrank  = Rrank;
	this->levels = foundationLength/2;
	this->bufferMode = bufferMode;
	executeFnc = 0;
	if(this->bufferMode == 1)
		reqs = (MPI_Request*) malloc(1 * sizeof(MPI_Request));
	else if(this->bufferMode == 2)
		reqs = (MPI_Request*) malloc(2 * sizeof(MPI_Request));		
	else
		reqs = (MPI_Request*) malloc(levels * sizeof(MPI_Request));

	recvReq  = (MPI_Request*) malloc(1 * sizeof(MPI_Request));	
	needReqFree = false;
}

void DiamondDiscretization1D::setInitialValue(double *initData)
{
	int index = 0;
	for(int i=0;i<foundationLength;i++)
	{
		foundation[index] = initData[i];
		if(i != foundationLength - 1)
		{
			for(int s=0;s<substeps;s++)
			{
				index++;			
				foundation[index] = 0;
			}
		}
		else
		{
			for(int s=0;s<substeps-1;s++)
			{
				index++;			
				foundation[index] = 0;
			}
		}

	}	
}

void DiamondDiscretization1D::printOutput()
{
	int index = 0;
	for(int i=0;i<foundationLength;i++)
	{
		printf("u[%d] = %f\n",i+(myrank*foundationLength),foundation[index]);
		if(i != foundationLength - 1)
		for(int s=0;s<substeps;s++)
		{
			index++;			
		}
	}	
}

void DiamondDiscretization1D::calculateInlined(int cycles)
{
	double *rightGhost,*leftGhost;	
	int fncIndex,executeFnc = 0;
	int Rindex,pointIndex = 0;
	PointStruct point;
	for(int cycle=1;cycle<=cycles;cycle++)
	{
		//DiamondTop Type 1:
		pointIndex = 0;
		fncIndex = executeFnc;
		Rindex = foundationLength*dataPointSize*(substeps)-(2*dataPointSize*(substeps));
		if(this->bufferMode == 1)
		{
			for(int i=0;i<2*dataPointSize*substeps;i++)
			{
				remoteGhost[i] = foundation[Rindex+i];
			}
		}
		else if(this->bufferMode == 2)
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
		//
		
		//Loop for each Level of this DiamondTop
		for(int i=1;i<levels;i++)	
		{
			point  = topPoints[pointIndex];		
			pointIndex++;

			for(int j=((levels-i)*2)-1;j>=0;j--)
			{
				memcpy(point.output+(1*dataPointSize),point.uInput,(substeps-1) * dataPointSize * sizeof(double));
				SubSteps::executeStepFnc(fncIndex,&point);
				if(j==foundationLength-(i*2)-2)
				{
					int index = i*2*dataPointSize*substeps;
					if(this->bufferMode == 1)
					{
						for(int c=0;c<2*dataPointSize*substeps;c++)
						{
							remoteGhost[index] = point.output[c];							
							index++;
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

						break;
					}							
				}
				else if(j == 0)
				{
					int index = i*2*dataPointSize*substeps;
					for(int c=0;c<2*dataPointSize*substeps;c++)
					{
						localGhost[index] = point.output[c];
						index++;
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
		
		if(this->bufferMode == 1)
		{
			MPI_Isend(remoteGhost,2*dataPointSize*substeps*levels,MPI_DOUBLE,Rrank,2,MPI_COMM_WORLD,&reqs[0]);			
		}
		else if(this->bufferMode == 2)
		{
			MPI_Isend(remoteGhost+2*dataPointSize*substeps*levels/2,2*dataPointSize*substeps*levels/2,MPI_DOUBLE,Rrank,2,MPI_COMM_WORLD,&reqs[1]);			
		}
		//DiamondTop Type 1 END
				
		//DiamondBottom Type 1:
		fncIndex = executeFnc;
		rightGhost = localGhost;
		leftGhost  = remoteGhost;			
		
		if(this->bufferMode == 1)
		{
			double start = MPI_Wtime();
			std::clock_t startTime = std::clock();
			MPI_Recv(leftGhost,2*dataPointSize*substeps*levels,MPI_DOUBLE,Lrank,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			std::clock_t endTime = std::clock();
			double totalTime = (endTime - startTime) / (double)CLOCKS_PER_SEC;
			double ellapsed = MPI_Wtime() - start;
			//printf("It took %f microseconds to get the ghost values!\n",ellapsed*1000000);
		}
		else if(this->bufferMode == 2)
		{
			double start = MPI_Wtime();
			std::clock_t startTime = std::clock();
			MPI_Recv(leftGhost,2*dataPointSize*substeps*levels/2,MPI_DOUBLE,Lrank,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			MPI_Irecv(leftGhost+2*dataPointSize*substeps*levels/2,2*dataPointSize*substeps*levels/2,MPI_DOUBLE,Lrank,2,MPI_COMM_WORLD,recvReq);
			std::clock_t endTime = std::clock();
			double totalTime = (endTime - startTime) / (double)CLOCKS_PER_SEC;
			double ellapsed = MPI_Wtime() - start;
		}
		pointIndex = 0;
		//PointStruct point;
	
		for(int i=1;i<this->levels;i++)
		{
			point  = bottomPoints[pointIndex];
			pointIndex++;
			if(this->bufferMode == 2)
			{
				if(i==levels/2)
				{
					double start = MPI_Wtime();		
					MPI_Wait(recvReq,MPI_STATUS_IGNORE);					
					double ellapsed = MPI_Wtime() - start;
				}
			}
			else if(this->bufferMode == 0)
			{
				double start = MPI_Wtime();
				//this->getGhostValues(i);
				int index = (i-1)*dataPointSize*2*substeps;
				MPI_Recv(&leftGhost[index],2*dataPointSize*substeps,MPI_DOUBLE,Lrank,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		
				double ellapsed = MPI_Wtime() - start;
			}			
			
			//adjustGhostValues(i);
			//
			int index = (i-1)*dataPointSize*2*substeps;
			double *destData;
			double *leftLocation;
			double *rightLocation;
		
			if(i%2 != 0)
			{
				destData   = foundation;
			}
			else
			{
				destData   = staging;
			}
			
			leftLocation  = destData;
			rightLocation = destData + (this->levels * dataPointSize * substeps);
			
			
			leftLocation  += (levels-i-1) * dataPointSize * substeps;		
			rightLocation += (i-1) * dataPointSize * substeps;
			
			memcpy(leftLocation ,&leftGhost[index] ,2 * dataPointSize * substeps * sizeof(double));
			memcpy(rightLocation,&rightGhost[index],2 * dataPointSize * substeps * sizeof(double));				
			
					
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
		if(this->bufferMode == 0)
		{
			double start = MPI_Wtime();
			//this->getGhostValues(i);
			int index = (levels-1)*dataPointSize*2*substeps;
			MPI_Recv(&leftGhost[index],2*dataPointSize*substeps,MPI_DOUBLE,Lrank,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		
			double ellapsed = MPI_Wtime() - start;
		}			
		
		int index = (levels-1)*dataPointSize*2*substeps;
		double *destData;
		double *leftLocation;
		double *rightLocation;
		
		if(levels%2 != 0)
		{
			destData   = foundation;
		}
		else
		{
			destData   = staging;
		}
		leftLocation  = destData;
		rightLocation = destData + (this->levels * dataPointSize * substeps);
		
		rightLocation += (levels-1) * dataPointSize * substeps;
		memcpy(leftLocation ,&leftGhost[(foundationLength * dataPointSize * substeps) - (dataPointSize * substeps)] ,1 * dataPointSize * substeps * sizeof(double));
		memcpy(rightLocation,&rightGhost[(foundationLength * dataPointSize * substeps) - (2 * dataPointSize * substeps)],1 * dataPointSize * substeps * sizeof(double));
		
		point = bottomPoints[pointIndex];
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
		point = bottomPoints[pointIndex];
		pointIndex++;
		memcpy(point.output+(1*dataPointSize),point.uInput,(substeps-1) * dataPointSize * sizeof(double));
		SubSteps::executeStepFnc(fncIndex,&point);
		point = bottomPoints[pointIndex];
		memcpy(point.output+(1*dataPointSize),point.uInput,(substeps-1) * dataPointSize * sizeof(double));
		SubSteps::executeStepFnc(fncIndex,&point);		
		fncIndex++;
		if(fncIndex==substeps)fncIndex=0;			
		executeFnc = fncIndex;
		//
		if(this->bufferMode == 1)	
			MPI_Request_free(&reqs[0]);
		else if(this->bufferMode == 2)
		{
			MPI_Request_free(&reqs[0]);
			MPI_Request_free(&reqs[1]);
		}
		else
		{
			for(int i=0;i<levels;i++)
				MPI_Request_free(&reqs[i]);
		}
		//
		//DiamondBottom Type 1 END

		//DiamondTop Type 2:
		fncIndex = executeFnc;
		pointIndex = 0;
	
		Rindex = foundationLength*dataPointSize*(substeps)-(2*dataPointSize*(substeps));
		for(int i=0;i<2*dataPointSize*substeps;i++)
		{
			localGhost[i] = foundation[Rindex+i];
		}
		if(this->bufferMode == 1)
		{
			for(int i=0;i<2*dataPointSize*substeps;i++)
			{
				remoteGhost[i] = foundation[i];
			}
		}
		else if(this->bufferMode == 2)
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
		
		//Loop for each Level of this DiamondTop
		for(int i=1;i<levels;i++)	
		{
			point  = topPoints[pointIndex];		
			pointIndex++;

			point.output  -= dataPointSize*substeps*(((levels-i)*2)-1);
			point.ulInput -= dataPointSize*substeps*(((levels-i)*2)-1);
			point.uInput  -= dataPointSize*substeps*(((levels-i)*2)-1);
			point.urInput -= dataPointSize*substeps*(((levels-i)*2)-1);
			for(int j=0;j<((levels-i)*2);j++)
			{
				memcpy(point.output+(1*dataPointSize),point.uInput,(substeps-1) * dataPointSize * sizeof(double));
				SubSteps::executeStepFnc(fncIndex,&point);

				if(j == 1)
				{
					point.output  -= dataPointSize*substeps;
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
					if(i==levels-1)
					{
						int index = i*2*dataPointSize*substeps;
						for(int c=0;c<2*dataPointSize*substeps;c++)
						{
							localGhost[index] = point.output[c];
							index++;
						}						
						break;
					}				
					point.output  += dataPointSize*substeps;
				}
				else if(j == ((levels-i)*2)-1)
				{
					point.output  -= dataPointSize*substeps;
					int index = i*2*dataPointSize*substeps;
					for(int c=0;c<2*dataPointSize*substeps;c++)
					{
						localGhost[index] = point.output[c];
						index++;
					}
					point.output  += dataPointSize*substeps;
				}
				point.output  += dataPointSize*substeps;
				point.ulInput += dataPointSize*substeps;
				point.uInput  += dataPointSize*substeps;
				point.urInput += dataPointSize*substeps;

			}
			fncIndex++;
			if(fncIndex==substeps)fncIndex=0;
		}	
		
		if(this->bufferMode == 1)
		{
			MPI_Isend(remoteGhost,2*dataPointSize*substeps*levels,MPI_DOUBLE,Lrank,1,MPI_COMM_WORLD,&reqs[0]);
		}
		else if(this->bufferMode == 2)
		{
			MPI_Isend(remoteGhost+2*dataPointSize*substeps*levels/2,2*dataPointSize*substeps*levels/2,MPI_DOUBLE,Lrank,1,MPI_COMM_WORLD,&reqs[1]);
		}		
		//DiamondTop Type 2 END

		//DiamondBottom Type 2:
		fncIndex = executeFnc;
		rightGhost = remoteGhost;
		leftGhost  = localGhost;
		
		if(this->bufferMode == 1)
		{
			double start = MPI_Wtime();
			std::clock_t startTime = std::clock();
			MPI_Recv(rightGhost,2*dataPointSize*substeps*levels,MPI_DOUBLE,Rrank,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		
			std::clock_t endTime = std::clock();
			double totalTime = (endTime - startTime) / (double)CLOCKS_PER_SEC;
			double ellapsed = MPI_Wtime() - start;
			//printf("It took %f microseconds to get the ghost values!\n",ellapsed*1000000);
		}
		else if(this->bufferMode == 2)
		{
			double start = MPI_Wtime();
			MPI_Recv(rightGhost,2*dataPointSize*substeps*levels/2,MPI_DOUBLE,Rrank,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			MPI_Irecv(rightGhost+2*dataPointSize*substeps*levels/2,2*dataPointSize*substeps*levels/2,MPI_DOUBLE,Rrank,1,MPI_COMM_WORLD,recvReq);
			double ellapsed = MPI_Wtime() - start;
		}
		pointIndex = 0;
	
		for(int i=1;i<this->levels;i++)
		{
			point  = bottomPoints[pointIndex];
			pointIndex++;
		
			if(this->bufferMode == 2)
			{
				if(i==levels/2)
				{
					double start = MPI_Wtime();
					MPI_Wait(recvReq,MPI_STATUS_IGNORE);					
					double ellapsed = MPI_Wtime() - start;
				}
			}
			else if(this->bufferMode !=1)
			{
				double start = MPI_Wtime();
				int index = (i-1)*dataPointSize*2*substeps;
				MPI_Recv(&rightGhost[index],2*dataPointSize*substeps,MPI_DOUBLE,Rrank,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		
				double ellapsed = MPI_Wtime() - start;
			}

			int index = (i-1)*dataPointSize*2*substeps;
			double *destData;
			double *leftLocation;
			double *rightLocation;
		
			if(i%2 != 0)
			{
				destData   = foundation;
			}
			else
			{
				destData   = staging;
			}
			leftLocation  = destData;
			rightLocation = destData + (this->levels * dataPointSize * substeps);
			
			leftLocation  += (levels-i-1) * dataPointSize * substeps;		
			rightLocation += (i-1) * dataPointSize * substeps;
			memcpy(leftLocation ,&leftGhost[index] ,2 * dataPointSize * substeps * sizeof(double));
			memcpy(rightLocation,&rightGhost[index],2 * dataPointSize * substeps * sizeof(double));
						
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
		if(this->bufferMode == 0)
		{
			double start = MPI_Wtime();
			int index = (levels-1)*dataPointSize*2*substeps;
			MPI_Recv(&rightGhost[index],2*dataPointSize*substeps,MPI_DOUBLE,Rrank,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		
			double ellapsed = MPI_Wtime() - start;
		}
		
		index = (levels-1)*dataPointSize*2*substeps;		
		
		if(levels%2 != 0)
		{
			destData   = foundation;
		}
		else
		{
			destData   = staging;
		}
		leftLocation  = destData;
		rightLocation = destData + (this->levels * dataPointSize * substeps);
		
		rightLocation += (levels-1) * dataPointSize * substeps;
		memcpy(leftLocation ,&leftGhost[(foundationLength * dataPointSize * substeps) - (dataPointSize * substeps)] ,1 * dataPointSize * substeps * sizeof(double));
		memcpy(rightLocation,&rightGhost[(foundationLength * dataPointSize * substeps) - (2 * dataPointSize * substeps)],1 * dataPointSize * substeps * sizeof(double));
		
		point = bottomPoints[pointIndex];
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
		point = bottomPoints[pointIndex];
		//A fix for DiamondBottom Type 2:
		point.ulInput = &this->localGhost[(foundationLength * dataPointSize * substeps) - (2 * dataPointSize * substeps)];
		//Fix End
		pointIndex++;
		memcpy(point.output+(1*dataPointSize),point.uInput,(substeps-1) * dataPointSize * sizeof(double));
		SubSteps::executeStepFnc(fncIndex,&point);
		point = bottomPoints[pointIndex];
		//A fix for DiamondBottom Type 2:
		point.urInput = &this->remoteGhost[(foundationLength * dataPointSize * substeps) - (1 * dataPointSize * substeps)];
		//Fix End
		memcpy(point.output+(1*dataPointSize),point.uInput,(substeps-1) * dataPointSize * sizeof(double));
		SubSteps::executeStepFnc(fncIndex,&point);
		fncIndex++;
		if(fncIndex==substeps)fncIndex=0;
	    executeFnc= fncIndex;
		//
		if(this->bufferMode == 1)
		{
			MPI_Request_free(&reqs[0]);
		}
		else if(this->bufferMode == 2)
		{
			MPI_Request_free(&reqs[0]);
			MPI_Request_free(&reqs[1]);
		}
		else
		{
			for(int i=0;i<levels;i++)
				MPI_Request_free(&reqs[i]);
		}
		//
		//DiamondBottom Type 2 END
	}

}

void DiamondDiscretization1D::calculate(int cycles)
{
	double buffer[4];
	for(int i=0;i<cycles;i++)
	{
		dt1->calculate(executeFnc);
		executeFnc    = db1->calculate(executeFnc);
		dt2->calculate(executeFnc);
		executeFnc    = db2->calculate(executeFnc);
	}
}

DiamondDiscretization1D::~DiamondDiscretization1D()
{
	delete[] foundation;
	delete[] staging;
	delete[] localGhost;
	delete[] remoteGhost;
	free(topPoints);
	free(bottomPoints);
	delete dt1,dt2,db1,db2;	
}
