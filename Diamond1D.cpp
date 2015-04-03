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

#include <iostream>
#include "MpiGlobals.h"
#include "DiamondDiscretization1D.h"
#include "ClassicDiscretization1D.h"
#include "stdlib.h"
#include <ctime>
#include <cmath>

using namespace std;

const int substeps = 4;
const int dataPointSize = 1;
long step1Count = 0,step2Count = 0,step3Count = 0,step4Count = 0,count1,count2,count3,count4 ;
int totalElements,numGrids;
double *foundationData;

double init(double x) {
    const double PI = atan(1.0) * 4;
    return cos(x / 128. * 19 * PI) * 2.;
}

int main(int argc,char *argv[])
{
	initMPI(argc,argv);
	if(*argv[1] == 't')
	{
		totalElements = atoi(argv[2]);
		numGrids = totalElements/size;
	}		
	else if(*argv[1] == 'i')
	{
		numGrids = atoi(argv[2]);
		totalElements = numGrids*size;
	}
	else
	{
		if(myrank == 0)
		printf("Invalid command line options!\n");
	}
	int totalCycles = atoi(argv[3]);
	int bufferMode = atoi(argv[4]);
	if(bufferMode == 1)
	{
		if(myrank == 0)
		printf("Working in Buffered Ghost Mode\n");
	}
	else if(bufferMode == 2)
	{
		if(myrank == 0)
		 printf("Working in Half-Buffered Ghost Mode\n");
	}
	else
	{
		if(myrank == 0)
		printf("Non-Buffered Ghost Mode is still underdevelopment to replace MPI by C Socket!\n Exiting...");
		exit(-1);
	}
	
	//Initial values
	int foundationLength = numGrids;
    foundationData = new double[foundationLength];
    double startX = myrank*numGrids*dx;
    for(int i=1;i<=foundationLength;i++)
    {
		double x = startX+(i*dx);
		foundationData[i-1] = init(x);
    }
		
    DiamondDiscretization1D dis1d(foundationLength,dataPointSize,substeps,Lrank,myrank,Rrank,bufferMode);
    dis1d.setInitialValue(foundationData);        
		
	double start = MPI_Wtime();
    std::clock_t startTime = std::clock();
	dis1d.calculate(totalCycles);
	//dis1d.calculateInlined(totalCycles);
	double ellapsed = MPI_Wtime() - start;
	//dis1d.printOutput();
	
	MPI_Reduce(&step1Count,&count1,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&step2Count,&count2,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&step3Count,&count3,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&step4Count,&count4,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    if(myrank == 0)
    {
		std::clock_t endTime = std::clock();
        int totalSubSteps = totalCycles*foundationLength;
		
        printf("Diamond Job DONE - Completed %d substeps for %d spacial points! counts: %ld,%ld,%ld,%ld\n",totalSubSteps,numGrids*size,count1,count2,count3,count4);
        double totalTime = (endTime - startTime) / (double)CLOCKS_PER_SEC;
		double substepTime = (totalTime * 1000000)/(totalSubSteps);
		double substepTimePerElement = substepTime/totalElements;
		//printf("%f microseconds per substep\nDone!\n",substepTime);
		printf("%f microseconds per substep\nDone!\n\n",(ellapsed*1000000)/totalSubSteps);
    }

	//Classic Partitioning Part
	step1Count = 0 ; step2Count = 0 ; step3Count = 0; step4Count = 0;
	ClassicDiscretization1D classic1D(foundationLength,dataPointSize,substeps,Lrank,Rrank);
	for(int i=1;i<=foundationLength;i++)
    {
		double x = startX+(i*dx);
		//printf("myrank = %d , x = %f\n",myrank,x);
		foundationData[i-1] = init(x);
    }
	classic1D.setInitialValue(foundationData);
	start = MPI_Wtime();
	classic1D.calculate(totalCycles*foundationLength);
    ellapsed = MPI_Wtime() - start;
	//classic1D.printOutput();

	MPI_Reduce(&step1Count,&count1,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&step2Count,&count2,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&step3Count,&count3,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&step4Count,&count4,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    if(myrank == 0)
    {
		std::clock_t endTime = std::clock();
        int totalSubSteps = totalCycles*foundationLength;
        printf("Classic Job DONE - Completed %d substeps for %d spacial points! counts: %ld,%ld,%ld,%ld\n",totalSubSteps,numGrids*size,count1,count2,count3,count4);
        double totalTime = (endTime - startTime) / (double)CLOCKS_PER_SEC;
		double substepTime = (totalTime * 1000000)/(totalSubSteps);
		double substepTimePerElement = substepTime/totalElements;
		//printf("%f microseconds per substep\nDone!\n",substepTime);
		printf("%f microseconds per substep\nDone!\n",(ellapsed*1000000)/totalSubSteps);
    }	
	delete[] foundationData;
	finMPI();
	return 0;
}
