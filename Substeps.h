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

#ifndef H_SUBSTEPS
#define H_SUBSTEPS

#include "DiamondGlobals.h"

const double dx = .5;
const double DT = 0.005;
extern long step1Count,step2Count,step3Count,step4Count;
class SubSteps
{
public:
	
	static inline void executeStepFnc(int executeFnc,PointStruct *point)
	{
		switch(executeFnc)
		{
			case(0):
			{
				double u = point->uInput[0],
				uL = point->ulInput[0],
				uR = point->urInput[0];
				point->output[0] = (uL + uR - 2 * u) / (dx * dx);								
				step1Count++;
				break;
			}
			case(1):
			{
				double u = point->uInput[1],
				uL = point->ulInput[1],
				uR = point->urInput[1];
				double uxx = point->uInput[0],
				uxxL = point->ulInput[0],
				uxxR = point->urInput[0];;
				double conv = (uR*uR - uL*uL) / (4 * dx);
				double diff = ((uL + uxxL) + (uR + uxxR) - 2 * (u + uxx)) / (dx * dx);
				double dudt = -conv - diff;
				point->output[0] = u + 0.5 * DT * dudt;				
				step2Count++;
				break;
			}
			case(2):
			{
				double u0 = point->uInput[2],
				u  = point->uInput[0],
				uL = point->ulInput[0],
				uR = point->urInput[0];
				point->output[0] = (uL + uR - 2 * u) / (dx * dx);				
				step3Count++;
				break;
			}
			case(3):
			{
				double u0 = point->uInput[3];
				double u  = point->uInput[1],
				uL = point->ulInput[1],
				uR = point->urInput[1];;
				double uxx  = point->uInput[0],
				uxxL = point->ulInput[0],
				uxxR = point->urInput[0];
				double conv = (uR*uR - uL*uL) / (4 * dx);
				double diff = ((uL + uxxL) + (uR + uxxR) - 2 * (u + uxx)) / (dx * dx);
				double dudt = -conv - diff;
				point->output[0] = u0 + 0.5 * DT * dudt;				
				step4Count++;
				break;
			}
		}
	}	
};

#endif
