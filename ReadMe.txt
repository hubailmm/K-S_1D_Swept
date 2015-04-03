========================================================================
    CONSOLE APPLICATION : Diamond1D Project Overview
========================================================================
This Diamond1D Project implements the Swept-Rule decomposition for the 1D K-S
Equation.  Use the included compile script to compile the code.
The code is executed as follows:

mpirun -np NUM_PROC ./diamond.exe i NUM_SPACIAL_POINTS NUM_INTEGRATION_CYCLES
BUFFER_MODE

Where:
NUM_PROC: number of MPI processes to run
NUM_SPACIAL_POINTS: number of grid points in each MPI process
NUM_INTEGRATION_CYCLES: number of substeps to run, where each cycle has
NUM_SPACIAL_POINTS substeps
BUFFER_MODE: choose between 1 for full buffering or 2 for half-buffering.
Non-Buffered mode (0) is still under development

A typicla run is as follow:

mpirun -np 2 ./diamond.exe i 128 1000 1

