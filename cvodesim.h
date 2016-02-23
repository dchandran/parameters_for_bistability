#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef GA_CVODE_WRAPPER_FILE
#define GA_CVODE_WRAPPER_FILE

#include <cvode/cvode.h>             /* prototypes for CVODE fcts. and consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., and macros */
#include <sundials/sundials_band.h>  /* definitions of type BandMat and macros */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS and EXP */
//#include "mathfunc.h"   /*eigenvalue computation */

#define SUNDIALS_DOUBLE_PRECISION 1

/*
 * set the flags
 * @param: only positive values
*/
void ODEflags(int);

/*
 * specify number of variables in the ODE system
 * @param: relative error allowed
 * @param: absolute error allowed
*/
void ODEtolerance(double,double);

/*
 * The Simulate function using Cvode (double precision)
 * @param: number of variables
 * @param: array of initial values
 * @param: ode function pointer
 * @param: start time for the simulation
 * @param: ending time for the simulation
 * @param: time increments for the simulation
 * @param: user data type for storing other information
 * @ret: 2D array made into linear array -- use getValue(array,N,i,j)
 */
double* ODEsim(int N, double * initValues, void (*odefnc)(double,double*,double*,void*), double startTime, double endTime, double stepSize, void * params);


/*
 * Gets jacobian matrix of the system at the given point
 * @param: number of variables
 * @param: array of values (point where Jacobian will be calculated)
 * @param: ode function pointer
 * @param: additional parameters needed for ode function
 * @ret: 2D array made into linear array -- use getValue(array,N,i,j)
 */
double* jacobian(int N, double * point,  void (*odefnc)(double,double*,double*,void*), void * params);

/*
 * Bring a system to steady state
 * @param: number of variables
 * @param: array of initial values
 * @param: ode function pointer
 * @param: minimum allowed value
 * @param: maximum time for simulation
 * @param: the difference in time to use for estimating steady state
 * @ret: array of values
 */
double* steadyState(int N, double * initialValues, void (*odefnc)(double,double*,double*,void*), void * params, double minerr, double maxtime, double delta);

/*
 * Find the rates of change after simulating for the given amount of time
 * @param: number of variables
 * @param: array of initial values
 * @param: ode function pointer
 * @param: time for simulation
 * @ret: array of values
 */
double* getDerivatives(int N, double * initValues, void (*odefnc)(double,double*,double*,void*), double startTime, double endTime, double stepSize, void * params);

#define getValue(array, N, i, j) ( array[ (((i)*(N)) + (j)) ] )

/*
* print a linearized 2D table to a file
* @param: filename to write to
* @param: data to write
* @param: number of rows
* @param: number of columns
*/
void writeToFile(char* filename, double* data, int rows, int cols);

#endif
