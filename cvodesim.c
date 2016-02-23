#include "cvodesim.h"

/*
 * The differential equations function
 * @param: time
 * @param: current values of the variables
 * @param: derviatives array used as the return array
 * @param: any other data pointer that is needed for the simulation
*/
static void (*ODEfunc)(double, double*, double*, void*) = NULL;

/*
 * relative error tolerance
 * absolute error tolerance
*/
double RelTol = 0, AbsTol = 1.0e-5;

int ODE_POSITIVE_VALUES_ONLY = 0;

/*
 * set the flags
 * @param: only positive values
*/
void ODEflags(int i)
{
   ODE_POSITIVE_VALUES_ONLY = i;
}

/*
 * @param: relative error tolerance allowed
 * @param: absolute error tolerance allowed
*/
void ODEtolerance(double relerr,double abserr)
{
   RelTol = relerr;
   AbsTol = abserr;
}

/**/
typedef struct
{
  void (*ODEfunc)(double, double*, double*, void*);
  void *userData;
} UserFunction;

/* f routine. Compute f(t,u). */

static int f(realtype t, N_Vector u, N_Vector udot, void * userFunc)
{
  realtype *udata, *dudata;

  udata = NV_DATA_S(u);
  dudata = NV_DATA_S(udot);

  UserFunction * info = (UserFunction*) userFunc;

  if ((*info).ODEfunc != NULL)
      ((*info).ODEfunc)(t,udata,dudata,(*info).userData);

  return (0);
}


/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */

  if (opt == 0 && flagvalue == NULL) {
    //fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
    //        funcname);
    return(1); }

  /* Check if flag < 0 */

  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
       //fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
       //      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */

  else if (opt == 2 && flagvalue == NULL) {
      //fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
      //        funcname);
    return(1); }

  return(0);
}

/*
 * The Simulate function using Cvode (double precision)
 * @param: number of variables
 * @param: array of initial values
 * @param: ode function pointer
 * @param: start time for the simulation
 * @param: ending time for the simulation
 * @param: time increments for the simulation
 * @param: user data type for storing other information
 * @ret: 2D array with time in the first column and values in the rest 
 */
double* ODEsim(int N, double* initialValues, void (*odefnc)(double,double*,double*,void*), double startTime, double endTime, double stepSize, void * params)
{
  if (startTime < 0) startTime = 0;
  if (endTime < startTime) return 0;

  if ( (2*stepSize) > (endTime-startTime) ) stepSize = (endTime - startTime)/2.0;

  double reltol = 0.0, abstol = 1.0e-5;
  double t = 0.0, tout = 0.0;
  void * cvode_mem = 0;
  N_Vector u;
  int flag, i, j;

  /*setup tolerance*/

  reltol = RelTol; 
  abstol = AbsTol;

  /*setup ode func*/

  ODEfunc = odefnc; 

  if (N < 1) return (0);  /*no variables in the system*/

  u = N_VNew_Serial(N);  /* Allocate u vector */
  if(check_flag((void*)u, "N_VNew_Serial", 0)) return(0);

  /* Initialize u vector */

  realtype * udata = NV_DATA_S(u);

  if (initialValues != NULL)
     for (i=0; i < N; ++i)
        udata[i] = initialValues[i];

  /* allocate output matrix */

  int M = (endTime - startTime) / stepSize;
  double* data = malloc ((N+1) * (M+1)  * sizeof(double) );

  /* setup CVODE */

  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(0);

  flag = CVodeMalloc(cvode_mem, f, 0, u, CV_SS, reltol, &abstol);
  if (check_flag(&flag, "CVodeMalloc", 1))
  {
     CVodeFree(&cvode_mem);
     N_VDestroy_Serial(u);
     if (data) free(data);
     return(0);
  }

  UserFunction * funcData = malloc( sizeof(UserFunction) );
  (*funcData).ODEfunc = odefnc;
  (*funcData).userData = params;

  flag = CVodeSetFdata(cvode_mem, funcData);
  if(check_flag(&flag, "CVodeSetFdata", 1))
  {
     CVodeFree(&cvode_mem);
     N_VDestroy_Serial(u);
     free(funcData);
     if (data) free(data);
     return(0);
  }

  flag = CVBand(cvode_mem, N, 0, N-1);
  if (check_flag(&flag, "CVBand", 1))
  {
     CVodeFree(&cvode_mem);
     N_VDestroy_Serial(u);
     free(funcData);
     if (data) free(data);
     return(0);
  }

   /* setup for simulation */

  startTime = 0.0;
  t = startTime;
  tout = startTime;
  i = 0;

  /*main simulation loop*/

  while ((tout <= endTime) && (i <= M))
  {
    /*store data*/
    if (data)
    {
       getValue(data,N+1,i,0) = t;
       for (j=0; j < N; ++j)
       {
           if (ODE_POSITIVE_VALUES_ONLY && (NV_DATA_S(u))[j] < 0) //special for bio networks
           {
              CVodeFree(&cvode_mem);
              N_VDestroy_Serial(u);
              free(funcData);
              if (data) free(data);
              data = 0;
              return 0;
           }
           else
             getValue(data,N+1,i,j+1) = (NV_DATA_S(u))[j]; //normal case
       }
    }
    ++i;

    tout = t + stepSize;
    flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    if (check_flag(&flag, "CVode", 1))
    {
       CVodeFree(&cvode_mem);
       N_VDestroy_Serial(u);
       free(funcData);
       if (data) free(data);
       data = 0;
       return 0;
    }
  }

  CVodeFree(&cvode_mem);
  N_VDestroy_Serial(u);
  free(funcData);
  return(data);   /*return outptus*/
}

/*
 * Gets jacobian matrix of the system at the given point
 * @param: number of variables
 * @param: array of values (point where Jacobian will be calculated)
 * @param: ode function pointer
 * @param: additional parameters needed for ode function
 * @ret: 2D array made into linear array -- use getValue(array,N,i,j)
 */
double* jacobian(int N, double * point,  void (*odefnc)(double,double*,double*,void*), void * params)
{
   if (odefnc == 0 || point == 0) return (0);
   double * J = (double*) malloc( N*N*sizeof(double));

   double dx = 1.0e-5;
   double * dy0 = (double*) malloc( N*sizeof(double) ),
          * dy1 = (double*) malloc( N*sizeof(double) );
   int i,j,k;

   for (i=0; i < N; ++i)
   {
     point[i] -= dx;             //x = x0-h
     odefnc(1.0,point,dy0,params);  //dy0 = f(x-h)
     point[i] += 2*dx;           //x = x0+h
     odefnc(1.0,point,dy1,params);  //dy1 = f(x+h)
     point[i] -= dx;             //x = x0
     for (j=0; j < N; ++j)
     {
        getValue(J,N,j,i) = (dy1[j] - dy0[j])/(dx+dx);  // J[j,i] = f(x+h) - f(x-h) / 2h
     }
   }
   free (dy0);
   free (dy1);
   return (J);
}

/*
 * Bring a system to steady state
 * @param: number of variables
 * @param: array of initial values
 * @param: ode function pointer
 * @param: maximum allowed value
 * @ret: array of values
 */
double* steadyState(int N, double * initialValues, void (*odefnc)(double,double*,double*,void*), void * params, double maxerr, double maxtime, double delta)
{
  double startTime = 0;
  double endTime = maxtime;

  double stepSize = 0.1;

  double reltol = 0.0, abstol = 1.0e-5;
  double t = 0.0, tout = 0.0;
  void * cvode_mem = 0;
  N_Vector u;
  int flag, i, j;
  /*setup tolerance*/

  reltol = RelTol;
  abstol = AbsTol;

  /*setup ode func*/
  ODEfunc = odefnc;
  if (N < 1) return (0);  /*no variables in the system*/

  u = N_VNew_Serial(N);  /* Allocate u vector */
  if(check_flag((void*)u, "N_VNew_Serial", 0)) return(0);

  /* allocate output matrix */

  double* ss = malloc (N * sizeof(double) );

  /* Initialize u vector */

  realtype * udata = NV_DATA_S(u);
  realtype * u0 = malloc(N*sizeof(realtype));
  if (initialValues != NULL)
     for (i=0; i < N; ++i)
        udata[i] = u0[i] = initialValues[i];

  /* setup CVODE */

  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(0);

  flag = CVodeMalloc(cvode_mem, f, 0, u, CV_SS, reltol, &abstol);
  if (check_flag(&flag, "CVodeMalloc", 1))
  {
     CVodeFree(&cvode_mem);
     N_VDestroy_Serial(u);
     if (ss) free(ss);
     if (u0) free(u0);
     return(0);
  }

  UserFunction * funcData = malloc( sizeof(UserFunction) );
  (*funcData).ODEfunc = odefnc;
  (*funcData).userData = params;

  flag = CVodeSetFdata(cvode_mem, funcData);
  if(check_flag(&flag, "CVodeSetFdata", 1))
  {
     CVodeFree(&cvode_mem);
     N_VDestroy_Serial(u);
     free(funcData);
     if (ss) free(ss);
     if (u0) free(u0);
     return(0);
  }
  flag = CVBand(cvode_mem, N, 0, N-1);
  if (check_flag(&flag, "CVBand", 1))
  {
     CVodeFree(&cvode_mem);
     N_VDestroy_Serial(u);
     free(funcData);
     if (ss) free(ss);
     if (u0) free(u0);
     return(0);
  }
  /* setup for simulation */

  double t0 = 0.0;
  startTime = 0.0;
  t = startTime;
  tout = startTime;
  i = 0;
  double err = maxerr + 1, temp;

  /*main simulation loop*/
  while (tout <= endTime)
  {
    tout = t + stepSize;
    flag = CVode(cvode_mem, tout, u, &t, CV_NORMAL);
    if (check_flag(&flag, "CVode", 1))
    {
       CVodeFree(&cvode_mem);
       N_VDestroy_Serial(u);
       free(funcData);
       if (ss) free(ss);
       if (u0) free(u0);
       ss = 0;
       u0 = 0;
       return 0;
    }
    if (ss && (tout - t0) >= delta)  //measure difference between y[t] - y[t-delta]
    {
      t0 = tout;
      err = ( (NV_DATA_S(u))[0] - u0[0] )*( (NV_DATA_S(u))[0] - u0[0] );
      for (j=0; j < N; ++j)
      {
         temp = ( (NV_DATA_S(u))[j] - u0[j] )*( (NV_DATA_S(u))[j] - u0[j] );
         if (temp > err) err = temp;               //max value from all dx/dt
         ss[j] = u0[j] = (NV_DATA_S(u))[j];  //next y points
         if (ODE_POSITIVE_VALUES_ONLY && (NV_DATA_S(u))[j] < 0)
         {
              CVodeFree(&cvode_mem);
              N_VDestroy_Serial(u);
              free(funcData);
              if (ss) free(ss);
              if (u0) free(u0);
              ss = 0; 
              u0 = 0;
              return 0;
         }
      }
    }
    /* check for steady state reached */
    if (err <= maxerr)
    {
       break;
       /*double * J = jacobian(N, ss, odefnc, params);  //get jacobian at steady state
       if (J)
       {
          double * wr = 0, * wi = 0;        //wr = real component, wi = imaginary component
          int k = eigenvalues(J,N,&wr,&wi); //get eigenvalues of J in wr and wi
          free(J);
          int stablePt = 1;
          if (k)  //if everything is ok (CLAPACK)
          {
              for (j=0; j<N; ++j) 
                 if (wr[j] >= 0)    //this is not a stable point
                 {
                     stablePt = 0;
                     break;
                 }
              free(wr);
              free(wi);
              if (stablePt)
              {
                 break;   //stable point found
              }
          }
        }*/
    }
  }
  if (tout >= endTime) //steady state not reached in the given amount of time
  {
      if (ss) free(ss);
      if (u0) free(u0);
      ss = 0;
      u0 = 0;
  }
  if (u0) free(u0);
  CVodeFree(&cvode_mem);  /* Free the integrator memory */
  N_VDestroy_Serial(u);
  free(funcData);
  return(ss);   /*return outptus*/
}


/*
 * Find the rates of change after simulating for the given amount of time
 * @param: number of variables
 * @param: array of initial values
 * @param: ode function pointer
 * @param: time for simulation
 * @ret: array of values
 */
double* getDerivatives(int N, double * initialValues, void (*odefnc)(double,double*,double*,void*), double startTime, double endTime, double stepSize, void * params)
{
  double * y = ODEsim(N,initialValues,odefnc,startTime,endTime,stepSize,params);
  if (y == 0) return 0;
  int sz = (int)((endTime-startTime)/stepSize);
  double * dy = malloc(N * sizeof(double));

  int i;
  for (i=0; i < N; ++i)
  {
      dy[i] = ( getValue(y,1+N,sz,1+i) - getValue(y,1+N,sz - 1,1+i) )/ stepSize;
  }
  free(y);
  return(dy);   /*return outptus*/
}


/*
* print a linearized 2D table to a file
*/
void writeToFile(char* filename, double* data, int rows, int cols)
{
   int i,j;
   FILE * out = fopen(filename,"w");
   for (i=0; i < rows; ++i)
   {
      fprintf(out, "%lf",getValue(data,cols,i,0));

      for (j=1; j < cols; ++j)
         fprintf(out, "\t%lf",getValue(data,cols,i,j));

      fprintf(out, "\n");
   }
   fclose(out);
}


