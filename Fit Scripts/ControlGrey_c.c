/*   Copyright 2005-2008 The MathWorks, Inc. */
/*   $Revision: 1.1.8.4 $ $Date: 2008/04/28 03:17:22 $ */
/*   Written by Peter Lindskog. */

/* Include libraries. */
#include "mex.h"
#include <math.h>

/* Specify the number of outputs here. */
#define NY 1

/* State equations. */
void compute_dx(double *dx, double t, double *x, double *u, double **p,
                const mxArray *auxvar)

                {
     double  *kon, *delg,*km,*kg,*Kdox,*Kgfp, *ksgfp, *delx, *del;  
    
    del = p[0];
    kon = p[1];    
    km = p[2];
    kg = p[3];
    Kdox = p[4];
    Kgfp = p[5];
    ksgfp = p[6];
    delx = p[7];
    delg = p[8];
    

    dx[0] = km[0]*pow(u[0],1.0)/(Kdox[0] + pow(u[0],1.0))  - del[0]*x[0] -delx[0]*x[0] -kon[0]*x[0]*(u[1] - x[1]) + Kgfp[0]*kon[0]*x[1];
    dx[1] = kon[0]*x[0]*(u[1] - x[1]) - Kgfp[0]*kon[0]*x[1];
    dx[2] = ksgfp[0]+kg[0]*pow(x[0],1.0)/(Kgfp[0] + pow(x[0],1.0)) -del[0]*x[2] - delg[0]*x[2];


}
                
                
/* Output equation. */
void compute_y(double *y, double t, double *x, double *u, double **p,
               const mxArray *auxvar)
{
  
    y[0] = x[2];
}



/*----------------------------------------------------------------------- *
   DO NOT MODIFY THE CODE BELOW UNLESS YOU NEED TO PASS ADDITIONAL
   INFORMATION TO COMPUTE_DX AND COMPUTE_Y
 
   To add extra arguments to compute_dx and compute_y (e.g., size
   information), modify the definitions above and calls below.
 *-----------------------------------------------------------------------*/

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* Declaration of input and output arguments. */
    double *x, *u, **p, *dx, *y, *t;
    int     i, np, nu, nx;
    const mxArray *auxvar = NULL; /* Cell array of additional data. */
    
    if (nrhs < 3) {
        mexErrMsgIdAndTxt("IDNLGREY:ODE_FILE:InvalidSyntax",
        "At least 3 inputs expected (t, u, x).");
    }
    
    /* Determine if auxiliary variables were passed as last input.  */
    if ((nrhs > 3) && (mxIsCell(prhs[nrhs-1]))) {
        /* Auxiliary variables were passed as input. */
        auxvar = prhs[nrhs-1];
        np = nrhs - 4; /* Number of parameters (could be 0). */
    } else {
        /* Auxiliary variables were not passed. */
        np = nrhs - 3; /* Number of parameters. */
    }
    
    /* Determine number of inputs and states. */
    nx = mxGetNumberOfElements(prhs[1]); /* Number of states. */
    nu = mxGetNumberOfElements(prhs[2]); /* Number of inputs. */
    
    /* Obtain double data pointers from mxArrays. */
    t = mxGetPr(prhs[0]);  /* Current time value (scalar). */
    x = mxGetPr(prhs[1]);  /* States at time t. */
    u = mxGetPr(prhs[2]);  /* Inputs at time t. */
    
    p = mxCalloc(np, sizeof(double*));
    for (i = 0; i < np; i++) {
        p[i] = mxGetPr(prhs[3+i]); /* Parameter arrays. */
    }
    
    /* Create matrix for the return arguments. */
    plhs[0] = mxCreateDoubleMatrix(nx, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(NY, 1, mxREAL);
    dx      = mxGetPr(plhs[0]); /* State derivative values. */
    y       = mxGetPr(plhs[1]); /* Output values. */
    
    /*
      Call the state and output update functions.
      
      Note: You may also pass other inputs that you might need,
      such as number of states (nx) and number of parameters (np).
      You may also omit unused inputs (such as auxvar).
      
      For example, you may want to use orders nx and nu, but not time (t)
      or auxiliary data (auxvar). You may write these functions as:
          compute_dx(dx, nx, nu, x, u, p);
          compute_y(y, nx, nu, x, u, p);
    */
    
    /* Call function for state derivative update. */
    compute_dx(dx, t[0], x, u, p, auxvar);
    
    /* Call function for output update. */
    compute_y(y, t[0], x, u, p, auxvar);
    
    /* Clean up. */
    mxFree(p);
}