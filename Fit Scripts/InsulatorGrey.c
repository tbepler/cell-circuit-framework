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
    /* Retrieve model parameters. */ 
    double *Kdox, *del, *kp, *kpp, *k1, *k2, *k3, *k4, *k5, *k6, *kon, *koff, *ksgfp, *Kdgfp, *kdgfp, *k7, *km, *kg, *n1, *n2,*k8,*k9,*k10,*k11,*k12,*kon2;
    Kdox = p[0];   
    del  = p[1];   
    kp = p[2];   
    kpp  = p[3];   
    k1 = p[4];   
    k2 = p[5];  
    k3 = p[6];
    k4 = p[7];
    k5 = p[8];
    k6 = p[9];
    kon = p[10];
    koff = p[11];
    ksgfp = p[12];
    Kdgfp = p[13];
    kdgfp = p[14];
    k7 = p[15];
    km = p[16];
    kg = p[17];
    n1 = p[18];
    n2 = p[19];
    k8 = p[20]; /* C* -> C */
    k9 = p[21]; /* C** -> C */
    k10 = p[22];/* C+W*->C* + W */
    k11 = p[23];/* C*+W -> C + W* */
    k12 = p[24]; /* Skn7 kon */
    kon2 = p[25]; /* Skn7* kon */
   

/*   x[0] = Z; */
/* x[1] = Zp; */
/* x[2] = Wp; */
/* x[3] = Xp; */
/* x[4] = Xpp; */
/* x[5] = Cp; */
/* x[6] = Cpp; */
/* x[7] = GFP; */

dx[0] = km[0]*pow(u[0],n1[0])/(Kdox[0] + pow(u[0],n1[0])) - del[0]*x[0] - k2[0]*x[2]*x[0] + k1[0]*x[1]*(u[3]-x[2]) - kp[0]*x[0] + kpp[0]*x[1];
dx[1] =-k1[0]*x[1]*(u[3]-x[2]) + k2[0]*x[2]*x[0] + kp[0]*x[0] - kpp[0]*x[1] - del[0]*x[1];
dx[2] = k1[0]*x[1]*(u[3]-x[2]) - k2[0]*x[2]*x[0] - k3[0]*(u[2] - x[3] - x[4] - x[5] - x[6])*x[2] + k4[0]*x[3]*(u[3]-x[2]) - k3[0]*x[3]*x[2] + k4[0]*x[4]*(u[3]-x[2]) -k7[0]*x[2] + k12[0]*x[5]*(u[3] - x[2]) -k9[0]*x[5]*x[2] + k10[0]*x[6]*(u[3] - x[2]);
dx[3] = k3[0]*(u[2] - x[3] - x[4] - x[5] - x[6])*x[2] - k4[0]*x[3]*(u[3]-x[2]) - k3[0]*x[3]*x[2] +k4[0]*x[4]*(u[3]-x[2]) - k5[0]*x[3] + k6[0]*x[4]-kon2[0]*x[3]*(u[1] - x[5] - x[6]) + koff[0]*x[5];
dx[4] = k3[0]*x[3]*x[2] - k4[0]*x[4]*(u[3] - x[2]) - k6[0]*x[4]-kon[0]*x[4]*(u[1] - x[5] - x[6]) + koff[0]*x[6];
dx[5] = kon2[0]*x[3]*(u[1] - x[5] - x[6]) - koff[0]*x[5] - k9[0]*x[5]*x[2] + k10[0]*x[6]*(u[3] - x[2]) + k8[0]*x[6] -k11[0]*x[5] - k12[0]*x[5]*(u[3]-x[2]);
dx[6] = kon[0]*x[4]*(u[1] - x[5] - x[6]) - koff[0]*x[6] + k9[0]*x[5]*x[2] - k10[0]*x[6]*(u[3] - x[2])- k8[0]*x[6];
dx[7] = ksgfp[0] + kg[0]*pow(x[4],n2[0])/(Kdgfp[0] + pow(x[4],n2[0])) - kdgfp[0]*x[7];
}

/* Output equation. */
void compute_y(double *y, double t, double *x, double *u, double **p,
               const mxArray *auxvar)
{
    y[0] = x[7];
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