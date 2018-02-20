#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mex.h"
#include </usr/include/cvodes/cvodes.h>           /* prototypes for CVODES fcts. and consts. */
#include </usr/include/cvodes/cvodes_dense.h>     /* prototype for CVDENSE fcts. and constants */
#include </usr/include/nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */
#include </usr/include/sundials/sundials_types.h> /* def. of type realtype */
#include </usr/include/sundials/sundials_math.h>  /* definition of ABS */
/* Path to CVode to be changed as required above */


/* Input Arguments
 * - T_IN - Simulation Time and Time Step input in the form [T_0:time_step:T1]
 * - Y_IN - Initial Conditions
 * - P_IN - Parameter values the model is to be simulated with */

#define	T_IN	prhs[0]
#define	Y_IN	prhs[1]
#define P_IN    prhs[2]
#define S_IN    prhs[3]

/* Output Arguments
 * - YP_OUT is a vector containing the probability of being in the open state
 * at each time point */

#define	YP_OUT	plhs[0]
#define S_OUT   plhs[1]

static int f(realtype t, N_Vector y, N_Vector ydot, void* pr0);
/*static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, 
               DlsMat J, void *pr0, 
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);*/
static realtype ReturnVoltage( realtype t, realtype* pr );

/* Mex Function*/
static void MexHHScaled(N_Vector ydot, realtype t, N_Vector y, double* pr, double* T, double* Y0, int M, int NS, double* yout, double* sout) {
    
    int N=(M-1);
    
    /*Work out T0 and Tfinal from T_IN*/
    realtype T0 =T[0];
    
    realtype Tfinal = ((M-1)*((T[1])-(T[0])));
    
    /* Define tolerances to be used*/
    N_Vector abstol = NULL;
    
    abstol = N_VNew_Serial(3);
    
    NV_Ith_S(abstol, 0) = 1e-8;
    NV_Ith_S(abstol, 1) = 1e-8;
    NV_Ith_S(abstol, 2) = 1e-8;
    realtype reltol = 1e-8;
    
    /* Set up CVode*/
    int flag, k, ks, idx;
    N_Vector y0 = NULL;
    void* cvode_mem = NULL;
    y0 = N_VNew_Serial(3);
    NV_Ith_S(y0, 0) = Y0[0];
    NV_Ith_S(y0, 1) = Y0[1];
    NV_Ith_S(y0, 2) = Y0[2];
    
    realtype NState = 3;
    
    /* Solver options*/
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    
    /* Initialize the memory */
    flag = CVodeSetUserData(cvode_mem, pr);
    
    flag = CVodeInit(cvode_mem, &f, T0, y0);
    
    /* Set tolerances and maximum step*/
    flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
    flag = CVDense(cvode_mem, NState);
    flag = CVodeSetMaxStep(cvode_mem, 0.1);
    
    /* Call CVode test for error and add result to output vector yout */
    for (k = 1; k < N; ++k) {
        double tout = (k*(Tfinal)/N);
        
        
        if (CVode(cvode_mem, tout, y0, &t, CV_NORMAL) < 0) {
            fprintf(stderr, "Error in CVode: %d\n", flag);
        }
        
        /*Probability of being in open state is equal to 1-probability of being in any other state*/
        /*All states are required to calculate fluxes*/
        yout[k]     = NV_Ith_S(y0, 0);
        yout[k+N]   = NV_Ith_S(y0, 1);
        yout[k+2*N] = NV_Ith_S(y0, 2); // Note: open probability is now the third column
        
        for ( ks = 0; ks < NS; ++ks ){
            idx = k + N * ks;
            sout[idx] = (double) k+ks; // for now, assign value k+ks to sout to check something happens on every loop.
            }
    }
    
    /*Free memory*/
    N_VDestroy_Serial(y0);
    CVodeFree(&cvode_mem);
    
}


static int f(realtype t, N_Vector y, N_Vector ydot, void* pr0) {
    
    double* pr = pr0;
    double v;
    /* The number corresponding to the protocol to be simulated is passed at the front of the parameter values vector*/
   
    v = ReturnVoltage( t, pr );
    /* Parameters from input parameters - as the first element of the vector corresponds to the protocol number, the parameter values start from the second element*/
    
    
    double P0 = pr[1];
    double P1 = pr[2];
    double P2 = pr[3];
    double P3 = pr[4];
    double P4 = pr[5];
    double P5 = pr[6];
    double P6 = pr[7];
    double P7 = pr[8];
    
    
    
    /*Ensures microscopic reversibility condition satisfied*/
    
    const double y1 = NV_Ith_S(y, 0);
    const double y2 = NV_Ith_S(y, 1);
    const double y3 = NV_Ith_S(y, 2);
    
    const double y4 = (1.0-y1-y2-y3);
    /* Model equations*/
    
    
    
    const double k32 = P4*exp(P5*v);
    const double k23 = P6*exp(-P7*v);
    
    const double k43 = P0*exp(P1*v);
    const double k34 = P2*exp(-P3*v);
    
    const double k12 = k43;
    const double k21 = k34;
    
    const double k41 = k32;
    const double k14 = k23;
    
    
    NV_Ith_S(ydot, 0) = -k12*y1 + k21*y2 + k41*y4 - k14*y1;
    NV_Ith_S(ydot, 1) = -k23*y2 + k32*y3 + k12*y1 - k21*y2;
    NV_Ith_S(ydot, 2) = -k34*y3 + k43*y4 + k23*y2 - k32*y3;
    
    
    return 0;
}


/* 
 * Jacobian routine. Compute J(t,y). 
 */
// static int Jac(long int N, realtype t,
//                N_Vector y, N_Vector fy, 
//                DlsMat J, void *pr0, 
//                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
// {
//     realtype* pr = pr0;
//     realtype v;
//     /* The number corresponding to the protocol to be simulated is passed at the front of the parameter values vector*/
//     v = ReturnVoltage( t, pr );
//     /* Parameters from input parameters - as the first element of the vector corresponds to the protocol number, the parameter values start from the second element*/
//         
//     realtype P0 = pr[1];
//     realtype P1 = pr[2];
//     realtype P2 = pr[3];
//     realtype P3 = pr[4];
//     realtype P4 = pr[5];
//     realtype P5 = pr[6];
//     realtype P6 = pr[7];
//     realtype P7 = pr[8];
//        
//     /* Model equations*/
//     realtype k32 = P4*exp(P5*v);
//     realtype k23 = P6*exp(-P7*v);
//     
//     realtype k43 = P0*exp(P1*v);
//     realtype k34 = P2*exp(-P3*v);
//     
//     realtype k12 = k43;
//     realtype k21 = k34;
//     
//     realtype k41 = k32;
//     realtype k14 = k23;
//     
//     //realtype input;
//     //input = 0;
//     DENSE_ELEM(J,0,0) = -(k12+k14)-k41;   DENSE_ELEM(J,0,1) = k21-k41;      DENSE_ELEM(J,0,2) = -k41;
//     DENSE_ELEM(J,1,0) =   k12;            DENSE_ELEM(J,1,1) = -(k23+k21);   DENSE_ELEM(J,1,2) = k12;
//     DENSE_ELEM(J,2,0) =  -k43;            DENSE_ELEM(J,2,1) = k23;          DENSE_ELEM(J,2,2) = -(k34+k32)-k43;
// 
//     return(0);
// }

/* Mex function definition */
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray*prhs[] )
        
{
    N_Vector ydot;
    realtype t;
    N_Vector y;
    int M;
    int NS;
    
    /* Pointer to input variables/parameters*/
    
    double* pr;
    pr = mxGetPr(prhs[2]);
    
    M = mxGetN(T_IN);
    double* T;
    T = mxGetPr(prhs[0]);
    
    double* Y0;
    Y0 = mxGetPr(prhs[1]);
    
    double *yout;
    double *sout;
    
    int Numvar;
    Numvar = mxGetN(Y_IN);
    NS = ( mxGetN(S_IN) - 1 ) * Numvar; // Number of sensitivity parameters is num parameters * num Variables  
                                        // Note: first parameter is protocol definition
    
    /* Create a matrix for the return vector of state occupancy probabilities */
    YP_OUT = mxCreateDoubleMatrix(M-1, Numvar, mxREAL);
    /* Create a matrix for the return vector of sensitivties of the open vector */
    S_OUT  = mxCreateDoubleMatrix(M-1, NS, mxREAL);
    
    /* Assign pointer output */
    yout = mxGetPr(YP_OUT);
    /* Assign pointer output */
    sout = mxGetPr(S_OUT);
    
    /* Mex function call */
    MexHHScaled(ydot, t, y, pr, T, Y0, M, NS, yout, sout);
    return;
    
}

/*Calculate voltage at current time*/
static double ReturnVoltage( realtype t, double* pr ){
    
    /* Define voltage protocol to be used*/
    double protocol_number = pr[0];
    double v;
    /* This shift is needed for simulated protocol to match the protocol recorded in experiment, which is shifted by 0.1ms as compared to the original input protocol. Consequently, each step is held for 0.1ms longer in this version of the protocol as compared to the input.*/
    double shift = 0.1;
    
    /* sine wave*/
    if (protocol_number==1) {
        double C[6] = {54, 26, 10, 0.007/(2*M_PI), 0.037/(2*M_PI), 0.19/(2*M_PI)};
        
        if (t>=0 && t<250+shift)
        {
            v = -80;
        }
        
        if (t>=250+shift && t<300+shift)
        {
            v = -120;
        }
        
        if (t>=300+shift && t<500+shift)
        {
            v = -80;
        }
        
        if (t>=500+shift && t<1500+shift)
        {
            v = 40;
        }
        
        
        if (t>=1500+shift && t<2000+shift)
        {
            v = -120;
        }
        if (t>=2000+shift && t<3000+shift)
        {
            v = -80;
        }
        
        
        if (t>=3000+shift && t<6500+shift)
        {
            v=-30+C[0]*(sin(2*M_PI*C[3]*(t-2500-shift))) + C[1]*(sin(2*M_PI*C[4]*(t-2500-shift))) + C[2]*(sin(2*M_PI*C[5]*(t-2500-shift)));
            
        }
        if (t>=6500+shift && t<7000+shift)
        {
            v=-120;
        }
        if (t>= 7000+shift && t<8000+shift)
        {
            v = -80;
        }
    }
    
    
    if (protocol_number==7)
    {
        
        int l = floor(10*t);
        v=pr[l+9];
        
    }
    
    if (protocol_number==9){
        
        int l = floor(10*t);
        v=pr[l+9];
        
    }
    
    if (protocol_number==11){
        
        int l = floor(10*t);
        v=pr[l+9];
        
    }
    if (protocol_number==12){
        
        int l = floor(10*t);
        v=pr[l+9];
        
    }
    if (protocol_number==13){
        
        int l = floor(10*t);
        v=pr[l+9];
        
    }     if (protocol_number==14){
        
        int l = floor(10*t);
        v=pr[l+9];
        
    }
    if (protocol_number==18){
        
        int l = floor(10*t);
        v=pr[l+9];
        
    }     if (protocol_number==19){
        
        int l = floor(10*t);
        v=pr[l+9];
        
    }
    if (protocol_number==20){
        
        int l = floor(10*t);
        v=pr[l+9];
        
    }
    
    return v;
}
