#include "mex.h"
#include <time.h>
#include <math.h>
#include <string.h>

/*Round function.*/
double round(double x) { return (x-floor(x))>0.5 ? ceil(x) : floor(x); }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Input/output variables. */
    double *data;
    double *xk;
    
    double *idx; 
  
    /* Intermediate variables.*/
    int datam,datan;      
    int i, j;
    
    /* Check for proper number of arguments. */    
    if (nrhs != 2)
    {
        mexErrMsgTxt("Two inputs required.");
    }
    else if (nlhs > 1)
    {
        mexErrMsgTxt("Wrong number of output arguments.");
    }
  
    /* Assign pointers to inputs. */
    data = mxGetPr(prhs[0]);
    xk   = mxGetPr(prhs[1]);
    
    /* Get sizes of input matrices (images, transformations, etc.).*/
    datam = mxGetM(prhs[0]);
    datan = mxGetN(prhs[0]);
    
    /* Create matrix for the return arguments. */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);  
         
    /* Assign pointers to output canvas (warped image2). */
    idx = mxGetPr(plhs[0]);
    
    /* Start computations. */    
        
    /* For each pixel in the target image... */
    idx[0] = 0;
    for(i=0;i<datam;i++)
    {
        if(data[i] == xk[0]){
            for(j=1;j<datan && data[(datam*j)+i] == xk[j];j++);
            if(j==datan)
            {
                idx[0] = i+1;
                return;
            }
        }
    }
    
    /* Bye bye.*/
    return;
}
