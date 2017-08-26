/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * C MEX-file for performing WLS control allocation.
 *
 * To build this C MEX file, different commands should be used for
 * different platforms.
 *
 * On SGI, LINUX, Solaris, Alpha, and Macintosh platforms, LAPACK and
 * BLAS function names end with an underscore. Therefore use
 *
 *  >> mex -DUSCORE wlsc_alloc.c
 *
 * On a PC Windows platform, you need to explicitly specify a LAPACK
 * library file to link with. Use this command if you are using the
 * Lcc compiler that ships with MATLAB:
 *
 *  >> mex wlsc_alloc.c -L'C:/Program Files/MATLAB/R2012b/extern/lib/win64/microsoft' libmwlapack.lib  libmwblas.lib 
 *  >> mex wlsc_alloc.c libmwlapack.lib  libmwblas.lib  

 *
 *
 * where <matlab> is the Matlab root directory.
 *
 * For more info, see the section on using LAPACK in the Matlab
 * documentation (External Interfaces/API: Creating C Language
 * MEX-Files: Using LAPACK and BLAS Functions).
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "mex.h"
#include "math.h" /* sqrt */
#define PREC 1e-8

#ifdef USCORE
#define DGEMM dgemm_
#define DGEMV dgemv_
#define DCOPY dcopy_
#define DGELS dgels_
#define DSCAL dscal_
#else
#define DGEMM dgemm
#define DGEMV dgemv
#define DCOPY dcopy
#define DGELS dgels
#define DSCAL dscal
#endif

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * ~~~          Subroutine prototypes          ~~~
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void insert_Wu(double *Wu, double *A, int k, int m);
int copy_columns(double *A, double *B, double *W, int m, int n);
void add_vector_cond(double *a, double *b, double *W, int m);
int is_feasible(double *u, double *umin, double *umax, double *W, int m);
int lagrange_mult(double *lambda_ns, double *W, int m);
float max_step(double *u, double *p_free, double *umin, double *umax,
	       double *W, int m, int *i_block, int *if_block);
void print_vector(double *b, int n);
void print_matrix(double *A, int m, int n);

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * ~~~              Core function              ~~~
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

int wls_solve(double *B, double *v, double *umin, double *umax,
	      double *Wv, double *Wu, double *ud, double gam,
	      double *u, double *W, int imax, int k, int m)
{
  /* Variables */
  double mu, *A, *d, *A_free, *p_free, *u_opt, *lambda_ns, alpha;
  int i_block, if_block, i_neg, iter, m_free, verbose = 0;

  /* Variables used in BLAS calls */
  double done = 1.0, dzero = 0.0, dmone = -1.0, *work, *A_tmp;
  int km = k+m, ione = 1, info, lwork, tmp;
  char *chn = "N", *chl = "L", *cht = "T";
  
  if (imax<0) {
    /* Turn on verbosity (output info to command window). */
    imax = -imax;
    verbose = 1; /* verbose = 2 for even more details */
  }

  /* Allocate variables. */
  A = mxCalloc((k+m)*m, sizeof(double));
  d = mxCalloc(k+m, sizeof(double));
  A_free = mxCalloc((k+m)*m, sizeof(double));
  A_tmp = mxCalloc((k+m)*m, sizeof(double));
  p_free = mxCalloc(k+m, sizeof(double));
  u_opt = mxCalloc(m, sizeof(double));
  lambda_ns = mxCalloc(m, sizeof(double));

  /* Work area for least squares solves DGELS */
  tmp = 64; /* Generous block size, see
	       http://www.netlib.org/lapack/lug/node120.html#secblocksize*/
  lwork = m+km*tmp;
  work = mxCalloc(lwork,sizeof(double));

  /* Reformulate problem:
   *    min ||Wu(u-ud)||^2 + gam ||Wv(Bu-v)||^2 <=>
   *    min ||mu Wv(Bu-v)|| = min ||Au-b||
   *        ||   Wu(u-ud)||
   * where mu^2 = gam, A = [mu Wv*B ; Wu], b = [mu Wv*v ; Wu*ud] 
   *
   * Change of variables: u = u0+p -> ||Au-b|| = ||Ap-d||
   * where d = b - A*u0
   */
  
  mu = sqrt(gam);
  
  /* Create A = [mu*Wv*B ; Wu] */
  DGEMM(chn,chn,&k,&m,&k,&mu,Wv,&k,B,&k,&dzero,A,&km);
  insert_Wu(Wu,A,k,m);

  /* Create b = [mu Wv*v ; Wu*ud], store in d  */
  DGEMV(chn,&k,&k,&mu,Wv,&k,v,&ione,&dzero,d,&ione);
  DGEMV(chn,&m,&m,&done,Wu,&m,ud,&ione,&dzero,&d[k],&ione);
  
  /* Create d = b - A*u0 */
  DGEMV(chn,&km,&m,&dmone,A,&km,u,&ione,&done,d,&ione);

  /* ----------------------------------------------------------
   *  Iterate until optimum is found, or the maximum number of
   *  iterations is reached.
   * ---------------------------------------------------------- */
  
  for (iter=1 ; iter<=imax ; iter++) {

    if (verbose) {
      /* Output iteration info */
      mexPrintf("------------------------------------------\n");
      mexPrintf("Iteration: %d\n",iter);
      mexPrintf(" u = ");
      print_vector(u,m);
      mexPrintf(" W = ");
      print_vector(W,m);
    }
    
    /* ----------------------------------------
     *  Compute optimal perturbation vector p.
     * ---------------------------------------- */

    /* Eliminate saturated variables (create A_free). */
    m_free = copy_columns(A,A_free,W,km,m);
    if (verbose > 1) {
      mexPrintf(" A_free =\n");
      print_matrix(A_free,km,m_free);
      mexPrintf(" d = ");
      print_vector(d,km);
    }

    /* Solve the reduced optimization problem for the free variables, 
       i.e., compute least squares solution p_free = A_free\d. */
    DCOPY(&km,d,&ione,p_free,&ione);      /* p_free = d */
    tmp = km*m_free;
    DCOPY(&tmp,A_free,&ione,A_tmp,&ione); /* A_tmp = A_free */
                                           /* p_free = A_tmp\p_free */
    DGELS(chn,&km,&m_free,&ione,A_tmp,&km,p_free,&km,work,&lwork,&info);

    if (verbose) {  
      mexPrintf(" p_free = ");
      print_vector(p_free,m_free);
    }
    
    /* ----------------------------
     *  Is the new point feasible?
     * ---------------------------- */

    /* Create u_opt = u+p */
    DCOPY(&m,u,&ione,u_opt,&ione);    /* u_opt = u */
    add_vector_cond(p_free,u_opt,W,m); /* u_opt += p */
    if (verbose) {
      mexPrintf(" u + p = ");
      print_vector(u_opt,m);
    }

    if (is_feasible(u_opt,umin,umax,W,m)) {
      
       /* ----------------------------
       *  Yes, check for optimality.
       * ---------------------------- */

      if (verbose) {
	mexPrintf("Feasible w.r.t. constraints.\n");
      }

      /* Update point and residuals. */
      DCOPY(&m,u_opt,&ione,u,&ione); /* u = u_opt */
				      /* d = d - A_free*p_free */
      DGEMV(chn,&km,&m_free,&dmone,A_free,&km,p_free,&ione,&done,d,&ione);

      /* Compute Lagrange multipliers, not considering the sign due to
	 upper/lower bound. lambda_ns = A'*d */
      DGEMV(cht,&km,&m,&done,A,&km,d,&ione,&dzero,lambda_ns,&ione);
      /* i_neg is the index of the most negative lambda = W.*lambda_ns */
      i_neg = lagrange_mult(lambda_ns,W,m);
      
      if (i_neg<0) { /* All Lagrange multipliers positive */
	/*------------------------ \
	| Optimum found, bail out. |
	\ ------------------------*/
	if (verbose) {
	  mexPrintf("Optimum reached. Optimal solution:\n u = ");
	  print_vector(u,m);
	  mexPrintf("------------------------------------------\n");
	}

	/* Clean up */
	mxFree(A);
	mxFree(d);
	mxFree(A_free);
	mxFree(A_tmp);
	mxFree(p_free);
	mxFree(u_opt);
	mxFree(lambda_ns);
	mxFree(work);
	
	/* Return nr of iterations */
	return iter;
    }
    
      /* --------------------------------------------------
       *  Optimum not found, remove one active constraint.
       * -------------------------------------------------- */

      /* Remove constraint with most negative Lagrange multiplier
	 lambda from working set. */
      W[i_neg] = 0;
      if (verbose) {
	mexPrintf("Removing variable %d from working set.\n",i_neg+1);
      }
      
    } else /* is_feasible */ {
      
      /* ---------------------------------------
       *  No, find primary blocking constraint.
       * --------------------------------------- */

      if (verbose) {
	mexPrintf("Infeasible w.r.t. constraints.\n");
      }
      
      /* Compute primary blocking constraint (i_block, if_block) and
	 the maximum feasible step length along p (alpha). */
      alpha = max_step(u,p_free,umin,umax,W,m,&i_block,&if_block);
      
      /* Update point and residual. */
      DSCAL(&m_free,&alpha,p_free,&ione); /* p_free = alpha*p_free */
      add_vector_cond(p_free,u,W,m);       /* u += p */
      DGEMV(chn,&km,&m_free,&dmone,A_free,&km,
	     p_free,&ione,&done,d,&ione);  /* d -= A_free*p_free */
      
      /* Include the blocking constraint in the working set. 
	 Lower bound: -1, upper bound: +1 */
      W[i_block] = (p_free[if_block]>0) ? 1 : -1;
      
      if (verbose) {
	if (verbose>1) {
	  mexPrintf("Scaled perturbation: p_free = ");
	  print_vector(p_free,m_free);
	}
	mexPrintf("Adding variable %d (%s) to working set.\n",i_block+1,
	       (W[i_block]==1 ? "max" : "min"));
	mexPrintf("Reduced step length: alpha = %.3g\n",alpha);
      }

    } /* is_feasible */

  } /* for-loop */

  if (verbose) {
    mexPrintf("Max nr of iterations reached. Suboptimal solution:\n u = ");
    print_vector(u,m);
    mexPrintf("------------------------------------------\n");
  }

  /* Free allocated variables */
  mxFree(A);
  mxFree(d);
  mxFree(A_free);
  mxFree(A_tmp);
  mxFree(p_free);
  mxFree(u_opt);
  mxFree(lambda_ns);
  mxFree(work);

  /* Return nr of iterations. */
  return imax;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * ~~~             Gateway routine             ~~~
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  /* Input variables */
  double *B,*v,*umin,*umax,*Wv,*Wu,*ud,gam,*u0,*W0;
  int imax, ione = 1;

  /* Output variables */
  double *u,*W,*iter;

  /* Other variables */
  int i,k,m;

  /* Check nr of inputs */
  if (nrhs != 11) {
    mexErrMsgTxt("Wrong number of input arguments.");
  }

  /* Create a pointers to input variables */
  B	  = mxGetPr(prhs[0]);
  v	  = mxGetPr(prhs[1]);
  umin	  = mxGetPr(prhs[2]);
  umax	  = mxGetPr(prhs[3]);
  Wv	  = mxGetPr(prhs[4]);
  Wu	  = mxGetPr(prhs[5]);
  ud	  = mxGetPr(prhs[6]);
  gam	  = mxGetScalar(prhs[7]);
  u0	  = mxGetPr(prhs[8]);
  W0	  = mxGetPr(prhs[9]);
  imax	  = mxGetScalar(prhs[10]);
  
  /* Number of variables */
  k = mxGetM(prhs[0]);
  m = mxGetN(prhs[0]);

  /* Allocate space for output variables */
  plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(m, 1, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);

  /* Create pointers to output variables */
  u    = mxGetPr(plhs[0]);
  W    = mxGetPr(plhs[1]);
  iter = mxGetPr(plhs[2]);

  /* Copy initial values of u and W */
  DCOPY(&m,u0,&ione,u,&ione); /* u = u0 */
  DCOPY(&m,W0,&ione,W,&ione); /* W = W0 */

  /* Call wls_solve */
  *iter = wls_solve(B,v,umin,umax,Wv,Wu,ud,gam,u,W,imax,k,m);
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * ~~~               Subroutines               ~~~
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void insert_Wu(double *Wu, double *A, int k, int m)
/* Insert Wu into A */ 
{
  int i, j, idx_A, idx_Wu;
  
  idx_Wu = 0;
  idx_A = k;

  for (j = 0 ; j < m ; j++) {
    for (i = 0 ; i < m ; i++) {
      A[idx_A] = Wu[idx_Wu];
      idx_A++;
      idx_Wu++;
    }
    idx_A += k;
  }
}  

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

int copy_columns(double *A, double *B, double *W, int m, int n)
/* Copy columns with zero entry in W from A (m x n) to B.
   Return number of columns copied = # zero elements in W */
{
  int i, idx_A = 0, idx_B = 0, ione = 1, cols = 0;
  
  for (i=0 ; i<n ; i++) {
    /* For each column of A */
    if (!W[i]) {
      /* Copy i:th column of A to B */
      DCOPY(&m,&A[idx_A],&ione,&B[idx_B],&ione);
      idx_B += m;
      cols++;
    }
    idx_A += m;
  }
  return cols;
}
 
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void add_vector_cond(double *a, double *b, double *W, int m)
/* Set b(not(W)) += a */
{
  int i, j=-1;

  for (i=0 ; i<m ; i++) {
    if (!W[i]) {
      j++;
      b[i] += a[j];
    }
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

int is_feasible(double *u, double *umin, double *umax, double *W, int m)
/* Check if umin <= u <= umax for free components of u */
{
  int i;
  
  for (i=0 ; i<m ; i++) {
    if (W[i]==0) /* Free variable */
      if ((u[i] < umin[i]) || (u[i] > umax[i]))
	return 0;
  }
  return 1;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

int lagrange_mult(double *lambda_ns, double *W, int m)
/* Return the index of the most negative multiplier. Return -1 if all
   multipliers are positive. */
{
  double lambda, lambda_min;
  int i, i_min;
  
  /* Check if all lambda corresponding to active constraints are all
     positive. If not, find the most negative one. */
  lambda_min = -PREC;
  i_min = -1;
  for (i=0 ; i<m ; i++) {
    if (W[i]) { /* active constraint */
      lambda = W[i]*lambda_ns[i];
      if (lambda < lambda_min) {
	lambda_min = lambda;
	i_min = i;
      }
    }
  }

  /* Return index of constraint to be removed from working set. */
  return i_min;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

float max_step(double *u, double *p_free, double *umin, double *umax,
	       double *W, int m, int *i_block, int *if_block)
/* Compute the maximum step length, a, such that u+a*p_free is
   feasible w.r.t. umin and umax. Store index of blocking constraint
   in i_block and if_block. */
{
  int i, j = -1;
  double a, a_min = 1;

  /* Determine step length a_min and the corresponding variable index. */
  for (i=0 ; i<m ; i++) {
    if (W[i]==0) { /* free variable */
      j++; /* index in p_free */

      /* Check the direction of the i:th perturbation component. */
      if (p_free[j]>0) {
	a = (umax[i]-u[i])/p_free[j]; /* step length to upper bound */
      } else {
	if(p_free[j]<0) {
	  a = (umin[i]-u[i])/p_free[j]; /* step length to lower bound */
	} else { /* p_free[j]=0, degenerate case */
	  a = 2; /* arbitrary number > 1 */
	}
      }

      /* Most restrictive step length so far? */
      if (a<a_min) {
	/* Yes, store step length and indeces. */
	a_min = a;
	*i_block = i;
	*if_block = j;
      }
    }
  }
  
  /* Return maximum feasible step length */
  return a_min;
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void print_vector(double *b, int n)
/* Display contents of an n vector. */
{
  int i;
  
  mexPrintf("[ ");
  for (i=0 ; i<n ; i++)
    mexPrintf("%.3g ",b[i]);
  mexPrintf("]\n");
}

void print_matrix(double *A, int m, int n)
/* Display contents of an m x n matrix. */
{
  int i, j;

  for (i=0 ; i<m ; i++) {
    if (i==0) mexPrintf("[ ");
    else mexPrintf("  ");

    for (j=0 ; j<n ; j++)
      mexPrintf("%8.3g ",A[i+j*m]);
    if (i<m-1) mexPrintf("\n");
    else mexPrintf("]\n");
  }
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

