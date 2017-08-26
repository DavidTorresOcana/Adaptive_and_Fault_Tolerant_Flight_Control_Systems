#include"mex.h"
#include"mexndinterp.c"


double	*getALPHA1(){
FILE *fp = fopen("aerodata/ALPHA1.dat","r");
int i;
double *alpha1,data;

if(fp==NULL)
	mexErrMsgTxt("Can't find file ALPHA1.dat");

alpha1 = doubleVector(20);

for(i=0;i<20;i++){
	fscanf(fp,"%lf",&data);
	alpha1[i] = data;
	}
fclose(fp);
return(alpha1);
}



double	*getALPHA2(){
FILE *fp = fopen("aerodata/ALPHA2.dat","r");
int i;
double *alpha2,data;

if(fp==NULL)
	mexErrMsgTxt("Can't find file ALPHA2.dat");

alpha2 = doubleVector(14);

for(i=0;i<14;i++){
	fscanf(fp,"%lf",&data);
	alpha2[i] = data;
	}
fclose(fp);
return(alpha2);
}



double	*getBETA1(){
FILE *fp = fopen("aerodata/BETA1.dat","r");
int i;
double *beta1,data;

if(fp==NULL)
	mexErrMsgTxt("Can't find file BETA1.dat");

beta1 = doubleVector(19);

for(i=0;i<19;i++){
	fscanf(fp,"%lf",&data);
	beta1[i] = data;
	}
fclose(fp);
return(beta1);
}



double	*getDH1(){
FILE *fp = fopen("aerodata/DH1.dat","r");
int i;
double *dh1,data;

if(fp==NULL)
	mexErrMsgTxt("Can't find file DH1.dat");

dh1 = doubleVector(5);

for(i=0;i<5;i++){
	fscanf(fp,"%lf",&data);
	dh1[i] = data;
	}
fclose(fp);
return(dh1);
}

	

double	*getDH2(){
FILE *fp = fopen("aerodata/DH2.dat","r");
int i;
double *dh2,data;

if(fp==NULL)
	mexErrMsgTxt("Can't find file DH2.dat");

dh2 = doubleVector(3);

for(i=0;i<3;i++){
	fscanf(fp,"%lf",&data);
	dh2[i] = data;
	}
fclose(fp);
return(dh2);
}



double _Cx(double alpha,double beta,double dele){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 3; 
	double x[3];	
	FILESIZE = 1900;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19; 
		ndinfo.nPoints[2] = 5; 
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		X[2] = getDH1();
		fp = fopen("aerodata/CX0120_ALPHA1_BETA1_DH1_201.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CX0120_ALPHA1_BETA1_DH1_201.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
	x[1] = beta;
	x[2] = dele;

    return interpn(X,DATA,x,ndinfo);
}/* End of function(...) */


double _Cz(double alpha,double beta, double dele){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 3; /* alpha,beta,dele */
	double x[3];	/* Number of dimension */

	FILESIZE = 1900;	/* There are 1900 elements in the 20x19x5 3D array */

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); /* There are 1900 elements */
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	/* Alpha npoints */
		ndinfo.nPoints[1] = 19; /* Beta npoints  */
		ndinfo.nPoints[2] = 5;  /* dele npoints  */
		X = (double **) malloc(nDimension*sizeof(double*));
		X[0] = getALPHA1();
		X[1] = getBETA1();
		X[2] = getDH1();

		fp = fopen("aerodata/CZ0120_ALPHA1_BETA1_DH1_301.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CZ0120_ALPHA1_BETA1_DH1_301.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}
	x[0] = alpha;
	x[1] = beta;
	x[2] = dele;
    return interpn(X,DATA,x,ndinfo);
}/* End of function(...) */


double _Cm(double alpha,double beta,double dele){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 3; 
	double x[3];	
	FILESIZE = 1900;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19; 
		ndinfo.nPoints[2] = 5; 
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		X[2] = getDH1();
		fp = fopen("aerodata/CM0120_ALPHA1_BETA1_DH1_101.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CM0120_ALPHA1_BETA1_DH1_101.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
	x[1] = beta;
	x[2] = dele;
    return	interpn(X,DATA,x,ndinfo);
}/* End of function(...) */


double _Cy(double alpha,double beta){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 380;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19; 
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		fp = fopen("aerodata/CY0320_ALPHA1_BETA1_401.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CY0320_ALPHA1_BETA1_401.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
	x[1] = beta;
    return	interpn(X,DATA,x,ndinfo);
}/* End of function(...) */


double _Cn(double alpha, double beta, double dele){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 3; 
	double x[3];	
	FILESIZE = 1140;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19;	
		ndinfo.nPoints[2] = 3;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		X[2] = getDH2();
		fp = fopen("aerodata/CN0120_ALPHA1_BETA1_DH2_501.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CN0120_ALPHA1_BETA1_DH2_501.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
	x[1] = beta;
	x[2] = dele;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cl(double alpha, double beta,double dele){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 3; 
	double x[3];	
	FILESIZE = 1140;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19;	
		ndinfo.nPoints[2] = 3;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		X[2] = getDH2();
		fp = fopen("aerodata/CL0120_ALPHA1_BETA1_DH2_601.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CL0120_ALPHA1_BETA1_DH2_601.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
	x[1] = beta;
	x[2] = dele;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cx_lef(double alpha,double beta){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 266;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		ndinfo.nPoints[1] = 19; 
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		X[1] = getBETA1();
		fp = fopen("aerodata/CX0820_ALPHA2_BETA1_202.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CX0820_ALPHA2_BETA1_202.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
	x[1] = beta;
    return	interpn(X,DATA,x,ndinfo);
}/* End of function(...) */


double _Cz_lef(double alpha,double beta){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 266;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		ndinfo.nPoints[1] = 19; 
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		X[1] = getBETA1();
		fp = fopen("aerodata/CZ0820_ALPHA2_BETA1_302.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CZ0820_ALPHA2_BETA1_302.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
	x[1] = beta;
    return	interpn(X,DATA,x,ndinfo);
}/* End of function(...) */


double _Cm_lef(double alpha,double beta){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 266;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		ndinfo.nPoints[1] = 19; 
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		X[1] = getBETA1();
		fp = fopen("aerodata/CM0820_ALPHA2_BETA1_102.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CM0820_ALPHA2_BETA1_102.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
	x[1] = beta;
    return	interpn(X,DATA,x,ndinfo);
}/* End of function(...) */


double _Cy_lef(double alpha,double beta){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 266;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		ndinfo.nPoints[1] = 19; 
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		X[1] = getBETA1();
		fp = fopen("aerodata/CY0820_ALPHA2_BETA1_402.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CY0820_ALPHA2_BETA1_402.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
	x[1] = beta;
    return	interpn(X,DATA,x,ndinfo);
}/* End of function(...) */


double _Cn_lef(double alpha,double beta){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 266;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		ndinfo.nPoints[1] = 19; 
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		X[1] = getBETA1();
		fp = fopen("aerodata/CN0820_ALPHA2_BETA1_502.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CN0820_ALPHA2_BETA1_502.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
	x[1] = beta;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cl_lef(double alpha,double beta){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; /* alpha,beta*/
	double x[2];	/* Number of dimension */
	FILESIZE = 266;	/* There are 266 elements in the 14x19 2D array */

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	/* Alpha npoints */
		ndinfo.nPoints[1] = 19; /* Beta npoints  */
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		X[1] = getBETA1();
		fp = fopen("aerodata/CL0820_ALPHA2_BETA1_602.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CL0820_ALPHA2_BETA1_602.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
	x[1] = beta;
    return interpn(X,DATA,x,ndinfo);
}/* End of function(...) */


double _CXq(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/CX1120_ALPHA1_204.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CX1120_ALPHA1_204.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _CZq(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/CZ1120_ALPHA1_304.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CZ1120_ALPHA1_304.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _CMq(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/CM1120_ALPHA1_104.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CM1120_ALPHA1_104.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _CYp(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/CY1220_ALPHA1_408.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CY1220_ALPHA1_408.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _CYr(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/CY1320_ALPHA1_406.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CY1320_ALPHA1_406.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _CNr(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/CN1320_ALPHA1_506.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CN1320_ALPHA1_506.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _CNp(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/CN1220_ALPHA1_508.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CN1220_ALPHA1_508.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _CLp(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/CL1220_ALPHA1_608.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CL1220_ALPHA1_608.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _CLr(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/CL1320_ALPHA1_606.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CL1320_ALPHA1_606.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CXq_lef(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 14;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		fp = fopen("aerodata/CX1420_ALPHA2_205.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CX1420_ALPHA2_205.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CYr_lef(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 14;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		fp = fopen("aerodata/CY1620_ALPHA2_407.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CY1620_ALPHA2_407.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CYp_lef(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 14;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		fp = fopen("aerodata/CY1520_ALPHA2_409.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CY1520_ALPHA2_409.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CZq_lef(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 14;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		fp = fopen("aerodata/CZ1420_ALPHA2_305.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CZ1420_ALPHA2_305.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CLr_lef(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 14;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		fp = fopen("aerodata/CL1620_ALPHA2_607.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CL1620_ALPHA2_607.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CLp_lef(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 14;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		fp = fopen("aerodata/CL1520_ALPHA2_609.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CL1520_ALPHA2_609.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CMq_lef(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 14;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		fp = fopen("aerodata/CM1420_ALPHA2_105.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CM1420_ALPHA2_105.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CNr_lef(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 14;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		fp = fopen("aerodata/CN1620_ALPHA2_507.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CN1620_ALPHA2_507.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CNp_lef(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 14;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		fp = fopen("aerodata/CN1520_ALPHA2_509.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CN1520_ALPHA2_509.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cy_r30(double alpha, double beta){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 380;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		fp = fopen("aerodata/CY0720_ALPHA1_BETA1_405.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CY0720_ALPHA1_BETA1_405.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
	x[1] = beta;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cn_r30(double alpha, double beta){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 380;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		fp = fopen("aerodata/CN0720_ALPHA1_BETA1_503.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CN0720_ALPHA1_BETA1_503.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
	x[1] = beta;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cl_r30(double alpha, double beta){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 380;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		fp = fopen("aerodata/CL0720_ALPHA1_BETA1_603.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CL0720_ALPHA1_BETA1_603.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
	x[1] = beta;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cy_a20(double alpha, double beta){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 380;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		fp = fopen("aerodata/CY0620_ALPHA1_BETA1_403.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CY0620_ALPHA1_BETA1_403.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
	x[1] = beta;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cy_a20_lef(double alpha, double beta){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 266;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		ndinfo.nPoints[1] = 19;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		X[1] = getBETA1();
		fp = fopen("aerodata/CY0920_ALPHA2_BETA1_404.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CY0920_ALPHA2_BETA1_404.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
	x[1] = beta;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cn_a20(double alpha, double beta){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 380;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		fp = fopen("aerodata/CN0620_ALPHA1_BETA1_504.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CN0620_ALPHA1_BETA1_504.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
	x[1] = beta;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cn_a20_lef(double alpha, double beta){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 266;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		ndinfo.nPoints[1] = 19;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		X[1] = getBETA1();
		fp = fopen("aerodata/CN0920_ALPHA2_BETA1_505.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CN0920_ALPHA2_BETA1_505.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
	x[1] = beta;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cl_a20(double alpha, double beta){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 380;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		ndinfo.nPoints[1] = 19;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		X[1] = getBETA1();
		fp = fopen("aerodata/CL0620_ALPHA1_BETA1_604.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CL0620_ALPHA1_BETA1_604.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}
	x[0] = alpha;
	x[1] = beta;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _Cl_a20_lef(double alpha, double beta){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 2; 
	double x[2];	
	FILESIZE = 266;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 14;	
		ndinfo.nPoints[1] = 19;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA2();
		X[1] = getBETA1();
		fp = fopen("aerodata/CL0920_ALPHA2_BETA1_605.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CL0920_ALPHA2_BETA1_605.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
	x[1] = beta;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CNbeta(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/CN9999_ALPHA1_brett.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CN9999_ALPHA1_brett.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_CLbeta(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/CL9999_ALPHA1_brett.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CL9999_ALPHA1_brett.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _delta_Cm(double alpha){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 20;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 20;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getALPHA1();
		fp = fopen("aerodata/CM9999_ALPHA1_brett.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file CM9999_ALPHA1_brett.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = alpha;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


double _eta_el(double el){
	static int flag = 0;
	static double *DATA = (double*) NULL;
	static double **X;
	static ND_INFO ndinfo ;
	
	FILE *fp;
	double data;
	int i,FILESIZE;
	int nDimension = 1; 
	double x[1];	
	FILESIZE = 5;	

	/* Initialise everything when this function is called for the first time */
	if(flag==0){
		flag = 1;	/* Set to FILE_READ_TAG */
		DATA = (double*) malloc(FILESIZE*sizeof(double)); 
		ndinfo.nDimension = nDimension;
		ndinfo.nPoints = intVector(nDimension);
		ndinfo.nPoints[0] = 5;	
		X = (double **) malloc(nDimension*sizeof(double*));

		X[0] = getDH1();
		fp = fopen("aerodata/ETA_DH1_brett.dat","r");
		if(fp==(FILE*) NULL)
			mexErrMsgTxt("Cannot find file ETA_DH1_brett.dat in current directory");

		for(i=0;i<FILESIZE;i++){
			fscanf(fp,"%lf",&data);
			DATA[i] = data;
			}
		fclose(fp);
		}

	x[0] = el;
    return (interpn(X,DATA,x,ndinfo));
}/* End of function(...) */


/*
double _delta_Cm_ds(double alpha, double el){
...............
...............
} End of function(...) */


void hifi_C(double alpha,double beta,double el,double *retVal){
        retVal[0] = _Cx(alpha,beta,el);
        retVal[1] = _Cz(alpha,beta,el);
        retVal[2] = _Cm(alpha,beta,el);
        retVal[3] = _Cy(alpha,beta);
        retVal[4] = _Cn(alpha,beta,el);
        retVal[5] = _Cl(alpha,beta,el);
}

void hifi_damping(double alpha, double *retVal){
	retVal[0] = _CXq(alpha);
	retVal[1] = _CYr(alpha);
	retVal[2] = _CYp(alpha);
	retVal[3] = _CZq(alpha);
	retVal[4] = _CLr(alpha);
	retVal[5] = _CLp(alpha);
	retVal[6] = _CMq(alpha);
	retVal[7] = _CNr(alpha);
	retVal[8] = _CNp(alpha);
}

void hifi_C_lef(double alpha,double beta, double *retVal){
        retVal[0] = _Cx_lef(alpha,beta) - _Cx(alpha,beta,0);
        retVal[1] = _Cz_lef(alpha,beta) - _Cz(alpha,beta,0);
        retVal[2] = _Cm_lef(alpha,beta) - _Cm(alpha,beta,0);
        retVal[3] = _Cy_lef(alpha,beta) - _Cy(alpha,beta);
        retVal[4] = _Cn_lef(alpha,beta) - _Cn(alpha,beta,0);
        retVal[5] = _Cl_lef(alpha,beta) - _Cl(alpha,beta,0);
}

void hifi_damping_lef(double alpha, double *retVal){
        retVal[0] = _delta_CXq_lef(alpha);
        retVal[1] = _delta_CYr_lef(alpha);
        retVal[2] = _delta_CYp_lef(alpha);
        retVal[3] = _delta_CZq_lef(alpha);
        retVal[4] = _delta_CLr_lef(alpha);
        retVal[5] = _delta_CLp_lef(alpha);
        retVal[6] = _delta_CMq_lef(alpha);
        retVal[7] = _delta_CNr_lef(alpha);
        retVal[8] = _delta_CNp_lef(alpha);
}

void hifi_rudder(double alpha, double beta, double *retVal){
        retVal[0] = _Cy_r30(alpha,beta) - _Cy(alpha,beta);
        retVal[1] = _Cn_r30(alpha,beta) - _Cn(alpha,beta,0);
        retVal[2] = _Cl_r30(alpha,beta) - _Cl(alpha,beta,0);
}

void hifi_ailerons(double alpha, double beta, double *retVal){
        retVal[0] = _Cy_a20(alpha,beta) - _Cy(alpha,beta);
        retVal[1] = _Cy_a20_lef(alpha,beta) - _Cy_lef(alpha,beta) - retVal[0];
        retVal[2] = _Cn_a20(alpha,beta) - _Cn(alpha,beta,0);
        retVal[3] = _Cn_a20_lef(alpha,beta) - _Cn_lef(alpha,beta) - retVal[2];
        retVal[4] = _Cl_a20(alpha,beta) - _Cl(alpha,beta,0);
        retVal[5] = _Cl_a20_lef(alpha,beta) - _Cl_lef(alpha,beta) - retVal[4];
}

void hifi_other_coeffs(double alpha, double el, double *retVal){
        retVal[0] = _delta_CNbeta(alpha);
        retVal[1] = _delta_CLbeta(alpha);
        retVal[2] = _delta_Cm(alpha);
        retVal[3] = _eta_el(el);
        retVal[4] = 0;       /* ignore deep-stall regime, delta_Cm_ds = 0 */
}


