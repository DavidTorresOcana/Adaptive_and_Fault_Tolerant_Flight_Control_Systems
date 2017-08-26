#include "math.h"

/*  Merging the nlplant.c (lofi) and nlplant_hifi.c to use
    same equations of motion, navigation equations and use 
    own look-up tables decided by a flag.                   */
    
void atmos(double,double,double*);          /* Used by both */
void accels(double*,double*,double*);       /* Used by both */

#include "aerodata/lofi_F16_AeroData.c"              /* LOFI Look-up header file*/
#include "aerodata/hifi_F16_AeroData.c"              /* HIFI Look-up header file*/
#include "aerodata/engine_model.c"          /*engine model */

void nlplant(double*,double*);
void nlplant_Grad(double*,double*);

/*########################################*/
/*### Added for mex function in matlab ###*/
/*########################################*/

int fix(double);
int sign(double);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

#define XU prhs[0]
#define XDOTY plhs[0]
#define INPUTSIZE 15 
#define OUTPUTSIZE 21 

int i;
double *xup, *Moms_Coef_Grad;

if (mxGetM(XU)==INPUTSIZE && mxGetN(XU)==1){ 

      /* Calling Program */
      xup = mxGetPr(XU);
      XDOTY = mxCreateDoubleMatrix(OUTPUTSIZE, 1, mxREAL);
//       XDOTY = mxCreateDoubleMatrix(3, 7, mxREAL);
      Moms_Coef_Grad = mxGetPr(XDOTY);
      

      nlplant_Grad(xup,Moms_Coef_Grad);


} /* End if */
else{ 
      mexErrMsgTxt("Input and/or output is wrong size.");
} /* End else */

} /* end mexFunction */

/*########################################*/
/*### Evaluation of gradient            ###*/
/*########################################*/

void  nlplant_Grad(double *INPUTS, double *B) {
    double *XU_pert = INPUTS, *Moms_Coef = (double *)malloc(3*sizeof(double));
    double *Moms_Coef_Left = (double *)malloc(3*sizeof(double));
    double *Moms_Coef_Right = (double *)malloc(4*sizeof(double)) ;
    double  epsilon = 2;
    int i,j;

//     Perturbate inputs from index 8 to 14
    for(i=8;i<=14;i++){
        
        // Left
        XU_pert[i] = XU_pert[i]-epsilon;
        // Reduce perturbations to limits??
        nlplant(XU_pert,Moms_Coef_Left);
        XU_pert[i] = XU_pert[i]+epsilon; // Undo what done

        // Right
        XU_pert[i] = XU_pert[i]+epsilon;
        nlplant(XU_pert,Moms_Coef_Right);
        XU_pert[i] = XU_pert[i]-epsilon; // Undo what done

        
        // Put into matrix(Array in this case)
        B[(int) (3*(i-8))]=(Moms_Coef_Right[0]-Moms_Coef_Left[0])*0.5/(epsilon);
        B[(int) (3*(i-8)+1)]=(Moms_Coef_Right[1]-Moms_Coef_Left[1])*0.5/(epsilon);
        B[(int) (3*(i-8)+2)]=(Moms_Coef_Right[2]-Moms_Coef_Left[2])*0.5/(epsilon);

    }

//     printf("\n\n\n");
// 
//     for(j=0;j<=OUTPUTSIZE;j++){
//         printf("%f \n", B[j]);
//     }


}


/*########################################*/
/*### Actual Non-Linear plant          ###*/
/*########################################*/

void nlplant(double *INPUTS, double *OUTPUTS){
//     double j;     
//     printf("\n");
//         for(j=0;j<14;j++){
//             printf("%f  ", INPUTS[j]);
//         }
  /* #include f16_constants */
  double g    = 9.81;          /* gravity, m/s^2 */
  double m    = 9295.44;         /* mass, kg */
  double B    = 9.144;             /* span, m */
  double S    = 27.87;            /* planform area, m^2 */
  double cbar = 3.45;          /* mean aero chord, m */
  double xcgr = 0.35;      /* reference center of gravity as a fraction of cbar */
  double xcg  = 0.30;      /* center of gravity as a fraction of cbar. */

  double Heng = 216.9; /* engine angular momentum, assumed fixed */
  double pi   = acos(-1);
  double r2d  = 180.0/pi;     /* radians to degrees */
  
  double l_LEF = 2.54;    /*m distance from center line to LEF AC*/
  double l_a  = 3.82;     /*m distance from center line to aileron AC*/
  double l_e  =  1.69  ;  /*m distance from center line to elevator AC*/

/*NasaData        %translated via eq. 2.4-6 on pg 80 of Stevens and Lewis*/

  double Jy  = 75673.6;           /* kg-m^2 */ 
  double Jxz = 1331.4;             /* kg-m^2 */   
  double Jz  = 85552.1;           /* kg-m^2 */
  double Jx  = 12874.8;            /* kg-m^2 */

  double *temp;
  int i, j;
  
  double npos, epos, alt, phi, theta, psi, vt, alpha, beta, P, Q, R, pow, cpow;
  double sa, ca, sb, cb, tb, st, ct, tt, sphi, cphi, spsi, cpsi;
  double Throttle, el_L, el_R, ail_L, ail_R, rud,lef_L, lef_R, Ael, Sel, Adail, Sdail, dail_L, dail_R, drud, dlef_L,dlef_R, Adlef, Sdlef;
  double qbar, mach, ps;
  double S_LEF_0 , Thrust;
  double L_tot, M_tot, N_tot, denom;
  
  

  double Cx_tot, Cx_D, Cx_X, Cx_L_D, Cx_R_D, delta_Cx_lef, dXdQ, Cxq, delta_Cxq_lef;
  double Cz_tot, Cz_X, Cz_D, Cz_L_D, Cz_R_D, delta_Cz_lef, dZdQ, Czq, delta_Czq_lef;
  double Cm_tot, Cm_X, Cm, Cm_L, Cm_R, eta_el, eta_el_X, eta_el_L, eta_el_R, delta_Cm_lef, dMdQ, Cmq, delta_Cmq_lef, delta_Cm, delta_Cm_X,  delta_Cm_L, delta_Cm_R, delta_Cm_ds;
  double Cy_tot, Cy_X, Cy_D, Cy_L_D, Cy_R_D, delta_Cy_lef, dYdail, delta_Cy_r30, dYdR, dYdP;
  double delta_Cy_a20, delta_Cy_a20_lef, Cyr, delta_Cyr_lef, Cyp, delta_Cyp_lef;
  double Cn_tot, Cn_X, Cn_D, Cn_L_D, Cn_R_D, delta_Cn_lef, dNdail, delta_Cn_r30, dNdR, dNdP , delta_Cnbeta, delta_Cnbeta_X,delta_Cnbeta_L, delta_Cnbeta_R;
  double delta_Cn_a20, delta_Cn_a20_lef, Cnr, delta_Cnr_lef, Cnp, delta_Cnp_lef;
  double Cl_tot, Cl_X, Cl_D, Cl_L_D, Cl_R_D, delta_Cl_lef, dLdail, delta_Cl_r30, dLdR, dLdP, delta_Clbeta, delta_Clbeta_X, delta_Clbeta_L, delta_Clbeta_R;
  double delta_Cl_a20, delta_Cl_a20_lef, Clr, delta_Clr_lef, Clp, delta_Clp_lef;
  double Cx_cross, Cy_cross, Cz_cross, Cm_cross, Cn_cross, Cl_cross, Cydail, Cndail, Cldail;
  
  temp = (double *)malloc(9*sizeof(double));  /*size of 9.1 array*/

  /* %%%%%%%%%%%%%%%%%%%
           States
     %%%%%%%%%%%%%%%%%%% */


  alt   = INPUTS[0];   /* altitude */
  
  vt    = INPUTS[1];     /* total velocity */
  alpha = INPUTS[2]*r2d; /* angle of attack in degrees */
  beta  = INPUTS[3]*r2d; /* sideslip angle in degrees */
  
  P     = INPUTS[4];     /* Roll Rate --- rolling  moment is Lbar */
  Q     = INPUTS[5];    /* Pitch Rate--- pitching moment is M */
  R     = INPUTS[6];    /* Yaw Rate  --- yawing   moment is N */


  sa    = sin(INPUTS[2]); /* sin(alpha) */
  ca    = cos(INPUTS[2]); /* cos(alpha) */
  sb    = sin(INPUTS[3]); /* sin(beta)  */
  cb    = cos(INPUTS[3]); /* cos(beta)  */
  tb    = tan(INPUTS[3]); /* tan(beta)  */


  if (vt <= 0.01) {vt = 0.01;}

  /* %%%%%%%%%%%%%%%%%%%
     Control inputs
     %%%%%%%%%%%%%%%%%%% */

  Throttle  = INPUTS[7];   /* Throttle [0-1] */
  el_L    = INPUTS[8];   /* Left Elevator setting in degrees. */
  el_R    = INPUTS[9];   /* Right Elevator setting in degrees. */
  ail_L   = INPUTS[10];   /* Left Ailerons mex setting in degrees. */
  ail_R   = INPUTS[11];   /* Right Ailerons mex setting in degrees. */
  rud     = INPUTS[12];   /* Rudder setting in degrees. */
  lef_L   = INPUTS[13];   /* Left Leading edge flap setting in degrees */
  lef_R   = INPUTS[14];   /* Right Leading edge flap setting in degrees */
  
  Ael=(el_L-el_R)/2;        /*Asymetrical elevator def */
  Sel=(el_L+el_R)/2;         /*Symetrical elevator def */
  
  Adail  = (ail_L-ail_R)/(2*21.5);     /*Asymetrical aileron normalized against max angle */
  Sdail  = (ail_L+ail_R)/(2*21.5);     /*Symetrical aileron normalized against max angle */
  dail_L =  ail_L/21.5;
  dail_R =  ail_R/21.5;
  
  drud  = rud/30.0;  /* rudder normalized against max angle */
  
    /* %%%%%%%%%%%%%%%%%%
     Atmospheric effects
     sets dynamic pressure, mach number and gravity
     %%%%%%%%%%%%%%%%%% */

    atmos(alt,vt,temp);
       mach = temp[0];
       qbar = temp[1];
       ps   = temp[2];
       g    = temp[3];
       
    /* %%%%%%%%%%%%%%%%%%%
     	LEF inputs. Including the Symetric defletion for stall control!!!!!!!!!!!
     %%%%%%%%%%%%%%%%%%% */
  S_LEF_0 = 1.45 -9.05*qbar/ps +1.38* alpha;
  lef_L = lef_L + S_LEF_0;
  lef_R = lef_R + S_LEF_0;
  
  dlef_L  = (1 - lef_L/25.0);  /* Left leading edge flap normalized against max angle */
  dlef_R  = (1 - lef_R/25.0);  /* Right leading edge flap normalized against max angle */
  Adlef   = (dlef_L-dlef_R)/2; /* Asymetrical leading edge flap normalized against max angle */
  Sdlef   = (dlef_L+dlef_R)/2; /* Symetrical leading edge flap normalized against max angle */

// // //     /* Debug*/
// // //   printf("\nNorm:Aile_sym=%f  &  Norm LEF_asym=%f\n",Sdail,Adlef);

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
  /* Engine model, based on Ying Huo's m-files */
    cpow = tgear (Throttle);
    Thrust = thrst ( cpow, alt , mach ); 
    
/* ASYMETRIC MODEL*/
    /* HIFI Table */

                        /* Just X*/
        hifi_C(alpha,beta,0,temp);
            Cx_X = temp[0];
            Cz_X = temp[1];
            Cm_X = temp[2];
            Cy_X = temp[3];
            Cn_X = temp[4];
            Cl_X = temp[5];
                        /* Only deltas*/
        hifi_C(alpha,beta,el_L,temp);
            Cx_L_D = temp[0]-Cx_X;
            Cz_L_D = temp[1]-Cz_X;
            Cm_L = temp[2];
            Cy_L_D = temp[3]-Cy_X;
            Cn_L_D = temp[4]-Cn_X;
            Cl_L_D = temp[5]-Cl_X;
        hifi_C(alpha,beta,el_R,temp);
            Cx_R_D = temp[0]-Cx_X;
            Cz_R_D = temp[1]-Cz_X;
            Cm_R = temp[2];
            Cy_R_D = temp[3]-Cy_X;
            Cn_R_D = temp[4]-Cn_X;
            Cl_R_D = temp[5]-Cl_X;

        hifi_C_lef(alpha,beta,temp);
            delta_Cx_lef = temp[0];
            delta_Cz_lef = temp[1];
            delta_Cm_lef = temp[2];
            delta_Cy_lef = temp[3];
            delta_Cn_lef = temp[4];
            delta_Cl_lef = temp[5];

        hifi_damping_lef(alpha,temp);
            delta_Cxq_lef = temp[0];
            delta_Cyr_lef = temp[1];
            delta_Cyp_lef = temp[2];
            delta_Czq_lef = temp[3];
            delta_Clr_lef = temp[4];
            delta_Clp_lef = temp[5];
            delta_Cmq_lef = temp[6];
            delta_Cnr_lef = temp[7];
            delta_Cnp_lef = temp[8];

        hifi_rudder(alpha,beta,temp);
            delta_Cy_r30 = temp[0];
            delta_Cn_r30 = temp[1];
            delta_Cl_r30 = temp[2];

        hifi_ailerons(alpha,beta,temp);
            delta_Cy_a20     = temp[0];
            delta_Cy_a20_lef = temp[1];
            delta_Cn_a20     = temp[2];
            delta_Cn_a20_lef = temp[3];
            delta_Cl_a20     = temp[4];
            delta_Cl_a20_lef = temp[5];

        hifi_other_coeffs(alpha,el_L,temp);
            delta_Cnbeta_L = temp[0];
            delta_Clbeta_L = temp[1];
            delta_Cm_L     = temp[2];
            eta_el_L       = temp[3];
            delta_Cm_ds  = 0;        /* ignore deep-stall effect */
        hifi_other_coeffs(alpha,el_R,temp);
            delta_Cnbeta_R = temp[0];
            delta_Clbeta_R = temp[1];
            delta_Cm_R     = temp[2];
            eta_el_R       = temp[3];
            delta_Cm_ds  = 0;        /* ignore deep-stall effect */ 
            
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    compute Cx_tot, Cz_tot, Cm_tot, Cy_tot, Cn_tot, and Cl_tot
    (as on NASA report p37-40)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
            
    /* XXXXXXXX Cx_tot XXXXXXXX */
    dXdQ = (cbar/(2*vt))*(Cxq + delta_Cxq_lef*Sdlef);
    Cx_cross = - abs( 0.5*(delta_Cn_a20 + delta_Cn_a20_lef*dlef_L)*dail_L+0.5*(delta_Cn_a20 + delta_Cn_a20_lef*dlef_R)*dail_R ) * B/l_e;

    Cx_tot = Cx_X + 0.5*Cx_L_D+0.5*Cx_R_D + delta_Cx_lef*Sdlef + dXdQ*Q + Cx_cross ;

       /* ZZZZZZZZ Cz_tot ZZZZZZZZ */ 
    dZdQ = (cbar/(2*vt))*(Czq + delta_Cz_lef*Sdlef);
    Cz_cross = -B/l_a*(  0.5*(delta_Cl_a20 + delta_Cl_a20_lef*dlef_L)*dail_L + 0.5*(delta_Cl_a20 + delta_Cl_a20_lef*dlef_R)*dail_R  );

    Cz_tot = Cz_X + 0.5*Cz_L_D+0.5*Cz_R_D + delta_Cz_lef*Sdlef + dZdQ*Q + Cz_cross;

       /* YYYYYYYY Cy_tot YYYYYYYY */
    Cydail = 0.5*(delta_Cy_a20 + delta_Cy_a20_lef*dlef_L)*dail_L - 0.5*(delta_Cy_a20 + delta_Cy_a20_lef*dlef_R)*dail_R; 
// //     Cydail = (delta_Cy_a20 + delta_Cy_a20_lef*Sdlef)*Adail; 

    dYdR = (B/(2*vt))*(Cyr + delta_Cyr_lef*Sdlef);

    dYdP = (B/(2*vt))*(Cyp + delta_Cyp_lef*Sdlef);
    Cy_cross =  0;

    Cy_tot = Cy_X + 0.5*Cy_L_D+0.5*Cy_R_D + delta_Cy_lef*Sdlef + Cydail + delta_Cy_r30*drud + dYdR*R + dYdP*P + Cy_cross;

       /* MMMMMMMM Cm_tot MMMMMMMM */ 
    dMdQ = (cbar/(2*vt))*(Cmq + delta_Cmq_lef*Sdlef);
    Cm_cross =  0;

    Cm_tot = Cm_X*eta_el_X + 0.5*( ( Cm_L*eta_el_L- Cm_X*eta_el_X) + (Cm_R*eta_el_R- Cm_X*eta_el_X) )  + Cz_tot*(xcgr-xcg) + delta_Cm_lef*Sdlef + dMdQ*Q + 0.5*(delta_Cm_L+delta_Cm_R) + delta_Cm_ds + Cm_cross;

        /* NNNNNNNN Cn_tot NNNNNNNN */ 
    Cndail = 0.5*(delta_Cn_a20 + delta_Cn_a20_lef*dlef_L)*dail_L-0.5*(delta_Cn_a20 + delta_Cn_a20_lef*dlef_R)*dail_R;

    dNdR = (B/(2*vt))*(Cnr + delta_Cnr_lef*Sdlef);

    dNdP = (B/(2*vt))*(Cnp + delta_Cnp_lef*Sdlef);
    Cn_cross = l_e/B*(0.5*Cx_L_D-0.5*Cx_R_D) +l_LEF/B*delta_Cx_lef*Adlef;

    Cn_tot = Cn_X + 0.5*Cn_L_D+0.5*Cn_R_D + delta_Cn_lef*Sdlef - Cy_tot*(xcgr-xcg)*(cbar/B) + Cndail + delta_Cn_r30*drud + dNdR*R + dNdP*P + 0.5*(delta_Cnbeta_L+delta_Cnbeta_R)*beta + Cn_cross;

        /* LLLLLLLL Cl_tot LLLLLLLL */
    Cldail = 0.5*(delta_Cl_a20 + delta_Cl_a20_lef*dlef_L)*dail_L-0.5*(delta_Cl_a20 + delta_Cl_a20_lef*dlef_R)*dail_R;

    dLdR = (B/(2*vt))*(Clr + delta_Clr_lef*Sdlef);

    dLdP = (B/(2*vt))*(Clp + delta_Clp_lef*Sdlef);
    Cl_cross = - l_e/B*( 0.5*Cz_L_D - 0.5*Cz_R_D )  + l_LEF/B*delta_Cz_lef*Adlef;

    Cl_tot = Cl_X + 0.5*Cl_L_D + 0.5*Cl_R_D + delta_Cl_lef*Sdlef + Cldail + delta_Cl_r30*drud + dLdR*R + dLdP*P + 0.5*(delta_Clbeta_L+delta_Clbeta_R)*beta + Cl_cross;

    /* End of aerodynamic models*/

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     compute Pdot, Qdot, and Rdot (as in Stevens and Lewis p32)
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

// // // // L_tot = Cl_tot*qbar*S*B;       /* get moments from coefficients */
// // // // M_tot = Cm_tot*qbar*S*cbar;
// // // // N_tot = Cn_tot*qbar*S*B;

OUTPUTS[0]= Cl_tot;       /* get moments from coefficients */
OUTPUTS[1]= Cm_tot;
OUTPUTS[2]= Cn_tot;


//         printf("\n");
//         for(j=0;j<3;j++){
//             printf("%f  \n ", OUTPUTS[j]);
//             
//         }
/*########################################*/
/*########################################*/

free(temp);

}; /*##### END of nlplant() ####*/








/*########################################*/
/*### Called Sub-Functions  ##############*/
/*########################################*/

/*########################################*/
/* Function for mach and qbar */
/*########################################*/

void atmos(double alt, double vt, double *coeff ){

    double rho0 = 1.225;
    double Re = 6371000;
    double R = 287.05;
    double T0 = 288.15;
    double g0 = 9.81;
    double gamma = 1.4;
    double temp, rho, mach, qbar, grav, ps;

    temp = T0 - 0.0065 * alt;
    rho = rho0 * pow(  temp/T0 ,g0/(R*0.0065)-1  );
    if (alt >= 11000.0) {
       temp = 216.65;
       rho  = 0.3629;
    }


//     printf(" %f",rho);
    mach = vt / sqrt(gamma * R * temp);
    qbar = .5 * rho * vt * vt;
    grav = g0*(Re*Re/((Re+alt)*(Re+alt)));
    ps = R*rho*temp;
    
    coeff[0] = mach;
    coeff[1] = qbar;
    coeff[2] = ps;
    coeff[3] = grav;
}

/*########################################*/
/*########################################*/


/*########################################*/
/*### Port from matlab fix() function ####*/
/*########################################*/

int fix(double in){
    int out;

    if (in >= 0.0){
       out = (int)floor(in);
    }
    else if (in < 0.0){
       out = (int)ceil(in);
    }

    return out;
}

/* port from matlab sign() function */
int sign(double in){
    int out;

    if (in > 0.0){
       out = 1;
    }
    else if (in < 0.0){
       out = -1;
    }
    else if (in == 0.0){
       out = 0;
    }
    return out;
}

/*########################################*/
/*########################################*/



/*########################################*/
/*########################################*/

/*########################################*/
/*########################################*/

