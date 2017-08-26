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
#define INPUTSIZE 34 
#define OUTPUTSIZE 2 

int i;
double *xup, *Moms_Coef;

//         printf("\n %f \n",mxGetN(XU));
if (mxGetM(XU)==INPUTSIZE && mxGetN(XU)==1){ 

      /* Calling Program */
      xup = mxGetPr(XU);
      XDOTY = mxCreateDoubleMatrix(OUTPUTSIZE, 1, mxREAL);
      Moms_Coef = mxGetPr(XDOTY);

      nlplant(xup,Moms_Coef);

      /* debug
      for (i=0;i<=14;i++){
        printf("xdotp(%d) = %e\n",i+1,xdotp[i]);
      }
      end debug */

} /* End if */
else{ 
      mexErrMsgTxt("Input and/or output is wrong size.");
} /* End else */

} /* end mexFunction */

/*########################################*/
/*########################################*/


void nlplant(double *INPUTS, double *OUTPUTS){
    
    int j ;

  /*################################################*/
  /*### Reconfiguration definitions and settings ###*/
  /*################################################*/
  /* 
  * Faults
      - Effectors Faults
            Blockades: Done in actuators dynamics block   DONE OUTSIDE!!
            Floating or loss: Eliminate all effects of the effector in aerodynamic coefficients. Can be understood as 100% of surface loss COMPUTED INSIDE!!
      - Structure Faults
            Loss of lifting surface(differenciating between wings!): Eliminate their entire effect in coeefficients and add the assymetry COMPUTED INSIDE!!
            Change in draging properties: Due to new configuration, external tank or loss of efficiency:
                        Multiply/add C_x_X+Delta_C_Lift_X*sin(alpha)-Delta_C_Drag_X*cos(alpha) ,  C_z_X and C_m_X  COMPUTED INSIDE!!
  
  * Deviations
      - Effectors Deviations
            Loss part of a surface: Multiply each surface contribution to coefficients times % of remaining surface: 100% is loss of surface COMPUTED INSIDE!!
            Aerodynamics interferece, icing can be seen as reduction of effectivenes: Same procedure as previous case COMPUTED INSIDE!!
            Reduction of linearity or faulty behaviour of actuators to be done at actuators block    DONE OUTSIDE!
      - Structure Deviations
            Loss of lifting surface (wing or fin): Multiply each surface contribution to coefficients times % of remaining surface: 100% is loss of surface
                Wing: Can be seen as a change in C_x_X*(0.15 +0.85*(S_current/S_w)) , C_z_X*(0.05 +0.95*(S_current/S_w)) and C_m_X*(0.4 +0.6*(S_current/S_w))?? COMPUTED INSIDE!!
                            OR!!!  C_L=-C_z_X*cos(alpha)+C_x_X*sin(alpha), C_L_p=C_L*(0.05 +0.95*(S_current/S_w)); C_x_X =-C_D_p*cos(alpha)+C_L_p*sin(alpha);
                Fin: Change some how Cn_X*(0.2 +0.8*(S_fin_current/S_fin)) and Cl_X*(0.25 +0.75*(S_fin_current/S_fin)) COMPUTED INSIDE!!
            Change in draging properties: Due to new configuration, external tank or loss of efficiency. Multiply/add C_x_X+Delta_C_Lift_X*sin(alpha)-Delta_C_Drag_X*cos(alpha)  ,  C_z_X and C_m_X     COMPUTED INSIDE!!
      - Global changes:
            Changes due to Mach and h: They are already considered NA
            Mass properties variations: Change Moments of inertia, CG and mass according with variations. These should be inputs to this script  COMING FROM OUTSIDE!!
                 */

  /* RECONFIGURABLE PARAMETERS */
  double m    = INPUTS[15]; //9295.44; /*TBC*/        /* mass, kg */
  double S_current_L  = INPUTS[16];// 27.87/2; /*TBC*/           /* Left planform area, m^2 */
  double S_current_R    = INPUTS[17]; //27.87/2; /*TBC*/           /* Right planform area, m^2 */
  double S_current = S_current_L + S_current_R;
  double S_fin_current    = INPUTS[18]; //6.578; /*TBC*/           /* Fin area, m^2 */
  double xcg  =  INPUTS[19]; //0.30; /*TBC*/     /* center of gravity as a fraction of cbar. */
  double Jy  =  INPUTS[20]; //75673.6; /*TBC*/          /* kg-m^2 */ 
  double Jxz =  INPUTS[21]; //1331.4;  /*TBC*/           /* kg-m^2 */   
  double Jz  = INPUTS[22]; //85552.1; /*TBC*/          /* kg-m^2 */
  double Jx  = INPUTS[23]; //12874.8; /*TBC*/           /* kg-m^2 */
  double EFF_el_L = INPUTS[24] , EFF_el_R = INPUTS[25], EFF_ail_L = INPUTS[26] , EFF_ail_R = INPUTS[27], EFF_rud = INPUTS[28], EFF_lef_L = INPUTS[29], EFF_lef_R= INPUTS[30]; /* [0-1] Dealing with surface Efficientcy: EFF_ represent a effectiveness or remaining surface compared with the normal surface*/
  double Delta_C_Lift_X = INPUTS[31], Delta_C_Drag_X = INPUTS[32], Delta_C_m_X = INPUTS[33]; /* Dealing with Change in draging properties: Treat as increments*/
  /*########################################*/
  /*########################################*/


  
  /* #include f16_constants */
  double g    = 9.81;          /* gravity, m/s^2 */
  double S    = 27.87;           /* planform area, m^2 */
  double S_L    = 27.87/2;        /* Left planform area, m^2 */
  double S_R    = 27.87/2;        /* Right  planform area, m^2 */
  double lambda    = 0.21;        /* Wing taper ratio */
  double S_fin    = 6.578;         /*  Fin area, Including both ventral fins! m^2 */

  double B    = 9.144;             /* span, m */
  double cbar = 3.45;          /* mean aero chord, m */
  double c_root = 5.04;          /* Root chord, m */

  double B_current_L  = 2*S_current_L/(c_root*(1+lambda));     /* Current left semi-span, m */
  double B_current_R  = 2*S_current_R/(c_root*(1+lambda));     /* Current right semi-span, m */
  double y_CMA_L = B_current_L/3*(1+2*lambda)/(1+lambda) , y_CMA_R = B_current_R/3*(1+2*lambda)/(1+lambda) ; /* Positions of the CMAs, m */
  double y_AC = 0.5*(y_CMA_R - y_CMA_L); /* Position of the AC*/
  double xcgr = 0.35 + 2.11*(y_AC/(B/2))/cbar ; /* reference center of gravity as a fraction of cbar */

  double Heng = 216.9; /* engine angular momentum, assumed fixed */
  double pi   = acos(-1);
  double r2d  = 180.0/pi;     /* radians to degrees */
  
  double l_LEF = 2.54;    /*m distance from center line to LEF AC*/
  double l_a  = 3.82;     /*m distance from center line to aileron AC*/
  double l_e  =  1.69  ;  /*m distance from center line to elevator AC*/

/*NasaData        %translated via eq. 2.4-6 on pg 80 of Stevens and Lewis*/

  double *temp;  
  
  double npos, epos, alt, phi, theta, psi, vt, alpha, beta, P, Q, R, pow, cpow;
  double sa, ca, sb, cb, tb, st, ct, tt, sphi, cphi, spsi, cpsi;
  double Throttle, el_L, el_R, ail_L, ail_R, rud,lef_L, lef_R, Ael, Sel, Adail, Sdail, dail_L, dail_R, drud, dlef_L,dlef_R, Adlef, Sdlef;
  double qbar, mach, ps;
  double S_LEF_0 , Thrust;
  double L_tot, M_tot, N_tot, denom;
  
  double C_Lift, C_Lift_p, C_Drag , C_Drag_p;
  double Cx_tot, Cx_D, Cx_X, Cx_L_D, Cx_R_D, delta_Cx_lef, dXdQ, Cxq, delta_Cxq_lef;
  double Cz_tot, Cz_X, Cz_D, Cz_L_D, Cz_R_D, delta_Cz_lef, dZdQ, Czq, delta_Czq_lef;
  double Cm_tot, Cm_X, Cm, Cm_L, Cm_R, eta_el, eta_el_X, eta_el_L, eta_el_R, delta_Cm_lef, dMdQ, Cmq, delta_Cmq_lef, delta_Cm, delta_Cm_X,  delta_Cm_L_D, delta_Cm_R_D, delta_Cm_ds;
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
  
  Adail  = (EFF_ail_L*ail_L-EFF_ail_R*ail_R)/(2*21.5);     /*Asymetrical aileron normalized against max angle */
  Sdail  = (EFF_ail_L*ail_L+EFF_ail_R*ail_R)/(2*21.5);     /*Symetrical aileron normalized against max angle */
  dail_L =  EFF_ail_L*ail_L/21.5;
  dail_R =  EFF_ail_R*ail_R/21.5;
  
  drud  = EFF_rud*rud/30.0;  /* rudder normalized against max angle */
  
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
    
  dlef_L  = EFF_lef_L*(1 - lef_L/25.0);  /* Left leading edge flap normalized against max angle */
  dlef_R  = EFF_lef_R*(1 - lef_R/25.0);  /* Right leading edge flap normalized against max angle */
  Adlef   = (dlef_L-dlef_R)/2; /* Asymetrical leading edge flap normalized against max angle */
  Sdlef   = (dlef_L+dlef_R)/2; /* Symetrical leading edge flap normalized against max angle */

  
  
  
//     printf("\n");
//     /* Debug*/
//     for(j=0;j<=INPUTSIZE;j++){
//         printf("%f ",INPUTS[j]);
//     }
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
            
        hifi_damping(alpha,temp);
            Cxq = temp[0];
            Cyr = temp[1];
            Cyp = temp[2];
            Czq = temp[3];
            Clr = temp[4];
            Clp = temp[5];
            Cmq = temp[6];
            Cnr = temp[7];
            Cnp = temp[8];

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

        hifi_other_coeffs(alpha,0,temp);
            delta_Cnbeta_X = temp[0];
            delta_Clbeta_X = temp[1];
            delta_Cm_X     = temp[2];
            eta_el_X       = temp[3];
        hifi_other_coeffs(alpha,el_L,temp);
            delta_Cnbeta_L = temp[0];
            delta_Clbeta_L = temp[1];
            delta_Cm_L_D     = temp[2]-delta_Cm_X;
            eta_el_L       = temp[3];
        hifi_other_coeffs(alpha,el_R,temp);
            delta_Cnbeta_R = temp[0];
            delta_Clbeta_R = temp[1];
            delta_Cm_R_D     = temp[2]-delta_Cm_X;
            eta_el_R       = temp[3];
            delta_Cm_ds  = 0;        /* ignore deep-stall effect */ 
            
    /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    compute Cx_tot, Cz_tot, Cm_tot, Cy_tot, Cn_tot, and Cl_tot
    (as on NASA report p37-40)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
            
    /* XXXXXX  Computing the effect of WING durface loss   XXXXXX */

	C_Lift=-Cz_X*ca+Cx_X*sa;
    C_Lift_p=C_Lift*(0.05 +0.95*(S_current/S));
	C_Drag=-Cx_X*ca-Cz_X*sa;
    C_Drag_p=C_Drag*(0.15 +0.85*(S_current/S));
    
    Cx_X = -C_Drag_p*ca+C_Lift_p*sa;
    Cz_X = -C_Drag_p*sa-C_Lift_p*ca;
    Cm_X = Cm_X*(0.4 +0.6*(S_current/S));
// // // // //         printf("\n");
// // // // //     printf("C_Lift %f  C_Drag %f  C_M %f",C_Lift,C_Drag,Cm_X);
    /* XXXXXX  Computing the effect of FIN durface loss   XXXXXX */

    Cn_X = Cn_X*(0.2 +0.8*(S_fin_current/S_fin));
    Cl_X = Cl_X*(0.25 +0.75*(S_fin_current/S_fin));
    
    /* XXXXXX  Change in Draging and Lifting characteristics: there are deltas!!   XXXXXX */
    Cx_X = Cx_X + Delta_C_Lift_X*sa-Delta_C_Drag_X*ca;
    Cz_X = Cz_X - Delta_C_Lift_X*ca-Delta_C_Drag_X*sa;
    Cm_X = Cm_X + Delta_C_m_X;
// // // //     printf("\n");
// // // //     printf("C_x %f  C_z %f  C_m %f",Cx_X,Cz_X,Cm_X);
    /* XXXXXXXX Cx_tot XXXXXXXX */
    dXdQ = (cbar/(2*vt))*(Cxq + delta_Cxq_lef*Sdlef);
    Cx_cross = - abs( 0.5*(delta_Cn_a20 + delta_Cn_a20_lef*dlef_L)*dail_L+0.5*(delta_Cn_a20 + delta_Cn_a20_lef*dlef_R)*dail_R ) * B/l_e;

    Cx_tot = Cx_X + 0.5*EFF_el_L*Cx_L_D+0.5*EFF_el_R*Cx_R_D + delta_Cx_lef*Sdlef + dXdQ*Q + Cx_cross ;

       /* ZZZZZZZZ Cz_tot ZZZZZZZZ */ 
    dZdQ = (cbar/(2*vt))*(Czq + delta_Cz_lef*Sdlef);
    Cz_cross = -B/l_a*(  0.5*(delta_Cl_a20 + delta_Cl_a20_lef*dlef_L)*dail_L + 0.5*(delta_Cl_a20 + delta_Cl_a20_lef*dlef_R)*dail_R  );

    Cz_tot = Cz_X + 0.5*EFF_el_L*Cz_L_D+0.5*EFF_el_R*Cz_R_D + delta_Cz_lef*Sdlef + dZdQ*Q + Cz_cross;

       /* MMMMMMMM Cm_tot MMMMMMMM */ 
    dMdQ = (cbar/(2*vt))*(Cmq + delta_Cmq_lef*Sdlef);
    Cm_cross =  0;

    Cm_tot = Cm_X*eta_el_X + 0.5*( EFF_el_L*( Cm_L*eta_el_L- Cm_X*eta_el_X) + EFF_el_R*(Cm_R*eta_el_R- Cm_X*eta_el_X) )  + Cz_tot*(xcgr-xcg) + delta_Cm_lef*Sdlef + dMdQ*Q + delta_Cm_X + 0.5*(EFF_el_L*delta_Cm_L_D + EFF_el_R*delta_Cm_R_D ) + delta_Cm_ds + Cm_cross;

       /* YYYYYYYY Cy_tot YYYYYYYY */
    Cydail = 0.5*(delta_Cy_a20 + delta_Cy_a20_lef*dlef_L)*dail_L - 0.5*(delta_Cy_a20 + delta_Cy_a20_lef*dlef_R)*dail_R; 
// //     Cydail = (delta_Cy_a20 + delta_Cy_a20_lef*Sdlef)*Adail; 

    dYdR = (B/(2*vt))*(Cyr + delta_Cyr_lef*Sdlef);

    dYdP = (B/(2*vt))*(Cyp + delta_Cyp_lef*Sdlef);
    Cy_cross =  0;

    Cy_tot = Cy_X + 0.5*EFF_el_L*Cy_L_D+0.5*EFF_el_R*Cy_R_D + delta_Cy_lef*Sdlef + Cydail + delta_Cy_r30*drud + dYdR*R + dYdP*P + Cy_cross;

       /* NNNNNNNN Cn_tot NNNNNNNN */ 
    Cndail = 0.5*(delta_Cn_a20 + delta_Cn_a20_lef*dlef_L)*dail_L-0.5*(delta_Cn_a20 + delta_Cn_a20_lef*dlef_R)*dail_R;

    dNdR = (B/(2*vt))*(Cnr + delta_Cnr_lef*Sdlef);

    dNdP = (B/(2*vt))*(Cnp + delta_Cnp_lef*Sdlef);
    Cn_cross = l_e/B*(0.5*EFF_el_L*Cx_L_D-0.5*EFF_el_R*Cx_R_D) +l_LEF/B*delta_Cx_lef*Adlef;

    Cn_tot = Cn_X + 0.5*EFF_el_L*Cn_L_D+0.5*EFF_el_R*Cn_R_D + delta_Cn_lef*Sdlef - Cy_tot*(xcgr-xcg)*(cbar/B) + Cndail + delta_Cn_r30*drud + dNdR*R + dNdP*P + 0.5*(delta_Cnbeta_L+delta_Cnbeta_R)*beta + Cn_cross;

    Cn_tot = Cn_tot - y_AC/B*Cx_X;/* Add asymetric loss of wing surface cross effect */

       /* LLLLLLLL Cl_tot LLLLLLLL */
    Cldail = 0.5*(delta_Cl_a20 + delta_Cl_a20_lef*dlef_L)*dail_L-0.5*(delta_Cl_a20 + delta_Cl_a20_lef*dlef_R)*dail_R;

    dLdR = (B/(2*vt))*(Clr + delta_Clr_lef*Sdlef);

    dLdP = (B/(2*vt))*(Clp + delta_Clp_lef*Sdlef);
    Cl_cross = - l_e/B*( 0.5*EFF_el_L*Cz_L_D - 0.5*EFF_el_R*Cz_R_D )  + l_LEF/B*delta_Cz_lef*Adlef;

    Cl_tot = Cl_X + 0.5*EFF_el_L*Cl_L_D + 0.5*EFF_el_R*Cl_R_D + delta_Cl_lef*Sdlef + Cldail + delta_Cl_r30*drud + dLdR*R + dLdP*P + 0.5*(delta_Clbeta_L+delta_Clbeta_R)*beta + Cl_cross;
    
    Cl_tot = Cl_tot + y_AC/B*Cz_X; /* Add asymetric loss of wing surface cross effect */
    
    
    /* End of aerodynamic models*/

  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     compute Pdot, Qdot, and Rdot (as in Stevens and Lewis p32)
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

// // // // L_tot = Cl_tot*qbar*S*B;       /* get moments from coefficients */
// // // // M_tot = Cm_tot*qbar*S*cbar;
// // // // N_tot = Cn_tot*qbar*S*B;

OUTPUTS[0]= C_Lift_p;       /* get moments from coefficients */
OUTPUTS[1]= C_Drag_p;

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

