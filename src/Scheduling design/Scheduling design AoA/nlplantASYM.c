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
#define INPUTSIZE 23 
#define OUTPUTSIZE 19 

int i;
double *xup, *xdotp;

//         printf("\n %f \n",mxGetN(XU));
if (mxGetM(XU)==INPUTSIZE && mxGetN(XU)==1){ 

      /* Calling Program */
      xup = mxGetPr(XU);
      XDOTY = mxCreateDoubleMatrix(OUTPUTSIZE, 1, mxREAL);
      xdotp = mxGetPr(XDOTY);

      nlplant(xup,xdotp);

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


void nlplant(double *xu, double *xdot){

  int fi_flag;
  int A_SMod;   /* Used symetric or asymetric model */

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

  double npos, epos, alt, phi, theta, psi, vt, alpha, beta, P, Q, R, pow, cpow;
  double sa, ca, sb, cb, tb, st, ct, tt, sphi, cphi, spsi, cpsi;
  double Throttle, el_L, el_R, ail_L, ail_R, rud,lef_L, lef_R, Ael, Sel, Adail, Sdail, dail_L, dail_R, drud, dlef_L,dlef_R, Adlef, Sdlef;
  double qbar, mach, ps;
  double U, V, W, Udot,Vdot,Wdot, Thrust;
  double L_tot, M_tot, N_tot, denom;
  
  
  double Cx_tot, Cx, Cx_L, Cx_R, delta_Cx_lef, dXdQ, Cxq, delta_Cxq_lef;
  double Cz_tot, Cz, Cz_L, Cz_R, delta_Cz_lef, dZdQ, Czq, delta_Czq_lef;
  double Cm_tot, Cm, Cm_L, Cm_R, eta_el, eta_el_L, eta_el_R, delta_Cm_lef, dMdQ, Cmq, delta_Cmq_lef, delta_Cm, delta_Cm_L, delta_Cm_R, delta_Cm_ds;
  double Cy_tot, Cy, Cy_L, Cy_R, delta_Cy_lef, dYdail, delta_Cy_r30, dYdR, dYdP;
  double delta_Cy_a20, delta_Cy_a20_lef, Cyr, delta_Cyr_lef, Cyp, delta_Cyp_lef;
  double Cn_tot, Cn, Cn_L, Cn_R, delta_Cn_lef, dNdail, delta_Cn_r30, dNdR, dNdP, delta_Cnbeta,delta_Cnbeta_L, delta_Cnbeta_R;
  double delta_Cn_a20, delta_Cn_a20_lef, Cnr, delta_Cnr_lef, Cnp, delta_Cnp_lef;
  double Cl_tot, Cl, Cl_L, Cl_R, delta_Cl_lef, dLdail, delta_Cl_r30, dLdR, dLdP, delta_Clbeta, delta_Clbeta_L, delta_Clbeta_R;
  double delta_Cl_a20, delta_Cl_a20_lef, Clr, delta_Clr_lef, Clp, delta_Clp_lef;
  double Cx_cross, Cy_cross, Cz_cross, Cm_cross, Cn_cross, Cl_cross, Cydail, Cndail, Cldail;
  
  temp = (double *)malloc(9*sizeof(double));  /*size of 9.1 array*/
  

  /* %%%%%%%%%%%%%%%%%%%
           States
     %%%%%%%%%%%%%%%%%%% */

  npos  = xu[0];   /* north position */
  epos  = xu[1];   /* east position */
  alt   = xu[2];   /* altitude */
  phi   = xu[3];   /* orientation angles in rad. */
  theta = xu[4];
  psi   = xu[5];

  vt    = xu[6];     /* total velocity */
  alpha = xu[7]*r2d; /* angle of attack in degrees */
  beta  = xu[8]*r2d; /* sideslip angle in degrees */
  P     = xu[9];     /* Roll Rate --- rolling  moment is Lbar */
  Q     = xu[10];    /* Pitch Rate--- pitching moment is M */
  R     = xu[11];    /* Yaw Rate  --- yawing   moment is N */

  pow     = xu[12];    /* Power level 0-100% */

  sa    = sin(xu[7]); /* sin(alpha) */
  ca    = cos(xu[7]); /* cos(alpha) */
  sb    = sin(xu[8]); /* sin(beta)  */
  cb    = cos(xu[8]); /* cos(beta)  */
  tb    = tan(xu[8]); /* tan(beta)  */

  st    = sin(theta);
  ct    = cos(theta);
  tt    = tan(theta);
  sphi  = sin(phi);
  cphi  = cos(phi);
  spsi  = sin(psi);
  cpsi  = cos(psi);

  if (vt <= 0.01) {vt = 0.01;}

  /* %%%%%%%%%%%%%%%%%%%
     Control inputs
     %%%%%%%%%%%%%%%%%%% */

  Throttle  = xu[13];   /* Throttle [0-1] */
  el_L    = xu[14];   /* Left Elevator setting in degrees. */
  el_R    = xu[15];   /* Right Elevator setting in degrees. */
  ail_L   = xu[16];   /* Left Ailerons mex setting in degrees. */
  ail_R   = xu[17];   /* Right Ailerons mex setting in degrees. */
  rud     = xu[18];   /* Rudder setting in degrees. */
  lef_L   = xu[19];   /* Left Leading edge flap setting in degrees */
  lef_R   = xu[20];   /* Right Leading edge flap setting in degrees */
  
  fi_flag = xu[21]/1;       /* fi_flag */
  A_SMod  = xu[22];       /* Flag of type of model: Use symetric or asymetric model */   
  
  Ael=(el_L-el_R)/2;        /*Asymetrical elevator def */
  Sel=(el_L+el_R)/2;         /*Symetrical elevator def */
  
  Adail  = (ail_L-ail_R)/(2*21.5);     /*Asymetrical aileron normalized against max angle */
  Sdail  = (ail_L+ail_R)/(2*21.5);     /*Symetrical aileron normalized against max angle */
  dail_L =  ail_L/21.5;
  dail_R =  ail_R/21.5;
  
  drud  = rud/30.0;  /* rudder normalized against max angle */
  
  dlef_L  = (1 - lef_L/25.0);  /* Left leading edge flap normalized against max angle */
  dlef_R  = (1 - lef_R/25.0);  /* Right leading edge flap normalized against max angle */
  Adlef   = (dlef_L-dlef_R)/2; /* Asymetrical leading edge flap normalized against max angle */
  Sdlef   = (dlef_L+dlef_R)/2; /* Symetrical leading edge flap normalized against max angle */

// // //     /* Debug*/
// // //   printf("\nNorm:Aile_sym=%f  &  Norm LEF_asym=%f\n",Sdail,Adlef);

  
  /* %%%%%%%%%%%%%%%%%%
     Atmospheric effects
     sets dynamic pressure, mach number and gravity
     %%%%%%%%%%%%%%%%%% */

atmos(alt,vt,temp);
   mach = temp[0];
   qbar = temp[1];
   ps   = temp[2];
   g    = temp[3];

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Dynamics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
  /* Engine model, based on Ying Huo's m-files */
    cpow = tgear (Throttle);
    xdot[12] = pdot ( pow, cpow );
    Thrust = thrst ( pow, alt , mach ); 
    
  /* %%%%%%%%%%%%%%%%%%
     Navigation Equations
     %%%%%%%%%%%%%%%%%% */

   U = vt*ca*cb;  /* directional velocities. */
   V = vt*sb;
   W = vt*sa*cb;

/* nposdot */
xdot[0] = U*(ct*cpsi) + 
            V*(sphi*cpsi*st - cphi*spsi) + 
            W*(cphi*st*cpsi + sphi*spsi);

/* eposdot */  
xdot[1] = U*(ct*spsi) + 
            V*(sphi*spsi*st + cphi*cpsi) + 
            W*(cphi*st*spsi - sphi*cpsi);

/* altdot */
xdot[2] = U*st - V*(sphi*ct) - W*(cphi*ct);

  /* %%%%%%%%%%%%%%%%%%%
     Kinematic equations
     %%%%%%%%%%%%%%%%%%% */
/* phidot */
xdot[3] = P + tt*(Q*sphi + R*cphi);

/* theta dot */
xdot[4] = Q*cphi - R*sphi;

/* psidot */
xdot[5] = (Q*sphi + R*cphi)/ct;

//  Select symetric or asymetric model for AERODYNAMICS
if(A_SMod==1){          /* ASYMETRIC MODEL*/
    /* HIFI Table */

        hifi_C(alpha,beta,el_L,temp);
            Cx_L = temp[0];
            Cz_L = temp[1];
            Cm_L = temp[2];
            Cy_L = temp[3];
            Cn_L = temp[4];
            Cl_L = temp[5];
        hifi_C(alpha,beta,el_R,temp);
            Cx_R = temp[0];
            Cz_R = temp[1];
            Cm_R = temp[2];
            Cy_R = temp[3];
            Cn_R = temp[4];
            Cl_R = temp[5];
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

    Cx_tot = 0.5*Cx_L+0.5*Cx_R + delta_Cx_lef*Sdlef + dXdQ*Q + Cx_cross ;

       /* ZZZZZZZZ Cz_tot ZZZZZZZZ */ 
    dZdQ = (cbar/(2*vt))*(Czq + delta_Cz_lef*Sdlef);
    Cz_cross = -B/l_a*(  0.5*(delta_Cl_a20 + delta_Cl_a20_lef*dlef_L)*dail_L + 0.5*(delta_Cl_a20 + delta_Cl_a20_lef*dlef_R)*dail_R  );

    Cz_tot = 0.5*Cz_L+0.5*Cz_R + delta_Cz_lef*Sdlef + dZdQ*Q + Cz_cross;

       /* MMMMMMMM Cm_tot MMMMMMMM */ 
    dMdQ = (cbar/(2*vt))*(Cmq + delta_Cmq_lef*Sdlef);
    Cm_cross =  0;

    Cm_tot = 0.5*(Cm_L*eta_el_L+Cm_R*eta_el_R)  + Cz_tot*(xcgr-xcg) + delta_Cm_lef*Sdlef + dMdQ*Q + 0.5*(delta_Cm_L+delta_Cm_R) + delta_Cm_ds + Cm_cross;

       /* YYYYYYYY Cy_tot YYYYYYYY */
    Cydail = 0.5*(delta_Cy_a20 + delta_Cy_a20_lef*dlef_L)*dail_L - 0.5*(delta_Cy_a20 + delta_Cy_a20_lef*dlef_R)*dail_R; 
// //     Cydail = (delta_Cy_a20 + delta_Cy_a20_lef*Sdlef)*Adail; 

    dYdR = (B/(2*vt))*(Cyr + delta_Cyr_lef*Sdlef);

    dYdP = (B/(2*vt))*(Cyp + delta_Cyp_lef*Sdlef);
    Cy_cross =  0;

    Cy_tot = 0.5*Cy_L+0.5*Cy_R + delta_Cy_lef*Sdlef + Cydail + delta_Cy_r30*drud + dYdR*R + dYdP*P + Cy_cross;

       /* NNNNNNNN Cn_tot NNNNNNNN */ 
    Cndail = 0.5*(delta_Cn_a20 + delta_Cn_a20_lef*dlef_L)*dail_L-0.5*(delta_Cn_a20 + delta_Cn_a20_lef*dlef_R)*dail_R;
// //     Cndail = (delta_Cn_a20 + delta_Cn_a20_lef*Sdlef)*Adail;

    dNdR = (B/(2*vt))*(Cnr + delta_Cnr_lef*Sdlef);

    dNdP = (B/(2*vt))*(Cnp + delta_Cnp_lef*Sdlef);
    Cn_cross = l_e/B*(0.5*Cx_L-0.5*Cx_R) +l_LEF/B*delta_Cx_lef*Adlef;

    Cn_tot = 0.5*Cn_L+0.5*Cn_R + delta_Cn_lef*Sdlef - Cy_tot*(xcgr-xcg)*(cbar/B) + Cndail + delta_Cn_r30*drud + dNdR*R + dNdP*P + 0.5*(delta_Cnbeta_L+delta_Cnbeta_R)*beta + Cn_cross;

       /* LLLLLLLL Cl_tot LLLLLLLL */
    Cldail = 0.5*(delta_Cl_a20 + delta_Cl_a20_lef*dlef_L)*dail_L-0.5*(delta_Cl_a20 + delta_Cl_a20_lef*dlef_R)*dail_R;
// //     Cldail = (delta_Cl_a20 + delta_Cl_a20_lef*Sdlef)*Adail;

    dLdR = (B/(2*vt))*(Clr + delta_Clr_lef*Sdlef);

    dLdP = (B/(2*vt))*(Clp + delta_Clp_lef*Sdlef);
    Cl_cross = - l_e/B*( 0.5*Cz_L - 0.5*Cz_R )  + l_LEF/B*delta_Cz_lef*Adlef;

    Cl_tot = 0.5*Cl_L+0.5*Cl_R + delta_Cl_lef*Sdlef + Cldail + delta_Cl_r30*drud + dLdR*R + dLdP*P + 0.5*(delta_Clbeta_L+delta_Clbeta_R)*beta + Cl_cross;

}
else{               /* SYMETRIC MODEL*/
        /* HIFI Table */

        hifi_C(alpha,beta,Sel,temp);
            Cx = temp[0];
            Cz = temp[1];
            Cm = temp[2];
            Cy = temp[3];
            Cn = temp[4];
            Cl = temp[5];

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

        hifi_other_coeffs(alpha,Sel,temp);
            delta_Cnbeta = temp[0];
            delta_Clbeta = temp[1];
            delta_Cm     = temp[2];
            eta_el       = temp[3];
            delta_Cm_ds  = 0;        /* ignore deep-stall effect */

    /* XXXXXXXX Cx_tot XXXXXXXX */

    dXdQ = (cbar/(2*vt))*(Cxq + delta_Cxq_lef*Sdlef);

    Cx_tot = Cx + delta_Cx_lef*Sdlef + dXdQ*Q;

       /* ZZZZZZZZ Cz_tot ZZZZZZZZ */ 

    dZdQ = (cbar/(2*vt))*(Czq + delta_Cz_lef*Sdlef);

    Cz_tot = Cz + delta_Cz_lef*Sdlef + dZdQ*Q;

       /* MMMMMMMM Cm_tot MMMMMMMM */ 

    dMdQ = (cbar/(2*vt))*(Cmq + delta_Cmq_lef*Sdlef);

    Cm_tot = Cm*eta_el + Cz_tot*(xcgr-xcg) + delta_Cm_lef*Sdlef + dMdQ*Q + delta_Cm + delta_Cm_ds;

       /* YYYYYYYY Cy_tot YYYYYYYY */

    dYdail = delta_Cy_a20 + delta_Cy_a20_lef*Sdlef;

    dYdR = (B/(2*vt))*(Cyr + delta_Cyr_lef*Sdlef);

    dYdP = (B/(2*vt))*(Cyp + delta_Cyp_lef*Sdlef);

    Cy_tot = Cy + delta_Cy_lef*Sdlef + dYdail*Adail + delta_Cy_r30*drud + dYdR*R + dYdP*P;

       /* NNNNNNNN Cn_tot NNNNNNNN */ 

    dNdail = delta_Cn_a20 + delta_Cn_a20_lef*Sdlef;

    dNdR = (B/(2*vt))*(Cnr + delta_Cnr_lef*Sdlef);

    dNdP = (B/(2*vt))*(Cnp + delta_Cnp_lef*Sdlef);

    Cn_tot = Cn + delta_Cn_lef*Sdlef - Cy_tot*(xcgr-xcg)*(cbar/B) + dNdail*Adail + delta_Cn_r30*drud + dNdR*R + dNdP*P + delta_Cnbeta*beta;

       /* LLLLLLLL Cl_tot LLLLLLLL */

    dLdail = delta_Cl_a20 + delta_Cl_a20_lef*Sdlef;

    dLdR = (B/(2*vt))*(Clr + delta_Clr_lef*Sdlef);

    dLdP = (B/(2*vt))*(Clp + delta_Clp_lef*Sdlef);

    Cl_tot = Cl + delta_Cl_lef*Sdlef + dLdail*Adail + delta_Cl_r30*drud + dLdR*R + dLdP*P + delta_Clbeta*beta;

        
}  /* End of aerodynamic models*/

   /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      compute Udot,Vdot, Wdot,(as on NASA report p36)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

Udot = R*V - Q*W - g*st + qbar*S*Cx_tot/m + Thrust/m;

Vdot = P*W - R*U + g*ct*sphi + qbar*S*Cy_tot/m;

Wdot = Q*U - P*V + g*ct*cphi + qbar*S*Cz_tot/m;

   /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      vt_dot equation (from S&L, p82)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

xdot[6] = (U*Udot + V*Vdot + W*Wdot)/vt;

   /* %%%%%%%%%%%%%%%%%%
      alpha_dot equation
      %%%%%%%%%%%%%%%%%% */

xdot[7] = (U*Wdot - W*Udot)/(U*U + W*W);

  /* %%%%%%%%%%%%%%%%%
     beta_dot equation
     %%%%%%%%%%%%%%%%% */

xdot[8] = (Vdot*vt - V*xdot[6])/(vt*vt*cb);



  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     compute Pdot, Qdot, and Rdot (as in Stevens and Lewis p32)
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

L_tot = Cl_tot*qbar*S*B;       /* get moments from coefficients */
M_tot = Cm_tot*qbar*S*cbar;
N_tot = Cn_tot*qbar*S*B;

denom = Jx*Jz - Jxz*Jxz;

  /* %%%%%%%%%%%%%%%%%%%%%%%
     Pdot
     %%%%%%%%%%%%%%%%%%%%%%% */

xdot[9] =  (Jz*L_tot + Jxz*N_tot - (Jz*(Jz-Jy)+Jxz*Jxz)*Q*R + Jxz*(Jx-Jy+Jz)*P*Q + Jxz*Q*Heng)/denom;

  /* %%%%%%%%%%%%%%%%%%%%%%%
     Qdot
     %%%%%%%%%%%%%%%%%%%%%%% */

xdot[10] = (M_tot + (Jz-Jx)*P*R - Jxz*(P*P-R*R) - R*Heng)/Jy;

  /* %%%%%%%%%%%%%%%%%%%%%%%
     Rdot
     %%%%%%%%%%%%%%%%%%%%%%% */

xdot[11] = (Jx*N_tot + Jxz*L_tot + (Jx*(Jx-Jy)+Jxz*Jxz)*P*Q - Jxz*(Jx-Jy+Jz)*Q*R +  Jx*Q*Heng)/denom;



/*########################################*/
/*### Create accelerations anx_cg, any_cg */
/*### ans anz_cg as outputs ##############*/
/*########################################*/

accels(xu,xdot,temp);

xdot[13]  = temp[0];	/* anx_cg */
xdot[14]  = temp[1];	/* any_cg */
xdot[15]  = temp[2];	/* anz_cg */
xdot[16]  = mach;
xdot[17]  = qbar;
xdot[18]  = ps;

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
/*### Calculate accelerations from states */
/*### and state derivatives. ############ */
/*########################################*/

void accels(double *state,
            double *xdot,
            double *y)
{


#define grav 9.81 

double sina, cosa, sinb, cosb ;
double vel_u, vel_v, vel_w ;
double u_dot, v_dot, w_dot ;
double nx_cg, ny_cg, nz_cg ;

sina = sin(state[7]) ;
cosa = cos(state[7]) ;
sinb = sin(state[8]) ;
cosb = cos(state[8]) ;
vel_u = state[6]*cosb*cosa ;
vel_v = state[6]*sinb ;
vel_w = state[6]*cosb*sina ;
u_dot =          cosb*cosa*xdot[6]
      - state[6]*sinb*cosa*xdot[8] 
      - state[6]*cosb*sina*xdot[7] ;
v_dot =          sinb*xdot[6] 
      + state[6]*cosb*xdot[8] ;
w_dot =          cosb*sina*xdot[6]
      - state[6]*sinb*sina*xdot[8] 
      + state[6]*cosb*cosa*xdot[7] ;

nx_cg = -1.0/grav*(u_dot + state[10]*vel_w - state[11]*vel_v)
      - sin(state[4]) ;
ny_cg = -1.0/grav*(v_dot + state[11]*vel_u - state[9]*vel_w)
      + cos(state[4])*sin(state[3]) ;
nz_cg = -1.0/grav*(w_dot + state[9]*vel_v - state[10]*vel_u)
      + cos(state[4])*cos(state[3]) ;

y[0] = nx_cg ;
y[1] = ny_cg ;
y[2] = nz_cg ;


} 

/*########################################*/
/*########################################*/

/*########################################*/
/*########################################*/

