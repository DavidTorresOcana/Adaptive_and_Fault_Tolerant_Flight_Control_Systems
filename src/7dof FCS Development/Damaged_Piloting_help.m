%% This script is intended to help the pilot when damage or big deviation occurs, 
% since the way of flying has to be changed. For instantce, when a % of
% wing is loss, the aircraft cannot fly with zero sideslip

close all
clc;
addpath(genpath('aerodata'),genpath('Used Functions'))
%% Inputs:
%           altitude (m)
%           V_TAS (m/s)
% 
% Inputs:   alpha (rad)
%           beta (rad)
%           P (rad/s)
%           Q (rad/s)
%           R (rad/s)
%           
%           Throttle [0-1]
% 
%           elevator_left (deg)
%           elevator_right (deg)
%           aileron_left (deg)
%           aileron_right (deg)
%           rudder (deg)
%           LEF_left (deg)  (additionally to the symetric defflection)
%           LEF_right (deg) (additionally to the symetric defflection)
%
%           mass  Estimated mass
%           S_current_L
%           S_current_R
%           S_fin_current
%           x_cg
%           Jy
%           Jxz
%           Jz
%           Jx
%           EFF_el_L
%           EFF_el_R
%           EFF_ail_L
%           EFF_ail_R
%           EFF_rud
%           EFF_lef_L
%           EFF_lef_R
%           Delta_C_Lift_X
%           Delta_C_Drag_X
%           Delta_C_m_X
%
% Outputs   C_L
%           C_M 
%           C_N
mex Moms_Coeff_Reconf_estimator.c
% Test
Inputs=[4400,  160,deg2rad(2),0  ,0,0,0,  0.26,   -2,-2,  0,0,  0,  2,2]';

% Inputs fail 
Delta_S_R =0;
Delta_S_L =0.4;
Health_state= [9295.44,27.87/2*0.3 + 0.7*27.87/2*(1-Delta_S_L),27.87/2*0.3 + 0.7*27.87/2*(1-Delta_S_R),6.578,0.3,75673.6,1331.4,85552.1,12874.8, [1,1,0,1,1,0,1],zeros(1,3)]; %TBC/
Inputs = [Inputs;Health_state'];
% Test
tic
    Moms_Coeff_Reconf_estimator(Inputs)
toc
return
%% Trim with any configuration. Alpha not constant
% The variables to be trimmed are: alpha, beta ,elevator_left, elevator_right, Ail_L,Ail_r, rudder,LEF_L,LEF_r

Cost = @(VARS) 1000*norm( Moms_Coeff_Reconf_estimator( [ trim_state(3),trim_state(7),VARS(1),VARS(2),trim_state(10:12)' , trim_throttle ,VARS(3),VARS(4),VARS(5),VARS(6),VARS(7),VARS(8),VARS(9),Health_state]' ) ) ;

lims=[ -deg2rad(20),-deg2rad(20),-24,-24,-21,-21,-29,-5,-5 ;...
    deg2rad(20),deg2rad(20),24,24,21,21,29,25,25 ];
[EQ_vars,cost]=fmincon(Cost,[trim_state(8:9)',trim_control(1),trim_control(1),trim_control(2),-trim_control(2),trim_control(3),0,0],[],[],[],[],lims(1,:),lims(2,:),[]);

cost
AoA_eq=EQ_vars(1)
Beta_eq=EQ_vars(2)
Constrols_eq= EQ_vars(3:end)

%% Trim with any configuration. Alpha to be constant=trim AoA
% The variables to be trimmed are:  beta ,elevator_left, elevator_right, Ail_L,Ail_r, rudder,LEF_L,LEF_r

% % Cost = @(VARS) 1000*norm( Moms_Coeff_Reconf_estimator( [ trim_state(3),trim_state(7),0.01,VARS(1),trim_state(10:12)' , trim_throttle ,VARS(2),VARS(3),VARS(4),VARS(5),VARS(6),VARS(7),VARS(8),Health_state]' ) ) ;
Cost = @(VARS) 1000*norm( Moms_Coeff_Reconf_estimator( [ trim_state(3),trim_state(7),trim_state(8),VARS(1),trim_state(10:12)' , trim_throttle ,VARS(2),VARS(3),VARS(4),VARS(5),VARS(6),VARS(7),VARS(8),Health_state]' ) ) ;

lims=[ -deg2rad(40),-24,-24,-21,-21,-29,-5,-5 ;...
    deg2rad(40),24,24,21,21,29,25,25 ];
[EQ_vars,cost]=fmincon(Cost,[20*trim_state(9)',trim_control(1),trim_control(1),trim_control(2),-trim_control(2),trim_control(3),0,0],[],[],[],[],lims(1,:),lims(2,:),[]);

cost
Beta_eq=EQ_vars(1)
Constrols_eq= EQ_vars(3:end)











