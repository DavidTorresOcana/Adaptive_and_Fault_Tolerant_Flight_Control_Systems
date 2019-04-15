%% This script trim the Symetric model (3dof) in steady flight, generates 
%  a linear model or linearized model for CA and start simulation for asymetric(7dof) model
close all
clear all
clc;
addpath(genpath('aerodata'),genpath('Used Functions'),addpath( genpath('../../lib')))
% mex nlplantASYM.c
% mex nlplantASYM_Reconf.c
% mex Moms_Coeff_Reconf_estimator.c
% mex Moms_Coeff_Grad_Reconf_estimator.c
% mex Ang_accels_Reconf_estimator.c
%% General params
T_sim=0.02; % Sample time of the simulation
% T_sim=0.04643990929705216 ; % Sample time of stall warning

%% Trim aircraft (Symetric model)
% This rutine trim the aircraft at any flight condition
runF16Sim;

% pause;
%% Settings of the CA
% Can choose either Linearized B around trim flight or on-line linearized
CA_settings

% F16ASYM_Controlled

%% NN settings and initialization
Setup_ANN

%% Initialize Control GUI

% Control_GUI

%% Stall protection: in load factor and AoA
Stall_protection

close all
%% Auto-trim aid
load('Autotrim_maps.mat')
% Auto_trim_aid

%% Scheduling in AoA chanel
load('Scheduling_map_AoA.mat')
tau_q_extrap_AoA = tau_q_extrap; % Madly important
K_b_q_extrap_AoA  = K_b_q_extrap;
surf(M_eval_sche_unc,H_eval_sche_unc,tau_AoA_extrap);
surf(M_eval_sche_unc,H_eval_sche_unc,tau_q_extrap_AoA);
surf(M_eval_sche_unc,H_eval_sche_unc,K_b_AoA_extrap);
surf(M_eval_sche_unc,H_eval_sche_unc,K_b_q_extrap_AoA);


%Plot
% plot(F_envelope.M_1g ,F_envelope.H_1g/1000,'--r','LineWidth',2)
% h=scatter3(reshape(M_eval_sche_unc,1,210),reshape(H_eval_sche_unc/1000,1,210),reshape(tau_AoA_extrap,1,210),[],reshape(tau_AoA_extrap,1,210),'o','filled');
% hChildren = get(h, 'Children');
% set(hChildren, 'Markersize', 15)
hold on

surf(M_eval_sche_unc,H_eval_sche_unc/1000,K_b_q_extrap_AoA);
title('\fontsize{16} K_{q} Map')
% xlabel('\fontsize{12} Mach');ylabel('\fontsize{12} Altitude (km)');zlabel('\fontsize{12} K_{q} ')
view(0,90)
%% Scheduling in ACC chanel
load('Scheduling_map_ACC.mat')
tau_q_extrap_ACC = tau_q_extrap; % Madly important
K_b_q_extrap_ACC  = K_b_q_extrap;

surf(M_eval_sche,H_eval_sche,tau_ACC_extrap);
surf(M_eval_sche,H_eval_sche,tau_q_extrap_ACC);
surf(M_eval_sche,H_eval_sche,K_b_ACC_extrap);
surf(M_eval_sche,H_eval_sche,K_b_q_extrap_ACC);
close all


%% IMU params load
IMU_param;


%% Run simulation and development environment
open_system( 'F16ASYM_Controlled' )
