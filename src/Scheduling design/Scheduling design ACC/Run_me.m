%% This script trim the Symetric model (3dof) in steady flight, generates 
%  a linear model or linearized model for CA and start simulation for asymetric(7dof) model
close all
clc;
addpath(genpath('aerodata'),genpath('Used Functions'))
mex Moms_Coeff_Reconf_estimator.c
mex Moms_Coeff_Grad_Reconf_estimator.c
mex Ang_accels_Reconf_estimator.c
%% General params
T_sim=0.04; % Sample time of the simulation
% T_sim=0.04643990929705216 ; % Sample time of the simulation

%% Trim aircraft (Symetric model)
% This rutine trim the aircraft at any flight condition
runF16Sim;

pause

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
