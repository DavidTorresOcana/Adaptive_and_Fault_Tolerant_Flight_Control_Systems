%% Main script to check specs compliance
close all


%% Partition of the FE
FE_partition


%% SPPO
% Execute for the different cases in lfight partitions
% variables used are: AoA_comm, AoA ,Pitch_rate, load factor
% The method is to inject a pulse in elevators
SPPO_analysis

%% DR analysis
% Execute for the different cases in flight partitions
% variables used are: p and q
% The method is to inject a pulse in rudder
DR_compl

%% Roll analysis
% Execute for the different cases in flight partitions
% variables used are: p 
Roll_compl

%% Roll bank angles
% Speeds
% Cat A


Bank_FE








