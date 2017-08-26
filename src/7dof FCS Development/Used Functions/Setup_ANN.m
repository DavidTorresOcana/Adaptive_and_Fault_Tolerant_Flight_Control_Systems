%% This script makes the set-up of the Neural networks
%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Isolated chanels Neural networks paramrters %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%# ROLL Chanel
    %# Roll 1 Hidden layer. Sizes:
    %   inputs: 13 [1,altitude, Mach,Mach^2,alpha,alpha^2,beta,r, p_err,p_sf,p_sf_dot,v_adp]
    %   Hidden layer neurons=30. Activation fun: sigmoidal
    %   Output:1
% %     
% %     V_roll_last = 0.1*rand(30,12);
% %     W_roll_last = 0.1*rand(1,30+1);

                % % % % %     V_roll_lastMOD = 0.1*rand(30,8);
                % % % % %     W_roll_lastMOD = 0.1*rand(1,30+1);

    %# Roll 2 Hidden layers. Sizes:
    %   inputs: 13 [1,altitude, Mach,Mach^2,alpha,alpha^2,beta,r, p_err,p_sf,p_sf_dot,v_adp,]
    %   Hidden layer 1 neurons=30.  Activation fun: sigmoidal
    %   Hidden layer 1 neurons=20.  Activation fun: lineal
    %   Output: 1 
% 
% %     Theta_1_roll_last = 0.1*rand(30,12);
% %     Theta_2_roll_last = 0.1*rand(20,30+1);
% %     Theta_3_roll_last = 0.1*rand(1,20+1);

%# LONG Chanel
    %# Long in AoA 1 Hidden layer. Sizes:
    %   inputs: 11 [1,altitude, Mach , Mach^2 , alpha_err, q_err , alpha_sf , alpha_sf_dot, q_sf_dot, v_ad_alpha]
    %   Hidden layer neurons=30. Activation fun: sigmoidal
    %   Output:1
% %     
% %     V_long_lastAoA = 0.1*rand(30,10);
% %     W_long_lastAoA = 0.1*rand(1,30+1);

    %# Long in A_w_z 1 Hidden layer. Sizes:
    %   inputs: 11 [1,altitude, Mach , Mach^2 , alpha_err, q_err , alpha_sf , alpha_sf_dot, q_sf_dot, v_ad_alpha]
    %   Hidden layer neurons=30. Activation fun: sigmoidal
    %   Output:1
% %     
% %     V_long_lastA_w_z = 2*rand(30,10)-1;
% %     W_long_lastA_w_z = 2*rand(1,30+1)-1;


%# YAW Chanel
    %# YAW 1 Hidden layer. Sizes:
    %   inputs: 13 [1,altitude, Mach , Mach^2 , alpha ,p , beta_err, r_err, beta_sf , beta_sf_dot, r_sf_dot , v_ad_beta ]
    %   Hidden layer neurons=30. Activation fun: sigmoidal
    %  Output:1
% %     
% %     V_yaw_last = 0.1*rand(30,12);
% %     W_yaw_last = 0.1*rand(1,30+1);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Crosssed Chanels Neural networks paramrters %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All chanles are mixed here. Only a single NN
%   inputs: 25
%   [1,altitude, Mach,Mach^2,alpha,alpha^2,beta,p,r, [p_err,p_sf,p_sf_dot], [[alpha_err, q_err , alpha_sf , alpha_sf_dot, q_sf_dot] or [A_err, q_err , A_sf , A_sf_dot, q_sf_dot] ], [beta_err, r_err , beta_sf , beta_sf_dot, r_sf_dot],  [v_ad_p,v_ad_q,v_ad_r]]
%   Hidden layer neurons=50. Activation fun: sigmoidal
%   Output: 3: [v_ad_p,v_ad_q,v_ad_r]

% % % % V_Comp_last = 0.05*(2*rand(50,25)-1);
% % % % W_Comp_last = 0.05*(2*rand(3,51)-1);
% % % % 
% % % % 
% % % % 
% % % % V_Comp_A_w_z_last = 0.05*(2*rand(50,25)-1);
% % % % W_Comp_A_w_z_last = 0.05*(2*rand(3,51)-1);


% % V_Comp_last_ONE = 0.05*(2*rand(3,25)-1)