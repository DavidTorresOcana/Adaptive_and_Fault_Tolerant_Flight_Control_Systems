%% save the values of the neural networks
% W_roll_last
% V_roll_last
% 
% Theta_1_roll_last
% Theta_2_roll_last
% Theta_3_roll_last
% 
% V_long_lastAoA
% W_long_lastAoA
% V_long_lastA_w_z
% W_long_lastA_w_z
% 
% V_yaw_last
% W_yaw_last


% % save('NN_i_weights','W_roll_last','V_roll_last','Theta_3_roll_last','Theta_2_roll_last','Theta_1_roll_last','V_long_lastAoA','W_long_lastAoA',...
% %     'V_long_lastA_w_z','W_long_lastA_w_z','V_yaw_last','W_yaw_last');

%% Cross NN
V_Comp_last
W_Comp_last
V_Comp_A_w_z_last
W_Comp_A_w_z_last

save('NN_c_weights','V_Comp_last','W_Comp_last','V_Comp_A_w_z_last','W_Comp_A_w_z_last');






%% Dont touch this!! 
% % % save('NN_weights_LAST_STABLE','W_roll_last','V_roll_last','Theta_3_roll_last','Theta_2_roll_last','Theta_1_roll_last','V_long_lastAoA','W_long_lastAoA',...
% % %     'V_long_lastA_w_z','W_long_lastA_w_z','V_yaw_last','W_yaw_last');
