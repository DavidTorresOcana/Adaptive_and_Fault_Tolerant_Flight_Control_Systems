%% This script plot the results of the simulation when executed.
close all
Responses_out

%% Plot of the responses 
Acc_response=Responses_out.Data(:,1:2);
AoA_response=Responses_out.Data(:,3:4);
Roll_response=Responses_out.Data(:,5:6);
Beta_response=Responses_out.Data(:,7:8);
TAS_response=Responses_out.Data(:,9:10);

figure
subplot(4,1,1)
plot(Responses_out.Time,Acc_response)
legend('\fontsize{12} Acc comm (m/s^2)','\fontsize{12} Acc actual (m/s^2)')
title('\fontsize{16} Responses and commands');ylabel('\fontsize{12} ACC (m/s^2)')
subplot(4,1,2)
plot(Responses_out.Time,AoA_response)
legend('\fontsize{12} AoA comm (º)','\fontsize{12} AoA actual (º)')
ylabel('\fontsize{12} AoA (º)')
subplot(4,1,3)
plot(Responses_out.Time,Roll_response)
legend('\fontsize{12} Roll rate comm (º/s)','\fontsize{12} Roll rate actual (º/s)')
ylabel('\fontsize{12} p (º/s)')
subplot(4,1,4)
plot(Responses_out.Time,Beta_response)
legend('\fontsize{12} SS comm (º)','\fontsize{12} SS actual (º)')
ylabel('\fontsize{12} \beta (º)')
% subplot(5,1,5)
% plot(Responses_out.Time,TAS_response)
% legend('\fontsize{12} TAS des (m/s)','\fontsize{12} TAS actual (m/s)')
xlabel('\fontsize{12} Time (s)')

%% PErformance adaptive!
%% Mass changes
% W_online
% V_online
% NN_output
% % Comm_ang_acce
% % 
% % % plot(V_online) % not useful
% % figure
% % subplot(2,1,1)
% % plot(NN_output.Time,NN_output.Data)
% % title('\fontsize{14} NN outputs' )
% % legend('\fontsize{12} Roll channel','\fontsize{12} Long channel','\fontsize{12} Yaw channel')
% % ylabel('\fontsize{12} Ang. acc. comm (rad/s^2)')
% % 
% % subplot(2,1,2)
% % plot(Comm_ang_acce.Time,Comm_ang_acce.Data)
% % title('\fontsize{14} Total Comm. Ang. accelerations' )
% % legend('\fontsize{12} Roll channel','\fontsize{12} Long channel','\fontsize{12} Yaw channel')
% % xlabel('\fontsize{12} Time (s)');ylabel('\fontsize{12} Ang. acc. comm (rad/s^2)')

%% Loss of aileron


% plot(V_online) % not useful
figure
subplot(3,1,1)
plot(Responses_out.Time,Roll_response)
legend('\fontsize{12} Roll rate comm (º/s)','\fontsize{12} Roll rate actual (º/s)')
ylabel('\fontsize{12} p (º/s)')
title('\fontsize{14} Roll rates' )

subplot(3,1,2)
try
    plot(NN_output.Time,NN_output.Data)
catch
    disp('Could nto load the data')
end
title('\fontsize{14} NN outputs: Only for Adaptive Configs' )
legend('\fontsize{12} Roll channel','\fontsize{12} Long channel','\fontsize{12} Yaw channel')
ylabel('\fontsize{12} Ang. acc. comm (rad/s^2)')

subplot(3,1,3)
try
    plot(Comm_ang_acce.Time,Comm_ang_acce.Data)
catch
    disp('Could nto load the data')
end
title('\fontsize{14} Total Comm. Ang. accelerations Only for Adaptive Configs' )
legend('\fontsize{12} Roll channel','\fontsize{12} Long channel','\fontsize{12} Yaw channel')
xlabel('\fontsize{12} Time (s)');ylabel('\fontsize{12} Ang. acc. comm (rad/s^2)')
return

%% PErformance Fault tolerant!
% Loss of both ailerons!!

% plot(V_online) % not useful
figure
subplot(2,2,1)
plot(Responses_out.Time,Roll_response)
legend('\fontsize{12} Roll rate comm (º/s)','\fontsize{12} Roll rate actual (º/s)')
ylabel('\fontsize{12} p (º/s)')
title('\fontsize{14} Roll rates' )
grid
subplot(2,2,2)
plot(NN_output.Time,NN_output.Data)
title('\fontsize{14} NN outputs' )
legend('\fontsize{12} Roll channel','\fontsize{12} Long channel','\fontsize{12} Yaw channel')
ylabel('\fontsize{12} Ang. acc. comm (rad/s^2)')
grid
subplot(2,2,4)
plot(Comm_ang_acce.Time,Comm_ang_acce.Data)
title('\fontsize{14} Total Comm. Ang. accelerations' )
legend('\fontsize{12} Roll channel','\fontsize{12} Long channel','\fontsize{12} Yaw channel')
xlabel('\fontsize{12} Time (s)');ylabel('\fontsize{12} Ang. acc. comm (rad/s^2)')
grid

subplot(2,2,3)
plot(Control_defl.Time,Control_defl.Data(:,[2,4,6,8,10]))
title('\fontsize{14} Control surface deflections ' )
legend('\fontsize{12} Left Elev','\fontsize{12} Right Elev','\fontsize{12} Left Ail','\fontsize{12} Right Ail','\fontsize{12} Rudder','\fontsize{12} Left LEF','\fontsize{12} Right LEF')
xlabel('\fontsize{12} Time (s)');ylabel('\fontsize{12} Deflection (º)')
grid
