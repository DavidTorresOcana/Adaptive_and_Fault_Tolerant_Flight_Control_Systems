
%% Comparison of the classical controller with basic
load('Classic_Comparison.mat')

figure
subplot(2,1,1)
Color=lines(4);
plot(Responses_out_Classic_DP.Time,Responses_out_Classic_DP.Data(:,5),'Color',Color(1,:))
hold on
plot(Responses_out_Classic_DP.Time,Responses_out_Classic_DP.Data(:,6),'Color',Color(2,:))
plot(Responses_out_Classic_DP.Time,Responses_out_Classic_Back_GS.Data(:,6),'Color',Color(3,:))
plot(Responses_out_Classic_DP.Time,Responses_out_Classic_Back_DPout.Data(:,6),'Color',Color(4,:))

legend('\fontsize{12} Commands ','\fontsize{12} Design Point','\fontsize{12} Gain Scheduling ','\fontsize{12} Off-design Point ')
ylabel('\fontsize{12} p (º/s)')
title('\fontsize{14} Roll rates' )
grid

% % subplot(2,1,2)
% % plot(Control_defl_Classic_DP.Time,[Control_defl_Classic_DP.Data(:,[6,8]),Control_defl_Classic_Back_GS.Data(:,[6,8]),Control_defl_Classic_Back_DPout.Data(:,[6,8])])
% % title('\fontsize{14} Control surface deflections ' )
% % legend('\fontsize{12} Left Elev','\fontsize{12} Right Elev','\fontsize{12} Left Ail','\fontsize{12} Right Ail','\fontsize{12} Rudder','\fontsize{12} Left LEF','\fontsize{12} Right LEF')
% % xlabel('\fontsize{12} Time (s)');ylabel('\fontsize{12} Deflection (º)')
% % grid
