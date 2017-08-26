%% Comparison between Basic and Advanced

% Responses_Ad_BFAULT = Responses_out
% ALL_Responses_Ad_BFAULT = ALL_Responses_out
% Control_defl_Ad_BFAULT = Control_defl
% NN_output_Ad_BFAULT = NN_output

% load('Fault_damage_Compar1.mat')
load('Fault_damage_Compar2.mat')

close all

fi_po=min( Responses_Ad_BFAULT.Length,Responses_FT_BFAULT.Length);

figure

subplot(3,1,1)
plot(Responses_Ad_BFAULT.Time(1:fi_po),[Responses_Ad_BFAULT.Data(1:fi_po,3),Responses_Ad_BFAULT.Data(1:fi_po,4),Responses_FT_BFAULT.Data(1:fi_po,4)])
title('\fontsize{14} AoA responses ')
axis([0 Responses_Ad_BFAULT.Time(fi_po) -5 10])
legend('\fontsize{12} Comm','\fontsize{12} Basic','\fontsize{12} Advanced')
ylabel('\fontsize{12} AoA (º)')
grid

subplot(3,1,2)
Color = lines(4);
plot(Responses_Ad_BFAULT.Time(1:fi_po),Responses_Ad_BFAULT.Data(1:fi_po,5),'Color',Color(1,:))
hold on
plot(Responses_Ad_BFAULT.Time(1:fi_po),Responses_Ad_BFAULT.Data(1:fi_po,6),'Color',Color(2,:))

title('\fontsize{14} Roll rate response Basic controller ')
axis([0 Responses_Ad_BFAULT.Time(fi_po) -30 30])
legend('\fontsize{12} Comm','\fontsize{12} Resp')
ylabel('\fontsize{12} Roll rate p (º/s)')
grid


subplot(3,1,3)
Color = lines(4);
hold on;
plot(Responses_Ad_BFAULT.Time(1:fi_po),Responses_FT_BFAULT.Data(1:fi_po,5),'Color',Color(1,:))
plot(Responses_Ad_BFAULT.Time(1:fi_po),Responses_FT_BFAULT.Data(1:fi_po,6),'Color',Color(3,:))
axis([0 Responses_Ad_BFAULT.Time(fi_po) -30 30])
title('\fontsize{14} Roll rate response Advanced controller ')
legend('\fontsize{12} Comm','\fontsize{12} Resp')
ylabel('\fontsize{12} Roll rate p (º/s)')
grid

xlabel('\fontsize{12} Time (s)')

return
%% PErformance Fault tolerant!
% Loss of both ailerons!!

% plot(V_online) % not useful
figure

subplot(2,2,1)
plot(Responses_out.Time,Responses_out.Data(:,3:4))
legend('\fontsize{12} Roll rate comm (º/s)','\fontsize{12} Roll rate actual (º/s)')
ylabel('\fontsize{12} p (º/s)')
title('\fontsize{14} Roll rates' )
grid

subplot(2,2,3)
plot(Responses_out.Time,Responses_out.Data(:,5:6))
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










