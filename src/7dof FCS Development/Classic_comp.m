%% Classical controller comparison with Basic
close all
start=505;
%% Classical
load('Classic_Comparison.mat')

figure
subplot(2,1,1)
Color=lines(5);
plot(Responses_out_Classic_DP.Time(start:end),Responses_out_Classic_DP.Data(start:end,5),'Color',Color(1,:))
hold on
plot(Responses_out_Classic_DP.Time(start:end),Responses_out_Classic_DP.Data(start:end,6),'Color',Color(2,:))
plot(Responses_out_Classic_DP.Time(start:end),Responses_out_Classic_Back_GS.Data(start:end,6),'Color',Color(3,:))
plot(Responses_out_Classic_DP.Time(start:end),Responses_out_Classic_Back_DPout.Data(start:end,6),'Color','k')

legend('\fontsize{12} Commands ','\fontsize{12} Design Point','\fontsize{12} Gain scheduling ','\fontsize{12} Off-design Point ')
ylabel('\fontsize{12} p (º/s)')
title('\fontsize{14} Roll rates Classical approach' )
grid

% Response params
start2=550;
fin=1000;
plot(Responses_out_Classic_DP.Time(start2:fin),Responses_out_Classic_DP.Data(start2:fin,5),'Color',Color(1,:))
hold on
plot(Responses_out_Classic_DP.Time(start2:fin),Responses_out_Classic_DP.Data(start2:fin,6),'Color',Color(2,:))
plot(Responses_out_Classic_DP.Time(start2:fin),Responses_out_Classic_Back_GS.Data(start2:fin,6),'Color',Color(3,:))
plot(Responses_out_Classic_DP.Time(start2:fin),Responses_out_Classic_Back_DPout.Data(start2:fin,6),'Color','k')

Classic_response_DP=stepinfo(Responses_out_Classic_DP.Data(start2:fin,6), Responses_out_Classic_DP.Time(start2:fin) , Responses_out_Classic_DP.Data(fin,5) )
Classic_response_GS=stepinfo(Responses_out_Classic_Back_GS.Data(start2:fin,6), Responses_out_Classic_Back_GS.Time(start2:fin) , Responses_out_Classic_DP.Data(fin,5) )
Classic_response_DPout=stepinfo(Responses_out_Classic_Back_DPout.Data(start2:fin,6), Responses_out_Classic_Back_DPout.Time(start2:fin) , Responses_out_Classic_DP.Data(fin,5) )





% subplot(2,1,2)
% plot(Control_defl_Classic_DP.Time,[Control_defl_Classic_DP.Data(:,[6,8]),Control_defl_Classic_Back_GS.Data(:,[6,8]),Control_defl_Classic_Back_DPout.Data(:,[6,8])])
% title('\fontsize{14} Control surface deflections ' )
% legend('\fontsize{12} Left Elev','\fontsize{12} Right Elev','\fontsize{12} Left Ail','\fontsize{12} Right Ail','\fontsize{12} Rudder','\fontsize{12} Left LEF','\fontsize{12} Right LEF')
% xlabel('\fontsize{12} Time (s)');ylabel('\fontsize{12} Deflection (º)')
% grid



%% Basic concept
% Responses_out_off_DP=Responses_out;

% Obtain desired response
tau=0.25;
K_ref=1/tau;	
TF = tf([K_ref],[1 K_ref]);

Response_ideal = lsim(TF,Responses_out.Data(:,5),Responses_out.Time );


subplot(2,1,2)
Color=lines(5);
plot(Responses_out.Time(start:end),Responses_out.Data(start:end,5),'Color',Color(1,:))
hold on
plot(Responses_out.Time(start:end),Response_ideal(start:end),'Color',Color(2,:))
plot(Responses_out.Time(start:end),Responses_out_DP.Data(start:end,6),'Color',Color(3,:))
plot(Responses_out_off_DP.Time(start:end),Responses_out_off_DP.Data(start:end,6),'Color','k')

legend('\fontsize{12} Commands ','\fontsize{12} Desired response','\fontsize{12} Design point ','\fontsize{12} Off-design Point ')
ylabel('\fontsize{12} p (º/s)')
xlabel('\fontsize{12} Time (s)')

title('\fontsize{14} Roll rates Basic concept' )
grid

Basic_response_DP=stepinfo(Responses_out_DP.Data(start2:fin,6), Responses_out_DP.Time(start2:fin) , Responses_out.Data(fin,5) )
Basic_response_off_DP=stepinfo(Responses_out_off_DP.Data(start2:fin,6), Responses_out_off_DP.Time(start2:fin) , Responses_out.Data(fin,5) )

