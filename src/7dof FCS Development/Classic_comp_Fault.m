%% Classical controller comparison with Basic
close all
clc
start=505;
%% Classical
% Responses_out_Classic_DP_fault=Responses_out
load('Classic_Comparison_Fault.mat')
start3=900;
fin3=1350;

figure
subplot(2,1,1)
Color=lines(8);
plot(Responses_out_Classic_DP.Time(start:end),Responses_out_Classic_DP.Data(start:end,5),'Color',Color(1,:))
hold on
plot(Responses_out_Classic_DP.Time(start:end),Responses_out_Classic_DP.Data(start:end,6),'Color',Color(2,:))
plot(Responses_out_Classic_DP.Time(start:end),Responses_out_Classic_DP_fault.Data(start:end,6),'Color',Color(3,:))

legend('\fontsize{12} Commands ','\fontsize{12} Normal','\fontsize{12} Fault ')
ylabel('\fontsize{12} p (º/s)')
title('\fontsize{14} Roll rates Classical approach Faults' )
grid



Response_spec_Classic_Normal=stepinfo(-Responses_out_Classic_DP.Data(start3:fin3,6)+Responses_out_Classic_DP.Data(start3,6), Responses_out_Classic_DP.Time(start3:fin3) , -Responses_out_Classic_DP.Data(fin3,5)+Responses_out_Classic_DP.Data(start3,5)  )
Response_spec_Classic_Fault=stepinfo(-Responses_out_Classic_DP_fault.Data(start3:fin3,6)+Responses_out_Classic_DP_fault.Data(start3,6), Responses_out_Classic_DP_fault.Time(start3:fin3) , -Responses_out_Classic_DP_fault.Data(fin3,5)+Responses_out_Classic_DP_fault.Data(start3,5) )

% return
%% Advandes and basic concepts
% Responses_out_Bas_ad_fault=Responses_out

% Obtain desired response
tau=0.25;
K_ref=1/tau;	
TF = tf(K_ref,[1 K_ref]);

Response_ideal = lsim(TF,Responses_out.Data(:,5),Responses_out.Time );


subplot(2,1,2)
Color=lines(8);
plot(Responses_out.Time(start:end),Responses_out.Data(start:end,5),'Color',Color(1,:))
hold on
plot(Responses_out.Time(start:end),Response_ideal(start:end),'Color',Color(2,:))
plot(Responses_out.Time(start:end),Responses_out_DP.Data(start:end,6),'Color',Color(3,:))
plot(Responses_out.Time(start:end),Responses_out_Ad_fault.Data(start:end,6),'Color',Color(4,:))
plot(Responses_out.Time(start:end),Responses_out_Bas_ad_fault.Data(start:end,6),'Color',Color(5,:))
plot(Responses_out.Time(start:end),Responses_out_Bas_noad_fault.Data(start:end,6),'Color','k')

legend('\fontsize{12} Commands ','\fontsize{12} Desired response','\fontsize{12} Normal response ',...
        '\fontsize{12} Fault: Advanced','\fontsize{12} Fault: Basic ','\fontsize{12} Fault: Basic No adapt ')
ylabel('\fontsize{12} p (º/s)')
xlabel('\fontsize{12} Time (s)')
title('\fontsize{14} Roll rates Basic&Advanced Faults' )
grid


Response_spec_Normal=stepinfo(-Responses_out_DP.Data(start3:fin3,6)+Responses_out_DP.Data(start3,6), Responses_out_DP.Time(start3:fin3) , -Responses_out.Data(fin3,5)+Responses_out.Data(start3,5)  )

Response_spec_Ad_fault=stepinfo(-Responses_out_Ad_fault.Data(start3:fin3,6)+Responses_out_Ad_fault.Data(start3,6), Responses_out_Ad_fault.Time(start3:fin3) , -Responses_out.Data(fin3,5)+Responses_out.Data(start3,5) )
Response_spec_Bas_ad_fault=stepinfo(-Responses_out_Bas_ad_fault.Data(start3:fin3,6)+Responses_out_Bas_ad_fault.Data(start3,6), Responses_out_Bas_ad_fault.Time(start3:fin3) ,  -Responses_out.Data(fin3,5)+Responses_out.Data(start3,5) )
Response_spec_Bas_noad_fault=stepinfo(-Responses_out_Bas_noad_fault.Data(start3:fin3,6)+Responses_out_Bas_noad_fault.Data(start3,6), Responses_out_Bas_noad_fault.Time(start3:fin3) , -Responses_out.Data(fin3,5)+Responses_out.Data(start3,5) )


