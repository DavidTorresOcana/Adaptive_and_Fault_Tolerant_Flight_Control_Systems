%% This script explain the way used to protect from stall
% The idea is to generate a parameter that multiplies the command in
% Acceleration or AoA and contemplate the AoA and closeness to stall
% clear all
close all
clc
%% Load factor protection
% % F_envelope.M_1g=M;
% % F_envelope.H_1g=Hm;
% % F_envelope.M_3g=M_2(1:end-4);
% % F_envelope.H_3g=Hm_2(1:end-4);
% % F_envelope.M_5g=M_3(1:end-3);
% % F_envelope.H_5g=Hm_3(1:end-3);
% % F_envelope.M_7g=M_4;
% % F_envelope.H_7g=Hm_4;
% % F_envelope.M_9g=M_8;
% % F_envelope.H_9g=FlightEnvelope

% F_envelope.M_9g=[F_envelope.M_9g(1);F_envelope.M_9g;F_envelope.M_9g(end)]
% F_envelope.H_9g=[0;F_envelope.H_9g;0]

load('FlightEnvelope.mat')
plot(F_envelope.M_1g ,F_envelope.H_1g )
hold on
plot(F_envelope.M_3g ,F_envelope.H_3g,'r' )
plot(F_envelope.M_5g ,F_envelope.H_5g,'g' )
plot(F_envelope.M_7g ,F_envelope.H_7g,'m' )
plot(F_envelope.M_9g ,F_envelope.H_9g,'k' )

F_envelope3d.M =[F_envelope.M_1g;F_envelope.M_3g;F_envelope.M_5g;F_envelope.M_7g;F_envelope.M_9g];
F_envelope3d.H =[F_envelope.H_1g;F_envelope.H_3g;F_envelope.H_5g;F_envelope.H_7g;F_envelope.H_9g];
F_envelope3d.g = [1*ones(size(F_envelope.M_1g));3*ones(size(F_envelope.M_3g));5*ones(size(F_envelope.M_5g));7*ones(size(F_envelope.M_7g));9*ones(size(F_envelope.M_9g))];

[fitresult, gof]=Interpolate_FE(F_envelope3d.M,F_envelope3d.H,F_envelope3d.g);
view(0,90)

% Generate a grid of points
M_sq_prev=[];
H_sq_prev=[];
% M_sq_prev = 0:0.2:max(F_envelope.M_1g);
% H_sq_prev = min( F_envelope.H_1g):2000:max(F_envelope.H_1g);
% % M_sq = [M_sq,(F_envelope3d.M)'];
% % H_sq = [H_sq,(F_envelope3d.H)'];
M_sq_prev = [M_sq_prev,(F_envelope.M_1g)'];
H_sq_prev = [H_sq_prev,(F_envelope.H_1g)'];

M_sq_prev = unique(sort(M_sq_prev)); H_sq = unique(sort(H_sq_prev));
    % Make data evenly distributed
    k=0;
    for i=2:size(M_sq_prev,2)
        if M_sq_prev(i)>M_sq_prev(i-1)+0.01 
            M_sq(i-k)=M_sq_prev(i);
        else
            k=k+1;
        end
    end
    k=0;
    for i=2:size(M_sq_prev,2)
        if H_sq_prev(i)>H_sq_prev(i-1)+50
            H_sq(i-k)=H_sq_prev(i);
        else
            k=k+1;
        end
    end
    M_sq =unique(sort(M_sq));
    H_sq =unique(sort(H_sq));
    
% Eval inside?
M_eval=NaN*ones(size(M_sq,2),size(H_sq,2));
H_eval=NaN*ones(size(M_sq,2),size(H_sq,2));
G_eval = NaN*ones(size(M_sq,2),size(H_sq,2));
for I=1:size(M_sq,2)
    for J=1:size(H_sq,2)
        if inpolygon(M_sq(I),H_sq(J),F_envelope.M_1g,F_envelope.H_1g) 
            M_eval(I,J)=M_sq(I);
            H_eval(I,J)=H_sq(J);
            G_eval(I,J) = fitresult(M_sq(I),H_sq(J));
        end
    end
end


hold on
surf(M_eval,H_eval,G_eval)
view(0,90)
close all
            
%% AoA  protection

Stall_curve.H = [0 3000 5000 9000  12000 15000];
Stall_curve.M = [0.15  0.18 0.23 0.3 0.41 0.59];
% % plot(Stall_curve.M ,Stall_curve.H )
% hold on
Stall_curve.H_2 = F_envelope.H_1g(1:14);
Stall_curve.M_2 = F_envelope.M_1g(1:14);

Stall_curve.M_2(1)=0.99*Stall_curve.M_2(1);
Stall_curve.M_2=unique(sort(Stall_curve.M_2));
Stall_curve.H_2=unique(sort(Stall_curve.H_2));


%% Switching between AoA and load factor

plot(Stall_curve.M_2 ,(Stall_curve.H_2) ,'r')
hold on
 plot(F_envelope.M_1g ,F_envelope.H_1g)
plot(Stall_curve.M_2*2 ,(Stall_curve.H_2) ,'g')


% Estimation lift curve
mex Lift_Drag_curves_reconf.c

Inputs = [  trim_state(3),trim_state(7),trim_state(8),trim_state(9),trim_state(10:12)'  ,  trim_throttle  ,trim_control(1),trim_control(1),trim_control(2),-trim_control(2),trim_control(3),0,0]';

Delta_S_R =0;
Delta_S_L =0;
Health_state= [9295.44,27.87/2*0.3 + 0.7*27.87/2*(1-Delta_S_L),27.87/2*0.3 + 0.7*27.87/2*(1-Delta_S_R),6.578,0.3,75673.6,1331.4,85552.1,12874.8, [0,0,0,0,0,0,0],zeros(1,3)]; %TBC/
Inputs = [Inputs;Health_state'];

tic
    Lift_Drag_curves_reconf(Inputs)
toc

AoA= deg2rad( [-10,-7,-5,-2,0,2,4,6,8,10,12,14,16,18,20] );

for i=1:size(AoA,2)
    Inputs = [  trim_state(3),trim_state(7),AoA(i),trim_state(9),trim_state(10:12)'  ,  trim_throttle  ,trim_control(1),trim_control(1),trim_control(2),-trim_control(2),trim_control(3),0,0]';
    Inputs = [Inputs;Health_state'];
    outp=Lift_Drag_curves_reconf(Inputs);
    C_L(i)=outp(1);
end

plot(rad2deg(AoA),C_L)

[fitresult, gof] = Fit_lift_curve(AoA, C_L);
C_L_0 = fitresult.p2;
a1 = fitresult.p1;










% % % % % % % AoA_min=-deg2rad(20); % rad
% % % % % % % AoA_max = deg2rad(20); % rad
% % % % % % % x = 1.01*AoA_min:0.001:1.01*AoA_max;
% % % % % % % 
% % % % % % % C =1;
% % % % % % % a=1;
% % % % % % % 
% % % % % % % C=[0.01 0.01];
% % % % % % % a=[0.5 1 2];
% % % % % % % hold on
% % % % % % % color= lines(size(C,2)*size(a,2));
% % % % % % % str={};
% % % % % % % k=0;
% % % % % % % for i=1:size(C,2)
% % % % % % %     for j=1:size(a,2)
% % % % % % %         k=k+1;
% % % % % % %         if mod(a(j),2)
% % % % % % %             y= 1 - C(i)./(x-AoA_min).^a(j)+C(i)./(x-AoA_max).^a(j);
% % % % % % %         else
% % % % % % %             y= 1 - C(i)./(x-AoA_min).^a(j)-C(i)./(x-AoA_max).^a(j);
% % % % % % %         end
% % % % % % %         str=[str;['C=',num2str(C(i)),' a=',num2str(a(j))]];
% % % % % % %         plot(x,y,'Color',color(k,:))
% % % % % % %     end
% % % % % % % end
% % % % % % % legend(str)
% % % % % % % axis([AoA_min AoA_max ,-1 2])
% % % % % % % 
% % % % % % % 
% % % % % % % y= 1 - C(1)./(x-AoA_min).^a(2) + C(1)./(x-AoA_max).^a(2);
% % % % % % % plot(x,y)
