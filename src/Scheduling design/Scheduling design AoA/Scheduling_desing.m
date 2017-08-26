%% This scrip is intended to iteratively find the best parameters set for 
% each setpoint
% clear all
close all
clc
clear Model_data

global altitude fi_type velocity fi_flag_Simulink 

addpath(genpath('aerodata'),genpath('Used Functions'))
mex Moms_Coeff_Reconf_estimator.c
mex Moms_Coeff_Grad_Reconf_estimator.c
mex Ang_accels_Reconf_estimator.c
mex nlplantASYM.c
mex Moms_Coeff_Reconf_estimator.c
mex Moms_Coeff_Grad_Reconf_estimator.c
mex Moms_Coeff_Grad_estimator.c
mex Moms_Coeff_estimator.c
T_sim=0.02; % Sample time of the simulation

Model_used_flag = 1;
%% Flight envelope
% Using plottodata
load('FlightEnvelope.mat')

plot(F_envelope.M_1g,F_envelope.H_1g)
title('\fontsize{14} Flight Envelope and available data')
xlabel('Mach')
ylabel('H (m)')
hold on
%% Avialable model data 
% Speed from 90 to 270 m/s
% Altitude form 1550 to 12200
v= 90:10:270;
alt = 1500:200:12200;
for I=1:size(alt,2)
    [T,a]=atmosisa(alt(I));
    if inpolygon(v(1)/a,alt(I),F_envelope.M_1g,F_envelope.H_1g)
        Model_data.M(I)=v(1)/a;
        Model_data.H(I)=alt(I);
    else
        Model_data.M(I)=NaN;
        Model_data.H(I)=NaN;
    end
end
for I=1:size(v,2)
    [T,a]=atmosisa(alt(end));
    if inpolygon(v(I)/a,alt(end),F_envelope.M_1g,F_envelope.H_1g)
        Model_data.M(end+1)=v(I)/a;
        Model_data.H(end+1)=alt(end);
    end
end
for I=size(alt,2):-1:1
    [T,a]=atmosisa(alt(I));
    Model_data.M(end+1)=v(end)/a;
    Model_data.H(end+1)=alt(I);
end
for I=size(v,2):-1:1
    [T,a]=atmosisa(alt(1));
    Model_data.M(end+1)=v(I)/a;
    Model_data.H(end+1)=alt(1);
end

plot(Model_data.M,Model_data.H,'r')
%% Generation evauation points inside flight envelope
M_sq_sche = 0:0.1:max(F_envelope.M_1g);
H_sq_sche = min( F_envelope.H_1g):2000:max(F_envelope.H_1g);

M_eval_sche=NaN*ones(size(M_sq_sche,2),size(H_sq_sche,2));
H_eval_sche=NaN*ones(size(M_sq_sche,2),size(H_sq_sche,2));

for I=1:size(M_sq_sche,2)
    for J=1:size(H_sq_sche,2)
        if inpolygon(M_sq_sche(I),H_sq_sche(J),F_envelope.M_1g,F_envelope.H_1g)
            M_eval_sche(I,J)=M_sq_sche(I);
            H_eval_sche(I,J)=H_sq_sche(J);
        end
    end
end

% ADD ANY DESIRED POINT??
I=4;J=5;
M_eval_sche(I,J)=1.1*M_sq_sche(I);
H_eval_sche(I,J)=H_sq_sche(J);

I=3;J=2;
M_eval_sche(I,J)=1.2*M_sq_sche(I);
H_eval_sche(I,J)=H_sq_sche(J);

I=3;J=3;
M_eval_sche(I,J)=0.219;
H_eval_sche(I,J)=H_sq_sche(J);

I=7;J=8;
M_eval_sche(I,J)=M_sq_sche(I);
H_eval_sche(I,J)=13400;

I=9;J=9;
M_eval_sche(I,J)=M_sq_sche(I);
H_eval_sche(I,J)=15500;

plot(M_eval_sche,H_eval_sche,'*r')
return
%% Iterative evaluation of setpoints

% ANN parameters. To be optimized to give best performance in all flight
% envelope
learning = 20;
reg_param =0.5;

% Run firstly a test point
runF16Sim;
Stall_protection
close all
Setup_ANN
load('Autotrim_maps.mat')
Scheduling_map_AoA.tau_AoA = NaN*ones(size(M_sq_sche,2),size(H_sq_sche,2));
Scheduling_map_AoA.tau_q   = NaN*ones(size(M_sq_sche,2),size(H_sq_sche,2));
Scheduling_map_AoA.K_b_AoA = NaN*ones(size(M_sq_sche,2),size(H_sq_sche,2));
Scheduling_map_AoA.K_b_q   = NaN*ones(size(M_sq_sche,2),size(H_sq_sche,2));
%% Go interatively

for I=9:14%size(M_sq_sche,2)
    for J=1:2%size(H_sq_sche,2)
        if isnan(M_eval_sche(I,J))==0 && M_eval_sche(I,J) >= 1.1*min( M_eval_sche(:,J) )  &&isnan( Scheduling_map_AoA.tau_AoA(I,J)  )
            altitude = H_eval_sche(I,J);
            [T,a,~,rho]=atmosisa(altitude);
            velocity = M_eval_sche(I,J)*a;
            % Trim here
            q=1/2*rho*velocity^2;
            alpha = rad2deg((9800*9.81/(q*27)-0.015)/4.22);
            [trim_state, itrim_throttle, trim_control, dLEF, UX,cost] = trim_F16_fast(throttle, elevator,beta,  alpha, aileron, rudder, velocity, altitude);
            % CA here
            CA_settings
            
            I
            J
            pause(0.5)
            % Run optimization
            Find_opt
            % Retrieve the optimal
            Scheduling_map_AoA.tau_AoA(I,J) = params_opt(1);
            Scheduling_map_AoA.tau_q(I,J) = params_opt(2);
            Scheduling_map_AoA.K_b_AoA(I,J) = params_opt(3);
            Scheduling_map_AoA.K_b_q(I,J) = params_opt(4);
            close all
        end
    end
end


%% Evaluation and shape
surf(M_eval_sche,H_eval_sche,Scheduling_map_AoA.tau_AoA)
% tau_AoA
[fitresult, gof] = createFit1(M_eval_sche, H_eval_sche, Scheduling_map_AoA.tau_AoA);
tau_AoA_extrap = NaN*ones(size(M_sq_sche,2),size(H_sq_sche,2));
for I=1:size(M_sq_sche,2)
    for J=1:size(H_sq_sche,2)
%         min( M_eval_sche(I,:) )
        if isnan(M_eval_sche(I,J))==0 %&& M_eval_sche(I,J) <= 2*min( M_eval_sche(I,:) )
%             if M_eval_sche(I,J)<=1 
                tau_AoA_extrap(I,J) = fitresult( M_eval_sche(I,J) , H_eval_sche(I,J));
%             elseif M_eval_sche(I,J)>1 && J<=9
%                 tau_AoA_extrap(I,J) = fitresult( M_eval_sche(11,J) , H_eval_sche(11,J));
%             elseif J>9
%                 tau_AoA_extrap(I,J) = tau_AoA_extrap(I,9);
%             end
        end
    end
end
figure
surf(M_eval_sche,H_eval_sche,tau_AoA_extrap);
view(0,90)
% axis([0.2 2 0 18000 0 2])
pause
close all
% tau_q
[fitresult, gof] = createFit1(M_eval_sche, H_eval_sche, Scheduling_map_AoA.tau_q);
tau_q_extrap = NaN*ones(size(M_sq_sche,2),size(H_sq_sche,2));
for I=1:size(M_sq_sche,2)
    for J=1:size(H_sq_sche,2)
        if isnan(M_eval_sche(I,J))==0 %&& M_eval_sche(I,J) <= 2*min( M_eval_sche(I,:) )
%             if M_eval_sche(I,J)<=1 
                tau_q_extrap(I,J) = fitresult( M_eval_sche(I,J) , H_eval_sche(I,J));
%             elseif M_eval_sche(I,J)>1 && J<=9
%                 tau_q_extrap(I,J) = fitresult( M_eval_sche(11,J) , H_eval_sche(11,J));
%             elseif J>9
%                 tau_q_extrap(I,J) = tau_q_extrap(I,9);
%             end
        end
    end
end
figure
surf(M_eval_sche,H_eval_sche,tau_q_extrap);
view(0,90)
% axis([0.2 2 0 18000 0 1])

pause
close all


% K_b_AoA
[fitresult, gof] = createFit1(M_eval_sche, H_eval_sche, Scheduling_map_AoA.K_b_AoA);
K_b_AoA_extrap = NaN*ones(size(M_sq_sche,2),size(H_sq_sche,2));
for I=1:size(M_sq_sche,2)
    for J=1:size(H_sq_sche,2)
        if isnan(M_eval_sche(I,J))==0 %&& M_eval_sche(I,J) <= 2*min( M_eval_sche(I,:) )
%             if M_eval_sche(I,J)<=1 
                K_b_AoA_extrap(I,J) = fitresult( M_eval_sche(I,J) , H_eval_sche(I,J));
%             elseif M_eval_sche(I,J)>1 && J<=9
%                 K_b_AoA_extrap(I,J) = fitresult( M_eval_sche(11,J) , H_eval_sche(11,J));
%             elseif J>9
%                 K_b_AoA_extrap(I,J) = K_b_AoA_extrap(I,9);
%             end
        end
    end
end
figure
surf(M_eval_sche,H_eval_sche,K_b_AoA_extrap);
view(0,90)
axis([0.2 2 0 18000 0 10])

pause
close all
% K_b_q
[fitresult, gof] = createFit1(M_eval_sche, H_eval_sche, Scheduling_map_AoA.K_b_q);
K_b_q_extrap = NaN*ones(size(M_sq_sche,2),size(H_sq_sche,2));
for I=1:size(M_sq_sche,2)
    for J=1:size(H_sq_sche,2)
        if isnan(M_eval_sche(I,J))==0 %&& M_eval_sche(I,J) <= 2*min( M_eval_sche(I,:) )
%             if M_eval_sche(I,J)<=1 
                K_b_q_extrap(I,J) = fitresult( M_eval_sche(I,J) , H_eval_sche(I,J));
%             elseif M_eval_sche(I,J)>1 && J<=9
%                 K_b_q_extrap(I,J) = fitresult( M_eval_sche(11,J) , H_eval_sche(11,J));
%             elseif J>9
%                 K_b_q_extrap(I,J) = K_b_q_extrap(I,9);
%             end
        end
    end
end
figure
surf(M_eval_sche,H_eval_sche,K_b_q_extrap);
view(0,90)
axis([0.2 2 0 18000 0 10])

%%
% We can consider the gains to be constant and the tau's to be larger when
% aproaching stall
nanmean( nanmean(Scheduling_map_AoA.K_b_AoA) )
nanmean( nanmean(Scheduling_map_AoA.K_b_q) )
% We have to smooth the surfaces
surf(M_eval_sche,H_eval_sche,tau_AoA_extrap);
surf(M_eval_sche,H_eval_sche,tau_q_extrap);
pause
close all
    % AoA
    [fitresult, gof] = Poly_fit(M_eval_sche, H_eval_sche, tau_AoA_extrap);
    tau_AoA_extrap_smooth = NaN*ones(size(M_sq_sche,2),size(H_sq_sche,2));
    for I=1:size(M_sq_sche,2)
        for J=1:size(H_sq_sche,2)
            tau_AoA_extrap_smooth(I,J)= fitresult( M_sq_sche(I) , H_sq_sche(J));
            M_eval_sche_unc(I,J) = M_sq_sche(I);
            H_eval_sche_unc(I,J) = H_sq_sche(J);
        end
    end
    surf(M_eval_sche_unc,H_eval_sche_unc,tau_AoA_extrap_smooth);
    pause
    close all
    % q
    [fitresult, gof] = Poly_fit(M_eval_sche, H_eval_sche, tau_q_extrap);
    tau_q_extrap_smooth = NaN*ones(size(M_sq_sche,2),size(H_sq_sche,2));
    for I=1:size(M_sq_sche,2)
        for J=1:size(H_sq_sche,2)
            tau_q_extrap_smooth(I,J)= fitresult( M_sq_sche(I) , H_sq_sche(J));
            M_eval_sche_unc(I,J) = M_sq_sche(I);
            H_eval_sche_unc(I,J) = H_sq_sche(J);
        end
    end
    surf(M_eval_sche_unc,H_eval_sche_unc,tau_q_extrap_smooth);

    