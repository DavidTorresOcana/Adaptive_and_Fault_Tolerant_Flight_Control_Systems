function [trim_state, trim_throttle, trim_control, dLEF, xu,cost] = trim_F16_fast(throttle, elevator,beta, alpha, ail, rud, vel, alt)
%================================================
%     F16 nonlinear model trimming routine
%  for longitudinal motion, steady level flight
%
% Author: T. Keviczky
%
%
%      Added addtional functionality.
%      This trim function can now trim at three 
%      additional flight conditions
%         -  Steady Turning Flight given turn rate
%         -  Steady Pull-up flight - given pull-up rate
%         -  Steady Roll - given roll rate
%
% Coauthor: Richard S. Russell
% Date:     18/04/2014
% Revised: David Torres

%
%================================================

global altitude velocity fi_flag_Simulink
global phi psi p q r phi_weight theta_weight psi_weight pow

altitude = alt;
velocity = vel;
alpha = alpha*pi/180;  %convert to radians
beta = beta*pi/180;  %convert to radians

% OUTPUTS: trimmed values for states and controls
% INPUTS:  guess values for throttle, elevator, alpha,beta  (assuming steady level flight)

% Initial Guess for free parameters
UX0 = [throttle; elevator; beta;alpha; ail; rud];  % free parameters: 3 control values & angle of attack

% Initialize some varibles
%
phi = 0; psi = 0;
p = 0; q = 0; r = 0;
phi_weight = 10; theta_weight = 10; psi_weight = 10;

FC_flag = 1;


% Initializing optimization options and running optimization:
OPTIONS = optimset('TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',5e+10,'MaxIter',1e+10);

iter = 1;
while iter == 1
   
    [UX,FVAL,EXITFLAG,OUTPUT] = fminsearch('trimfun',UX0,OPTIONS);
   
    [cost, Xdot, xu] = trimfun(UX);
    
    fprintf('\n');
    disp('Trim Values and Cost:');
    disp(['cost   = ' num2str(cost)])
    disp(['throttle = ' num2str(xu(14)*100) ' %'])
    disp(['elev   = ' num2str(xu(15)) ' deg'])
    disp(['ail    = ' num2str(xu(16)) ' deg'])
    disp(['rud    = ' num2str(xu(17)) ' deg'])
    disp(['alpha  = ' num2str(xu(8)*180/pi) ' deg'])
    disp(['beta   = ' num2str(xu(9)*180/pi) ' deg'])
    disp(['dLEF   = ' num2str(xu(18)) ' deg'])
    disp(['Vel.   = ' num2str(velocity) 'm/s']) 
% %     flag = input('Continue trim rountine iterations? (y/n):  ','s'); 
        flag ='n';
    if flag == 'n'
        iter = 0;
    end
    UX0 = UX;
end

% % UX0 = [throttle; elevator; beta;alpha; ail; rud];  
% For simulink:
trim_state=xu(1:13);
trim_throttle=UX(1);
trim_ele=UX(2);
trim_ail=UX(5);
trim_rud=UX(6);
trim_control=[UX(2);UX(5);UX(6)];
dLEF = xu(18);