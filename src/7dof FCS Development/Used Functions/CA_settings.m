%% Eval Moments Coeff Estimator (Non-linear Plant): Moms_Coeff_estimator 
% This estimator is based on the actual plant including the actual dynamics
% and coefficients
%           altitude (m)
%           V_TAS (m/s)
% 
% Inputs:   alpha (rad)
%           beta (rad)
%           P (rad/s)
%           Q (rad/s)
%           R (rad/s)
%           
%           Throttle [0-1]
% 
%           elevator_left (deg)
%           elevator_right (deg)
%           aileron_left (deg)
%           aileron_right (deg)
%           rudder (deg)
%           LEF_left (deg)  (additionally to the symetric defflection)
%           LEF_right (deg) (additionally to the symetric defflection)
%
%
% Outputs   C_L
%           C_M 
%           C_N
fprintf('\n \n\n Now aerodynamic coefficients model is compiled and evaluated\n')
% mex Moms_Coeff_estimator.c
% Test
Inputs_test=[4400,  160,deg2rad(2),0  ,0,0,0,  0.26,   -2,-2,  0,0,  0,  2,2]';
Inputs_test = [  trim_state(3),trim_state(7:12)'  ,  trim_throttle  ,trim_control(1),trim_control(1),trim_control(2),-trim_control(2),trim_control(3),0,0]';
tic
    Moms_Coeff_estimator(Inputs_test)
toc

% Trim aircraft (Asymetric model and Only Moments)
% This rutine trim the aircraft at any flight condition
fprintf('\n Now the actual Asymetric model is trim from the trim point obtained ')
fprintf('  with the Symetric Model (should be same result)\n \t Only Rotational dynamics is taken into account and symetric deflections are supposed\n');
AsymTrimCF = @(deltas) norm( Moms_Coeff_estimator( [ trim_state(3),trim_state(7:12)' , trim_throttle ,deltas(1),deltas(1),deltas(2),-deltas(2),deltas(3),0,0]' ) ) ;
A= [1,0,0; -1,0,0;0,1,0;0,-1,0;0,0,1;0,0,-1];
b=[24,24,21,21,29,29];

options = optimset('Display','off','Algorithm','active-set','TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',5e+04,'MaxIter',1e+04);
tic
[trim_control,cost] = fmincon(AsymTrimCF ,trim_control,A,b,[],[],[],[],[],options);
toc
    disp('Trim state once again')
    disp(['cost   = ' num2str(cost)])
    disp(['throttle = ' num2str(trim_throttle*100) ' %'])
    disp(['Left elev   = ' num2str(trim_control(1)) ' deg'])
    disp(['Right elev   = ' num2str(trim_control(1)) ' deg'])
    disp(['Left ail    = ' num2str(trim_control(2)) ' deg'])
    disp(['Right ail    = ' num2str(-trim_control(2)) ' deg'])
    disp(['rud    = ' num2str(trim_control(3)) ' deg'])
    disp(['alpha  = ' num2str(trim_state(8)*180/pi) ' deg'])
    disp(['beta   = ' num2str(trim_state(9)*180/pi) ' deg'])
    disp(['dLEF   = ' num2str(dLEF) ' deg'])
    disp(['Vel.   = ' num2str(trim_state(7)) 'm/s']) 
% pause
%% Linearization 
% In this section, the plant (just moments) is linearize at the triming points chosen
%           elevator_left (deg)
%           elevator_right (deg)
%           aileron_left (deg)
%           aileron_right (deg)
%           rudder (deg)
%           LEF_left (deg)  (additionally to the symetric defflection)
%           LEF_right (deg) (additionally to the symetric defflection)
fprintf('  \n Now Rotational dynamics is linearized  about the triming point\n ')
fprintf('  The aim is to find linear relationships between {C_L,C_M,C_N} and control deflections:\n ')
fprintf('   {C_L,C_M,C_N} = B*(deflections-deflections_trim) \n ')
fprintf('  \n The controls effectors taken into account are:\n ')
fprintf('   \t  elevator_left (deg)\n ')
fprintf('   \t  elevator_right (deg)\n ')
fprintf('   \t  aileron_left (deg)\n ')
fprintf('   \t  aileron_right (deg)\n ')
fprintf('   \t  rudder (deg)\n ')
fprintf('   \t  LEF_left (deg)  (additionally to the symetric defflection)\n ')
fprintf('   \t  LEF_right (deg) (additionally to the symetric defflection)\n ')
% pause

Inputs_test = [  trim_state(3),trim_state(7:12)'  ,  trim_throttle  ,trim_control(1),trim_control(1),trim_control(2),-trim_control(2),trim_control(3),0,0]';
Moms_Coeff_estimator(Inputs_test);
% Linearization
epsilon=2;
Moments = @(deltas) Moms_Coeff_estimator( [ trim_state(3),trim_state(7:12)' , trim_throttle  ,deltas(1),deltas(2),deltas(3),deltas(4),deltas(5),deltas(6),deltas(7)]' );
B = zeros(3,7);
tic
for i =1:7
    deltas = [trim_control(1),trim_control(1),trim_control(2),-trim_control(2),trim_control(3),0,0];
    deltas(i) = deltas(i) - epsilon;
    B_left = feval(Moments, deltas );
    
    deltas = [trim_control(1),trim_control(1),trim_control(2),-trim_control(2),trim_control(3),0,0];
    deltas(i) = deltas(i) + epsilon;
    B_right = feval(Moments, deltas );
    B(:,i) = (B_right- B_left )./(2*epsilon);
end
toc
fprintf('  \n The found effectiveness matrix is:\n ');

B
[k_B,m_B]=size(B);
% qcatdoc

% Checking for correctness of B
deltas = [trim_control(1),trim_control(1),trim_control(2),-trim_control(2),trim_control(3),0,0];
deltas(2) = deltas(2) - epsilon;
CM = feval(Moments, deltas )
B*(deltas-[trim_control(1),trim_control(1),trim_control(2),-trim_control(2),trim_control(3),0,0])'
% pause
%% Script for on-line evaluation of B
fprintf('  \n Now the effectivenes matrix is computed with the On-line function\n ')
% mex Moms_Coeff_Grad_estimator.c
Inputs_test = [  trim_state(3),trim_state(7:12)'  ,  trim_throttle  ,trim_control(1),trim_control(1),trim_control(2),-trim_control(2),trim_control(3),0,0]';
% B_OL=zeros(B,3*7,1);

tic 
    B_OL=Moms_Coeff_Grad_estimator(Inputs_test);
toc


B_OL=reshape(B_OL,3,7)

B;

B./B_OL

% Checking for correctness of B_OL
deltas = [trim_control(1),trim_control(1),trim_control(2),-trim_control(2),trim_control(3),0,0];
deltas(7) = deltas(7) - epsilon;
CM = feval(Moments, deltas )
B_OL*(deltas-[trim_control(1),trim_control(1),trim_control(2),-trim_control(2),trim_control(3),0,0])'
% pause
%% Performance of C online-linearized

% % % % for i=1:100
% % % %     Inputs_test = [  trim_state(3),trim_state(7:12)'  ,  trim_throttle  ,trim_control(1),rand(1)*trim_control(1),trim_control(2),-rand(1)*trim_control(2),trim_control(3),0,0]';
% % % %     tic
% % % %     B_OL=Moms_Coeff_Grad_estimator(Inputs_test);
% % % %     elapsedTime(i)=toc;
% % % % end
% % % % 
% % % % mean(elapsedTime)

%% Compile reconfigurated models

% mex Moms_Coeff_Reconf_estimator.c

% mex Moms_Coeff_Grad_Reconf_estimator.c
