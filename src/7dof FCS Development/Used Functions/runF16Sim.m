%================================================
%     Matlab Script File used to run the 
%     non-linear F-16 Simulation.  The results
%     will also be saved to a file and plotted.
%
% Author: Richard S. Russell
% Revised: David Torres
%================================================

% clear;
% clc;
% close all
global altitude fi_type velocity fi_flag_Simulink;
global surface1 surface2 surface3;
global ElevatorDis AileronDis RudderDis;

surface1 = 'ele_';
surface2 = 'ail_';
surface3 = 'rud_';


newline = sprintf('\n');

disp('This is an F-16 Simulation.');
disp('The simulation will begin by asking you for the flight ');
disp('conditions for which the simulation will be performed.');
disp(newline);
disp('Accpetable values for flight condition parameters are:');
disp(newline);
disp('                                  Model');
disp('  Variable                LOFI            HIFI');
disp('              Units   Min     Max     Min     Max');
disp('  Altitude:   m      1550     12200   1550    12200');
disp('  AOA         deg    -10      45     -10      90');
disp('  Thrust      [0-1]  0        1       0        1');
disp('  Elevator    deg    -25.0    25.0   -25.0    25.0');
disp('  Aileron     deg    -21.5    21.5   -21.5    21.5');
disp('  Rudder      deg    -30      30     -30      30');
disp('  Velocity    m/s    90       270     90     270');
disp(newline);
disp(newline);
disp('The flight condition you choose will be used to trim the F16.');
disp('Note:  The trim routine will trim to the desired');
disp('altitude and velocity.  All other parameters');
disp('will be varied until level flight is achieved.  ');
disp('You may need to view the results of the simulation');
disp(' and retrim accordingly.');
disp(newline);
disp(newline);

%% Ask user which simulation to run.
%%
disp('Which model would you like to use to trim the aircraft:')
disp('  1. Low Fidelity F-16 Trim')
disp('  2. High Fidelity F-16 Trim')
% % fi_flag = input('Your Selection:  ');
disp(newline);
disp(newline);
fi_flag =2;
disp('  High Fidelity Model CHOOSEN')

% pause
%% Determine from flag the correct simulation.
%%
if fi_flag == 1;
  fi_type = 'lofi';
  fi_flag_Simulink = 0;
elseif fi_flag == 2;
  fi_type = 'hifi';
  fi_flag_Simulink = 1;
else
  disp('Invalid selection');
  break;
end

%% Trim aircraft to desired altitude and velocity
%%
altitude = input('Enter the altitude for the simulation (m)  :  ');
velocity = input('Enter the velocity for the simulation (m/s):  ');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize some varibles used to create disturbances. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DisEle_1 = 0;    DisEle_2 = 0;    DisEle_3 = 0;
DisAil_1 = 0;    DisAil_2 = 0;    DisAil_3 = 0;
DisRud_1 = 0;    DisRud_2 = 0;    DisRud_3 = 0;
ElevatorDis = 0; AileronDis = 0;  RudderDis = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find out which surface to create a disturbance on.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% % dis_flag = input('Would you like to create a disturbance on a surface (y/n):  ', 's');
dis_flag='n';
if dis_flag == 'y' 
  ElevatorDis = input('Enter the elevator distrubance deflection           (deg) :  ');
  DisEle_1 = ElevatorDis;    DisEle_2 = -2*ElevatorDis;    DisEle_3 = ElevatorDis; surfacedef = 'elevator';
  
  AileronDis = input('Enter the aileron distrubance deflection            (deg) :  ');
  DisAil_1 = AileronDis;    DisAil_2 = -2*AileronDis;    DisAil_3 = AileronDis; surfacedef = 'aileron';
  
  RudderDis = input('Enter the rudder distrubance deflection             (deg) :  ');
  DisRud_1 = RudderDis;    DisRud_2 = -2*RudderDis;    DisRud_3 = RudderDis; surfacedef = 'rudder';
  
elseif dis_flag == 'n'
  surfacedef = 'none';  %do nothing
else
  disp('Invalid Selection');
  break;
end
disp(newline);
disp(newline);


delta_T = 0.001;

TStart = 0; TFinal = 30;

%% SYMETRIC TRIMMING: Initial Conditions for trim routine.
%% FlightGear Visualization

% Cranfield Airport
% airport_ID : EGTC
% runway_ID : 1
latitude  = 52.073105 ;
longitude = -000.616697;

% San Fransisco Airport
% airport_ID : KSFO
% runway_ID : 10L
% latitude  = 37.76   
% longitude = -122.4
%% The following values seem to trim to most
%% flight condition.  If the F16 does not trim    Change these values.
throttle = 0.26;         % throttle, [0-1]
elevator = -0.09;       % elevator, degrees
beta = 0;               % sideslip, degrees
alpha = 8.49;           % AOA, degrees
rudder = -0.01;         % rudder angle, degrees
aileron = 0.01;         % aileron, degrees

[trim_state, trim_throttle, trim_control, dLEF, UX] = trim_F16(throttle, elevator,beta,  alpha, aileron, rudder, velocity, altitude);
trim_LEF_abs=dLEF;

%%

% % % % % % % % % % % % % % % % % % %% Check the goodness of the Symetric tim
% % % % % % % % % % % % % % % % % % sim( 'F16BlockSYM' ,[TStart TFinal]);
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % trim_file = sprintf('%s%.3f%s%.3f%s%.3f_%smodel_alt%0.f_vel%.0f.txt', surface1, ElevatorDis, surface2, AileronDis, surface3, RudderDis, fi_type, altitude, velocity);
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % fid_trim = fopen(trim_file, 'w');
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % heading1 = sprintf('%% \n\t\t  %s DATA Trim-Doublet on %s: Alt %.0f, Alpha %.0f\n\n', fi_type, surfacedef, altitude, alpha);
% % % % % % % % % % % % % % % % % % heading2 = sprintf('\ntime,npos,epos,alt,phi,theta,psi,vel,alpha,beta,p,q,r,pow,nx,ny,nz,mach,qbar,ps,\n\n');
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % fprintf(fid_trim,heading1);
% % % % % % % % % % % % % % % % % % fprintf(fid_trim,heading2);
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % fid_trim = fopen(trim_file, 'a');
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % for row = 1 : 1 : length(y_sim(:,1))
% % % % % % % % % % % % % % % % % %   fprintf(fid_trim,'%8.5f,',T(row,:));
% % % % % % % % % % % % % % % % % %   for column = 1 : 1 : length(y_sim(1,:))
% % % % % % % % % % % % % % % % % %     fprintf(fid_trim,'%8.5f,',y_sim(row,column));
% % % % % % % % % % % % % % % % % %   end
% % % % % % % % % % % % % % % % % %   for column = 1:1:length(surfaces(1,:))
% % % % % % % % % % % % % % % % % %     fprintf(fid_trim,'%8.5f,',surfaces(row,column));
% % % % % % % % % % % % % % % % % %   end
% % % % % % % % % % % % % % % % % %   fprintf(fid_trim,'\n');
% % % % % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % fclose(fid_trim);
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % plot_flag = input('Plot results (y/n):  ', 's');
% % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % if plot_flag == 'n'
% % % % % % % % % % % % % % % % % %   break;
% % % % % % % % % % % % % % % % % % else
% % % % % % % % % % % % % % % % % %   graphF16;
% % % % % % % % % % % % % % % % % % end

%% ASYMMETRIC PLANT
% This plant has aymmetric deflections of ailerons, elevator, LEFs. 
% As the LEFs are inteded to mantain high AoA flights and. hence, they are
% determined by the state: Here we suppose linear aerodynamics and take
% them also as control surfaces, adding the control signal to the
% determined by the AoA
fprintf('\n\n\n\n\n\n\n Now the Asymetric model will be compiled\n\n\n\n\n\n\n')

%  Un-comment if ANY change made to nlplantASYM.c
% mex nlplantASYM.c

Model_used_flag = 1;



