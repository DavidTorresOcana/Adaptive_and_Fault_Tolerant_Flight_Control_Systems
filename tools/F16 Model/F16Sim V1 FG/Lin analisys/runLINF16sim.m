%================================================
%     Matlab Script File used to linearize the 
%     non-linear F-16 model. The program will 
%     also run the linearized F-16 Simulation,
%     save the LTI state-space matrices to a file,
%     save and plot the simulation results.
%
% Author: Richard S. Russell
% 
%================================================

clear;
clc;

global fi_flag_Simulink;
global ElevatorDis AileronDis RudderDis;

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
disp('  Altitude:   ft      5000    40000   5000    40000');
disp('  AOA         deg    -10      45     -10      90');
disp('  Thrust      lbs     1000    19000   1000    19000');
disp('  Elevator    deg    -25.0    25.0   -25.0    25.0');
disp('  Aileron     deg    -21.5    21.5   -21.5    21.5');
disp('  Rudder      deg    -30      30     -30      30');
disp('  Velocity    ft/s    300     900     300     900');
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

%% Trim aircraft to desired altitude and velocity
%%
altitude = input('Enter the altitude for the simulation (ft)  :  ');
velocity = input('Enter the velocity for the simulation (ft/s):  ');

%% Initial guess for trim
%%
thrust = 5000;          % thrust, lbs
elevator = -0.09;       % elevator, degrees
alpha = 8.49;              % AOA, degrees
rudder = -0.01;             % rudder angle, degrees
aileron = 0.01;            % aileron, degrees

%% Find trim for Hifi model at desired altitude and velocity
%%
disp('Trimming High Fidelity Model:');
fi_flag_Simulink = 1;
[trim_state_hi, trim_thrust_hi, trim_control_hi, dLEF, xu_hi] = trim_F16(thrust, elevator, alpha, aileron, rudder, velocity, altitude);
trim_state_lin = trim_state_hi; trim_thrust_lin = trim_thrust_hi; trim_control_lin = trim_control_hi;

%% Find the state space model for the hifi model at the desired alt and vel.
%%
[A_hi,B_hi,C_hi,D_hi] = linmod('LIN_F16Block', [trim_state_lin; trim_thrust_lin; trim_control_lin(1); trim_control_lin(2); trim_control_lin(3); dLEF; -trim_state_lin(8)*180/pi], [trim_thrust_lin; trim_control_lin(1); trim_control_lin(2); trim_control_lin(3)]);

%% Find trim for lofi model at desired altitude and velocity
%%
disp('Trimming Low Fidelity Model:');
fi_flag_Simulink = 0;
[trim_state_lo, trim_thrust_lo, trim_control_lo, dLEF, xu_lo] = trim_F16(thrust, elevator, alpha, aileron, rudder, velocity, altitude);
trim_state_lin = trim_state_lo; trim_thrust_lin = trim_thrust_lo; trim_control_lin = trim_control_lo;

%% Find the state space model for the hifi model at the desired alt and vel.
%%
[A_lo,B_lo,C_lo,D_lo] = linmod('LIN_F16Block', [trim_state_lin; trim_thrust_lin; trim_control_lin(1); trim_control_lin(2); trim_control_lin(3); dLEF; -trim_state_lin(8)*180/pi], [trim_thrust_lin; trim_control_lin(1); trim_control_lin(2); trim_control_lin(3)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save State Space and eigenvalues to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trim_file = sprintf('StateSpace_alt%.0f_vel%.0f.txt', altitude, velocity);
fid_trim = fopen(trim_file, 'w');
%% For Hifi
%% Print A
%%
fprintf(fid_trim,'A_hi = \n');
for i = 1:1:length(A_hi(:,1))
    for j = 1:1:length(A_hi(1,:))
        fprintf(fid_trim, '%8.5f,', A_hi(i,j));
    end
    fprintf(fid_trim, '\n');
end
fprintf(fid_trim, '\n\n');

%% Print B
%%
fprintf(fid_trim,'B_hi = \n');
for i = 1:1:length(B_hi(:,1))
    for j = 1:1:length(B_hi(1,:))
        fprintf(fid_trim, '%8.5f,', B_hi(i,j));
    end
    fprintf(fid_trim, '\n');
end
fprintf(fid_trim, '\n\n');

%% Print C
%%
fprintf(fid_trim,'C_hi = \n');
for i = 1:1:length(C_hi(:,1))
    for j = 1:1:length(C_hi(1,:))
        fprintf(fid_trim, '%8.5f,', C_hi(i,j));
    end
    fprintf(fid_trim, '\n');
end
fprintf(fid_trim, '\n\n');


%% Print D
%%
fprintf(fid_trim,'D_hi = \n');
for i = 1:1:length(D_hi(:,1))
    for j = 1:1:length(D_hi(1,:))
        fprintf(fid_trim, '%8.5f,', D_hi(i,j));
    end
    fprintf(fid_trim, '\n');
end
fprintf(fid_trim, '\n\n');

%% For Lofi
%% Print A
%%
fprintf(fid_trim,'A_lo = \n');
for i = 1:1:length(A_lo(:,1))
    for j = 1:1:length(A_lo(1,:))
        fprintf(fid_trim, '%8.5f,', A_lo(i,j));
    end
    fprintf(fid_trim, '\n');
end
fprintf(fid_trim, '\n\n');

%% Print B
%%
fprintf(fid_trim,'B_lo = \n');
for i = 1:1:length(B_lo(:,1))
    for j = 1:1:length(B_lo(1,:))
        fprintf(fid_trim, '%8.5f,', B_lo(i,j));
    end
    fprintf(fid_trim, '\n');
end
fprintf(fid_trim, '\n\n');

%% Print C
%%
fprintf(fid_trim,'C_lo = \n');
for i = 1:1:length(C_lo(:,1))
    for j = 1:1:length(C_lo(1,:))
        fprintf(fid_trim, '%8.5f,', C_lo(i,j));
    end
    fprintf(fid_trim, '\n');
end
fprintf(fid_trim, '\n\n');


%% Print D
%%
fprintf(fid_trim,'D_lo = \n');
for i = 1:1:length(D_lo(:,1))
    for j = 1:1:length(D_lo(1,:))
        fprintf(fid_trim, '%8.5f,', D_lo(i,j));
    end
    fprintf(fid_trim, '\n');
end
fprintf(fid_trim, '\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize some varibles used to create disturbances. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DisEle_1 = 0;    DisEle_2 = 0;    DisEle_3 = 0;
DisAil_1 = 0;    DisAil_2 = 0;    DisAil_3 = 0;
DisRud_1 = 0;    DisRud_2 = 0;    DisRud_3 = 0;
ElevatorDis = 0; AileronDis = 0;  RudderDis = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find out which surface to creat a disturbance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dis_flag = input('Would you like to create a disturbance on a surface (y/n):  ', 's');

if dis_flag == 'y'
    ElevatorDis = input('Enter the elevator distrubance deflection           (deg) :  ');
    DisEle_1 = ElevatorDis;    DisEle_2 = -2*ElevatorDis;    DisEle_3 = ElevatorDis;
    
    AileronDis = input('Enter the aileron distrubance deflection            (deg) :  ');
    DisAil_1 = AileronDis;    DisAil_2 = -2*AileronDis;    DisAil_3 = AileronDis;
    
    RudderDis = input('Enter the rudder distrubance deflection             (deg) :  ');
    DisRud_1 = RudderDis;    DisRud_2 = -2*RudderDis;    DisRud_3 = RudderDis;
elseif dis_flag == 'n'
    %do nothing
else
    disp('Invalid Selection');
    break;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Conditions for model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thrust = 0;  % Since this a linear model

deltaT = 0.001;
TStart = 0;
TFinal = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run and save hifi then lofi linearized simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for fi_flag_Simulink = 0:1:1
    if fi_flag_Simulink == 0
        fi_model = 'hifi';
        A = A_hi; B = B_hi; C = C_hi; D = D_hi;
        trim_state = xu_hi;
        trim_thrust = trim_thrust_hi;
        trim_control = trim_control_hi;
    else
        fi_model = 'lofi';
        A = A_lo; B = B_lo; C = C_lo; D = D_lo;
        trim_state = xu_lo;
        trim_thrust = trim_thrust_lo;
        trim_control = trim_control_lo;
    end
    
     
    sim( 'SS_F16_Block' ,[TStart TFinal]);
    
    trim_file = sprintf('%s%.3f%s%.3f%s%.3f_%smodel_alt%0.f_vel%.0f_LTI.txt', 'ele_', ElevatorDis, 'ail_', AileronDis, 'rud_', RudderDis, fi_model, altitude, velocity);
    fid_trim = fopen(trim_file, 'w');
    
    heading = sprintf('\ntime,npos,epos,alt,phi,theta,psi,vel,alpha,beta,p,q,r,nx,ny,nz,mach,qbar,ps,thrust,ele,ail,rud\n\n');
    
    fprintf(fid_trim,heading);
    
    fid_trim = fopen(trim_file, 'a');
    
    for row = 1:1:length(simout(:,1))
        fprintf(fid_trim,'%8.5f,',T(row,:));
        for column = 1:1:length(simout(1,:))
            fprintf(fid_trim,'%8.5f,',simout(row,column));
        end
        for column = 1:1:length(controls(1,:))
            fprintf(fid_trim,'%8.5f,',controls(row,column));
        end
        fprintf(fid_trim,'\n');
    end
    
    fclose(fid_trim);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


plot_flag = input('Plot results (y/n):  ', 's');

if plot_flag == 'n'
    break;
else
    graphF16_all;
end
