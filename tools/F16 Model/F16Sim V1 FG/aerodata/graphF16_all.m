function graphF16_all()

%% Function to Graph the results of runsim
%%

global altitude velocity fi_type;
global ElevatorDis AileronDis RudderDis;

surface1 = 'ele_';
surface2 = 'ail_';
surface3 = 'rud_';

lofi_data_file = sprintf('%s%.3f%s%.3f%s%.3f_%smodel_alt%0.f_vel%.0f.txt', surface1, ElevatorDis, surface2, AileronDis, surface3, RudderDis, 'lofi', altitude, velocity);
hifi_data_file = sprintf('%s%.3f%s%.3f%s%.3f_%smodel_alt%0.f_vel%.0f.txt', surface1, ElevatorDis, surface2, AileronDis, surface3, RudderDis, 'hifi', altitude, velocity);
lofi_LTI_file = sprintf('%s%.3f%s%.3f%s%.3f_%smodel_alt%0.f_vel%.0f_LTI.txt', surface1, ElevatorDis, surface2, AileronDis, surface3, RudderDis, 'lofi', altitude, velocity);
hifi_LTI_file = sprintf('%s%.3f%s%.3f%s%.3f_%smodel_alt%0.f_vel%.0f_LTI.txt', surface1, ElevatorDis, surface2, AileronDis, surface3, RudderDis, 'hifi', altitude, velocity);


lofiID = fopen(lofi_data_file);
hifiID = fopen(hifi_data_file);
lofiLTIID = fopen(lofi_LTI_file);
hifiLTIID = fopen(hifi_LTI_file);

if (lofiID > 0)
    [time_lo, npos_lo, epos_lo, alt_lo, phi_lo, theta_lo, psi_lo, vel_lo, alpha_lo, sideslip_lo, roll_lo, pitch_lo, yaw_lo, nx_lo, ny_lo, nz_lo, mach_lo, qbar_lo, ps_lo, thrust_lo, ele_lo, ail_lo, rud_lo] = textread(lofi_data_file,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'delimiter', ',' , 'headerlines',3);
else
    time_lo = 0; npos_lo = 0; epos_lo = 0; alt_lo = 0; phi_lo= 0; theta_lo= 0; psi_lo= 0; vel_lo= 0; alpha_lo= 0; sideslip_lo= 0; roll_lo= 0; pitch_lo= 0; yaw_lo= 0; nx_lo= 0; ny_lo= 0; nz_lo= 0; mach_lo= 0; qbar_lo= 0; ps_lo= 0; thrust_lo= 0; ele_lo= 0; ail_lo= 0; rud_lo= 0; 
end

if (hifiID > 0)
    [time_hi, npos_hi, epos_hi, alt_hi, phi_hi, theta_hi, psi_hi, vel_hi, alpha_hi, sideslip_hi, roll_hi, pitch_hi, yaw_hi, nx_hi, ny_hi, nz_hi, mach_hi, qbar_hi, ps_hi, thrust_hi, ele_hi, ail_hi, rud_hi] = textread(hifi_data_file,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'delimiter', ',' , 'headerlines',3);
else
    time_hi = 0; npos_hi = 0; epos_hi = 0; alt_hi = 0; phi_hi= 0; theta_hi= 0; psi_hi= 0; vel_hi= 0; alpha_hi= 0; sideslip_hi= 0; roll_hi= 0; pitch_hi= 0; yaw_hi= 0; nx_hi= 0; ny_hi= 0; nz_hi= 0; mach_hi= 0; qbar_hi= 0; ps_hi= 0; thrust_hi= 0; ele_hi= 0; ail_hi= 0; rud_hi= 0; 
end

if (lofiLTIID > 0)
    [time_lo_LTI, npos_lo_LTI, epos_lo_LTI, alt_lo_LTI, phi_lo_LTI, theta_lo_LTI, psi_lo_LTI, vel_lo_LTI, alpha_lo_LTI, sideslip_lo_LTI, roll_lo_LTI, pitch_lo_LTI, yaw_lo_LTI, nx_lo_LTI, ny_lo_LTI, nz_lo_LTI, mach_lo_LTI, qbar_lo_LTI, ps_lo_LTI, thrust_lo_LTI, ele_lo_LTI, ail_lo_LTI, rud_lo_LTI] = textread(lofi_LTI_file,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'delimiter', ',' , 'headerlines',3);
else
    time_lo_LTI = 0; npos_lo_LTI = 0; epos_lo_LTI = 0; alt_lo_LTI = 0; phi_lo_LTI= 0; theta_lo_LTI= 0; psi_lo_LTI= 0; vel_lo_LTI= 0; alpha_lo_LTI= 0; sideslip_lo_LTI= 0; roll_lo_LTI= 0; pitch_lo_LTI= 0; yaw_lo_LTI= 0; nx_lo_LTI= 0; ny_lo_LTI= 0; nz_lo_LTI= 0; mach_lo_LTI= 0; qbar_lo_LTI= 0; ps_lo_LTI= 0; thrust_lo_LTI= 0; ele_lo_LTI= 0; ail_lo_LTI= 0; rud_lo_LTI= 0; 
end

if (hifiLTIID > 0)
    [time_hi_LTI, npos_hi_LTI, epos_hi_LTI, alt_hi_LTI, phi_hi_LTI, theta_hi_LTI, psi_hi_LTI, vel_hi_LTI, alpha_hi_LTI, sideslip_hi_LTI, roll_hi_LTI, pitch_hi_LTI, yaw_hi_LTI, nx_hi_LTI, ny_hi_LTI, nz_hi_LTI, mach_hi_LTI, qbar_hi_LTI, ps_hi_LTI, thrust_hi_LTI, ele_hi_LTI, ail_hi_LTI, rud_hi_LTI] = textread(hifi_LTI_file,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'delimiter', ',' , 'headerlines',3);
else
    time_hi_LTI = 0; npos_hi_LTI = 0; epos_hi_LTI = 0; alt_hi_LTI = 0; phi_hi_LTI= 0; theta_hi_LTI= 0; psi_hi_LTI= 0; vel_hi_LTI= 0; alpha_hi_LTI= 0; sideslip_hi_LTI= 0; roll_hi_LTI= 0; pitch_hi_LTI= 0; yaw_hi_LTI= 0; nx_hi_LTI= 0; ny_hi_LTI= 0; nz_hi_LTI= 0; mach_hi_LTI= 0; qbar_hi_LTI= 0; ps_hi_LTI= 0; thrust_hi_LTI= 0; ele_hi_LTI= 0; ail_hi_LTI= 0; rud_hi_LTI= 0; 
end

title_string = sprintf('Trim: Vel. = %.1f \n Alt. = %.1f', velocity, altitude);

figure(1);

subplot(231)
plot(time_lo, npos_lo,'-', time_hi, npos_hi, '-.',  time_lo_LTI, npos_lo_LTI, ':', time_hi_LTI, npos_hi_LTI, '--');
ylabel('North Pos.');
xlabel('Time (sec)');
legend('LOFI', 'HIFI', 'LOFI-LTI', 'HIFI-LTI');
title(title_string);

subplot(232)
plot(time_lo, epos_lo,'-', time_hi, epos_hi, '-.', time_lo_LTI, epos_lo_LTI, ':',time_hi_LTI, epos_hi_LTI,  '--');     
ylabel('East Pos.');
xlabel('Time (sec)');
title(title_string);

subplot(233)
plot(time_lo, alt_lo, '-', time_hi , alt_hi, '-.', time_lo_LTI, alt_lo_LTI ,':',time_hi_LTI , alt_hi_LTI, '--') ;
ylabel('Altitude ft');
xlabel('Time (sec)');
title(title_string);

subplot(234)
plot(time_lo, phi_lo,'-', time_hi, phi_hi, '-.', time_lo_LTI, phi_lo_LTI, ':',time_hi_LTI, phi_hi_LTI, '--');
ylabel('PHI (degrees)');
xlabel('Time (sec)');
title(title_string);

subplot(235)
plot(time_lo, theta_lo, '-',time_hi, theta_hi,'-.', time_lo_LTI, theta_lo_LTI,':', time_hi_LTI, theta_hi_LTI,'--');     
ylabel('THETA (degrees)');
xlabel('Time (sec)');
title(title_string);

subplot(236)
plot(time_lo ,psi_lo,'-', time_hi, psi_hi, '-.', time_lo_LTI ,psi_lo_LTI, ':',time_hi_LTI, psi_hi_LTI,'--');
ylabel('PSI (degrees)');
xlabel('Time (sec)');
title(title_string);

%% Figure 2
%%
figure(2);
title(title_string);
subplot(231)
plot(time_lo,vel_lo,'-', time_hi, vel_hi, '-.', time_lo_LTI,vel_lo_LTI, ':',time_hi_LTI, vel_hi_LTI,'--');
ylabel('Velocity (ft/s)');
xlabel('Time (sec)');
legend('LOFI', 'HIFI', 'LOFI-LTI', 'HIFI-LTI');
title(title_string);

subplot(232)
plot(time_lo, alpha_lo,'-', time_hi, alpha_hi, '-.', time_lo_LTI, alpha_lo_LTI, ':',time_hi_LTI, alpha_hi_LTI,'--');     
ylabel('Angle of Attack (degrees)');
xlabel('Time (sec)');
title(title_string);


subplot(233)
plot(time_lo , sideslip_lo,'-', time_hi, sideslip_hi, '-.',  time_lo_LTI , sideslip_lo_LTI, ':',time_hi_LTI, sideslip_hi_LTI, '--');
ylabel('Side Slip (degrees)');
xlabel('Time (sec)');
title(title_string);


subplot(234)
plot(time_lo,roll_lo,'-', time_hi, roll_hi, '-.', time_lo_LTI,roll_lo_LTI,   ':', time_hi_LTI, roll_hi_LTI, '--');
ylabel('Roll Rate (deg/s)');
xlabel('Time (sec)');
title(title_string);


subplot(235)
plot(time_lo, pitch_lo, '-', time_hi, pitch_hi, '-.', time_lo_LTI, pitch_lo_LTI ,  ':',time_hi_LTI, pitch_hi_LTI,'--');
ylabel('Pitch Rate (deg/s)');
xlabel('Time (sec)');
title(title_string);


subplot(236)
plot(time_lo, yaw_lo,'-', time_hi, yaw_hi, '-.', time_lo_LTI, yaw_lo_LTI, ':', time_hi_LTI, yaw_hi_LTI,'--');
ylabel('Yaw Rate (deg/s)');
xlabel('Time (sec)');
title(title_string);


%% Figure 3
%%
figure(3);
title(title_string);
subplot(231)
plot(time_lo, nx_lo, '-',time_hi, nx_hi, '-.', time_lo_LTI, nx_lo_LTI, ':', time_hi_LTI, nx_hi_LTI,'--');
ylabel('acc x');
xlabel('Time (sec)');
legend('LOFI', 'HIFI', 'LOFI-LTI', 'HIFI-LTI');
title(title_string);

subplot(232)
plot(time_lo, ny_lo,'-', time_hi, ny_hi, '-.', time_lo_LTI, ny_lo_LTI,  ':',time_hi_LTI, ny_hi_LTI,'--');     
ylabel('acc y');
xlabel('Time (sec)');
title(title_string);


subplot(233)
plot(time_lo , nz_lo,'-', time_hi, nz_hi, '-.', time_lo_LTI , nz_lo_LTI, ':', time_hi_LTI, nz_hi_LTI, '--');
ylabel('acc z');
xlabel('Time (sec)');
title(title_string);


subplot(234)
plot(time_lo, mach_lo,'-', time_hi, mach_hi,'-.', time_lo_LTI, mach_lo_LTI, ':', time_hi_LTI, mach_hi_LTI,'--');
ylabel('Mach');
xlabel('Time (sec)');
title(title_string);


subplot(235)
plot(time_lo, qbar_lo,'-', time_hi, qbar_hi, '-.', time_lo_LTI, qbar_lo_LTI,  ':',time_hi_LTI, qbar_hi_LTI, '--');
ylabel('q bar)');
xlabel('Time (sec)');
title(title_string);


subplot(236)
plot(time_lo, ps_lo, '-', time_hi, ps_hi, '-.', time_lo_LTI, ps_lo_LTI, ':', time_hi_LTI, ps_hi_LTI, '--');
ylabel('ps');
xlabel('Time (sec)');
title(title_string);

% Figure 4
%
figure(4);
title(title_string);
subplot(221)
plot(time_lo, thrust_lo, '-',time_hi, thrust_hi, '-.', time_lo_LTI, thrust_lo_LTI, ':', time_hi_LTI, thrust_hi_LTI,'--');
ylabel('Thrust');
xlabel('Time (sec)');
legend('LOFI', 'HIFI', 'LOFI-LTI', 'HIFI-LTI');
title(title_string);

subplot(222)
plot(time_lo, ele_lo,'-', time_hi, ele_hi, '-.', time_lo_LTI, ele_lo_LTI,  ':',time_hi_LTI, ele_hi_LTI,'--');     
ylabel('Del Elevator');
xlabel('Time (sec)');
title(title_string);


subplot(223)
plot(time_lo , ail_lo,'-', time_hi, ail_hi, '-.', time_lo_LTI , ail_lo_LTI, ':', time_hi_LTI, ail_hi_LTI, '--');
ylabel('Del Aileron');
xlabel('Time (sec)');
title(title_string);


subplot(224)
plot(time_lo, rud_lo,'-', time_hi, rud_hi,'-.', time_lo_LTI, rud_lo_LTI, ':', time_hi_LTI, rud_hi_LTI,'--');
ylabel('Del Rudder');
xlabel('Time (sec)');
title(title_string);

% alt_min_hi = min(alt_hi) - 2000;
% alt_min(1:length(alt_lo),1) = alt_min_hi;
% figure(5)
% plot3(epos_lo,npos_lo,alt_lo, 'b-' , epos_hi,npos_hi,alt_hi,'g-.', epos_lo,npos_lo,alt_min, 'b:',  epos_hi,npos_hi,alt_min,'g:')
% grid;
% legend('LOFI', 'HIFI');
% xlabel('East Position')
% ylabel('North Position')
% zlabel('Altitude')











