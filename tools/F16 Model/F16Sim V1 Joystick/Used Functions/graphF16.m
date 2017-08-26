function graphF16()

%% Function to Graph the results of runsim
%%

global altitude velocity fi_type;
global surface1 surface2 surface3;
global ElevatorDis AileronDis RudderDis;

lofi_data_file = sprintf('%s%.3f%s%.3f%s%.3f_%smodel_alt%0.f_vel%.0f.txt', surface1, ElevatorDis, surface2, AileronDis, surface3, RudderDis, 'lofi', altitude, velocity);
hifi_data_file = sprintf('%s%.3f%s%.3f%s%.3f_%smodel_alt%0.f_vel%.0f.txt', surface1, ElevatorDis, surface2, AileronDis, surface3, RudderDis, 'hifi', altitude, velocity);

lofiID = fopen(lofi_data_file);
hifiID = fopen(hifi_data_file);

title_string = sprintf('Trimed at Velocity = %.1f \n Alt. = %.1f', velocity, altitude);
%        title_string = sprintf('');
if (lofiID > 0 & hifiID > 0)
    
    fclose(lofiID);
    fclose(hifiID);
    
    [time_lo, npos_lo, epos_lo, alt_lo, phi_lo, theta_lo, psi_lo, vel_lo, alpha_lo, sideslip_lo, roll_lo, pitch_lo, yaw_lo, nx_lo, ny_lo, nz_lo, mach_lo, qbar_lo, ps_lo, thrust_lo, ele_lo, ail_lo, rud_lo] = textread(lofi_data_file,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'delimiter', ',' , 'headerlines',3);
    [time_hi, npos_hi, epos_hi, alt_hi, phi_hi, theta_hi, psi_hi, vel_hi, alpha_hi, sideslip_hi, roll_hi, pitch_hi, yaw_hi, nx_hi, ny_hi, nz_hi, mach_hi, qbar_hi, ps_hi, thrust_hi, ele_hi, ail_hi, rud_hi] = textread(hifi_data_file,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'delimiter', ',' , 'headerlines',3);
    
    figure(1);
    title(title_string);
    subplot(231)
    plot(time_lo,npos_lo, time_hi, npos_hi, '--');
    ylabel('North Pos.');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(232)
    plot(time_lo, epos_lo, time_hi, epos_hi, '--');     
    ylabel('East Pos.');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(233)
    plot(time_lo, alt_lo,time_hi , alt_hi, '--') ;
    ylabel('Altitude ft');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(234)
    plot(time_lo, phi_lo, time_hi, phi_hi, '--');
    ylabel('PHI (degrees)');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(235)
    plot(time_lo, theta_lo, time_hi, theta_hi, '--');     
    ylabel('THETA (degrees)');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(236)
    plot(time_lo ,psi_lo, time_hi, psi_hi, '--');
    ylabel('PSI (degrees)');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    %% Figure 2
    %%
    
    figure(2)
    title(title_string);
    subplot(231)
    plot(time_lo,vel_lo, time_hi, vel_hi, '--');
    ylabel('Velocity (ft/s)');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(232)
    plot(time_lo, alpha_lo, time_hi, alpha_hi, '--');     
    ylabel('Angle of Attack (degrees)');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(233)
    plot(time_lo , sideslip_lo, time_hi, sideslip_hi, '--');
    ylabel('Side Slip (degrees)');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(234)
    plot(time_lo,roll_lo, time_hi, roll_hi, '--');
    ylabel('Roll Rate (deg/s)');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(235)
    plot(time_lo, pitch_lo , time_hi, pitch_hi, '--');
    ylabel('Pitch Rate (deg/s)');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    subplot(236)
    plot(time_lo, yaw_lo, time_hi, yaw_hi, '--');
    ylabel('Yaw Rate (deg/s)');
    xlabel('Time (sec)');
    title(title_string);
    legend('LOFI', 'HIFI');
    
    %% Figure 3
    %%
    figure(3);
    title(title_string);
    subplot(231)
    plot(time_lo, nx_lo, time_hi, nx_hi, '--');
    ylabel('acc x');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(232)
    plot(time_lo, ny_lo, time_hi, ny_hi, '--');     
    ylabel('acc y');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(233)
    plot(time_lo , nz_lo, time_hi, nz_hi, '--');
    ylabel('acc z');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(234)
    plot(time_lo, mach_lo, time_hi, mach_hi, '--');
    ylabel('Mach');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(235)
    plot(time_lo, qbar_lo, time_hi, qbar_hi, '--');
    ylabel('q bar)');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(236)
    plot(time_lo, ps_lo, time_hi, ps_hi, '--');
    ylabel('ps');
    xlabel('Time (sec)');
    title(title_string);
    
    %% Figure 4
    %%
    figure(4)
    title(title_string);
    subplot(221)
    plot(time_lo, thrust_lo, time_hi, thrust_hi, '--');
    ylabel('del Thrust');
    xlabel('Time (sec)');
    legend('LOFI', 'HIFI');
    title(title_string);
    
    subplot(222)
    plot(time_lo, ele_lo, time_hi, ele_hi, '--');     
    ylabel('del Elevator');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(223)
    plot(time_lo , ail_lo, time_hi, ail_hi, '--');
    ylabel('del Aileorn');
    xlabel('Time (sec)');
    title(title_string);
    
    
    subplot(224)
    plot(time_lo, rud_lo, time_hi, rud_hi,'--');
    ylabel('del Rudder');
    xlabel('Time (sec)');
    title(title_string);

    alt_min_hi = min(alt_hi) - 2000;
    alt_min(1:length(alt_lo),1) = alt_min_hi;
    figure(5)
    plot3(epos_lo,npos_lo,alt_lo, 'b-' , epos_hi,npos_hi,alt_hi,'g-.', epos_lo,npos_lo,alt_min, 'b:',  epos_hi,npos_hi,alt_min,'g:')
    grid;
    legend('LOFI', 'HIFI');
    xlabel('East Position')
    ylabel('North Position')
    zlabel('Altitude')
    
%     figure(6)
%     for j = 1:10:length(epos_lo)
%         plot3(epos_lo,npos_lo,alt_lo, 'b:',epos_hi,npos_hi,alt_hi, 'g:', epos_lo(j),npos_lo(j),alt_lo(j), 'ro', epos_hi(j),npos_hi(j),alt_hi(j), 'ro')
%         grid;
%         xlabel('East Position')
%         ylabel('North Position')
%         zlabel('Altitude')
%         legend('LOFI', 'HIFI')
%         F(j) = getframe;
%     end
     
else
    
    new_data_file = sprintf('%s%.3f%s%.3f%s%.3f_%smodel_alt%0.f_vel%.0f.txt', surface1, ElevatorDis, surface2, AileronDis, surface3, RudderDis, fi_type, altitude, velocity);
    [time_new, npos_new, epos_new, alt_new, phi_new, theta_new, psi_new, vel_new, alpha_new, sideslip_new, roll_new, pitch_new, yaw_new, nx_new, ny_new, nz_new, mach_new, qbar_new, ps_new, thrust_new, ele_new, ail_new, rud_new] = textread(new_data_file,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'delimiter', ',' , 'headerlines',3);
    
    % Plot angle of attack v. time
    figure(1);
    title(title_string);
    subplot(231)
    plot(time_new,npos_new);
    ylabel('North Pos.');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(232)
    plot(time_new, epos_new);     
    ylabel('East Pos.');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(233)
    plot(time_new , alt_new);
    ylabel('Altitude ft');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(234)
    plot(time_new, phi_new);
    ylabel('PHI (degrees)');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(235)
    plot(time_new, theta_new);     
    ylabel('THETA (degrees)');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(236)
    plot(time_new ,psi_new);
    ylabel('PSI (degrees)');
    xlabel('Time (sec)');
    title(title_string);
    
    %% Figure 2
    %%
    figure(2);
    title(title_string);
    subplot(231)
    plot(time_new,vel_new);
    ylabel('Velocity (ft/s)');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(232)
    plot(time_new, alpha_new);     
    ylabel('Angle of Attack (degrees)');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(233)
    plot(time_new , sideslip_new);
    ylabel('Side Slip (degrees)');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(234)
    plot(time_new,roll_new);
    ylabel('Roll Rate (deg/s)');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(235)
    plot(time_new,pitch_new);
    ylabel('Pitch Rate (deg/s)');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(236)
    plot(time_new, yaw_new);
    ylabel('Yaw Rate (deg/s)');
    xlabel('Time (sec)');
    title(title_string);
    
    %% Figure 3
    %%
    figure(3);
    title(title_string);
    subplot(231)
    plot(time_new, nx_new);
    ylabel('acc x');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(232)
    plot(time_new, ny_new);     
    ylabel('acc y');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(233)
    plot(time_new , nz_new);
    ylabel('acc z');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(234)
    plot(time_new, mach_new);
    ylabel('Mach');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(235)
    plot(time_new, qbar_new);
    ylabel('q bar)');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(236)
    plot(time_new, ps_new);
    ylabel('ps');
    xlabel('Time (sec)');
    title(title_string);
    
    %% Figure 4
    %%
    figure(4);
    title(title_string);
    subplot(221)
    plot(time_new, thrust_new);
    ylabel('Thrust');
    xlabel('Time (sec)');
    legend('LOFI', 'HIFI');
    title(title_string);
    
    subplot(222)
    plot(time_new, ele_new);     
    ylabel('del Elevator');
    xlabel('Time (sec)');
    title(title_string);
    
    subplot(223)
    plot(time_new , ail_new);
    ylabel('del Aileron');
    xlabel('Time (sec)');
    title(title_string);
    
    
    
    subplot(224)
    plot(time_new, rud_new);
    ylabel('del Rudder');
    xlabel('Time (sec)');
    title(title_string);
    
    
    alt_min_new = min(alt_new) - 2000;
    alt_min(1:length(alt_new),1) = alt_min_new;
    figure(5)
    plot3(epos_new, npos_new, alt_new, 'b',epos_new,npos_new,alt_min, 'b:')
    grid;
    legend('LOFI', 'HIFI');
    xlabel('East Position')
    ylabel('North Position')
    zlabel('Altitude')
%     figure(6)
%     for j = 1:10:length(epos_new)
%         plot3(epos_new, npos_new, alt_new, 'b:', epos_new(j),npos_new(j),alt_new(j), 'ro')
%         grid;
%         xlabel('East Position')
%         ylabel('North Position')
%         zlabel('Altitude')
%         legend('LOFI', 'HIFI')
%         F(j) = getframe;
%     end
end  % if end

%movie(F)