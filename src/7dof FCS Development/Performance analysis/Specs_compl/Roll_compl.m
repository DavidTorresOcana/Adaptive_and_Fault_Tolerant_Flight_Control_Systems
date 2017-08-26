%% Roll rates compl
close all
ROLL

start=200;
finih= 450;
time= ROLL.Time(start:finih);
p_roll=ROLL.Data(start:finih,1);

plot(time,p_roll)
title('\fontsize{16} Roll rate response set-point 4')
xlabel('\fontsize{12} Time (s)');ylabel('\fontsize{12} Roll rate(deg/s)');

Response_ROLL=stepinfo( p_roll,time,p_roll(end))

Response_ROLL.RiseTime/log(9)


%% PLots

subplot(2,1,1)
plot(Responses_out.Time,Roll_response)
legend('\fontsize{10} Roll rate comm (º/s)','\fontsize{10} Roll rate actual (º/s)')
title('\fontsize{16} Rol rate and bank angle responses Set-point 3')
subplot(2,1,2)
xlabel('\fontsize{10} Time (s)')
plot(EULER_PLOT.time,EULER_PLOT.signals(1).values)
legend('\fontsize{10} Bank angle (º)')
