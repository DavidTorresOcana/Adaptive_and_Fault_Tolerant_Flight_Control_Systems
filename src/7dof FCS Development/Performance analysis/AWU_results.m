%% This script is intended to show how AWU works
load('AWU_How.mat')
AWU
NN_action
Comm_ang_accel
Control_defl
Responses

for i=1:size(Control_defl.Data,1)
    if Control_defl.Data(i,5)<=-21
        Control_defl.Data(i,5)=-21;
    end
        if Control_defl.Data(i,7)<=-21
        Control_defl.Data(i,7)=-21;
    end
end

subplot(4,1,1)
plot(Responses.time,Responses.signals(3).values);axis([4 25 -200 0])
title(' \fontsize{14} Roll rate (deg/s)');
legend('\fontsize{14} Roll rate comm','\fontsize{14} Actual Roll rate defl');grid
subplot(4,1,2)
plot(Control_defl.Time,Control_defl.Data(:,5:8));axis([4 25 -22 22])
title(' \fontsize{14} Ailerons (deg/s)');
legend('\fontsize{14} Left ail comm','\fontsize{14} Left ail defl','\fontsize{14} Right ail comm','\fontsize{14} Right ail defl')
grid
subplot(4,1,3)
plot(Comm_ang_accel.Time,Comm_ang_accel.Data(:,1))
hold on;
plot(NN_action.Time,NN_action.Data(:,1),'r');axis([4 25 -20 12])
title(' \fontsize{14} Commands rad/s^2');legend('\fontsize{14} Total comm \nu','\fontsize{14} NN comm \nu_{ad}');grid
subplot(4,1,4)
plot(AWU.Time,AWU.Data(:,1));axis([4 25 -0.2 3])
title(' \fontsize{14} AWU signal rad/s^2')
xlabel('\fontsize{12} Time (s)')
grid

%% Eval physical limits

ratio = vview(B,[[24,24,21,21,29,25,25];[-24,-24,-21,-21,-29,0,0]]',pinv(B)); % Ok, same result (13.7%)
xlabel(' \fontsize{12} Cl')
ylabel('\fontsize{12}Cm')
zlabel('\fontsize{12}Cn')
title(' \fontsize{16} Achievable space')


  disp('---------------------------------------------------------------')
  disp(' ');
  disp(' Control effectiveness matrix: B ='),disp(' '),disp(B)
  disp(' Position limits: [umin umax]'' ='),disp(' '),disp(plim')
  disp(' Blue (outer) set: { v : v = B*u, umin < u < umax }')
  disp('  Feasible virtual control set with constrained allocation')
  disp(' ')
  disp(' Red (inner) set: { v : umin < P*v < umax } where P = pinv(B)')
  disp('  Feasible virtual control set with linear allocation, u = P*v')
  disp(' ')
  disp(sprintf(' Red to blue size ratio: %0.3g%%',ratio*100))
  disp(' ')
  disp('---------------------------------------------------------------')
  
