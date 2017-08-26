%% SI performance check
clear  DET Entropy
close all

SI_Surfaces_est
SI_P_matrix
SI_input

for i=1:size(SI_Surfaces_est.Time,1)
    if SI_Surfaces_est.Time(i)<5
        RWS(i)=27.8/2;
    else
        RWS(i)= 0.29*27.8/2 + 0.71*27.8/2*(0.5);
    end
   
end
subplot(6,1,[1:2])
plot(SI_input.Time,SI_input.Data)
title(' {\fontsize{16}  Estimated surfaces vs. time. Loss of 50% right wing}')
ylabel(' \fontsize{14} Ang acc (rad/s^2)')
legend('\fontsize{12} p_{dot}','\fontsize{12} q_{dot}','\fontsize{12} r_{dot}')
axis([0 80 -2 2])

subplot(6,1,[3:6])
plot(SI_Surfaces_est.Time,SI_Surfaces_est.Data,'Linewidth',1)
hold on
h=plot(SI_Surfaces_est.Time,RWS,'--m')
set(h,'Linewidth',1.5)
h=plot([5,5],[15,5],'k')
set(h,'Linewidth',1.5)
xlabel('\fontsize{14}  time (s)');ylabel(' \fontsize{14} Area (m^2)')
legend('\fontsize{12} Left wing Surface','\fontsize{12} Right wing Surface','\fontsize{12} Fin surface','\fontsize{12} Actual Right Surface','\fontsize{12} Damage injection')
axis([0 80 5 15])

return
%%
figure
for i=1:size(SI_P_matrix.Time,1)
   DET(i)=det(SI_P_matrix.Data(:,:,i)); 
%    Entropy(i)=1/2*log10(2*pi*exp(1))^2*DET(i);
   Entropy(i)=1/2*(log(DET(i)) +3+3*log(2*pi));

end
plot(SI_P_matrix.Time,DET)
hold on
plot(SI_P_matrix.Time,Entropy,'r')
title(' {\fontsize{16}  Estimations quality}')
xlabel('\fontsize{14}  time (s)');ylabel(' \fontsize{14} Uncertainty')
legend('\fontsize{12} det(P)','\fontsize{12} Entropy')
axis([0 80 1 50])

