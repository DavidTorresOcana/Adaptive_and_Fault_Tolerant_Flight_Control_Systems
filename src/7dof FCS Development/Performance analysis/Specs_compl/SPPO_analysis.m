%%SPPO analysis
close all
clear time AoA_comm AoA
% sub_0=220;
sub_0=452

time=SPPO4_ge.Time(sub_0:end);
AoA_comm=SPPO4_ge.Data(sub_0:end,3);
AoA=SPPO4_ge.Data(sub_0:end,4);
Pitch_rate=SPPO4_ge.Data(1:end,5);
ACC=SPPO4_ge.Data(sub_0:end,2);

plot(time,AoA)
hold on
plot(time,AoA_comm,'r')

%%
x_0=AoA(1) -AoA(end)  % 
t_0 = time(1) % 

% Find time that
x_1= x_0*0.7358 +AoA(end)  % 
x_2= x_0*0.4060 + AoA(end)  % 
x_3= x_0*0.19918 + AoA(end) % 

[ind,t_1p] = crossing(AoA,time,x_1)
[ind,t_2p] = crossing(AoA,time,x_2)
[ind,t_3p] = crossing(AoA,time,x_3)

t_1=t_1p(1)-t_0
t_2 = t_2p(1)-t_0
t_3 = t_3p(1)-t_0
plot(t_1p,x_1,'r*');plot(t_2p,x_2,'r*');plot(t_3p,x_3,'r*')

t_3/t_1
t_2/t_1
(t_3-t_2)/(t_2-t_1)

% Go to graphs with these values
relative_damping=mean([0.65,0.85,0.25])% aprox

omega_n= mean([0.85/t_1,1.48/t_2 ,1.8/t_3])
return
%% load factor
plot(SPPO4_ge.Time,SPPO4_ge.Data(:,4))
hold on
plot(SPPO4_ge.Time,SPPO4_ge.Data(:,2),'r')

(-2.43+3.25)/9.81/deg2rad(x_0) % Deltan/delta alpha   g's/rad





