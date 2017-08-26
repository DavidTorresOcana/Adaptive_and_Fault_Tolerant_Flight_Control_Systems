%% Specs 3.3.4.1 section
close all
% load('FlightEnvelope.mat')

h=plot(F_envelope.M_1g ,F_envelope.H_1g/1000,'k' )
hold on
load('Partition_FE.mat')
% [x_part,y_part] = ginput(8);
% plot(x_part,y_part,'*r' )
set(h,'Linewidth',1.5);
%% CO Operational flight envelope
for i=1:size(F_envelope.M_1g,1)
    if i<=13
        CO_op_FE.M(i)=F_envelope.M_1g(i)*1.4;
        CO_op_FE.H(i)=F_envelope.H_1g(i);
    elseif i>=18
        CO_op_FE.M(i-5)=F_envelope.M_1g(i);
        CO_op_FE.H(i-5)=F_envelope.H_1g(i);
    end
end

h=plot(CO_op_FE.M ,CO_op_FE.H/1000 )
set(h,'Linewidth',1.5);

%% CR Operational flight envelope
% Assumin C_D_0=0.013
% K=0.1638
% V_BR=3^(0.25)*C_D_min
C_L_opt=sqrt(0.013/0.1638);
for i=1:size(F_envelope.M_1g,1)
    if i<=13
        [~,a,~,rho]=atmosisa(F_envelope.H_1g(i));
        V_BR=3^0.25*1/(C_L_opt^0.5)*sqrt(2*9295*9.81/(rho*27.89));
        CR_op_FE.M(i)=V_BR/a;
        CR_op_FE.H(i)=F_envelope.H_1g(i);
    elseif i>=18
        CR_op_FE.M(i-5)=F_envelope.M_1g(i);
        CR_op_FE.H(i-5)=F_envelope.H_1g(i);
    end
end

h=plot(CR_op_FE.M ,CR_op_FE.H/1000,'r' )
set(h,'Linewidth',1.5);


Color=lines(8);
hold on
%% CO ADD points @ 4000m
[~,a_4000,~,rho]=atmosisa(6000);
CO_V_VL = 98;
CO_V_L = 120;
CO_V_M = 180;
CO_V_H = 420;
h=plot([CO_V_VL]/a_4000,4,'.')
set(h,'Color',Color(1,:),'Markersize',20)
h=plot([CO_V_L]/a_4000,4,'.')
set(h,'Color',Color(2,:),'Markersize',20)
h=plot([CO_V_M]/a_4000,4,'.')
set(h,'Color',Color(3,:),'Markersize',20)
h=plot([CO_V_H]/a_4000,4,'.')
set(h,'Color',Color(4,:),'Markersize',20)
%% CR ADD points @ 6000m
[~,a_6000,~,rho]=atmosisa(6000);
CR_V_VL = 0.85*a; % @ M=0.7 
CR_V_L = 0.88*a+20;
CR_V_M = 1.4*0.88*a;
CR_V_H = 480;
h=plot([CR_V_VL]/a_6000,6,'.')
set(h,'Color',Color(5,:),'Markersize',20)
h=plot([CR_V_L]/a_6000,6,'.')
set(h,'Color',Color(6,:),'Markersize',20)
h=plot([CR_V_M]/a_6000,6,'.')
set(h,'Color',Color(7,:),'Markersize',20)
h=plot([CR_V_H]/a_6000,6,'.')
set(h,'Color',Color(8,:),'Markersize',20)

%% Partition
Part_1.M=[F_envelope.M_1g(1:8);x_part(3:-1:1)];
Part_1.H=[F_envelope.H_1g(1:8)/1000;y_part(3:-1:1)];
plot(Part_1.M,Part_1.H,'r')
h=fill(Part_1.M,Part_1.H,'r')
set(h,'facealpha',.5)
set(h,'edgealpha',.5)

Part_2.M=[x_part(3);F_envelope.M_1g(9:17);x_part([5,4,2])];
Part_2.H=[y_part(3);F_envelope.H_1g(9:17)/1000;y_part([5,4,2])];
plot(Part_2.M,Part_2.H,'g')
h=fill(Part_2.M,Part_2.H,'g')
set(h,'facealpha',.5)
set(h,'edgealpha',.5)

Part_3.M=[x_part(5);F_envelope.M_1g(18:23);x_part([7,6,4,5])];
Part_3.H=[y_part(5);F_envelope.H_1g(18:23)/1000;y_part([7,6,4,5])];
plot(Part_3.M,Part_3.H,'c')
h=fill(Part_3.M,Part_3.H,'c')
set(h,'facealpha',.5)
set(h,'edgealpha',.5)

Part_4.M=[x_part(7);F_envelope.M_1g(24:end);x_part([8,6])];
Part_4.H=[y_part(7);F_envelope.H_1g(24:end)/1000;y_part([8,6])];
plot(Part_4.M,Part_4.H,'m')
h=fill(Part_4.M,Part_4.H,'m')
set(h,'facealpha',.5)
set(h,'edgealpha',.5)

Part_5.M=[x_part([1,2,4,6,8])];
Part_5.H=[y_part([1,2,4,6,8])];
plot(Part_5.M,Part_5.H,'b')
h=fill(Part_5.M,Part_5.H,'b')
set(h,'facealpha',.5)
set(h,'edgealpha',.5)

plot(F_envelope.M_1g ,F_envelope.H_1g/1000 ,'k')


% return

axis([0.1 2.5 0 20])
title(' \fontsize{14} 1g Flight envelope partitions and Op. FE')
xlabel(' \fontsize{12} M' )
ylabel(' \fontsize{12} Altitude (km)' )
legend('\fontsize{12} 1g Envelope','\fontsize{12} CO Op. F. envelope','\fontsize{12} CR Op. F. envelope',...
    '\fontsize{12} CO VL','\fontsize{12} CO L','\fontsize{12} CO M','\fontsize{12} CO H',...
    '\fontsize{12} CR VL','\fontsize{12} CR L','\fontsize{12} CR M','\fontsize{12} CR H')
