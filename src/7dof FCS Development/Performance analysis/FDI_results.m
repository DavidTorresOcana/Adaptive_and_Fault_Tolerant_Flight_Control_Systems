%% FDI performance check
close all
load('FDI_results.mat')
FDI_block_fake
FDI_block_real
FDI_control

FDI_effect_fake
FDI_effect_real
Control_defl
%% Effectiveness
Color=lines(4);
plot(FDI_effect_real.Time,FDI_effect_real.Data(:,3),'Linewidth',1.5,'Color',Color(1,:))
hold on
plot(FDI_effect_real.Time,FDI_effect_real.Data(:,4),'Linewidth',1.5,'Color',Color(2,:))

plot(FDI_effect_fake.Time,FDI_effect_fake.Data(:,3),'Linewidth',1,'Color',Color(3,:))
plot(FDI_effect_fake.Time,FDI_effect_fake.Data(:,4),'Linewidth',1,'Color',Color(4,:))
title(' \fontsize{14} Effectiveness FDI. Loss 50% Left aileron and 100% rigth aileron')
xlabel('\fontsize{12} Time (s)');ylabel('\fontsize{12}Effectivenes %')
legend('\fontsize{12} Actual eff. Left ail','\fontsize{12} Actual eff. Right ail','\fontsize{12} FDI eff. Left ail','\fontsize{12} FDI eff. right ail')
pause
%% Blocakges
figure
FDI_control_defl=reshape(FDI_control.Data(4,:,:),1,size(FDI_control.Data,3));
for i=1:size(FDI_block_fake.Data(:,4),1)
    if isnan(FDI_block_fake.Data(i,4))
        FDI_block_fake.Data(i,4)=FDI_control_defl(i);
    end
end

Color=lines(4);
plot(FDI_control.Time,FDI_control_defl,'Linewidth',1.45,'Color',Color(3,:))
hold on
plot(Control_defl.Time,Control_defl.Data(:,8),'Linewidth',1.5,'Color',Color(1,:))
plot(FDI_block_fake.Time,FDI_block_fake.Data(:,4),'Linewidth',1.5,'Color',Color(2,:))
axis([0 20 -20 20])
title(' \fontsize{14} FDI blockages estimation. Right aileron blockage case')
xlabel('\fontsize{12} Time (s)');ylabel('\fontsize{12} Deflection (deg)')
legend('\fontsize{12} Control command','\fontsize{12} Actual ','\fontsize{12} FDI estimated')


