function [ Cost ] = Cost_function( Long, params )
global  PERFORMANCE_AoA Response_AoA tau_AoA tau_q K_b_AoA K_b_q

tau_AoA = params(1);
tau_q = params(2);
K_b_AoA = params(3);
K_b_q= params(4);

try
    sim('F16ASYM_Controlled_AoA',[0 max(Long.time)])
    
    for i=2:size(PERFORMANCE_AoA.Data,1)
        if PERFORMANCE_AoA.Data(i,1)<PERFORMANCE_AoA.Data(i-1,1) % Falling flank
            flank_r =i;
            break
        end
    end
    
    if exist('flank_r') && max(PERFORMANCE_AoA.Time(:,1))>19.5
        
        
        % %     plot(PERFORMANCE_AoA.Time(flank_r:end,1)-PERFORMANCE_AoA.Time(flank_r,1),-PERFORMANCE_AoA.Data(flank_r:end,2)+PERFORMANCE_AoA.Data(flank_r-1,2))
        
        Response_AoA=stepinfo(-PERFORMANCE_AoA.Data(flank_r:end,2)+PERFORMANCE_AoA.Data(flank_r-1,2),PERFORMANCE_AoA.Time(flank_r:end,1)-PERFORMANCE_AoA.Time(flank_r,1),-PERFORMANCE_AoA.Data(end,2)+PERFORMANCE_AoA.Data(flank_r-1,2));
        % if Response_AoA.SettlingTime>3 || isnan(Response_AoA.SettlingTime)
        %     Response_AoA=stepinfo(PERFORMANCE_AoA.Data(flank_r:end,2)-PERFORMANCE_AoA.Data(flank_r-1,2),PERFORMANCE_AoA.Time(flank_r:end,1)-PERFORMANCE_AoA.Time(flank_r,1),-PERFORMANCE_AoA.Data(end,1)+PERFORMANCE_AoA.Data(flank_r-1,2));
        % end
        
        %% Cost
        % Cost = rms(  bsxfun(@minus,  -PERFORMANCE_AoA.Data(flank_r:end,2)+PERFORMANCE_AoA.Data(flank_r-1,2), -PERFORMANCE_AoA.Data(flank_r:end,1)+PERFORMANCE_AoA.Data(flank_r-1,1))   );
        
        Cost =  0*Response_AoA.RiseTime +  1*Response_AoA.SettlingTime +  0.4*Response_AoA.Overshoot;
        
    else
        Cost =  Inf;
        
    end
catch
    Cost =  Inf;
  
end

% var(-PERFORMANCE_AoA.Data(flank_r:end,2)+PERFORMANCE_AoA.Data(flank_r-1,2))




%% Plotting
% figure('Name','responses')
% plot(0,0)
% fig=findall(findobj('Name','responses'),'Type','axes');
% 
% plot(fig(end),PERFORMANCE_AoA.Time,PERFORMANCE_AoA.Data(:,2),'r')
% hold on
% plot(fig(end),PERFORMANCE_AoA.Time,PERFORMANCE_AoA.Data(:,1))


end

