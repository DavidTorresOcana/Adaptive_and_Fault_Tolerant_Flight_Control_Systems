function [ Cost ] = Cost_function( Long, params )
global  PERFORMANCE_ACC Response_ACC tau_ACC tau_q K_b_ACC K_b_q

tau_ACC = params(1);
tau_q = params(2);
K_b_ACC = params(3);
K_b_q= params(4);

try
    sim('F16ASYM_Controlled_ACC',[0 max(Long.time)])
    
    for i=2:size(PERFORMANCE_ACC.Data,1)
        if PERFORMANCE_ACC.Data(i,1)<PERFORMANCE_ACC.Data(i-1,1) % Falling flank
            flank_r =i;
            break
        end
    end
    
    if exist('flank_r') && max(PERFORMANCE_ACC.Time(:,1))>19.5
        
        
        % %     plot(PERFORMANCE_ACC.Time(flank_r:end,1)-PERFORMANCE_ACC.Time(flank_r,1),-PERFORMANCE_ACC.Data(flank_r:end,2)+PERFORMANCE_ACC.Data(flank_r-1,2))
        
        Response_ACC=stepinfo(-PERFORMANCE_ACC.Data(flank_r:end,2)+PERFORMANCE_ACC.Data(flank_r-1,2),PERFORMANCE_ACC.Time(flank_r:end,1)-PERFORMANCE_ACC.Time(flank_r,1),-PERFORMANCE_ACC.Data(end,2)+PERFORMANCE_ACC.Data(flank_r-1,2));
        % if Response_ACC.SettlingTime>3 || isnan(Response_ACC.SettlingTime)
        %     Response_ACC=stepinfo(PERFORMANCE_ACC.Data(flank_r:end,2)-PERFORMANCE_ACC.Data(flank_r-1,2),PERFORMANCE_ACC.Time(flank_r:end,1)-PERFORMANCE_ACC.Time(flank_r,1),-PERFORMANCE_ACC.Data(end,1)+PERFORMANCE_ACC.Data(flank_r-1,2));
        % end
        
        %% Cost
        % Cost = rms(  bsxfun(@minus,  -PERFORMANCE_ACC.Data(flank_r:end,2)+PERFORMANCE_ACC.Data(flank_r-1,2), -PERFORMANCE_ACC.Data(flank_r:end,1)+PERFORMANCE_ACC.Data(flank_r-1,1))   );
        
        Cost =  0.2*Response_ACC.RiseTime +  1*Response_ACC.SettlingTime +  1*Response_ACC.Overshoot+  0.1*Response_ACC.Undershoot;
        
    else
        Cost =  Inf;
        
    end
catch
    Cost =  Inf;
  
end

% var(-PERFORMANCE_ACC.Data(flank_r:end,2)+PERFORMANCE_ACC.Data(flank_r-1,2))




% % %% Plotting
% % figure('Name','responses')
% % plot(0,0)
% % fig=findall(findobj('Name','responses'),'Type','axes');
% % 
% % plot(fig(end),PERFORMANCE_ACC.Time,PERFORMANCE_ACC.Data(:,2),'r')
% % hold on
% % plot(fig(end),PERFORMANCE_ACC.Time,PERFORMANCE_ACC.Data(:,1))


end

