%% This script iteratively evaluate the performance of the model
% The input signal are several square signals
Long.time=(0:T_sim:20)';
Long.signals.values( 1:size(Long.time,1) ,1) = zeros(size(Long.time,1),1);
Long.signals.values( floor(size(Long.time,1)*5/8)+2:end ,1) =-0.05*ones(floor(size(Long.time,1)*3/8),1);

% plot(Long.time,Long.signals.values)


%% Test evaluation
global  PERFORMANCE_ACC Response_ACC tau_ACC tau_q K_b_ACC K_b_q

params_ini = [0.3, 0.07,2.5 , 3];
tic
Cost = Cost_function( Long ,params_ini ) % params in this order: tau_ACC tau_q K_b_ACC K_b_q
toc
%% Iterative searching

% % [params_opt,cost] = fmincon( @(params) Cost_function( Long ,params ), params_ini,[],[],[],[],[0.1, 0.01,1, 1],[2, 1,20 , 20] ,[])

options = pso;
options.PopulationSize	= 20;
options.Vectorized		= 'off';
options.BoundaryMethod	= 'penalize';
options.Display			= 'iter';
options.HybridFcn		= @fmincon;
options.Generations     = 3;

My_pop=[params_ini;params_opt];
% for i=1:options.PopulationSize-1
%     My_pop=[ My_pop; rand(1,4).*([2, 1,20 , 20]-[0.1, 0.01,1, 1])+[0.1, 0.01,1, 1] ];
% end
options.InitialPopulation = My_pop;

[params_opt,cost] = pso(@(params) Cost_function( Long ,params ),4,[-1 5 0 0],0,[],[],[0.1, 0.01,1, 1],[2, 0.2,5 , 4] ,[],options)


% Optimal at first set-point of four vars = 
%% Check result
Long.signals.values( floor(size(Long.time,1)*5/8)+2:end ,1) =-0.05*ones(floor(size(Long.time,1)*3/8),1);
% params_opt =[1,0.2,5,2]
Cost = Cost_function( Long ,params_opt ) % params in this order: tau_ACC tau_q K_b_ACC K_b_q

% Plotting
figure
plot(PERFORMANCE_ACC.Time,PERFORMANCE_ACC.Data(:,2),'r')
hold on
plot(PERFORMANCE_ACC.Time,PERFORMANCE_ACC.Data(:,1))

if Cost>100
    params_opt=params_opt*NaN;
end
params_opt
I
J
pause(10)

% 
% plot(F_envelope.M_1g,F_envelope.H_1g)
% title('\fontsize{14} Flight Envelope and available data')
% xlabel('Mach')
% ylabel('H (m)')
% hold on
% plot(M_eval_sche,H_eval_sche,'*r')
% plot(M_eval_sche(I,J),H_eval_sche(I,J),'*g')
