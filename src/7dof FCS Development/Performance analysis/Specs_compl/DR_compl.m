%% DR compl
close all
DR1
start=420;
time= DR1.Time(start:end);
p_DR=DR1.Data(start:end,1);
r_DR=DR1.Data(start:end,2);

plot(time,r_DR)
% return
[p_max,idx_max_p]=findpeaks(r_DR)
hold on
plot(time(idx_max_p),p_max,'r*')
[p_min,idx_min_p]=findpeaks(-r_DR)
plot(time(idx_min_p),-p_min,'r*')

% Damped freq
peaks = sort([idx_max_p;idx_min_p]);

for i=2:size(peaks,1)
    Delta(i-1)=time(peaks(i))-time(peaks(i-1));
end

omega_d = mean(2*pi./Delta);

% Damping ratio
plot(time(peaks),abs(r_DR(peaks)),'g*')

time_fit=time(peaks)-time(peaks(1));
p_fit=abs(r_DR(peaks));
[fitresult, gof] = Exp_fit(time_fit, p_fit);  % a + b*exp(-tau*x)

sigma = fitresult.tau;

xi= sqrt(1/(1+(omega_d/sigma)^2))
omega_n= omega_d/sqrt(1-xi^2)
return
%% step responses
plot(time(peaks:end),r_DR(peaks:end))
Response_DR=stepinfo( r_DR(peaks:end),time(peaks:end),r_DR(end))

