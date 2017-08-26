%% Auto_trim aid 
% This script perform a trim routine to find out a 3D map of the trim
% states of AoA and SS at level flight.
% These 3D maps would implemente later on as bias to the pilot inouts, so
% he does nt have to retrim or constantly give an input to perorm a level
% flight

clear beta_map  alpha_map

load('FlightEnvelope.mat')
plot(F_envelope.M_1g ,F_envelope.H_1g )
hold on

%% Generation evauation points inside flight envelope
M_sq_autotrim = 0:0.1:max(F_envelope.M_1g);
H_sq_autotrim = min( F_envelope.H_1g):2000:max(F_envelope.H_1g);

M_eval=NaN*ones(size(M_sq_autotrim,2),size(H_sq_autotrim,2));
H_eval=NaN*ones(size(M_sq_autotrim,2),size(H_sq_autotrim,2));

for I=1:size(M_sq_autotrim,2)
    for J=1:size(H_sq_autotrim,2)
        if inpolygon(M_sq_autotrim(I),H_sq_autotrim(J),F_envelope.M_1g,F_envelope.H_1g)
            M_eval(I,J)=M_sq_autotrim(I);
            H_eval(I,J)=H_sq_autotrim(J);
        end
    end
end

% ADD ANY DESIRED POINT??
I=4;J=5;
M_eval(I,J)=M_sq_autotrim(I);
H_eval(I,J)=H_sq_autotrim(J);

I=3;J=2;
M_eval(I,J)=1.1*M_sq_autotrim(I);
H_eval(I,J)=H_sq_autotrim(J);

I=3;J=3;
M_eval(I,J)=0.219;
H_eval(I,J)=H_sq_autotrim(J);

I=7;J=8;
M_eval(I,J)=M_sq_autotrim(I);
H_eval(I,J)=13400;


I=9;J=9;
M_eval(I,J)=M_sq_autotrim(I);
H_eval(I,J)=15500;

plot(M_eval,H_eval,'*r')
return
%% Iterative evaluation of setpoints
alpha_map = ones(size(M_eval))*NaN;
beta_map = ones(size(M_eval))*NaN;

for I=1:size(M_sq_autotrim,2)
    for J=1:size(H_sq_autotrim,2)
        if isnan(M_eval(I,J))==0
            altitude = H_eval(I,J);
            [T,a,~,rho]=atmosisa(altitude);
            velocity = M_eval(I,J)*a;
            
            % Run tests
            q=1/2*rho*velocity^2;
            alpha = rad2deg((9800*9.81/(q*27)-0.015)/4.22);
            [trim_state, itrim_throttle, trim_control, dLEF, UX,cost] = trim_F16_fast(throttle, elevator,beta,  alpha, aileron, rudder, velocity, altitude);
            if cost<0.1
%                 keyboard
                alpha_map(I,J) = trim_state(8);
                beta_map(I,J) = trim_state(9);
            else
                % Assign the nearest value
                dist=100000;
                for I2=1:size(M_sq_autotrim,2)
                    for J2=1:size(H_sq_autotrim,2)
                        if isnan(M_eval(I2,J2))==0 && dist>sqrt((M_eval(I2,J2)-M_eval(I,J))^2+(H_eval(I2,J2)-H_eval(I,J))^2) && I2~=I && J2~=J
                            dist = sqrt((M_eval(I2,J2)-M_eval(I,J))^2+(H_eval(I2,J2)-H_eval(I,J))^2);
                            alpha_map(I,J) = alpha_map(I2,J2);
                            beta_map(I,J) = beta_map(I2,J2);
                        end
                    end
                end
                
            end

            
%             pause
        end
    end
end

figure(1)
plot(F_envelope.M_1g ,F_envelope.H_1g )
hold on;
surf(M_eval,H_eval,rad2deg(alpha_map))

figure
plot(F_envelope.M_1g ,F_envelope.H_1g )
hold on;
surf(M_eval,H_eval,rad2deg(beta_map))
return
%% Interpolate the maps and evaluate
[fitresult_alpha, gof] = createFit(M_eval, H_eval, alpha_map);
[fitresult_beta, gof] = createFit(M_eval, H_eval, beta_map);

alpha_map_eval = ones(size(M_eval))*NaN;
beta_map_eval = ones(size(M_eval))*NaN;

for I=1:size(M_sq_autotrim,2)
    for J=1:size(H_sq_autotrim,2)
        if isnan(M_eval(I,J))==0
            
            alpha_map_eval(I,J) = fitresult_alpha(M_sq_autotrim(I),H_sq_autotrim(J));
            beta_map_eval(I,J) = fitresult_beta(M_sq_autotrim(I),H_sq_autotrim(J));
            
        end
    end
end

pause
close all

figure(1)
plot(F_envelope.M_1g ,F_envelope.H_1g )
hold on;
surf(M_eval,H_eval,rad2deg(alpha_map_eval))

figure
plot(F_envelope.M_1g ,F_envelope.H_1g )
hold on;
surf(M_eval,H_eval,rad2deg(beta_map_eval))

%% Interpolate, Extrapolate and smooth the maps and evaluate
[fitresult_alpha, gof] = createFit_2(M_eval, H_eval, alpha_map_eval);
[fitresult_beta, gof] = createFit_2(M_eval, H_eval, beta_map_eval);

alpha_map_eval = ones(size(M_eval))*NaN;
beta_map_eval = ones(size(M_eval))*NaN;

for I=1:size(M_sq_autotrim,2)
    for J=1:size(H_sq_autotrim,2)
%         if isnan(M_eval(I,J))==0
            
            alpha_map_eval_smooth(I,J) = fitresult_alpha(M_sq_autotrim(I),H_sq_autotrim(J));
            beta_map_eval_smooth(I,J) = fitresult_beta(M_sq_autotrim(I),H_sq_autotrim(J));
            M_eval_autotrim_unc(I,J) = M_sq_autotrim(I);
            H_eval_autotrim_unc(I,J) = H_sq_autotrim(J);

%         end
    end
end

pause
close all

figure(1)
plot(F_envelope.M_1g ,F_envelope.H_1g )
hold on;
surf(M_eval_autotrim_unc,H_eval_autotrim_unc,rad2deg(alpha_map_eval_smooth))

figure
plot(F_envelope.M_1g ,F_envelope.H_1g )
hold on;
surf(M_eval_autotrim_unc,H_eval_autotrim_unc,rad2deg(beta_map_eval_smooth))
