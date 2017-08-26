% Initialization parameters for IMU navigation model
%Sample time
TSample=T_sim;

%Gyro parameters
GYRO.bias_repeat= 0.003*(pi/180)/3600; % [rad/s] 0.003º/hr
GYRO.bias_stability=0.01*(pi/180)/3600; % [rad/s] 0.01º/hr
GYRO.noise_roll= 0.0002*(pi/180)/sqrt(3600); %ARW [rad/hr^0.5]
GYRO.noise_pith=0.0002*(pi/180)/sqrt(3600);%ARW [º/hr^0.5]
%Accelerometer parameter
ACCEL.bias_repeat= 25e-6; % [g] 
ACCEL.bias_stability=0.3e-6; % [g] 
ACCEL.noise= 0.01981/sqrt(3600); %ARW [m/s/hr^0.5]
