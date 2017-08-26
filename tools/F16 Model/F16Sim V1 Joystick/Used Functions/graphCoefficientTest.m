function graphCoefficientTest()

close all;

%% Change When Appropriate %%

Coef = 'Cy';
Beta = '0.0';
LEF = '';
%LEF = ' LEF = 0.0 deg' 


data_file = sprintf('Cx_file.txt');
[alpha, elevator, Cx_lo, Cx_hi] = textread(data_file, '%f %f %f %f', 'delimiter', ',');

lowerbound = -20;
upperbound = 20;
interval = 1;

k = 1;

number = (upperbound-lowerbound)/interval + 1;

X = zeros(number,number);
Y = zeros(number,number);
Z = zeros(number,number);

for i = 1:1:number,
    for j = 1:1:number,
        X(i,j) = alpha(k);
        Y(i,j) = elevator(k);
        Z_lo(i,j) = Cx_lo(k);
        Z_hi(i,j) = Cx_hi(k);
        k = k+1;
    end
end


figure
surf(X,Y,Z_lo)
title(sprintf('LOFI %s Beta = %s deg', Coef, Beta) )
%view(140, 40)
xlabel('alpha')
ylabel('elevator')
zlabel(sprintf('%s',Coef))
zlim([0.0 0.4])

figure
surf(X,Y,Z_hi)
title(sprintf('HIFI %s Beta = %s deg %s', Coef, Beta, LEF ))
view(140, 40)
xlabel('alpha')
ylabel('elevator')
zlabel(sprintf('%s',Coef))
%zlim([-0.5 0.5])