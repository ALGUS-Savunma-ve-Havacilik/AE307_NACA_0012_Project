%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% NACA_0012 (NACA_0012.m)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Nicholas Ngo Syuan Yaw (ERAU)
% AE307 01DB
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% NACA 0012
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear
clc
close

filename = 'NACA_0012_edited.txt';

data = dlmread(filename);               % data read

alpha = data(:,1);                      % angle of attack
Cl = data(:,2);                         % coefficient of lift

plot(alpha,Cl);
title('Cl vs Alpha (NACA 0012)');
xlabel('Angle of Attack');
ylabel('Coefficient of Lift');
grid on
hold on

ClTA = 2*pi*alpha;                      % Cl for thin airfoil theorem
C_ClTA = ClTA*(pi/180);                 % Data Conversion

plot(alpha,C_ClTA);
legend('Cl (XFLR5 v6.40)','Cl (Thin Airfoil Theorem)','Location','southeast')
hold off

%end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Nicholas Ngo (2018)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~