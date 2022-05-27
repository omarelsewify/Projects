 
clear; close all; clc;
 
h = 1;
n = 11;
T = ones(n,1); Told = T;
T(1) = 1; %Left boundary
T(n) = 10; %Right boundary
x = linspace(0,1,n);
dx = x(2)-x(1);
dt = 0.5*dx^2/3; %cfl condition 
 
error = 1;
TOL = 1e-6;
k = 0;
while error > TOL,
   Told = T;
   k = k+1;
   for i = 2:n-1
    T(i) = dt*(Told(i+1)-2*Told(i)+Told(i-1))/dx^2+Told(i);
   end
   error = max(abs(T-Told));
   if mod(k,5)==0, out(k,:) = T; end    
end
 
plot(x,out)
xlabel('x'),ylabel('Temperature'),
title(['Fourier Heat Conduction']),
%legend('Cooling Trend','Steady State')
