 clear all;
   close all;
   clc;
   
   %inputs
   L =1; % length of domain and it is a square
   %nx,ny no of grid points
   nx = 51;
   ny = nx;
   x =linspace(0,L,nx);
   y =linspace(0,L,ny);
   dx = x(2)-x(1);
   dy = y(2)-y(1);
   error = 9e9;
   tol = 1e-4;
   dt=0.1;
   alpha = 1e-3; %convective coeffecient
   t=5;
   
    %CFL number
  M = alpha*dt/(dx^2);
  
  % Define boundary conditions
  T_L = 400;
  T_T = 600;
  T_R = 800;
  T_B = 900;
  T = 298*ones(nx,ny);
  %Temperature at boundary
 T(:,1)=T_L;
 T(:,end)=T_R;
 T(end,:)=T_B;
 T(1,:)=T_T;
 
 %  Average Temperature at the corner matrix
 T(1,1) = (T_T+T_L)/2;
 T(1,end) = (T_T+T_R)/2;
 T(end,1) = (T_L+T_B)/2;
 T(end,end) = (T_R+T_B)/2;
 %Calculation of temperature distributon explicitly using iterative solvers
 k1 = (alpha*dt)/(dx^2);
 k2 = (alpha*dt)/(dy^2);
 
 %creating a copy of T
  [X,Y] = meshgrid(x,y);
  Told =T;
  nt=t/dt;
  % 2d transient unsteady state heat conduction equation in explicit method
  iterative_solver =1;
  tic;
  iterations =1;
  for v =1:nt
  if iterative_solver ==1
    
    error = 100;
    while(error >tol)
    for i=2:nx-1
     for j=2:ny-1
     term1= Told(i,j);
     term2=k1*(Told(i+1,j)-2*Told(i,j)+Told(i-1,j));
     term3=k2*(Told(i,j+1)-2*Told(i,j)+Told(i,j-1));
     T(i,j)=term1+term2+term3;
   end
 end
error = max(abs(Told-T));
    Told = T;
    iterations = iterations+1;
  end
end
time=toc;
%plotting 
figure(1);
[P,Q] = contourf(X,Y,T);
clabel(P,Q,"fontsize",12);
colormap(jet);
colorbar;
title_text=sprintf('Explicit method ,Iteration number=%d,computation time=%fs,M=%f ',iterations,time,M);
title(title_text);
xlabel('Xaxis');
ylabel('Yaxis');
pause(0.03);
  end