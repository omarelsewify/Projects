clear
%---------------------------------------------------------------
% Input parameters
%---------------------------------------------------------------
N = 6;              % Numbers of nodes in x-direction
L = 1;              % Length [m]
time = 0;           % Time
k = 40;             % Thermal conductivity [W/mK]
rho = 1000;         % Density [kg/m^3]
cp = 100;           % Specific heat [J/kgK]
tf = 100;           % Total time of simulation
dt = tf/(5*tf-1);   % dt value
alpha = k/(rho*cp); % Thermal diffusivity
dx = (L/(N-1));      % dx value
constant = 1/(alpha*dt) + 2/((L/(N-1))^2);
%---------------------------------------------------------------

%% Construct Matrices

% Looking to construct a linear system Ax=b
% Matrix A must be a square matrix with number of columns and rows being
% equal to the number of points to be simulated for (in this case 3000)

% Matrix x must be a vector with number of rows being equal to the
% number of data points that need to be found

% Matrix b is a vector the same size as x and it contains the solutions
% for all the equations 

A= zeros(5*N*tf-N+1);
x= zeros(size(A,1),1);
b= zeros(size(x));

while time<=tf
% Populate the A matrix by inserting the co-efffecients of the temperature
% equations in their respective positions

% Populate the b matrix by inserting the solutions of the temperature
% equations in their respective positions
for i= 1 : N
  for n = 1 : (tf/dt +1)
    % Get 1D index from i and n, use getID function which will assist in
    % the formation of the final results table
    in=getID(i, n, N);
    % Right Edge Boundary Condition
    if (i==N)
        % T_in = 100
        
        % Add coefficients into A and RHS of equation into b
         A(in ,in)=1;
         b(in)=100;
         
    % Left Edge Boundary Condition
    elseif (i==1)
        % T_in = 50
        
        % Add coefficients into A and RHS of equation into b
         A(in ,in)=1;
         b(in)=50;
         
    % Bottom Edge Boundary Condition
    elseif (n==1 && i~=1 && i~=N)
        % T_in = 20
        
        % Add coefficients into A and RHS of equation into b
         A(in ,in)=1;
         b(in)=20;
         
    % Domain Interior
    else
       % T_in = [(T_i+1n + T_i-1n)/dx^2 + (1/alpha)*(T_in-1/dt)]/constant
       
       % Get indices for the coefficients
       ip1n=getID(i+1,n,N);
       im1n=getID(i-1,n,N);
       inm1=getID(i,n-1,N);
       
       %Add coefficients into the A matrix
       A(in, in)= -1;
       A(in, ip1n)= (1/((dx^2)*constant));
       A(in, im1n)= (1/((dx^2)*constant));
       A(in, inm1)= (1/(alpha*dt*constant));
       
       % b values do not need to be changed as they are already set to zero
    end
  end
end
% March through code by changing the time by dt after every run
    time=time+dt;
end
%% Solve System

% Use direct solver in MATLAB
x= A\b;

%---------------------------------------------------------------
% Plotting temperatures along X at end of 100 seconds
%---------------------------------------------------------------
time=tf;
figure
axes('fontsize',16);
hold on
plot(linspace(0,L,N), x(end-(N-1):end),'k','linestyle', '-','linewidth',2);
hold off
xlabel('$x$ [$m$]','Interpreter','latex');
ylabel('$T$ [$^{\circ}{\rm C}$]','Interpreter','latex');
mytitle = ['t =\,',num2str(time),'\,s,\,','N =\,',num2str(N)];
title(mytitle,'Interpreter','latex');
grid on
box on
%---------------------------------------------------------------
% Plotting Temperature - Time graph at points 2 and 5
%---------------------------------------------------------------

% r2 T2 and r5 T5 combined create a data set containin the 2nd, 8th, 14th,...
% values in the x vector, these represent the progression throught time of
% the temperature at point 2

r2=(1:5*tf);
T2=x(2 + N*(r2-1));
r5=(1:5*tf);
T5=x(5 + N*(r5-1));
figure;
axes('fontsize',16);
hold on
plot(linspace(0,tf,5*tf), T2,'k','linestyle', '--','linewidth',2);
plot(linspace(0,tf,5*tf), T5,'k','linestyle', '-','linewidth',2);
hold off
xlabel('$t$ [$s$]','Interpreter','latex');
ylabel('$T$ [$^{\circ}{\rm C}$]','Interpreter','latex');
mytitle = 'N = 2 and N = 5';
title(mytitle,'Interpreter','latex');
legend({'N=2','N=5'},'Location','southeast')
grid on
box on