%This m-file solves the 2D transient conduction in a given rectangular...
%domain in x-y plane. In this case, boundary condition is constant...
%temperatures of T1...T4 for right, top, left and bottom edges...
%respectively. Heat generation can be added by setting non-zero positive...
%values for g.

clc;
clear all;
close all;

%----------------initial variables--------------------------------------
M = 41; %grid number for x axis
N = 121; %grid number for y axis

X = 1; %'rectangular width in x dimension',unit is m
Y = 3; %'rectangular width in Y dimension',unit is m

g = 100000; %'heat generation'
T1 = 400; %constant temperature at right side,unit is K
T2 = 400; %constant temperature at top side,unit is K
T3 = 400; %constant temperature at left side,unit is K
T4 = 400; %constant temperature at bttom side,unit is K

Ti = 298; %initial temperature of the plate
k = 237; %conduction coefficience,unit is W/m*K
alpha = 9.71e-5; %thermal diffusivity of the material, unit is m^2/s
EPS = 5e-9; %'convergence accuracy'
IfConverged = 0; %The 'IfConverged'parameter represents if the steady state has been reached.

dx = X/(M-1); %width of the grid
dy = Y/(N-1); %length of the grid

imove = 0;
percenterror=0;

%--------------------setup dt and the experiment time-------------
Fo = 0.1; % by making sure your Fo number is below 0.25 your solution will stay stable
experi_t = 5000; %set the experiment time, here the time is set up to 800 seconds
delta_t = Fo*(dx^2)/alpha; %set the basic dt. make sure to satisfy convergence condition.
counter_t = fix ( experi_t / delta_t); %set up the counter
minimum_t = experi_t; %initialize the minimum time equal to experiment time

%------------------initiate the temperature array of the plate-------
T = ones(M,N); %T array holds the temperature caculated from the last iteration
Tnew = ones(M,N) * Ti;
Tvstime = zeros(M,N,fix(counter_t/100)+1); % to save every 100 temp. profiles if your
% memory permits.
%------------------apply boundary condition--------------------
% since we have constant temperature boundary condition, we can apply it
% out of the main loop (boundary temperatures will not change inside the main loop).

%right boundary condition(constant temperature T1)
for j=1:N
i=M;
Tnew(i,j)=T1;
end
%top boundary condition (constant temperature T2)
for i=2:M-1
j=N;
Tnew(i,j)=T2;
end
%left boundary condition (constant temperature T3)
for j=1:N
i=1;
Tnew(i,j)=T3;
end
%bottom boundary condition(constant temperature T4)
for i=2:M-1
j=1;
Tnew(i,j)=T4;
end
%treat the corner points seperately in order to avoid singular gridpoints.
Tnew(1,N) = 0.5*(T3+T2);
Tnew(M,N) = 0.5*(T1+T2);
Tnew(M,1) = 0.5*(T4+T1);
Tnew(1,1) = 0.5*(T4+T3);
T = Tnew;

%---------------------main loop-------------------------------
%this loop keeps recaculating the temperature until the plate reach steady
%state, and caculate the time used to reach the steady state.

while counter_t > 0

%interior node
for i=2:M-1
for j=2:N-1
Tnew(i,j) = Fo*(T(i-1,j) +T(i+1,j)+ T(i,j-1)+T(i,j+1)) + (1-4*Fo)*T(i,j)+g*(dx)*(dx)/k*Fo ;
end
end

IfConverged = 1;
Tabs = abs(Tnew-T); %Tabs is used to check if steady state is reached.

for i=1:M %using absolute value of (Tnew-T)
for j=1:N
if Tabs(i,j)>EPS %if every element of Tabs is smaller than EPS, then the temperature reach the steady state, else not.
IfConverged=0;
end
end
end

T = Tnew; %copy the newly caculated temp value to the old array

if IfConverged == 1% if the steady state is reached
% the minimum_t is euqal to total experiment time minus time left,
% the value of counter_t*delta_t is the time left.
minimum_t = experi_t - counter_t * delta_t;
counter_t = 0; %clear the counter then exit the experiment
else
counter_t = counter_t - 1 %show on screen
end
if ( mod(counter_t,100)==0 )
imove = imove + 1;
Tvstime(:,:,imove) = Tnew(:,:); %keep Tnew in memory (every 100 steps), so that you can have intermediate temp profiles.
end

end

%---------compute heat conducted from edges----------------
% loss --> negative, gain--> positive

%q1 = -k * ( T(M-1,:) - T(M,:) ) * dy / dx; %1st order accurate RIGHT
q1 = -k * ( T(M-2,:) - 4 * T(M-1,:) + 3 * T(M,:)) * dy / (2*dx); %2nd order accurate
q1 (1,1) = q1 (1,1) / 2;
q1 (1,N) = q1 (1,N) / 2;

%q3 = -k * ( T(2,:) - T(1,:) ) * dy / dx; %1st order accurate TOP
q3 = -k * ( T(3,:) - 4 * T(2,:) + 3 * T(1,:) ) * dy / (2*dx); %2nd order accurate
q3 (1,1) = q3 (1,1) / 2;
q3 (1,N) = q3 (1,N) / 2;

%q2 = -k * ( T(:,N-1) - T(:,N) ) * dx / dy; %1st order accurate LEFT
q2 = -k * ( T(:,N-2) - 4 * T(:,N-1) + 3 * T(:,N)) * dx / (2*dy); %2nd order accurate
q2 (1,1) = q2 (1,1) / 2;
q2 (M,1) = q2 (M,1) / 2;

%q4 = -k * ( T(:,2) - T(:,1) ) * dx / dy; %1st order accurate  BOTTOM
q4 = -k * ( T(:,3) - 4 * T(:,2) + 3 * T(:,1)) * dx / (2*dy); %2nd order accurate
q4 (1,1) = q4 (1,1) / 2;
q4 (M,1) = q4 (M,1) / 2;

% qtot is the total amount of heat conducted from 4 edges which at steady state must be
% equal to qgeneration ideally
qgeneration = g * X * Y;
qtot = sum (q1) + sum(q2) + sum(q3) + sum(q4);
percenterror = abs((qtot-qgeneration)/qtot)*100;

%---------display the caculated temperature in the plate----------------

if IfConverged==1 %if the steady state is reached
steadytime = experi_t - minimum_t;
disp('--------------------------------------------------------');
disp('Steady state is reached. Time to reach steadys state in minutes is:');
disp(steadytime/60);
disp('seconds');
else %if the steady state is not reached
runtime = minimum_t;
disp('---------------------------------------------------');
disp('The steady state is not reached. Total run time is:');
disp(runtime);
disp('seconds');
end

disp('total heat rate leaving from the boundaries in Watts:');
disp(qtot);
disp('total energy being generated in the block in Watts:');
disp(qgeneration);
disp('Percent error in energy balance calculation over the entire block');
disp(percenterror);
subplot (1,2,1);
contourf(T',10); % plot 10 contours of the final temperature
colorbar; %show the colorbar of the graph, specific color shows the value of the specific temperature range.
title('Final temperature distribution');
xlabel('x(m)'); % anotation of the x axis
ylabel('y(m)'); %anotation of the y axis

subplot (1,2,2); % plot 10 contours of temperature at intermediate times
Tinter = Tvstime(:,:,floor(imove/10));
contourf(Tinter',10); % plot 10 contours of the temperature
colorbar; %show the colorbar of the graph, specific color shows the value of the specific temperature range.
title('temperature distribution at intermediate time');
disp(imove);
xlabel('x(m)'); % anotation of the x axis
ylabel('y(m)'); %anotation of the y axis

