clc,clear,close all
                            %%  FDM solver for Laplace's equation in 2D %%
%% Settings
% Resolution
res = 60;           
left1 = res/15;
left2 = res*7/15;
right1 = res*8/15;
right2 = res*14/15;
 
% Width of domain
W = 0.20;
 
% Height of domain
L = 0.20;
 
% Grid spacing
gs = W / (res - 1);
 
%% Battery and Thermal Characteristics
 
k = 68;                         % thermal conductivity (W/m.K)
cp = 1266;                      % specific heat (J/kg K)
rho = 2548.2;                   % density (kg/m^3)
alpha = k/(cp*rho);             % alpha (1/sec)
 
R             = 0.00154 ;               % Ohm Resistance of the battery
R_tabP        = 9.74e-5 ;               % Ohm Resistance of the positive tab
W_tabP        = 0.08 ;                  % Positive tab width
R_tabN        = 6.88e-5 ;               % Ohm Resistance of the negative tab
W_tabN        = 0.08 ;                  % Negative tab width
Bat_cap       = 53.0 ;                  % Battery capacity (Ah)
C_rate        = 3.0 ;                   % Charge/discharge rate
I             = Bat_cap * C_rate ;      % Current of the battery: 5C * 53Ah
 
%% Time parameters
 
tfinal = 1100;
time = 0;
dt = 0.1;

Fo = (alpha*dt)/(gs^2);             % Fourier Number

%% Setting boundary conditions
 
Twall = 25;
phi = Twall*ones(res);
 
%% Sensor 1
X_sensor1 = 30;
Y_sensor1 = 30;
Tt1 = [time phi(Y_sensor1, X_sensor1)] ;   % Tt=[0 25]
 
%% Sensor 2
X_sensor2 = 45;
Y_sensor2 = 12;
Tt2 = [time phi(Y_sensor2, X_sensor2)] ;   % Tt=[0 25]
 
%% Heat Flux Calculation
 
    qgen = ((I^2)*R)/(W*L);
    qtab_right = ((I^2)*R_tabN/W_tabN);
    qtab_left = ((I^2)*R_tabP/W_tabP);
    h = 0;
    
    Bi = h*0.2/k;
    
% Creating a copy of phi
  phi_old = phi;
  
while time < tfinal  
    
    % Internal Domain
    for i = 2:res-1
        for j = 2:res-1
            
            phi(i,j) = Fo*(phi_old(i+1,j)-4*phi_old(i,j)+phi_old(i-1,j)+phi_old(i,j+1)+phi_old(i,j-1))+...
                            phi_old(i,j)+qgen*alpha*dt;

        end
    end
    
    % Right Tab
    for i = right1:right2
        for j = res           
            phi(i,j) = Fo*(phi_old(i-1,j)+2*phi_old(i,j-1)+phi_old(i+1,j)+ 2*gs*qtab_right/398 +(1/Fo -4)*phi_old(i,j));       
        end
    end

    % Left Tab
    for i = left1:left2
        for j = res            
            phi(i,j) = Fo*(phi_old(i-1,j)+2*phi_old(i,j-1)+phi_old(i+1,j)+ 2*gs*qtab_left/398 +(1/Fo -4)*phi_old(i,j));
        end
    end
    
    % Bottom Edge
    for i = 2:res-1
        for j = 1                  
            phi(i,j) = Fo*(phi_old(i-1,j)+2*phi_old(i,j+1)+phi_old(i+1,j)+ 2*Bi*25 +(1/Fo -4 -2*Bi)*phi_old(i,j));
        end
    end
    
    % Right Edge
    for i = res
        for j = 2:res-1           
            phi(i,j) = Fo*(phi_old(i,j+1)+2*phi_old(i-1,j)+phi_old(i,j-1)+ 2*Bi*25 +(1/Fo -4 -2*Bi)*phi_old(i,j));      
        end
    end
    
    % Left Edge
    for i = 1
        for j = 2:res-1     
            phi(i,j) = Fo*(phi_old(i,j+1)+2*phi_old(i+1,j)+phi_old(i,j-1)+ 2*Bi*25 +(1/Fo -4 -2*Bi)*phi_old(i,j));
        end
    end
    
    % Top Edge Gaps
    for i = [2:left1-1  left2+1:right1-1  right2+1:res-1]
        for j = res                  
            phi(i,j) = Fo*(phi_old(i-1,j)+2*phi_old(i,j-1)+phi_old(i+1,j)+ 2*Bi*25 +(1/Fo -4 -2*Bi)*phi_old(i,j));
        end
    end
    
    % Bottom Left Corner
    for i = 1 
        for j = 1     
            phi(i,j) = 0.5*(phi_old(i,j+1)+phi_old(i+1,j));
        end
    end

    % Bottom Right Corner
    for i = res 
        for j = 1     
            phi(i,j) = 0.5*(phi_old(i,j+1)+phi_old(i-1,j));
        end
    end
    
    % Top Left Corner
    for i = 1 
        for j = res     
            phi(i,j) = 0.5*(phi_old(i,j-1)+phi_old(i+1,j));
        end
    end
    
    % Top Right Corner
    for i = res 
        for j = res     
            phi(i,j) = 0.5*(phi_old(i,j-1)+phi_old(i-1,j));
        end
    end
    
    phi_old = phi;
    time = time + dt;
    
    Tt1 = [Tt1 ; time phi(Y_sensor1, X_sensor1)] ; %% Tt1= [0 25 ; 0.1 25.1]
    Tt2 = [Tt2 ; time phi(Y_sensor2, X_sensor2)] ; %% Tt2= [0 25 ; 0.1 25.1]

end

%% Plot results

load colormap.mat
phi=transpose(phi);

% Create meshgrid of plotting points
[xplot, yplot] = meshgrid(linspace(0, W, res), linspace(0, L, res));
 
% Find Temp Gradient
[ux, uy] = gradient(phi);
ux = -ux;
uy = -uy;
 
mag = sqrt(ux.^2 + uy.^2);
uxn = ux ./ mag;
uyn = uy ./ mag;
 
figure
quiver(xplot, yplot, uxn, uyn);
title('Normalised Heat Flux Field','interpreter','latex','Fontsize',14)
xlabel('x','interpreter','latex','Fontsize',14)
ylabel('y','interpreter','latex','Fontsize',14)
axis equal
axis tight
 
figure
contourf(xplot, yplot, phi);
colorbar
colormap(map)
title('Solution to Laplace''s Equation','interpreter','latex','Fontsize',14)
xlabel('x','interpreter','latex','Fontsize',14)
ylabel('y','interpreter','latex','Fontsize',14)
hcb=colorbar;
caxis
title(hcb,'Temperature','interpreter','latex','Fontsize',10)
view(2)
axis equal
axis tight

% Sensor 1
figure 
plot(Tt1(1:end,1),Tt1(1:end,2))
title_text=sprintf('Temperature variation vs time at point (%d,%d) ',X_sensor1,Y_sensor1);
title(title_text,'interpreter','latex','Fontsize',14);
xlabel('Time (s)','interpreter','latex','Fontsize',14)
ylabel('Temperature (C)','interpreter','latex','Fontsize',14)

% Sensor 2
figure 
plot(Tt2(1:end,1),Tt2(1:end,2))
title_text=sprintf('Temperature variation vs time at point (%d,%d) ',X_sensor2,Y_sensor2);
title(title_text,'interpreter','latex','Fontsize',14);
xlabel('Time (s)','interpreter','latex','Fontsize',14)
ylabel('Temperature (C)','interpreter','latex','Fontsize',14)
