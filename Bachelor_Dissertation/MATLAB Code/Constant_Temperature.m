clc,clear,close all
                            %%  FDM solver for Laplace's equation in 2D %%
%% Settings
% Resolution
res = 50;           % ONLY  MULTIPLES OF 10
rdfive = res/5;
 
% Width of domain
W = 0.208;

% Height of domain
L = 0.184;

% Grid spacing
gs = W / (res - 1);
 
%% Battery and Thermal Characteristics

conductivity = 68;                                        % thermal conductivity (W/m.K)
specific_heat = 1266;                                     % specific heat (J/kg K)
density = 2548.2;                                         % density (kg/m^3)
alpha = conductivity/(specific_heat*density);             % alpha (1/sec)

%% Time parameters

tfinal = 500;
time = 0.1;
dt = 0.1;

%% Setting boundary conditions

Twall = 25;
Ttab = 60;
phi = 25*ones(50);

phi(:,1)=Twall;
phi(:,end)=Twall;
phi(end,:)=Twall;
phi(1,:)=Twall;

% Right Tab
phi(res, (3*rdfive)+1 : (4*rdfive))= Ttab;
% Left Tab
phi(res, (rdfive)+1 : (2*rdfive))= Ttab;


%% Sensor 1
X_sensor1 = 25;
Y_sensor1 = 25;
Tt1 = [time phi(Y_sensor1, X_sensor1)] ;   % Tt=[0 25]

%% Sensor 2
X_sensor2 = 12;
Y_sensor2 = 40;
Tt2 = [time phi(Y_sensor2, X_sensor2)] ;   % Tt=[0 25]

%% Current Schedule
for t=1:tfinal
        I = size(time);

        % Relax stage
        I(1:61) = 0;
        % Discharge stage
        I(61:120) = -20;
        % Relax stage
        I(121:240) = 0;
        % Charge stage
        I(241:tfinal) = 4;

        % Internal Resistance (Ohms)
        Rin = 50e-3;
        % Heat Flux Total
        qgen=size(time);
        
        for t=1:tfinal
            qgen(t)=(I(t)^2)*Rin;
        end
end
% Creating a copy of phi
  phi_old = phi;
  
while time < tfinal    
    for i = 2:res-1
        for j = 2:res-1
            
            term1= phi_old(i,j);
            term2=((alpha*dt)/gs^2)*(phi_old(i+1,j)-2*phi_old(i,j)+phi_old(i-1,j));
            term3=((alpha*dt)/gs^2)*(phi_old(i,j+1)-2*phi_old(i,j)+phi_old(i,j-1));
            term4=(qgen(70)*dt)/(specific_heat*density);
            
            phi(i,j) = term1+term2+term3+term4;

        end
    end

    phi_old = phi;
    time = time + dt;
    Tt1 = [Tt1 ; time phi(Y_sensor1, X_sensor1)] ; %% Tt1= [0 25 ; 0.1 25.1]
    Tt2 = [Tt2 ; time phi(Y_sensor2, X_sensor2)] ; %% Tt2= [0 25 ; 0.1 25.1]

end

%% Plot results
 
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
contourf(xplot, yplot, phi, 8);
colorbar
colormap(hot)
title('Solution to Laplace''s Equation','interpreter','latex','Fontsize',14)
xlabel('x','interpreter','latex','Fontsize',14)
ylabel('y','interpreter','latex','Fontsize',14)
hcb=colorbar;
caxis([25 70])
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
