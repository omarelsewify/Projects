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

% Convection heat transfer coefficient (W/Km^2)
h= 3;
% input('Enter the convection heat transfer coefficient (W/Km^2) : ');

% Conduction heat transfer coefficient (W/Km)
k = 30;

% Initial Ambient Temp (C)
T_amb = 25;
%input('Enter the initial ambient temperature (C) : ');

%% Construct matrices 
% Looking to construct  linear system Ax = b
 
A = zeros(res * res);
phi = zeros(size(A, 1), 1);
b = zeros(size(phi));

tfinal=540;
time=1:tfinal;
%prompt= 'What is the probe point in time? ';
    %tp=input(prompt);
 
% Populate the A matrix
for t=1:tfinal
        I = size(time);

        % Relax stage
        I(1:61) = 60;
        % Discharge stage
        I(61:120) = -20;
        % Relax stage
        I(121:240) = 0;
        % Charge stage
        I(241:tfinal) = 4;

        % Internal Resistance (Ohms)
        Rin = 0.00154;
        % Tab Resistance (Ohms)
        Rtab = 9.74e-5;
        % Heat Flux Total
        qgen=size(time);
        T_hg=size(time);
        
        for t=1:tfinal
            qgen(t)=(I(t)^2)*Rin;
        end
    % Initial qin from tabs
    qtab = 120; %(265^2)*(9.74e-5);   % W/m^2
    Ttab = qtab/h + T_amb;

    % qout from sides
    qout = 50;   % W/m^2
    Twall = qout/h + T_amb;
    
    for i = 1 : res
        for j = 1 : res
            
            % Get the position in the vector for unknown phi_ij
            ij = getId(i, j, res);

            % Top Right Corner
            if (i == res && j == res)

                % phi_ij - 0.5 * (phi_i-1j + phi_ij-1) = 0

                % Get indices for the coefficients
                im1j = getId(i-1, j, res);
                ijm1 = getId(i, j-1, res);

                % Add coefficients to the A matrix
                A(ij, ij) = 1;
                A(ij, im1j) = -0.5;
                A(ij, ijm1) = -0.5;

            % Top Left Corner
            elseif (i == 1 && j == res)

                % phi_ij - 0.5 * (phi_i+1j + phi_ij-1) = 0

                % Get indices for the coefficients
                ip1j = getId(i+1, j, res);
                ijm1 = getId(i, j-1, res);

                % Add coefficients to the A matrix
                A(ij, ij) = 1;
                A(ij, ip1j) = -0.5;
                A(ij, ijm1) = -0.5;

            % Right Edge and the Bottom Right Corner
            elseif (j ~= res && i == res)

                % phi_ij = T_amb

                % Add coefficients to the A matrix and RHS
                A(ij, ij) = 1;
                b(ij) = Twall;

            % Left Edge and the Bottom Left Corner
            elseif (j ~= res && i == 1)

                % phi_ij = T_amb

                % Add coefficients to the A matrix and RHS
                A(ij, ij) = 1;
                b(ij) = Twall;    

            % Bottom Edge
            elseif (j == 1 )

                % phi_ij = 0

                % Add coefficients to the A matrix and RHS
                A(ij, ij) = 1;
                b(ij) = Twall;

            % Top Edge
            elseif (j == res )

                % phi_ij = T_amb

                % Add coefficients to the A matrix and RHS
                A(ij, ij) = 1;
                b(ij) = Twall;

            % Domain (inner nodes)
            else

             % 1/h^2 (ph_i+1j + ph_i-1j + phi_ij+1 + phi_ij-1 - 4 * phi_ij) = 0

                % Get the Indices for the coefficients
                ip1j = getId(i+1, j, res);
                im1j = getId(i-1, j, res);
                ijp1 = getId(i, j+1, res);
                ijm1 = getId(i, j-1, res);

                % Add coefficients to the A matrix
                A(ij, ij) = -4 / gs^2;
                A(ij, ip1j) = 1 / gs^2;
                A(ij, im1j) = 1 / gs^2;
                A(ij, ijp1) = 1 / gs^2;
                A(ij, ijm1) = 1 / gs^2;   

                b(ij)=-T_hg;
            end

            % Left Tab 
            for LT = rdfive+1 : (2*rdfive)
            if (j == res && i == LT)

                % phi_ij = qtab

                % Add coefficients to the A matrix and RHS
                A(ij, ij) = 1;
                b(ij) = Ttab; 
            end
            end

            % Right Tab 
            for RT = (3*rdfive)+1 : (4*rdfive)
            if (j == res && i == RT)

                % phi_ij = qtab

                % Add coefficients to the A matrix and RHS
                A(ij, ij) = 1;
                b(ij) = Ttab;
            end
            end

            % Empty Top Edges
            for ETE = (2*rdfive)+1 : (3*rdfive) 
            if (j == res && i == ETE)
                % phi_ij = T_amb

                % Add coefficients to the A matrix and RHS
                A(ij, ij) = 1;
                b(ij) = Twall;

            end
            end
       end 
    end
end 
%% Solve system
    % Use direct solver in MATLAB
    phi = A \ b;

    % Reshape 1D results to 2D grid
    phi2 = reshape(phi, res, res);
    
%% Plot results
 
% Create meshgrid of plotting points
[xplot, yplot] = meshgrid(linspace(0, W, res), linspace(0, L, res));
 
% Find Temp Gradient
[ux, uy] = gradient(phi2);
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
contourf(xplot, yplot, phi2,8);
colorbar
colormap(hot)
title('Solution to Laplace''s Equation','interpreter','latex','Fontsize',14)
xlabel('x','interpreter','latex','Fontsize',14)
ylabel('y','interpreter','latex','Fontsize',14)
hcb=colorbar;
caxis([Twall 70])
title(hcb,'Temperature','interpreter','latex','Fontsize',10)
view(2)
axis equal
axis tight
