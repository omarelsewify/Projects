clear
close all

%  FDM solver for Laplace's equation in 2D

% BCs:      phi = 10        (right)
%           phi = 100       (top)
%           dphi/dx = 0     (left)
%           dphi/dy = 0     (bottom)


%% Settings
% Resoluton
res = 20;

% Length of domian
L = 1;

% Grid spacing
h = L / (res - 1);

%% Constrcut matrices

% Looking to construct  linear system Ax = b

A = zeros(res * res);
phi = zeros(size(A, 1), 1);
b = zeros(size(phi));

% Populate the A matrix
for i = 1 : res
    for j = 1 : res
        
        % Get the poition in the vector for unknown phi_ij
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
        
        % Top Edge and the Top Left Corner
        elseif (i ~= res && j == res)
            
            % phi_ij = 100
            
            % Add coefficients to the A matrix and RHS
            A(ij, ij) = 1;
            b(ij) = 100;
            
        % Right Edge and the Bottom Right Corner
        elseif (j ~= res && i == res)
            
            % phi_ij = 10
            
            % Add coefficients to the A matrix and RHS
            A(ij, ij) = 1;
            b(ij) = 10;
            
        % Bottom Left Corner
        elseif (i == 1 && j == 1)
            
            % 1/h^2 (2 * ph_i+1j + 2 * phi_ij+1 - 4 * phi_ij) = 0
            
            % Get the Indices for the coefficiens
            ip1j = getId(i+1, j, res);
            ijp1 = getId(i, j+1, res);
            
            % Add coefficients to the A matrix
            A(ij, ij) = -4 / h^2;
            A(ij, ip1j) = 2 / h^2;
            A(ij, ijp1) = 2 / h^2;
            
        % Bottom Edge
        elseif (j == 1)
            
            % 1/h^2 (ph_i+1j + ph_i-1j + 2 * phi_ij+1 - 4 * phi_ij) = 0
            
            % Get the Indices for the coefficiens
            ip1j = getId(i+1, j, res);
            im1j = getId(i-1, j, res);
            ijp1 = getId(i, j+1, res);
            
            % Add coefficients to the A matrix
            A(ij, ij) = -4 / h^2;
            A(ij, ip1j) = 1 / h^2;
            A(ij, im1j) = 1 / h^2;
            A(ij, ijp1) = 2 / h^2;
            
        % Left Edge
        elseif (i == 1)
            
            % 1/h^2 (2 * ph_i+1j + phi_ij+1 + phi_ij-1 - 4 * phi_ij) = 0
            
            % Get the Indices for the coefficiens
            ip1j = getId(i+1, j, res);
            ijp1 = getId(i, j+1, res);
            ijm1 = getId(i, j-1, res);
            
            % Add coefficients to the A matrix
            A(ij, ij) = -4 / h^2;
            A(ij, ip1j) = 2 / h^2;
            A(ij, ijm1) = 1 / h^2;
            A(ij, ijp1) = 1 / h^2;
            
        % Domain (inner nodes)
        else
            
            % 1/h^2 (ph_i+1j + ph_i-1j + phi_ij+1 + phi_ij-1 - 4 * phi_ij) = 0
            
            % Get the Indices for the coefficiens
            ip1j = getId(i+1, j, res);
            im1j = getId(i-1, j, res);
            ijp1 = getId(i, j+1, res);
            ijm1 = getId(i, j-1, res);
            
            % Add coefficients to the A matrix
            A(ij, ij) = -4 / h^2;
            A(ij, ip1j) = 1 / h^2;
            A(ij, im1j) = 1 / h^2;
            A(ij, ijp1) = 1 / h^2;
            A(ij, ijm1) = 1 / h^2;            
            
        end        
    end
end

%% Solve system

% Use direct solver in MATLAB
phi = A \ b;

%% Plot results

% Create meshgrid of plotting points
[xplot, yplot] = meshgrid(linspace(0, L, res), linspace(0, L, res));

% Reshae 1D results to 2D grid
phi2 = reshape(phi, res, res);

% Plot the surface
subplot(1,2,1)
surf(xplot, yplot, phi2);
colorbar
title('Solution to Laplace''s Equation','interpreter','latex','Fontsize',14)
xlabel('x','interpreter','latex','Fontsize',14)
ylabel('y','interpreter','latex','Fontsize',14)
view(2)
axis equal
axis tight

% Find velocity
[ux, uy] = gradient(phi2);
ux = -ux;
uy = -uy;

mag = sqrt(ux.^2 + uy.^2);
uxn = ux ./ mag;
uyn = uy ./ mag;

% Plot the velocity
subplot(1,2,2)
quiver(xplot, yplot, uxn, uyn);
title('Normalisd Velocity Field','interpreter','latex','Fontsize',14)
xlabel('x','interpreter','latex','Fontsize',14)
ylabel('y','interpreter','latex','Fontsize',14)
axis equal
axis tight


figure
surf(xplot, yplot, phi2);
colorbar
title('Solution to Laplace''s Equation','interpreter','latex','Fontsize',14)
xlabel('x','interpreter','latex','Fontsize',14)
ylabel('y','interpreter','latex','Fontsize',14)
view(2)
axis equal
axis tight

figure
quiver(xplot, yplot, uxn, uyn);
title('Normalisd Velocity Field','interpreter','latex','Fontsize',14)
xlabel('x','interpreter','latex','Fontsize',14)
ylabel('y','interpreter','latex','Fontsize',14)
axis equal
axis tight