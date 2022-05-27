
clear
close all
%  FDM solver for Laplace's equation in 2D
 
% BCs:      phi = qtab       (right tab)
%           phi = qtab       (left tab)
%           dphi/dx = 0     (left and  right)
%           dphi/dy = 0     (bottom and gaps on top)
 
%% Settings
% Resolution
res = 20;           % ONLY  MULTIPLES OF 10
rdfive = res/5;
 
% Length of domain
L = 1;
 
% Grid spacing
gs = L / (res - 1);
 
%% Battery and Thermal Characteristics
% Initial qin from tabs
qtab = 373;
 
% Initial Ambient Temp (K)
T_amb = 298;

% Initial Coolant Temp (K)
T_cool = 285;

% Internal Temp Generated
in_hg = 0;
 
% Convective heat transfer coefficient (W/Km^2)
h = 35;

% Conductive heat transfer coefficient (W/Km)
k = 0.7;
 
%% Construct matrices 
% Looking to construct  linear system Ax = b
 
A = zeros(res * res);
phi = zeros(size(A, 1), 1);
b = zeros(size(phi));
 
% Populate the A matrix
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
            
        % Right Edge 
        elseif (j ~= res && i == res && j~=1)
            
            % Get the Indices for the coefficients
            ip1j = getId(i+1, j, res);
            im1j = getId(i-1, j, res);
            ijp1 = getId(i, j+1, res);
            ijm1 = getId(i, j-1, res);
            
            % Add coefficients to the A matrix
            A(ij, ij) = (1 / gs^2) * (-2 * gs * (h/k) - 4);
            A(ij, im1j) = (1 / gs^2) * 2;
            A(ij, ijp1) = 1 / gs^2;
            A(ij, ijm1) = (1 / gs^2);   
            
            b(ij)=-in_hg - 2 * gs * (h/k)* T_cool;  
            
        % Left Edge 
        elseif (j ~= res && i == 1 && j ~= 1)
                        
            % Get the Indices for the coefficients
            ip1j = getId(i+1, j, res);
            im1j = getId(i-1, j, res);
            ijp1 = getId(i, j+1, res);
            ijm1 = getId(i, j-1, res);
            
            % Add coefficients to the A matrix
            A(ij, ij) = (1 / gs^2) * (-2 * gs * (h/k) - 4);
            A(ij, ip1j) = (1 / gs^2) * 2;
            A(ij, ijp1) = 1 / gs^2;
            A(ij, ijm1) = (1 / gs^2);   
            
            b(ij)=-in_hg - 2 * gs * (h/k)* T_cool;  
        
        % Bottom Edge
        elseif (j == 1 && i~=1 && i~=res)
            
            % Get the Indices for the coefficients
            ip1j = getId(i+1, j, res);
            im1j = getId(i-1, j, res);
            ijp1 = getId(i, j+1, res);
            ijm1 = getId(i, j-1, res);
            
            % Add coefficients to the A matrix
            A(ij, ij) = (1 / gs^2) * (-2 * gs * (h/k) - 4);
            A(ij, im1j) = (1 / gs^2);
            A(ij, ijp1) = (1 / gs^2)*2;
            A(ij, ip1j) = (1 / gs^2);   
            
            b(ij)=-in_hg - 2 * gs * (h/k) * T_cool;
            
            
       % Bottom Left Corner
        elseif (i == 1 && j == 1)

            % 1/h^2 (2 * ph_i+1j + 2 * phi_ij+1 - 4 * phi_ij) = 0

            % Get the Indices for the coefficiens
            ip1j = getId(i+1, j, res);
            ijp1 = getId(i, j+1, res);

            % Add coefficients to the A matrix

            A(ij, ij) = -4 / gs^2;
            A(ij, ip1j) = 2 / gs^2;
            A(ij, ijp1) = 2 / gs^2;
       
        % Bottom Right Corner
        elseif (i == res && j == 1)

            % 1/h^2 (2 * ph_i-1j + 2 * phi_ij+1 - 4 * phi_ij) = 0

            % Get the Indices for the coefficiens
            ip1j = getId(i+1, j, res);
            ijp1 = getId(i, j+1, res);

            % Add coefficients to the A matrix

            A(ij, ij) = -4 / gs^2;
            A(ij, im1j) = 2 / gs^2;
            A(ij, ijp1) = 2 / gs^2;
            
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
            
            b(ij)=-in_hg;
        end
        
        % Left Tab 
        for LT = rdfive+1 : (2*rdfive)
        if (j == res && i == LT)
            
            % phi_ij = qtab
            
            % Add coefficients to the A matrix and RHS
            A(ij, ij) = 1;
            b(ij) = qtab; 
        end
        end
        
        % Right Tab 
        for RT = (3*rdfive)+1 : (4*rdfive)
        if (j == res && i == RT)
       
            % phi_ij = qtab
            
            % Add coefficients to the A matrix and RHS
            A(ij, ij) = 1;
            b(ij) = qtab;
        end
        end
        
        % Empty Top Edges
        for ETE = (2*rdfive)+1 : (3*rdfive) 
        if (j == res && i == ETE)
            % phi_ij = T_amb
            
            % Add coefficients to the A matrix and RHS
            A(ij, ij) = 1;
            b(ij) = T_amb;
            
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
[xplot, yplot] = meshgrid(linspace(0, L, res), linspace(0, L, res));
 
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
contourf(xplot, yplot, phi2,10);
colorbar
colormap(hot)
title('Solution to Laplace''s Equation','interpreter','latex','Fontsize',14)
xlabel('x','interpreter','latex','Fontsize',14)
ylabel('y','interpreter','latex','Fontsize',14)
hcb=colorbar;
title(hcb,'Temperature','interpreter','latex','Fontsize',10)
view(2)
axis equal
axis tight