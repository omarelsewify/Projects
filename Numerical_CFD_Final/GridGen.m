function [ grid ] = GridGen( nx, ny )
%
% Grid generator used for the take-home final exam of AA214B
% Stanford University, Winter Quarter, 2015
% Written by Benjamin Blake
%
% Input parameters:
%   nx - the number of cells in the x direction
%   ny - the number of cells in the y direction
%
% Output parameters:
%   grid       - data object containing all necessary grid data
%   grid.nx    - number of cells in x
%   grid.ny    - number of cells in y
%   grid.x     - 2D array [nx+1,ny+1], contains x node locations of grid
%   grid.y     - 2D array [nx+1,ny+1], contains y node locations of grid
%   grid.xc    - 2D array [nx,ny], contains x cell center locations of grid
%   grid.yc    - 2D array [nx,ny], contains y cell center locations of grid
%   grid.area  - 2D array [nx,ny], contains areas of all cells
%   grid.edge  - 3D array [nx,ny,4], contains edge lengths of all cells
%   grid.normx - 3D array [nx,ny,4], contains x-component of outward normals of edges of all cells
%   grid.normy - 3D array [nx,ny,4], contains y-component of outward normals of edges of all cells
%

%% Set the number of cells
grid.nx = nx;
grid.ny = ny;

%% Construct uniform grid
x = linspace(0,4,nx+1);
y = linspace(0,1.5,ny+1);
y_ffb = 1.5;
grid.x = zeros(nx+1,ny+1);
grid.y = zeros(nx+1,ny+1);
for i=1:nx+1
    for j=1:ny+1
        grid.x(i,j) = x(i);
        grid.y(i,j) = y(j);
    end
end

%% Deform grid to account for ramp
for i=1:(nx+1)
    if (grid.x(i,1) >= 0.5)
        grid.y(i,1) = grid.y(i,1) + (grid.x(i,1) - 0.5) * tan(8*pi/180);
    end
end
for i=1:nx+1
    dy = grid.y(i,1);
    for j=2:ny+1
        ratio = 1 - (grid.y(i,j)) / y_ffb;
        deformation = dy * ratio;
        grid.y(i,j) = grid.y(i,j) + deformation;
    end
end

%% Generate array of cell-centered locations
grid.xc = zeros(nx,ny);
grid.yc = zeros(nx,ny);
for i=1:nx
    for j=1:ny
        grid.xc(i,j) = 0.25 * ( grid.x(i,j) + grid.x(i+1,j) + grid.x(i+1,j+1) + grid.x(i,j+1) );
        grid.yc(i,j) = 0.25 * ( grid.y(i,j) + grid.y(i+1,j) + grid.y(i+1,j+1) + grid.y(i,j+1) );
    end
end

%% Generate array of areas of all cells
grid.area = zeros(nx,ny);
for i=1:nx
    for j=1:ny
        grid.area(i,j) = 0.5 * abs( ...
            grid.x(i+1,j) * grid.y(i+1,j+1) + grid.x(i+1,j+1) * grid.y(i,j+1) ...
          + grid.x(i,j+1) * grid.y(i,j)     + grid.x(i,j) * grid.y(i+1,j) ...
          - grid.x(i+1,j+1) * grid.y(i+1,j) - grid.x(i,j+1) * grid.y(i+1,j+1) ...
          - grid.x(i,j) * grid.y(i,j+1)     - grid.x(i+1,j) * grid.y(i,j) );
    end
end

%% Generate array of cell edge lengths
grid.edge = zeros(nx,ny,4);
for i=1:nx
    for j=1:ny
        grid.edge(i,j,1) = sqrt( (grid.x(i+1,j+1) - grid.x(i+1,j))^2 ...
                               + (grid.y(i+1,j+1) - grid.y(i+1,j))^2 );
        grid.edge(i,j,2) = sqrt( (grid.x(i,j+1) - grid.x(i+1,j+1))^2 ...
                               + (grid.y(i,j+1) - grid.y(i+1,j+1))^2 );
        grid.edge(i,j,3) = sqrt( (grid.x(i,j) - grid.x(i,j+1))^2 ...
                               + (grid.y(i,j) - grid.y(i,j+1))^2 );
        grid.edge(i,j,4) = sqrt( (grid.x(i+1,j) - grid.x(i,j))^2 ...
                               + (grid.y(i+1,j) - grid.y(i,j))^2 );
    end
end

%% Generate array of outward normal vectors for all edges of every cell
grid.normx = zeros(nx,ny,4);
grid.normy = zeros(nx,ny,4);
rot = [0,1;-1,0];
for i=1:nx
    for j=1:ny
        v1 = [grid.x(i+1,j+1) - grid.x(i+1,j); ...
              grid.y(i+1,j+1) - grid.y(i+1,j)];
        v2 = [grid.x(i,j+1) - grid.x(i+1,j+1); ...
              grid.y(i,j+1) - grid.y(i+1,j+1)];
        v3 = [grid.x(i,j) - grid.x(i,j+1); ...
              grid.y(i,j) - grid.y(i,j+1)];
        v4 = [grid.x(i+1,j) - grid.x(i,j); ...
              grid.y(i+1,j) - grid.y(i,j)];
        
        v1 = rot * v1 / norm(v1);
        grid.normx(i,j,1) = v1(1);
        grid.normy(i,j,1) = v1(2);
        
        v2 = rot * v2 / norm(v2);
        grid.normx(i,j,2) = v2(1);
        grid.normy(i,j,2) = v2(2);
        
        v3 = rot*v3/norm(v3);
        grid.normx(i,j,3) = v3(1);
        grid.normy(i,j,3) = v3(2);
        
        v4 = rot*v4/norm(v4);
        grid.normx(i,j,4) = v4(1);
        grid.normy(i,j,4) = v4(2);
    end
end

end