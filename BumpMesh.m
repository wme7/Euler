function [vx,vy,EtoV,nE,BC] = BumpMesh()
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           Triangular Grid Generation for a Bump in a Channel
%         Katate Masatsuka, January 2013. http://www.cfdbooks.com
%           coded and modified by Manuel Diaz, NTU, 2015.05.22
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input 
%              nx = Number of nodes in x-direction.
%              ny = Number of nodes in y-direction.
%            xmin = Left end of the rectangular domain
%            xmax = Right end of the domain
%            ymin = Top end of the domain
%            ymax = Bottom end of the domain (except the bump)
%
% Output 
%           [vx,vy] = vertices of triangular elements
%            [EtoV] = connectivity matrix
%              [nE] = total number of elements built
%              [BC] = list of boundary conditions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear; clc; close all;

% 1. Input parameters
%
% ymax --------------------------------------
%      .                                    .
%      .                                    .
%      .                                    .
%      .                                    .
%      .               Bump                 .
%      .                ___                 .
% ymin ----------------/   \-----------------
%    xmin                                  xmax

nx   = 320/4;     % Number of nodes in x-direction
ny   =  40/4;     % Number of nodes in y-direction
xmin =-2.0;     % x-coordinate of the left end
xmax = 3.0;     % x-coordinate of the right end
ymin = 0.0;     % y-coordinates of the bottom
ymax = 2.0;     % y-coordinates of the top

% 2. Allocate grid coordiantes
x=zeros(ny,nx); y=zeros(ny,nx);

% 3. Build unifrom grids in x-directions
dx = (xmax-xmin)/(nx-1); % uniform spacing in x
x0 = xmin:dx:xmax;    

% 4. Build bump profile (http://http://turbmodels.larc.nasa.gov)
xb = x0(x0>0.3 & x0<1.2); 
yb = 0.05*(sin(pi*xb/0.9-(pi/3))).^4;
%plot(xb,yb); % verify!

% 5. Build domains y-values of the lower boundary
y0 = [0*x0(x0<=0.3), yb, 0*x0(x0>=1.2)];
dy = (ymax-y0)/(ny-1);

% 6. Build mesh
for j=1:ny 
    for i=1:nx
        y(j,i) = y0(i)+dy(i)*(j-1);
        x(j,i) = x0(i);
    end
end
%mesh(x,y,ones(size(x))) % verify!

% 7. Streaching nodes towards the bottom profile
stretch = 'yes';
switch stretch
    case 'yes'
        % set stretching factor
        sf = 3.5;
        % normalize y-dimension
        yn = (y-ymin)/(ymax-ymin);
        % using exponential stretching in [0,1]
        ys = (1-exp(sf*yn))/(1-exp(sf));
        % transform back to y's dimensions
        y = ys*(ymax-ymin) + ymin;
    otherwise
        % do nothing!
end
%mesh(x,y,ones(size(x))) % verify!

% 8. Perturb nodal coordinates to generate an irregular grid
for j=2:ny-1 
    for i=2:nx-1
        rn = rand; % random number
        ldx = (x(j,i+1)-x(j,i-1))/2;
        ldy = (y(j+1,i)-y(j-1,i))/2;
        y(j,i) = y(j,i)+(rn-0.5)*ldy/5;
        x(j,i) = x(j,i)-(rn-0.5)*ldx/5;
    end
end
%mesh(x,y,ones(size(x))); % verify!

% 9. Generate unstructured vertices
vx = zeros(ny*nx,1);
vy = zeros(ny*nx,1);
for i=1:nx 
    for j=1:ny
        inode = j+(i-1)*ny;
        vx(inode) = x(j,i);
        vy(inode) = y(j,i);
    end
end

% 10. Build elements (EtoV matrix)
xE=(nx-1); yE=(ny-1); nE=2*xE*yE; eNodes=3; EtoV=zeros(nE,eNodes); 
e=1; % element counter 
for i = 1:xE;
    for j = 1:yE;
        k = j+ny*(i-1); % dummy variable
        rn = rand;
        if rn < 0.5 
        EtoV( e ,:) = [k, ny+k,  ny+k+1];
        EtoV(e+1,:) = [k, ny+k+1,   k+1];
        else
        EtoV( e ,:) = [k,    ny+k,   k+1];
        EtoV(e+1,:) = [ny+k, ny+k+1, k+1];
        end
        e = e+2; % element counter
    end
end
%trimesh(EtoV,vx,vy); % verigy!

% 11. Build Boundary Node lists
BC.bottom = zeros(nx,1);
BC.right = zeros(ny,1);
BC.top = zeros(nx,1);
BC.left = zeros(ny,1);

% Boundary nodes are ordered counterclockwise 
j= 1; for i=1:nx; inode=j+(i-1)*ny; BC.bottom(i)=inode; end % Bottom
i=nx; for j=1:ny; inode=j+(i-1)*ny; BC.right(j)=inode; end  % Right 
j=ny; for i=1:nx; inode=j+(i-1)*ny; BC.top(i)=inode; end    % Top
i= 1; for j=1:ny; inode=j+(i-1)*ny; BC.left(j)=inode; end   % Left