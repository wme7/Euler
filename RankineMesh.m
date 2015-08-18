function [vx,vy,EtoV,nE,nN,BC] = RankineMesh(Etype)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           2D quad and triangular grids in a rectangular domain.
%          Katate Masatsuka, January 2012. http://www.cfdbooks.com
%           coded and modified by Manuel Diaz, NTU, 2015.05.25
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program generates 2D quad and triangular grids for a flow over
% Rankine's half body.
%                              Free Stream
%                 --------------------------------------
%                 .                                    . Outflow
%                 .                                 .  .
%                 .                           .      
%                 .                     .     
%                 .                  .
%                 .  inflow         .   (Wall)
%                 .  --->           .   Rankine's half body
%                 .                  .   
%                 .                     .
%                 .                           . 
%                 .                                 .  .
%                 .                                    . Outflow
%                 --------------------------------------
%                               Free Stream
% 3-Step Generation:
%
% 1. Generate a temporary structured grid data for nodes: xs(i,j) and ys(i,j)
% 2. Generate a 1D node array: x(1:nnodes), y(1:nnodes)
% 3. Generate element connectivity data: tria(1:ntria,3), quad(1:nquad,4)
%
%  Input: 
%     nxp = number of nodes over the half body
%      ny = number of nodes in the direction from the body to outer boundary
%
%        (NOTE: All input parameters are defined inside the program.)
%
% Output:
%     tria_rankine_grid = tecplot file of the triangular grid (with exact sol)
%     quad_rankine_grid = tecplot file of the quadrilateral grid (with exact sol)
%         project.bcmap = file that contains boundary condition info
%
%   Exact solution: We set sigma = Vinf, so that
%
%       Velocity potential: phi/Vinf = x + ln(x^2+y^2)/(4*pi)
%          Stream function: psi/Vinf = y + theta/(2*pi)
%               X-velocity:   u/Vinf = 1 + x/(2*pi)/(x^2+y^2)
%               Y-velocity:   u/Vinf =     y/(2*pi)/(x^2+y^2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear; close all; clc;

% 0. Define the grid size and allocate the structured grid data array.
%
%                             B5
%           ymax --------------------------------------  ^
%                .                                 B6 .  | ny nodes here
%                .                                    .  v 
%                .                              .     . theta=theta_min
%                .                       .     
%                .                   .
%                .                 .
%              B4.             B1  .  <-- Origin (0,0)    
%                .                  .   
%                .                     .
%                .                              .     . theta=theta_max
%                .                                    .  ^
%                .             B3                  B2 .  | ny nodes here
%           ymin --------------------------------------  v
%              xmin                                  xmax
%
% Boundaries are considered separately: 
% for the body(B1), lower outflow(B2), bottom outer boundary(B3), 
% left inflow boundary(B4), top outer boundary(B5) and top outflow(B6).

% Define the domain: here we define the coordinates that define the corners.
nnB3 = 20;	% Number of nodes The number of nodes along B3
nnB4 = 20;	% Number of nodes The number of nodes along B4
nnB5 = 20;	% Number of nodes The number of nodes along B5
 ny  = nnB3+nnB4+nnB5-2; % The number of nodes over the half-body(B1)
 nx  = 21;	% Number of nodes from B3,B4,B5 to B1
xmin =-1.0;	% x-coordinate of the left end
xmax = 1.0;	% x-coordinate of the right end
ymin =-1.0;	% y-coordinates of the bottom
ymax = 1.0;	% y-coordinates of the top
irregular = true;

% Allocate all nodes
xs = zeros(ny,nx);
ys = zeros(ny,nx);

%% 1. Generate nodes along the outer boundary (B3,B4,B5): i=1, j=1:nx
xB3=linspace(xmax,xmin,nnB3);  yB3=-ones(size(xB3));    % Bottom
yB4=linspace(ymin,ymax,nnB4);  xB4=-ones(size(yB4));    % Left
xB5=linspace(ymin,xmax,nnB5);  yB5= ones(size(xB5));    % Top
xs(:,1)=[xB3,xB4(2:nx-2),xB5]; ys(:,1)=[yB3,yB4(2:nx-2),yB5];

%% 2. Generate nodes along the Rankine half body (B1): i=nx, j=1:nx

% First, find the starting value of theta, which gives x=1.0
theta_min = 0;
for i = 1:100
    % Equation of theta, representing x=1.0.
    % This is basically streamfunction=sigma/2 with x=1.0.
    f = tan(theta_min) - 0.5*(1 - theta_min/pi);
    
    % If the equation is satisfied, we're done.
    if abs(f)<1E-18; disp(theta_min); break; end
    
    % Derivative (Jacobian) of the equation we are solving.
    df = 1/cos(theta_min)^2 + 0.5*pi;
    
    % Newton iteration:
    theta_min = theta_min - f/df;
end
% By symmetry, the maximum value of theta is just 2*pi-theta_min.
theta_max = 2*pi - theta_min;

% Build Half body
for i = 1:ny
    % i=1:nx corresponds to theta=[ theta_max,theta_min ]
    dtheta = (theta_max-theta_min)/(ny-1);
    theta = theta_max - dtheta*(i-1);
    
    % Apply some disturbance to 
    if irregular; rn=rand; theta=theta+dtheta*(rn-0.5); end
    
    % Apply stretching ibn theta to distribute nodes as uniformly as possible.
    sf = 4.0;  s = ( theta-theta_min ) / (theta_max-theta_min);	  % Transform in space s=[0,1]
    z1 = -tanh(sf*(-0.5));
    snew = 0.5*(z1+tanh(sf*(s-0.5))) / ( 0.5*(z1+tanh(sf*0.5)) ); % Apply stretching in s
    theta = (theta_max-theta_min)*snew + theta_min;               % Transform back in theta
    
    % Equation (7.11.19) on page 225
    ys(i,nx) = 0.5*(1-theta/pi);
    
    % From the definition: tan(theta)=y/x: 
    % Careful here, we choose nnB(even) such that we avoid theta=pi.
    xs(i,nx) = ys(i,nx)/tan(theta);
end

%% 3. Generate a structured 2D grid data, (i,j) data
%    (ex,ey) is the unit vector pointing from a node on (B3,B4,B5) to the
%    corresponding node at the outer boundary B1.
ex=xs(:,nx)-xs(:,1); ey=ys(:,nx)-ys(:,1);

% distance is the distance between the two nodes.
dist=sqrt(ex.^2+ey.^2); ex=ex./dist; ey=ey./dist; %<-- Unit vectors

% dr is the uniform spacing between the two nodes.
dr = dist/(nx-1);

% Generate interior nodes uniformly from the body(B1) to the outer(B3,B4,B5)
for j=2:nx-1; xs(:,j)=xs(:,1)+dr.*ex*(j-1); ys(:,j)=ys(:,1)+dr.*ey*(j-1); end

%% 4. Perturb nodal coordinates to generate an irregular grid
if irregular
    for j=2:ny-1
        for i=2:nx-1
            rn = rand; % random number
            ldx = (xs(j,i+1)-xs(j,i-1))/2;
            ldy = (ys(j+1,i)-ys(j-1,i))/2;
            ys(j,i) = ys(j,i)+(rn-0.5)*ldy/2;
            xs(j,i) = xs(j,i)-(rn-0.5)*ldx/2;
        end
    end
end
% verify!
figure; axis([-1,1,-1,1]); scatter(xs(:),ys(:));

%% 5. Generate unstructured data: 1D array to store the node information.
%  Total number of nodes
nN=nx*ny; vx=zeros(nN,1); vy=zeros(nN,1); bmark=zeros(ny,nx);
phi=zeros(ny,nx); psi=zeros(ny,nx); u=zeros(ny,nx); v=zeros(ny,nx);

% Node data: the nodes are ordered in 1D array.
for i = 1:ny   %Go up in y-direction.
    for j = 1:nx  %Go to the right in x-direction.
        inode = i+(j-1)*ny; % <- Node number in the lexcographic ordering
        vx(inode) = xs(i,j); 
        vy(inode) = ys(i,j);
        %Create boundary mark: =0 for interior nodes, =1 for boundary nodes
        bmark(inode) = 0;
        if (i== 1); bmark(inode) = 1; end
        if (i==ny); bmark(inode) = 1; end
        if (j== 1); bmark(inode) = 1; end
        if (j==nx); bmark(inode) = 1; end
        %Compute the exact solution
        distance = sqrt( vx(inode)^2 + vy(inode)^2 );
        if vy(inode)<0
            theta = 2*pi - acos( vx(inode)/distance );
        else
            theta = acos( vx(inode)/distance );
        end
        phi(inode) = vx(inode) + 1/(4*pi)*log(distance^2);
        psi(inode) = vy(inode) + theta/(2*pi);
        u(inode) = 1 + vx(inode)/(2*pi)/distance^2;
        v(inode) =     vy(inode)/(2*pi)/distance^2;
    end
end

% Verify exact solution
%subplot(221); surf(xs,ys,u);    subplot(222); surf(xs,ys,v)
%subplot(223); surf(xs,ys,phi);  subplot(224); surf(xs,ys,psi)

%% 6. Generate element connectivy matrix (EtoV matrix)
disp(Etype);
switch Etype
    case 'QUAD'
        xE=(nx-1); yE=(ny-1); nE=xE*yE; eNodes=4; EtoV=zeros(nE,eNodes);
        %
        %  inode+1    inode+ny+1    i4      i3
        %       o-------o           o-------o
        %       |       |           |       |
        %       |       |     or    |       |
        %       |       |           |       |
        %       o-------o           o-------o
        %     inode   inode+ny      i1      i2
        %
        % Quad is defined by the counterclockwise ordering of nodes.

        e=1; % element counter
        for i = 1:xE;
            for j = 1:yE;
                k = j+ny*(i-1); % dummy variable
                EtoV( e ,:) = [k, ny+k,  ny+k+1, k+1];
                e = e+1; % element counter
            end
        end
    case 'TRI'
        xE=(nx-1); yE=(ny-1); nE=2*xE*yE; eNodes=3; EtoV=zeros(nE,eNodes);
        % Trianguler grid with right-up diagonals (i.e., / ).
        %
        %  inode+1    inode+ny+1    i4      i3
        %       o-------o           o-------o
        %       |     . |           |     . |
        %       |   .   |     or    |   .   |
        %       | .     |           | .     |
        %       o-------o           o-------o
        %    inode    inode+ny      i1      i2
        %
        % Triangle is defined by the counterclockwise ordering of nodes.
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
    otherwise
        error('Wrong type of element')
end

%% 6. Build Boundary data

% Allocate
% BC.B1 = zeros(ny,1);
% BC.B2 = zeros(nx,1);
% BC.B3 = zeros(nnB3,1);
% BC.B4 = zeros(nnB4,1);
% BC.B5 = zeros(nnB5,1);
% BC.B6 = zeros(nx,1);

% Boundary nodes are ordered counterclockwise 
% i=1;  for j=1:ny; inode=j+(i-1)*ny; BC.B1(j)=inode; end % B1
% i=1;  for j=1:ny; inode=j+(i-1)*ny; BC.B2(j)=inode; end % B2
% i=1;  for j=1:ny; inode=j+(i-1)*ny; BC.B3(j)=inode; end % B3
% i=1;  for j=1:ny; inode=j+(i-1)*ny; BC.B4(j)=inode; end % B4
% i=1;  for j=1:ny; inode=j+(i-1)*ny; BC.B5(j)=inode; end % B5
% i=1;  for j=1:ny; inode=j+(i-1)*ny; BC.B6(j)=inode; end % B6

%% visual verification
% Plot Element as patches
figure(1); hold on;
for e = 1:nE
    patch(vx(EtoV(e,:))',vy(EtoV(e,:))',zeros(size(EtoV(e,:)')),'w');
end
%cell2table(EtoV) % cannot display all info :(

% Plot nodes and elements numbers
% for i = 1:nN
%     text(vx(i),vy(i),int2str(i),'fontsize',8,....
%         'fontweight','bold','Color','r');
% end
% for e = 1:nE
%     pos = [mean(vx(EtoV(e,:))),mean(vy(EtoV(e,:)))];
%     text(pos(1),pos(2),int2str(e),'fontsize',8,...
%         'fontweight','bold','color','b');
% end