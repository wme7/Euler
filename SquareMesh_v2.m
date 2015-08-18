function [vx,vy,EtoV,nE,nN,BC] = SquareMesh_v2(Etype)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           2D quad and triangular grids in a rectangular domain.
%          Katate Masatsuka, January 2012. http://www.cfdbooks.com
%           coded and modified by Manuel Diaz, NTU, 2015.05.22
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary information is set up for a channel flow problem
%
%                                Wall
%                         --------------------
%                         .                  .
%                         .                  .          
%                         .                  .
%                  Inflow .  -->             . Outflow
%                         .                  .
%                         .                  .
%                         .                  .
%                         --------------------
%                                Wall
%
%--------------------------------------------------------------------------
% Inuput: 
%        xmin, xmax = x-coordinates of the left and right ends
%        ymin, ymax = y-coordinates of the bottom and top ends
%                nx = number of nodes in x-direction
%                ny = number of nodes in y-direction
% Output:
%           [vx,vy] = vertices of triangular elements
%            [EtoV] = connectivity matrix
%              [nE] = total number of elements built
%              [BC] = list of boundary conditions
%
%**************************************************************************

%clear; clc; close all;

% 1. Define the grid size and allocate the structured grid data array.
%
% ymax --------------------------------------
%      .                                    .
%      .                                    .
%      .                                    .
%      .                                    .
%      .                                    .
%      .                                    .
%      .                                    .
%      .                                    .
%      .                                    .
%      .                                    .
%      .                                    .
%      .                                    .
%      .                                    .
%      .                                    .
%      .                                    .
% ymin --------------------------------------
%    xmin                                  xmax

nx   = 11;     % Number of nodes in x-direction
ny   = 11;     % Number of nodes in y-direction
xmin = 0;       % x-coordinate of the left end
xmax = 1;       % x-coordinate of the right end
ymin = 0;       % y-coordinates of the bottom
ymax = 1;       % y-coordinates of the top
%Etype = 'QUAD';  % TRI or QUAD
nN   = nx*ny;   % total number of nodesl

% 2. Generate a structured 2D grid data, (i,j) data: go up in y-direction%
%
% j=5 o--------o--------o--------o--------o
%     |        |        |        |        |
%     |        |        |        |        |   On the left is an example:
%     |        |        |        |        |         nx = 5
% j=4 o--------o--------o--------o--------o         ny = 5
%     |        |        |        |        |
%     |        |        |        |        |
%     |        |        |        |        |
% j=3 o--------o--------o--------o--------o
%     |        |        |        |        |
%     |        |        |        |        |
%     |        |        |        |        |
% j=2 o--------o--------o--------o--------o
%     |        |        |        |        |
%     |        |        |        |        |
%     |        |        |        |        |
% j=1 o--------o--------o--------o--------o
%     i=1      i=2      i=3      i=4      i=5

dx = (xmax-xmin)/(nx-1); % uniform spacing in x
dy = (ymax-ymin)/(ny-1); % uniform spacing in y
[x,y] = meshgrid(xmin:dx:xmax,ymin:dy:ymax);    

% 3. Generate unstructured data: elements vertices coordinates
%
%    - The so-called lexcographic ordering -
%
%    5       10       15       20       25
%     o--------o--------o--------o--------o
%     |        |        |        |        |
%     |        |        |        |        |   On the left is an example:
%    4|       9|      14|      19|      24|           nx = 5
%     o--------o--------o--------o--------o           ny = 5
%     |        |        |        |        |
%     |        |        |        |        |   nnodes = 5x5 = 25
%    3|       8|      13|      18|      23|
%     o--------o--------o--------o--------o
%     |        |        |        |        |
%     |        |        |        |        |
%    2|       7|      12|      17|      22|
%     o--------o--------o--------o--------o
%     |        |        |        |        |
%     |        |        |        |        |
%    1|       6|      11|      16|      21|
%     o--------o--------o--------o--------o
vx = zeros(nN,1);
vy = zeros(nN,1);
for i=1:nx 
    for j=1:ny
        inode = j+(i-1)*ny;
        vx(inode) = x(j,i);
        vy(inode) = y(j,i);
    end
end

% 4. Generate element connectivy matrix (EtoV matrix)
disp(Etype);
switch Etype
    case 'QUAD'
        xE=(nx-1); yE=(ny-1); nE=xE*yE; EtoV=cell(nE,1);
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
                EtoV{ e } = [k, ny+k,  ny+k+1, k+1];
                e = e+1; % element counter
            end
        end
    case 'TRI'
        xE=(nx-1); yE=(ny-1); nE=2*xE*yE; EtoV=cell(nE,1);
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
                    EtoV{ e } = [k, ny+k,  ny+k+1];
                    EtoV{e+1} = [k, ny+k+1,   k+1];
                else
                    EtoV{ e } = [k,    ny+k,   k+1];
                    EtoV{e+1} = [ny+k, ny+k+1, k+1];
                end
                e = e+2; % element counter
            end
        end
	case 'MIXED'
        xE=(nx-1); yE=(ny-1); nE=xE*yE; EtoV=cell(nE,1);
        %
        midxE = (xE+1)/2; % xE is assumed odd number
        midyE = (yE+1)/2; % xE is assumed odd number
        %
        e=1; % element counter
        for i = 1:xE;
            for j = 1:yE;
                k = j+ny*(i-1); % dummy variable
                if ((i<midxE) && (j<midyE)) || ((i>midxE) && (j>midyE))
                    rn = rand;
                    if rn < 0.5
                        EtoV{ e } = [k, ny+k,  ny+k+1];
                        EtoV{e+1} = [k, ny+k+1,   k+1];
                    else
                        EtoV{ e } = [k,    ny+k,   k+1];
                        EtoV{e+1} = [ny+k, ny+k+1, k+1];
                    end
                    e = e+2; % element counter    
                else
                   EtoV{ e } = [k, ny+k,  ny+k+1, k+1];
                    e = e+1; % element counter 
                end
            end
        end 
        nE = numel(EtoV); % true number of elements!
    otherwise
        error('Wrong type of element')
end

% 5. Build Boundary data
%
% NOTE: These boundary data are specific to the shock diffraction problem.
%
%  Example: nx=ny=7
%
%   in = inflow
%    w = wall
%    e = outflow
%    o = interior nodes
%  inw = this node belongs to both inflow and wall boundaries.
%   we = this node belongs to both wall and outflow boundaries.
%
%   inw----w----w----w----w----w----we
%     |    |    |    |    |    |    |
%    in----o----o----o----o----o----e
%     |    |    |    |    |    |    |
%    in----o----o----o----o----o----e
%     |    |    |    |    |    |    |
%    in----o----o----o----o----o----e
%     |    |    |    |    |    |    |
%    in----o----o----o----o----o----e
%     |    |    |    |    |    |    |
%    in----o----o----o----o----o----e
%     |    |    |    |    |    |    |
%   inw----w----w----w----w----w----we
%
% allocate
BC.inflow     = zeros(ny,1);
BC.bottomwall = zeros(nx,1);
BC.outflow    = zeros(ny,1);
BC.topwall    = zeros(nx,1);

% Boundary nodes are ordered counterclockwise 
i= 1; for j=1:ny; inode=j+(i-1)*ny; BC.inflow(j)=inode; end     % inflow
j= 1; for i=1:nx; inode=j+(i-1)*ny; BC.bottomwall(i)=inode; end    % bottom outflow
i=nx; for j=1:ny; inode=j+(i-1)*ny; BC.outflow(j)=inode; end	% right outflow
j=ny; for i=1:nx; inode=j+(i-1)*ny; BC.topwall(i)=inode; end    % top wall

%% visual verification
% Plot Element as patches
figure(1); hold on;
for e = 1:nE
    patch(vx(EtoV{e})',vy(EtoV{e})',zeros(size(EtoV{e}')),'w');
end
%cell2table(EtoV) % cannot display all info :(

% Plot nodes and elements numbers
for i = 1:nN
    text(vx(i),vy(i),int2str(i),'fontsize',8,....
        'fontweight','bold','Color','r');
end
for e = 1:nE
    pos = [mean(vx(EtoV{e})),mean(vy(EtoV{e}))];
    text(pos(1),pos(2),int2str(e),'fontsize',8,...
        'fontweight','bold','color','b');
end