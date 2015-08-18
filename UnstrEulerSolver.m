function [res] = UnstrEulerSolver(U,Winf,node,elem,edge,bound,gamma,fType,limitr)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computes the residual for a node-centered finite-volume method
% -------------------------------------------------------------------------
%  Input: the current solution
% Output: node(:)%res = the residual computed by the current solution.
% -------------------------------------------------------------------------
% Note: dU/dt + dF/dx + dG/dy = 0. Residuals are first computed as
%       the integral of (dF/dx + dG/dy), and at the end negative sign is
%       added so that we have dU/dt = Res at every node.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0. Load U and W states into the nodes
for i=1:nN; node(i).u=U(:,i); node(i).w=u2w(U(:,i),gamma); end

% 1/2. Eliminate the normal mass flux   
%--------------------------------------------------------------------------
% THIS IS A SPECIAL TREATMENT FOR THE SHOCK DIFFRACTION PROBLEM.
%
% NOTE: This is a corner point between the inflow boundary and
%       the lower-left wall. Enforce zero y-momentum, which is
%       not ensured by the standard BCs.
%       This special treatment is necessary because the domain
%       is rectangular (the left boundary is a straight ine) and
%       the midpoint node on the left boundary is actually a corner.
%
%       Our computational domain:
%
%                 ---------------
%          Inflow |             |
%                 |             |  o: Corner node
%          .......o             |
%            Wall |             |  This node is a corner.
%                 |             |
%                 ---------------
%
%       This is to simulate the actual domain shown below:
%      
%         -----------------------
% Inflow  |                     |
%         |                     |  o: Corner node
%         --------o             |
%            Wall |             |
%                 |             |
%                 ---------------
%      In effect, we're simulating this flow by a simplified
%      rectangular domain (easier to generate the grid).
%      So, an appropriate slip BC at the corner node needs to be applied,
%      which is "zero y-momentum", and that's all.
%
for i = 1:nbound
    if strcmp(bound(i).bcname,'slip_wall') % only slip wall BCs
        for j = 1:bound(i).nbnodes
            
            if (i==2 && j==1); % if is the corner node!
                inode = bound(i).bnode(j);             % get the node's id
                node(inode).u(3) = 0;                  % Make sure zero y-momentum.
                node(inode).w    = u2w(node(inode).u); % Update primitive variables
                % cycle bnodes_slip_wall % That's all we need. Go to the next.
                for  j = 1:bound(i).nbfaces
                    % Get Left and Right base data:
                    n1 = bound(i).bnode( j );       %Left node
                    n2 = bound(i).bnode(j+1);       %Right node
                    n12(1) = bound(i).bfnx(j);      %x-component of the unit face normal vector
                    n12(2) = bound(i).bfny(j);      %y-component of the unit face normal vector
                    mag_e12 = bound(i).bfn(j)*0.5;  %Half length of the boundary face, j.
                    
                    % 1. Left node
                    wL = node(n1).w;    wR = wL;
                    [num_flux,wsn] = fluxfunc2d(wL,wR,n12,fType);
                    bfluxL = num_flux;
                    node(n1).wsn = node(n1).wsn + wsn*mag_e12;

                    % 2. Right node
                    wL = node(n2).w;    wR = wL;
                    [num_flux,wsn] = fluxfunc2d(wL,wR,n12,fType);
                    bfluxR = num_flux;
                    node(n2).wsn = node(n2).wsn + wsn*mag_e12;
                    
                    % 3. Add contributions to the two nodes (See Nishikawa AIAA2010-5093)
                    switch elm( bound(i).belm(j) ).nV
                        case 3	% Triangle
                            node(n1).res = node(n1).res + (5*bfluxL+bfluxR)/6*mag_e12;
                            node(n2).res = node(n2).res + (5*bfluxR+bfluxL)/6*mag_e12;
                        case 4	% Quad
                            node(n1).res = node(n1).res + bfluxL.mag_e12;
                            node(n2).res = node(n2).res + bfluxR.mag_e12;
                        otherwise
                            error('Element is neither tria nor quad. Stop. ');
                    end
                    
                end
            end
            
            inode = bound(i).bnode(j);
            n12(1) = bound(i).bnx(j);
            n12(2) = bound(i).bny(j);
            
            normal_mass_flux = node(inode).u(2)*n12(1) + node(inode).u(3)*n12(2);
            
            node(inode).u(2) = node(inode).u(2) - normal_mass_flux.n12(1);
            node(inode).u(3) = node(inode).u(3) - normal_mass_flux.n12(2);
            
            node(inode).w    = u2w(node(inode).u); 
        end 
    end 
end
%--------------------------------------------------------------------------

% 1. Initialize residual and fluxes arrays in nodes
for i=1:nN; node(i).res=0; node(i).wsn=0; end 

% 2. Gradient Reconstruction
for i=1:nN; node(i).gradw=LSQgradients2d(node(i).x,node(i).y,node(i).w,...
        [node(node(i).nghbr).x],[node(node(i).nghbr).y],...
        [node(node(i).nghbr).w], node(i).invAtA ); 
end

% Flux computation across internal edges (to be accumulated in res(:))
%
%   node2              1. Extrapolate the solutions to the edge-midpoint
%       o                 from the nodes, n1 and n2.
%        \   face      2. Compute the numerical flux
%         \ -------c2  3. Add it to the residual for n1, and subtract it 
%        / \              from the residual for n2.
%   face/   \ edge
%      /     o         Directed area is the sum of the left and the right 
%    c1    node1       faces. Left/right face is defined by the edge-midpoint 
%                      and the centroid of the left/right element.
%                      Directed area is positive in n1 -> n2
%
% (c1, c2: element centroids)

% Left and right nodes of the i-th edge

    node1 = edge(i).n1;  % Left node of the edge
    node2 = edge(i).n2;  % Right node of the edge
      n12 = edge(i).dav; % This is the directed area vector (unit vector)
  mag_n12 = edge(i).da;  % Magnitude of the directed area vector
      e12 = edge(i).ev;  % This is the vector along the edge (uniti vector)
  mag_e12 = edge(i).e;   % Magnitude of the edge vector (Length of the edge)

% Solution gradient projected along the edge
%
%  NOTE: The gradient is multiplied by the distance.
%        So, it is equivalent to the solution difference.

  dwL = (node(node1).gradw(:,ix)*e12(ix) + node(node1).gradw(:,iy)*e12(iy) )*0.5*mag_e12;
  dwR = (node(node2).gradw(:,ix)*e12(ix) + node(node2).gradw(:,iy)*e12(iy) )*0.5*mag_e12;

%  It is now limiter time%

%  (1) No limiter (good for smooth solutions)

   limiter : if (trim(limiter_type) == "none") then

%      Simple linear extrapolation
       wL = node(node1).w + dwL;
       wR = node(node2).w - dwR;

%  (2) UMUSCL-type limiters: simple 1D limiting.

   elseif (trim(limiter_type) == "vanalbada") then

%--------------------------------------------------------------------------
% In 1D: dwp = w_{j+1}-w_j, dwm = w_j-w_{j-1} 
%               => limited_slope = limiter(dwm,dwp)
%
% We can do the same in 2D as follows.
%
% In 2D:    dwp = w_{neighbor}-w_j, dwm = 2*(grad_w_j*edge)-dwp
%               => limited_slope = limiter(dwm,dwp)
%
% NOTE: On a regular grid, grad_w_j*edge will be the central-difference,
%       so that the average (dwm+dwp)/2 will be the central-difference just
%       like in 1D. 
%--------------------------------------------------------------------------

% Edge derivative
dwij = 0.5*(node(node2).w - node(node1).w);

% Left face value (wL) with the Van Albada limiter
dwm  = 2*dwL-dwij;
dwp  = dwij;
wL  = node(node1).w + vanAlbada(dwm,dwp,mag_e12);

% Right face value (wR) with the Van Albada limiter
dwm  = -(2*dwR-dwij);
dwp  = -dwij;
wR  = node(node2).w + vanAlbada(dwm,dwp,mag_e12);

% Left and right nodes of the i-th edge
    node1 = edge(i).n1;  % Left node of the edge
    node2 = edge(i).n2;  % Right node of the edge
      n12 = edge(i).dav; % This is the directed area vector (unit vector)
  mag_n12 = edge(i).da;  % Magnitude of the directed area vector
      e12 = edge(i).ev;  % This is the vector along the edge (uniti vector)
  mag_e12 = edge(i).e;   % Magnitude of the edge vector (Length of the edge)

%--------------------------------------------------------------------------
% Solution gradient projected along the edge
%
%  NOTE: The gradient is multiplied by the distance.
%        So, it is equivalent to the solution difference.
%--------------------------------------------------------------------------

  dwL = (node(node1).gradw(:,ix)*e12(ix) + node(node1).gradw(:,iy)*e12(iy) )*0.5*mag_e12;
  dwR = (node(node2).gradw(:,ix)*e12(ix) + node(node2).gradw(:,iy)*e12(iy) )*0.5*mag_e12;

%  It is now limiter time%

%  (1) No limiter (good for smooth solutions)

   limiter : if (trim(limiter_type) == "none") then

%      Simple linear extrapolation
       wL = node(node1)%w + dwL
       wR = node(node2)%w - dwR

%  (2) UMUSCL-type limiters: simple 1D limiting.

   elseif (trim(limiter_type) == "vanalbada") then

%--------------------------------------------------------------------------
% In 1D: dwp = w_{j+1}-w_j, dwm = w_j-w_{j-1} 
%               => limited_slope=limiter(dwm,dwp)
%
% We can do the same in 2D as follows.
%
% In 2D: dwp = w_{neighbor}-w_j, dwm = 2*(grad_w_j*edge)-dwp
%               => limited_slope = limiter(dwm,dwp)
%
% NOTE: On a regular grid, grad_w_j*edge will be the central-difference,
%       so that the average (dwm+dwp)/2 will be the central-difference just
%       like in 1D. 
%--------------------------------------------------------------------------

%     Edge derivative
      dwij = 0.5*(node(node2).w - node(node1).w);

%     Left face value (wL) with the Van Albada limiter
      dwm  = 2*dwL-dwij;
      dwp  = dwij;
       wL  = node(node1).w + vanAlbada(dwm,dwp,mag_e12);

%     Right face value (wR) with the Van Albada limiter
      dwm  = -(2*dwR-dwij);
      dwp  = -dwij;
       wR  = node(node2).w + vanAlbada(dwm,dwp,mag_e12);

%-------------------------------------------------------------------------
% Close with the boundary flux using the element-based formula that is
% exact for linear fluxes (See Nishikawa AIAA2010-5093 for boundary weights
% that ensure the linear exactness for 2D/3D elements).
%
%      |  Interior Domain          |
%      |        .........          |
%      |        .       .          |
%      |        .       .          |
%      o--o--o-----o---------o--o--o  <- Boundary segment
%                  n1   |   n2
%                       v
%                     n12 (unit face normal vector)
%
% NOTE: We visit each boundary face, defined by the nodes n1 and n2,
%       and compute the flux across the boundary face: left half for node1,
%       and the right half for node2. In the above figure, the dots indicate
%       the control volume around the node n1. Clearly, the flux across the
%       left half of the face contributes to the node n1. Similarly for n2.
%
%--------------------------------------------------------------------------
%  BC: Upwind flux via freestream values
%
%      NOTE: If the final solution at the boundary node is far from
%            the freestream values, then the domain is probably is not
%            large enough. 
%--------------------------------------------------------------------------
%  BC: Solid body and Supersonic outflow
%
%      NOTE: Basically, simply compute the physical flux, which
%            can be done by calling the Roe flux with wR = wL.
%            It is equivalent to the interior-extrapolation condition.
%      NOTE: Tangency condition for solid body will be applied later.
%--------------------------------------------------------------------------
%  BC: Subsonic Outflow - Fixed Back Pressure
%
%      NOTE: Fix the pressure as freestream pressure
%            on the right side of the face (outside the domain).
%            Assumption is that the outflow boundary is far from the body.
%--------------------------------------------------------------------------

for i = 1:nbound
    for  j = 1:bound(i).nbfaces
        % Get Left and Right base data:
        n1 = bound(i).bnode( j );       %Left node
        n2 = bound(i).bnode(j+1);       %Right node
        n12(1) = bound(i).bfnx(j);      %x-component of the unit face normal vector
        n12(2) = bound(i).bfny(j);      %y-component of the unit face normal vector
        mag_e12 = bound(i).bfn(j)*0.5;  %Half length of the boundary face, j.
        switch bound(i).bcname
            case 'freestream' % bnodes_numerical_flux_via_freestream
                %   1. Left node
                wL = node(n1).w;    wR = Winf;
                [num_flux,wsn] = fluxfunc2d(wL,wR,n12,fType);
                bfluxL = num_flux;
                node(n1).wsn = node(n1).wsn + wsn*mag_e12;
                
                %   2. Right node
                wL = node(n2).w;    wR = Winf;
                [num_flux,wsn] = fluxfunc2d(wL,wR,n12,fType);
                bfluxR = num_flux;
                node(n2).wsn = node(n2).wsn + wsn*mag_e12;
                
            case {'slip_wall','outflow_supersonic'} % bnodes_slip_wall
                %   1. Left node
                wL = node(n1).w;    wR = wL;
                [num_flux,wsn] = fluxfunc2d(wL,wR,n12,fType);
                bfluxL = num_flux;
                node(n1).wsn = node(n1).wsn + wsn*mag_e12;
                
                %   2. Right node
                wL = node(n2).w;    wR = wL;
                [num_flux,wsn] = fluxfunc2d(wL,wR,n12,fType);
                bfluxR = num_flux;
                node(n2).wsn = node(n2).wsn + wsn*mag_e12;
                
            case 'outflow_back_pressure' % bnodes_outflow
                %   1. Left node
                wL = node(n1).w;    wR = wL;    wR(4) = Winf(4); %Fix the pressure
                [num_flux,wsn] = fluxfunc2d(wL,wR,n12,fType);
                bfluxL = num_flux;
                node(n1).wsn = node(n1).wsn + wsn*mag_e12;
                
                %   2. Right node
                wL = node(n2).w;    wR = wL;    wR(4) = Winf(4); %Fix the pressure
                [num_flux,wsn] = fluxfunc2d(wL,wR,n12,fType);
                bfluxR = num_flux;
                node(n2).wsn = node(n2).wsn + wsn*mag_e12;
        end
        %   3. Add contributions to the two nodes (See Nishikawa AIAA2010-5093)
        switch elm( bound(i).belm(j) ).nV
            case 3	% Triangle
                node(n1).res = node(n1).res + (5*bfluxL+bfluxR)/6*mag_e12;
                node(n2).res = node(n2).res + (5*bfluxR+bfluxL)/6*mag_e12;
            case 4	% Quad
                node(n1).res = node(n1).res + bfluxL.mag_e12;
                node(n2).res = node(n2).res + bfluxR.mag_e12;
            otherwise
                error('Element is neither tria nor quad. Stop. ');
        end
    end
end