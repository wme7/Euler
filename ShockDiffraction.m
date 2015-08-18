%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%         Shock-diffraction Problem using Nodal Centered FV method
%         Adapted from Dr. Katate Masatsuka (info[at]cfdbooks.com),
%               coded by Manuel A. Diaz, NTU, 2015.05.22
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                 Wall
%                         --------------------
%     Post-shock (inflow) |                  |
%                         |->Shock           |            o: Corner node
%                         |  Mach=5.09       |
%                  .......o                  | Outflow
%                    Wall |                  |
%                         |                  |
%                         |                  |
%                         --------------------
%                               Outflow
% Solver features:
% * Node-centered FV method for unstructured grids (quad/tri/mixed) used,
% * Roe flux with an entropy fix and Rotated-RHLL flux,
% * Gradient reconstruction by unweighted least-squares method,
% * Van Albada slope limiter to the primitive variable gradients,
% * 2-Stage Runge-Kutta global time-stepping towards the final time,
% * All quantities are nondimensionalized; velocity and pressure are
%   nondimensionalized based on the free stream speed of sound.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

% 1. Control Parameters
%M_inf  = 0;	% free stream Mach number (set in the IC)
n      = 5;     % Total DoF of gas 
cfl    = 0.95;  % CFL condition
tEnd   = 0.18;  % output time
flux   ='ROE';  % ROE, LLF
limitr ='MC';   % MC, MM, VA
eType  = 2;     % 1:TRI or 2:QUAD
gamma  = (n+2)/n;

% 2. Build grid data
switch eType
    case 1; % use Triangles
        [vx,vy,EtoV,nE,nN,BC] = SquareMesh('TRI');
    case 2; % use Quadrilaterals
        [vx,vy,EtoV,nE,nN,BC] = SquareMesh('QUAD');
	case 3; % use a Voronoi grid 
        [vx,vy,EtoV,nE,nN,BC] = SquareMesh('VORONOI'); % not yet!
end
[node,elem,edge,bound]= BuildUnstructuredMesh2d(vx,vy,EtoV,nE,nN,BC);

% 3. Test integrity of the mesh
verifyMesh(node,elem,edge,bound,vx,vy,EtoV,nE,nN)

% 4. Set initial condition
[w0,w_inf,M_inf] = ShockDiffractionIC(); %W0=repmat(w0,[1,nN]);
U0=repmat(w2u(w0,gamma),[1,nN]); Winf=w_inf; 

% 5. Solver Loop
%------------------------------------------------------
% Two-stage Runge-Kutta scheme: u^n is saved as u0(:,:)
%  1. u^*     = u^n - (dt/vol)*Res(u^n)
%  2. u^{n+1} = 1/2*u^n + 1/2*[u^* - (dt/vol)*Res(u^*)]
%------------------------------------------------------

% set initial time step
dt0 = 0.01;

% Load initial conditions
U=U0; t=0; it=0; dt=dt0; 

while t < tEnd
    % RK step 1
    Us = U - dt*UnstrEulerSolver(U,Winf,node,elem,edge,bound,gamma,flux,limitr);
    
    % RK step 2
    Un = 0.5*(U + (Us - dt*UnstrEulerSolver(Us,Winf,node,elem,edge,bound,gamma,flux,limitr)));
    
    % Update solution
    U = Un;
    
    % update time step (local and global)
    dt = dt0;
    t=t+dt; it=it+1;
    
    % Plot solution
    if rem(it,10)==0
        % plot solution
    end
end

% 6. Post-process
