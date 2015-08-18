%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Roe's Approximate 1-D Riemann solver for Euler equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Following the ideas of:
% 1. E.F. Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics 
% Manchester U.K., Springer Editorial, 2nd Ed. Chapert 11.
%
% This code solves the Sod's shock tube problem 
%
% t=0                                 t=tEnd
% Density                             Density
%   ****************|                 *********\
%                   |                           \
%                   |                            \
%                   |                             ****|
%                   |                                 |
%                   |                                 ****|
%                   ***************                       ***********
%
% coded by Manuel Diaz, 2012.12.25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc; 

%% Parameters
CFL      = 0.7;     % CFL number
tEnd     = 1.0;     % Final time
nx       = 200;     % number of cells
n        = 1;       % Number of degrees of freedom
IC       = 07;      % IC: {1}~{12}. See Euler_IC1d.m
etpfix   = 0.90;	% {#} Harten's sonic entropy fix value, {0} no entropy fix
plot_fig = 1;       % {1} plot figures, {0} do NOT plot figures
wrt_sol  = 1;       % {1} write solution, {0} do NOT write solution file

%% Domain & IC

% Ratio of specific heats for ideal di-atomic gas
gamma=(n+2)/n;

% Discretize spatial domain
a=0; b=1; x=linspace(a,b,nx); dx=(b-a)/nx;

% set IC
[r,u,p,~,~] = Euler_IC1d(x,IC);

% Specific internal energy of the IC
e_s = p./((gamma-1).*r);
% Internal energy of the IC : e = e_s*rho
e = p./(gamma-1);
% Total Enegy: E = e + 1/2 rho*u^2
E = e + 1/2*r.*u.^2;
% Total Enthalpy: H = g/(g-1)p/rho + 1/2 u^2
H = gamma/(gamma-1)*p./r + u.^2/2;

%% Exact Riemann Solution
[xe,rhoe,ue,pe,ee,te,Me,se] = ...
   EulerExact(r(1),u(1),p(1),r(nx),u(nx),p(nx),tEnd,n);

%% Compute primitive Variables @ cell center variables
% Build vector Q and dQgit , and F
    % Define:
    L = 1:nx-1; R = 2:nx; % inner nodes / middle points
    % Q = [Q1 Q2 Q3]' % defined at every cell center
    Q = [r; r.*u; E];
    % Q1 = Density, Q2 = Momentum, Q3 = Total Energy

% Initial time
t = 0; it = 0;
        
%% Main Loop
Q_next = zeros(3,nx);
%for time = t; 
while t < tEnd
    % Compute Roe Averages
    % Velovity 'u'
    u_bar = (sqrt (r(L)).*u(L) + sqrt (r(R)).*u(R)) ...
                ./ (sqrt(r(L))+sqrt(r(R)));
    % Total Entalpy 'H'
    H_bar = (sqrt (r(L)).*H(L) + sqrt (r(R)).*H(R)) ...
                ./ (sqrt(r(L))+sqrt(r(R)));
    % Sound Speed 'a'
    a_bar = sqrt((gamma-1)*(H_bar-0.5*u_bar.^2));

    % Compute Delta U's
    % dU = U(R)-U(L) at the cell boundaries { x_{i},x_{i+1}, ... }
    dQ = Q(:,R)-Q(:,L);

    % Compute Fluxes at the cell centers
    % F = [F1 F2 F3]
    F = [r.*u; r.*u.^2 + p; r.*u.*H]; 
    % Fluxes to the left and right of 'F_ {i+1/2}'
    FL = F(:,L); FR = F(:,R);
    
    % Scaled Right Eigenvectors are given by:
    k1_bar = [1*ones(1,nx-1); u_bar-a_bar; H_bar-u_bar.*a_bar];
    k2_bar = [1*ones(1,nx-1); u_bar      ; 1/2*u_bar.^2      ];
    k3_bar = [1*ones(1,nx-1); u_bar+a_bar; H_bar+u_bar.*a_bar];
    
    % compute Roe waves strength alpha_bar
    alpha2_bar = (gamma-1)./(a_bar.^2).*(dQ (1,:).*(H_bar-u_bar.^2) ...
        + u_bar.*dQ(2,:) - dQ(3,:));
    alpha1_bar = 1./(2*a_bar).*(dQ (1,:).*(u_bar+a_bar) ...
        - dQ(2,:)-a_bar.*alpha2_bar);
    alpha3_bar = dQ(1,:)-(alpha1_bar + alpha2_bar);
    
    % Eigenvalues of A_bar are (same as the original A mat)
    lambda_bar(1,:) = abs(u_bar - a_bar);
    lambda_bar(2,:) = abs(u_bar);
    lambda_bar(3,:) = abs(u_bar + a_bar);
    
    % Entropy Fix
    boolv1 = lambda_bar(1,:) < etpfix;
        lambda_bar(1,:) = 0.5*(etpfix + lambda_bar (1,:).^2/etpfix).*boolv1 ...
            + lambda_bar(1,:).*(1-boolv1);
    boolv2 = lambda_bar(3,:) < etpfix;
        lambda_bar(3,:) = 0.5*(etpfix + lambda_bar (3,:).^2/etpfix).*boolv2 ...
            + lambda_bar(3,:).*(1-boolv2);
        
    % Conditioning data
    alpha1_bar = repmat(alpha1_bar,3,1);
    alpha2_bar = repmat(alpha2_bar,3,1);
    alpha3_bar = repmat(alpha3_bar,3,1);
    lambda1_bar = repmat(lambda_bar(1,:),3,1);
    lambda2_bar = repmat(lambda_bar(2,:),3,1);
    lambda3_bar = repmat(lambda_bar(3,:),3,1);

    % Update time step
    dt=dx*CFL/max(abs(lambda_bar(3,:))); t=t+dt; dtdx = dt/dx;
    
    % Roe Fluxes
    Flux = 0.5*(FL+FR)-0.5*(alpha1_bar.*lambda1_bar.*k1_bar + ...
                            alpha2_bar.*lambda2_bar.*k2_bar + ...
                            alpha3_bar.*lambda3_bar.*k3_bar );
                        
    % Compute next time step
    for i = 2:nx-1
        Q_next(:,i) = Q(:,i) - dtdx*(Flux(:,i) - Flux(:,i-1));
    end
    
    % BCs
    Q_next(:,1)  = Q(:,2); % Neumann condition to the left
    Q_next(:,nx) = Q(:,nx-1); % Neumann condition to the right
    
    % Compute variables of the new time step
    r_next = Q_next(1,:);               % Density
    u_next = Q_next (2,:)./Q_next(1,:);  % Velocity
    E_next = Q_next(3,:);               % Total Energy
    p_next = (gamma-1).*(E_next-r_next.*u_next.^2/2);  % Pressure
    e_next = 1/(gamma-1)*(p_next./r_next);      % Internal Energy
    a_bar_next = sqrt(gamma*p_next./r_next);    % sound speed
    m_next = u_next./a_bar_next;        % Mach 
    s_next = log(p_next./r_next.^gamma);% Entropy
    H_next = (E_next + p_next)./r_next; % Enthalpy
    
    % Update info
    Q = Q_next;
    r = r_next;
    u = u_next;
    e = e_next;
    p = p_next;
    m = m_next;
    s = s_next;
    H = H_next;
        
    % Update Counter
    it=it+1;
    
    % Plot figure
    if rem(it,1) == 0
        if plot_fig == 1;
            subplot(2,3,1); plot(x,r,'.b');
            subplot(2,3,2); plot(x,u,'.m'); 
            subplot(2,3,3); plot(x,p,'.k'); 
            subplot(2,3,4); plot(x,s,'.c');
            subplot(2,3,5); plot(x,m,'.g');
            subplot(2,3,6); plot(x,e,'.r');
        end
    drawnow
    end
end

%% Write and plot final results
if wrt_sol == 1;
    s1=subplot(2,3,1); plot(x,r,'.b',xe,rhoe,'k'); xlabel('x'); ylabel('Density');
    s2=subplot(2,3,2); plot(x,u,'.m',xe,ue,'k'); xlabel('x'); ylabel('Velocity');
    s3=subplot(2,3,3); plot(x,p,'.k',xe,pe,'k'); xlabel('x'); ylabel('Pressure');
    s4=subplot(2,3,4); plot(x,s,'.c',xe,se,'k'); xlabel('x'); ylabel('Entropy/R gas');
    s5=subplot(2,3,5); plot(x,m,'.g',xe,Me,'k'); xlabel('x'); ylabel('Mach number');
    s6=subplot(2,3,6); plot(x,e,'.r',xe,ee,'k'); xlabel('x'); ylabel('Internal Energy');
    title(s1,['Roe Approximate Euler Solver for gas of','  \gamma:',num2str(rats(gamma))]);
end