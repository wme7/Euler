function [x,r,u,T,p,tau_nn,q_x] = StationaryShock(varargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Matlab Program for the shock structure calculations
%                  coded by Manuel Diaz, NTU, 2015.03.21
%
% Refs:
% [1] Xu, Kun. "A gas-kinetic BGK scheme for the Navier-Stokes equations
%     and its connection with artificial dissipation and Godunov method."
%     Journal of Computational Physics 171.1 (2001): 289-335.  
% [2] Gilbarg, Dj, and D. Paolucci. "The structure of shock waves in the
%     continuum theory of fluids." journal of Rational Mechanics and
%     Analysis 2.5 (1953): 617-642.  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    case 7
        % Controlling Parameters
        gamma=varargin{1};
        M  = varargin{2}; %mach number
        Pr = varargin{3};
        eps= varargin{4};
        x0 = varargin{5};
        dx = varargin{6};
        mu1= varargin{7};
        % upstream states:
        r1 = 1.0;
        u1 = 1.0;
        p1 = 1/gamma/M^2;
    case 10
        % Controlling Parameters
        gamma=varargin{1};
        M  = varargin{2}; %mach number
        Pr = varargin{3};
        eps= varargin{4};
        x0 = varargin{5};
        dx = varargin{6};
        mu1= varargin{7};
        % upstream states:
        r1 = varargin{8};
        u1 = varargin{9};
        p1 = varargin{10};
end
%Temperature on state 1
T1=2*p1/r1;

% Downstream states:
r2=(gamma+1)*M^2/(2+(gamma-1)*M^2)*r1;
u2=((gamma-1)/(gamma+1)+2/(gamma+1)/M^2)*u1;
p2=(2*gamma/(gamma+1)*M^2-(gamma-1)/(gamma+1))*p1;
T2=2*p2/r2;

% Display upstream and downstream conditions
fprintf('     1    2\n');
fprintf('r: %1.2f %1.2f\n',r1,r2);
fprintf('u: %1.2f %1.2f\n',u1,u2);
fprintf('t: %1.2f %1.2f\n',T1,T2);

% Conservations:
A=r1*u1;
B=r1*u1^2+p1;
C=(1/2*r1*u1^2+p1/(gamma-1)+p1)*u1;

% ODE solver:
y0=[T2; -u2*(1+eps)];       %downstream states [T2; -U2]
%tspan=-[(x0+dx),(x0-dx)];	%x=-tspan, solve from downstream to upstream
tspan=[-0.25,0.25];             %x=-tspan, solve from downstream to upstream
options=odeset('RelTol',1e-6,'AbsTol',[1e-9 1e-9]);
[t , y]=ode45(@Fb,tspan,y0,options);
x=-t; T=y(:,1); u=-y(:,2); mu=mu1*(T/T1).^0.8; 

% U_x and T_x are given by
u_x=-3/4./mu.*(B-A*u-1/2*A*T./u); 
T_x=4/5*Pr./mu.*(-1/2*A*u.^2+3/4*A*T-C+B*u); 

% Conserved properties
r=A./u;             % density
p=r.*T/2;           % pressure
c=sqrt(gamma*T);    % sound speed

% Normal stress and heat flux 
tau_nn = 4/3*mu.*u_x./(2*p);
q_x = -5/4/Pr*mu.*T_x./p./sqrt(gamma*T/2);

%% Extra function Fb
    function dy=Fb(t,y)
        %
        % Where:  y(1): T,  y(2): -U,  t: -x
        %
        Mu=mu1*(y(1)/T1)^0.8;
        %Mu=mu1;
        dy=[0;0];
        dy(1)=4/5*Pr/Mu*(1/2*A*y(2)^2-3/4*A*y(1)+C+B*y(2));
        dy(2)=3/4/Mu*(-B-A*y(2)-1/2*A*y(1)/y(2));
    end
end