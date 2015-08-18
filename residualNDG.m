function dF = residualNDG(q,K,Lift,Dr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Compute the Residual for 1d wave equation using Nodal DG 
%
%                       Residual = dF/dxi 
%                 where F = is our Flux function
%
%              coded by Manuel Diaz, NTU, 2013.10.29
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global gamma

% q values
q1=q(1:K+1,:); q2=q(K+2:2*(K+1),:); q3=q(2*K+3:3*(K+1),:);

% Compute Flux values
rho=q1; u=q2./q1; E=q3./q1;p=(gamma-1)*rho.*(E-0.5*u.^2);
f1=rho.*u; f2=rho.*u.^2+p; f3=u.*(rho.*E+p);

% Maximun velocity for LF flux
a = sqrt(gamma*p./rho); % Speed of sound
alpha = max(max(abs(u+a)));

% Interpolate q and flux values at the boundaries of Ij
% use only: {'LGL','ChebyshevMod'}
q1_lbd = q1(1,:);   q2_lbd = q2(1,:);	q3_lbd = q3(1,:);
q1_rbd = q1(end,:); q2_rbd = q2(end,:); q3_rbd = u(end,:);
f1_lbd = f1(1,:);   f2_lbd = f2(1,:);   f3_lbd = f3(1,:);
f1_rbd = f1(end,:); f2_rbd = f2(end,:); f3_rbd = f3(end,:);

% Build Numerical fluxes across faces
q1_pface = [q1_lbd,0]; f1_pface = [f1_lbd,0]; % + side 
q1_nface = [0,q1_rbd]; f1_nface = [0,f1_rbd]; % - side 
q2_pface = [q2_lbd,0]; f2_pface = [f2_lbd,0]; % + side 
q2_nface = [0,q2_rbd]; f2_nface = [0,f2_rbd]; % - side 
q3_pface = [q3_lbd,0]; f3_pface = [f3_lbd,0]; % + side 
q3_nface = [0,q3_rbd]; f3_nface = [0,f3_rbd]; % - side 

% Apply Neumann BCs
q1_nface(1)   = q1_pface(1);    f1_nface(1)   = f1_pface(1);    % left BD
q1_pface(end) = q1_nface(end);  f1_pface(end) = f1_nface(end);  % right BD
q2_nface(1)   = q2_pface(1);    f2_nface(1)   = f2_pface(1);    % left BD
q2_pface(end) = q2_nface(end);  f2_pface(end) = f2_nface(end);  % right BD
q3_nface(1)   = q3_pface(1);    f3_nface(1)   = f3_pface(1);    % left BD
q3_pface(end) = q3_nface(end);  f3_pface(end) = f3_nface(end);  % right BD

% LF numerical flux
nflux1 = 0.5*(f1_nface+f1_pface-alpha*(q1_pface-q1_nface));
nflux2 = 0.5*(f2_nface+f2_pface-alpha*(q2_pface-q2_nface));
nflux3 = 0.5*(f3_nface+f3_pface-alpha*(q3_pface-q3_nface));
nflux1L = nflux1(1:end-1); nflux1R = nflux1(2:end);
nflux2L = nflux2(1:end-1); nflux2R = nflux2(2:end);
nflux3L = nflux3(1:end-1); nflux3R = nflux3(2:end);

% Compute the derivate: F = f + lL*(nfluxL-f_bdL) - lR*(nfluxR-f_bdR)
dF1 = -Dr*f1 + Lift*[(nflux1L-f1_lbd);(f1_rbd-nflux1R)];
dF2 = -Dr*f2 + Lift*[(nflux2L-f2_lbd);(f2_rbd-nflux2R)];
dF3 = -Dr*f3 + Lift*[(nflux3L-f3_lbd);(f3_rbd-nflux3R)];
dF = -[dF1;dF2;dF3];