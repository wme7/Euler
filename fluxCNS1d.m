function [viscFlux] = fluxCNS1d(q,gamma,dx,nx,Pr,T_inf,Re_inf,M_inf,L_inf)

% Load conservative variables
r=q(1,:);       u=q(2,:)./r;        E=q(3,:)./r; 
p=(gamma-1)*r.*(E-0.5*u.^2);        T=gamma*p./r;

% compute transport coefs
[mu,k] = transportCoefs(q,gamma,Pr,T_inf,Re_inf,M_inf,L_inf);

% Initialize viscous fluxes array
cell(nx).vFlux = 0;
for j = 1:nx; cell(j).vFlux = zeros(3,1); end

% Compute viscous F_{i,j+1/2}
F = [0;0;0];
for j = 1:nx-1
        F(2)= (4*mu(j))/(3*dx)*(u(j+1)-u(j));
        F(3) = 0.5*(u(j+1)+u(j))*F(2) + (k(j)/dx)*(T(j+1)-T(j));
        cell( j ).vFlux = cell( j ).vFlux + F/dx;
        cell(j+1).vFlux = cell(j+1).vFlux - F/dx;
end

% Viscous flux
viscFlux = zeros(3,nx);
for j = 1:nx; viscFlux(:,j) = cell(j).vFlux; end