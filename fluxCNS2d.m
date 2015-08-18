function viscFlux = fluxCNS2d(q,gamma,dx,dy,nx,ny,mu,k)

% Load conservative variables
r=q(:,:,1);  u=q(:,:,2)./r;  v=q(:,:,3)./r;  E=q(:,:,4)./r; 
p=(gamma-1)*r.*(E-0.5*(u.^2+v.^2));  T=p./r;

% Initialize viscous fluxes array
cell(ny,nx).vFlux = 0;
parfor i = 1:ny
    for j = 1:nx
        cell(i,j).vFlux = zeros(4,1);
    end
end

% Compute F_{i,j+1/2}
F = [0;0;0;0];
for i = 2:ny-1
    for j = 1:nx-1
        F(2) = (4*mu)/(3*dx)*(u(i,j+1)-u(i,j)) - mu/(6*dy)*( (v(i+1,j+1)-v(i-1,j+1)) + (v(i+1,j)-v(i-1,j)) );
        F(3) = mu/(4*dy)*((u(i+1,j+1)-u(i-1,j+1))+(u(i+1,j)-u(i-1,j))) + (mu/dx)*(v(i,j+1)-v(i,j));
        F(4) = 0.5*(u(i,j+1)+u(i,j))*F(2) + 0.5*(v(i,j+1)+v(i,j))*F(3) + (k/dx)*(T(i,j+1)-T(i,j));
        cell( i,j ).vFlux = cell( i,j ).vFlux + F/dx;
        cell(i,j+1).vFlux = cell(i,j+1).vFlux - F/dx;
    end
end

% Compute G_{i+1/2,j} 
G = [0;0;0;0];
for i = 1:ny-1
    for j = 2:nx-1
        F(2) = (mu/dy)*(u(i+1,j)-u(i,j)) + mu/(4*dx)*((v(i+1,j+1)-v(i+1,j-1))+(v(i,j+1)-v(i,j-1)));
        F(3) = (4*mu)/(3*dx)*(v(i+1,j)-v(i,j)) - mu/(6*dy)*( (u(i+1,j+1)-u(i+1,j-1)) + (u(i,j+1)-u(i,j-1)) );
        F(4) = 0.5*(u(i+1,j)+u(i,j))*G(2) + 0.5*(v(i+1,j)+v(i,j))*G(3) + (k/dy)*(T(i+1,j)-T(i,j));
        cell( i,j ).vFlux = cell( i,j ).vFlux + G/dy;
        cell(i+1,j).vFlux = cell(i+1,j).vFlux - G/dy;
    end
end

% Viscous flux
viscFlux = zeros(ny,nx,4);
parfor i = 1:ny
    for j = 1:nx
        viscFlux(i,j,:) = cell(i,j).vFlux;
    end
end