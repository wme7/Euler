function Shock()
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Matlab Program for the shock structure calculations
%                  coded by Manuel Diaz, NTU, 2015.03.21
%
% Refs:
% [1] Xu, Kun. "A gas-kinetic BGK scheme for the Navier???Stokes equations
%     and its connection with artificial dissipation and Godunov method."
%     Journal of Computational Physics 171.1 (2001): 289-335.  
% [2] Gilbarg, Dj, and D. Paolucci. "The structure of shock waves in the
%     continuum theory of fluids." journal of Rational Mechanics and
%     Analysis 2.5 (1953): 617-642.  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; %clc; close all;

% parameters
M=1.5; %mach number
gamma=5/3;
Pr=2/3;
eps=1e-12;

%upstream states:
r1=1.0;
u1=1.0;
mu1=0.00025;
p1=1/gamma/M^2;
T1=2*p1/r1;

%downstream states:
r2=(gamma+1)*M^2/(2+(gamma-1)*M^2)*r1;
u2=((gamma-1)/(gamma+1)+2/(gamma+1)/M^2)*u1;
p2=(2*gamma/(gamma+1)*M^2-(gamma-1)/(gamma+1))*p1;
T2=2*p2/r2;

%conservations:
A=r1*u1;
B=r1*u1^2+p1;
C=(1/2*r1*u1^2+p1/(gamma-1)+p1)*u1;

%ODE solver:
y0=[T2; -u2*(1+eps)]; %downstream states [T2; -U2]
tspan=[-0.05,0.05]; %x=-tspan, solve from downstream to upstream
options = odeset('RelTol',1e-6,'AbsTol',[1e-8 1e-8]);
[t, y]=ode45(@Fb,tspan, y0, options);
x=-t;
T=y(:,1);
u=-y(:,2);
mu=mu1* (T/T1).^0.8; 
ux=-3/4./mu.*(B-A*u-1/2*A*T./u); %U_x
Tx=4/5*Pr./mu.*(-1/2*A*u.^2+3/4*A*T-C+B*u); %T_x
r=A./u; %density
p=r.*T/2; %pressure
c = sqrt(gamma*T); %sound speed

% normal stress and heat flux 
tau_nn = 4/3*mu.*ux./(2*p);
q_x = -5/4*(mu/Pr).*Tx./(2*p.*c);

%plot T and u:
figure(1)
subplot(221); plot(x,T); xlabel('x'); ylabel('T'); axis('square');
subplot(222); plot(x,u); xlabel('x'); ylabel('U'); axis('square');
subplot(223); plot(x,tau_nn); xlabel('x'); ylabel('\tau'); axis('square');
subplot(224); plot(x,q_x); xlabel('x'); ylabel('q_x'); axis('square');

end

function dy=Fb(t,y)
%y(1): T
%y(2): -U
%t: -x
M=1.5; %mach number
gamma=5/3;
Pr=2/3;
%upstream states:
r1=1.0;
u1=1.0;
mu1=0.0005;
p1=1/gamma/M^2;
T1=2*p1/r1;
A=r1*u1;
B=r1*u1^2+p1;
C=(1/2*r1*u1^2+p1/(gamma-1)+p1)*u1;
mu=mu1*(y(1)/T1)^0.8;
%mu=mu1;
dy=[0;0];
dy(1)=4/5*Pr/mu*(1/2*A*y(2)^2-3/4*A*y(1)+C+B*y(2)); % -T_x
dy(2)=3/4/mu*(-B-A*y(2)-1/2*A*y(1)/y(2)); % U_x
end

%% write output files:
% out(:,1)=x;
% out(:,2)=T;
% save t out -ascii %save [x T] in file 't'
% out(:,1)=x;
% out(:,2)=u;
% save u out -ascii %save [x U] in file 'u'
% out(:,1)=x;
% out(:,2)=p;
% save p out -ascii %save [x p] in file 'p'
% out(:,1)=x;
% out(:,2)=r;
% save r out -ascii %save [x rho] in file 'r'
% out(:,1)=x;
% out(:,2)=mu;
% save mu out -ascii %save [x mu] in file 'mu'
% out(:,1)=x;
% out(:,2)=ux;
% save ux out -ascii %save [x U_x] in file 'ux'
% out(:,1)=x;
% out(:,2)=Tx;
% save tx out -ascii %save [x T_x] in file 'tx'
% out(:,1)=x;
% out(:,2)=4/3*mu.*ux./(2*p);
% save v out -ascii %save in file 'v'
% out(:,1)=x;
% out(:,2)=5/4/Pr*mu.*Tx./p./sqrt(gamma*T/2);
% save h out -ascii %save in file 'h'
