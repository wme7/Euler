%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Build Exact solutions for BGK compressible N-S regimes 
%                  coded by Manuel Diaz, NTU, 2015.03.21
%
% Refs:
% [1] Xu, Kun. "A gas-kinetic BGK scheme for the Navier-Stokes equations
%     and its connection with artificial dissipation and Godunov method."
%     Journal of Computational Physics 171.1 (2001): 289-335.  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; %close all;

% Couette Flows ExactSolutions.
figure(1);

[x1,tR1] = CouetteFlow(40,0.72);
[x2,tR2] = CouetteFlow(40,1.00);
[x3,tR3] = CouetteFlow(40,2.50);

subplot(121); plot(x1,tR1,x2,tR2,x3,tR3);

[x1,tR1] = CouetteFlow(04,0.5);
[x2,tR2] = CouetteFlow(20,0.5);
[x3,tR3] = CouetteFlow(40,0.5);

subplot(122); plot(x1,tR1,x2,tR2,x3,tR3);

% Stationary Shock, Navier-Stokes Exact Solution.
figure(2);

x0=0.6; dx=0.1; gamma=5/3; M=10; Pr=2/3; mu1=0.0005; eps=1e-6;
[x,~,u,T,~,tau_nn,q_x] = StationaryShock(gamma,M,Pr,eps,x0,dx,mu1);

subplot(221); plot(x,T); ylabel('T'); xlabel('x');
subplot(222); plot(x,u); ylabel('u'); xlabel('x');
subplot(2,2,[3,4]); plot(u,tau_nn,u,q_x); xlabel('u/U_\infty');
legend('\tau_{nn}','q_x','location','southwest'); legend('boxoff'); 

% Stationary Shock - density profiles, Navier-Stokes Exac Solution.
figure(3)

x0=0.0; dx=0.5; gamma=5/3; M=3.0; Pr=3/4; mu1=0.0025; eps=1e-10;
r1 = 1.0; u1 = 1.0; p1 = 1/gamma/M^2;

[x,r,u,T] = StationaryShock(gamma,M,Pr,eps,x0,dx,mu1,r1,u1,p1);
plot(x,r);