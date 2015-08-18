function [flux,wsn] = fluxfunc2d(primL,primR,njk,fluxType)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Flux Functions form (Refs):
%
% [1] P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and
% Difference Schemes, Journal of Computational Physics, 43, pp. 357-372.
%
% [2] H. Nishikawa and K. Kitamura, Very Simple, Carbuncle-Free,
% Boundary-Layer Resolving, Rotated-Hybrid Riemann Solvers,
% Journal of Computational Physics, 227, pp. 2560-2581, 2008.
%
% -------------------------------------------------------------------------
%  Input:   primL(1:4) =  left state (rhoL, uL, vL, pL)
%           primR(1:4) = right state (rhoR, uR, vR, pR)
%               njk(2) = Face normal (L -> R). Must be a unit vector.
%
% Output:    flux(1:4) = numerical flux
%                  wsn = half the max wave speed
%                        (to be used for time step calculations)
% -------------------------------------------------------------------------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch fluxType
    case 'ROE' % -- Roe's Flux Function with entropy fix---

  nx = njk(1);
  ny = njk(2);

%Tangent vector (Do you like it? Actually, Roe flux can be implemented 
% without any tangent vector. See "I do like CFD, VOL.1" for details.)
  mx = -ny;
  my =  nx;

%Primitive and other variables.
%  Left state
    rhoL = primL(1);
      uL = primL(2);
      vL = primL(3);
     unL = uL*nx+vL*ny;
     umL = uL*mx+vL*my;
      pL = primL(4);
      aL = sqrt(gamma*pL/rhoL);
      HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL);
%  Right state
    rhoR = primR(1);
      uR = primR(2);
      vR = primR(3);
     unR = uR*nx+vR*ny;
     umR = uR*mx+vR*my;
      pR = primR(4);
      aR = sqrt(gamma*pR/rhoR);
      HR = aR*aR/(gamma-one) + half*(uR*uR+vR*vR);

%First compute the Roe Averages
    RT = sqrt(rhoR/rhoL);
   rho = RT*rhoL;
     u = (uL+RT*uR)/(one+RT);
     v = (vL+RT*vR)/(one+RT);
     H = (HL+RT* HR)/(one+RT);
     a = sqrt( (gamma-one)*(H-half*(u*u+v*v)) );
    un = u*nx+v*ny;
    um = u*mx+v*my;

%Wave Strengths
   drho = rhoR - rhoL; 
     dp =   pR - pL;
    dun =  unR - unL;
    dum =  umR - umL;

  LdU(1) = (dp - rho*a*dun )/(two*a*a);
  LdU(2) = rho*dum;
  LdU(3) =  drho - dp/(a*a);
  LdU(4) = (dp + rho*a*dun )/(two*a*a);

%Wave Speed
  ws(1) = abs(un-a);
  ws(2) = abs(un);
  ws(3) = abs(un);
  ws(4) = abs(un+a);

%Harten's Entropy Fix JCP(1983), 49, pp357-393:
% only for the nonlinear fields.
  dws(1)=1/5; if ws(1) < dws(1); ws(1) = 0.5*( ws(1)*ws(1)/dws(1)+dws(1) ); end
  dws(4)=1/5; if ws(4) < dws(4); ws(4) = 0.5*( ws(4)*ws(4)/dws(4)+dws(4) ); end

%Right Eigenvectors
  Rv(1,1) = 1;
  Rv(2,1) = u - a*nx;
  Rv(3,1) = v - a*ny;
  Rv(4,1) = H - un*a;

  Rv(1,2) = 0;
  Rv(2,2) = mx;
  Rv(3,2) = my;
  Rv(4,2) = um;

  Rv(1,3) = 1;
  Rv(2,3) = u;
  Rv(3,3) = v;
  Rv(4,3) = 0.5*(u*u+v*v);

  Rv(1,4) = 1;
  Rv(2,4) = u + a*nx;
  Rv(3,4) = v + a*ny;
  Rv(4,4) = H + un*a;

%Dissipation Term
  diss = zeros(1,4);
  for i=1:4;
   for j=1:4;
    diss(i) = diss(i) + ws(j)*LdU(j)*Rv(i,j);
   end 
  end

%Compute the flux.
  fL(1) = rhoL*unL;
  fL(2) = rhoL*unL * uL + pL*nx;
  fL(3) = rhoL*unL * vL + pL*ny;
  fL(4) = rhoL*unL * HL;

  fR(1) = rhoR*unR;
  fR(2) = rhoR*unR * uR + pR*nx;
  fR(3) = rhoR*unR * vR + pR*ny;
  fR(4) = rhoR*unR * HR;

  flux = 0.5 * (fL + fR - diss);
  wsn = 0.5*(abs(un) + a);  %Normal max wave speed times half

    case 'RHLL' % -- Rotated-RHLL Flux Function ---

  nx = njk(1);
  ny = njk(2);

%Tangent vector (Do you like it? Actually, Roe flux can be implemented 
% without any tangent vector. See "I do like CFD, VOL.1" for details.)
  %mx = -ny;
  %my =  nx;

%Primitive and other variables.
%  Left state
    rhoL = primL(1);
      uL = primL(2);
      vL = primL(3);
     unL = uL*nx+vL*ny;
     %umL = uL*mx+vL*my;
      pL = primL(4);
      aL = sqrt(gamma*pL/rhoL);
      HL = aL*aL/(gamma-1) + 0.5*(uL*uL+vL*vL);
%  Right state
    rhoR = primR(1);
      uR = primR(2);
      vR = primR(3);
     unR = uR*nx+vR*ny;
     %umR = uR*mx+vR*my;
      pR = primR(4);
      aR = sqrt(gamma*pR/rhoR);
      HR = aR*aR/(gamma-1) + 0.5*(uR*uR+vR*vR);

%Primitive and oth
%Compute the flux.
   fL(1) = rhoL*unL;
   fL(2) = rhoL*unL * uL + pL*nx;
   fL(3) = rhoL*unL * vL + pL*ny;
   fL(4) = rhoL*unL * HL;

   fR(1) = rhoR*unR;
   fR(2) = rhoR*unR * uR + pR*nx;
   fR(3) = rhoR*unR * vR + pR*ny;
   fR(4) = rhoR*unR * HR;

%Define n1 and n2, and compute alpha1 and alpha2: (4.2) in the original paper.
%(NB: n1 and n2 may need to be frozen at some point during 
%     a steady calculation to fully make it converge. For time-accurate 
%     calculation, this is fine.)
% NB: For a boundary face, set (nx2,ny2)=(nx,ny), (nx1,ny1)=(-ny,nx).

    eps = 1.0e-12 * M_inf;
    abs_dq = sqrt( (uR-uL)^2 + (vR-vL)^2 );

  if ( abs_dq > eps)
    nx1 = (uR-uL)/abs_dq;
    ny1 = (vR-vL)/abs_dq;
  else
    nx1 = -ny; 
    ny1 =  nx;
  end

%  Rey = 1000.0_p2
%  temp = ( tanh(Rey*(abs_dq-eps)) - tanh(-Rey*eps) ) &
%        /( tanh(Rey*(   one-eps)) - tanh(-Rey*eps) )
%  nx1 = temp*(uR-uL)/(abs_dq + eps) + (one-temp)*(-ny)
%  ny1 = temp*(vR-vL)/(abs_dq + eps) + (one-temp)*( nx)

    alpha1 = nx * nx1 + ny * ny1;
%   To make alpha1 always positive.
      temp = sign(one,alpha1);
       nx1 = temp * nx1;
       ny1 = temp * ny1;
    alpha1 = temp * alpha1;

% Take n2 as perpendicular to n1.
       nx2 = -ny1;
       ny2 =  nx1;
    alpha2 = nx * nx2 + ny * ny2;
%   To make alpha2 always positive.
      temp = sign(one,alpha2);
       nx2 = temp * nx2;
       ny2 = temp * ny2;
    alpha2 = temp * alpha2;

%Now we are going to compute the Roe flux with n2 as the normal
%and n1 as the tagent vector, with modified wave speeds (5.12)

     RT = sqrt(rhoR/rhoL);
    rho = RT*rhoL;
      u = (uL + RT*uR)/(one + RT);
      v = (vL + RT*vR)/(one + RT);
      H = (HL + RT*HR)/(one + RT);
      a = sqrt( (gamma-1)*(H-0.5*(u*u+v*v)) );
     un = u*nx2+v*ny2;
     um = u*nx1+v*ny1;

%Wave Strengths (remember that n2 is the normal and n1 is the tangent.)
    unL = uL*nx2 + vL*ny2;
    unR = uR*nx2 + vR*ny2;
    umL = uL*nx1 + vL*ny1;
    umR = uR*nx1 + vR*ny1;

   drho = rhoR - rhoL;
     dp =   pR - pL;
    dun =  unR - unL;
    dum =  umR - umL;

  LdU(1) = (dp - rho*a*dun )/(2*a*a);
  LdU(2) =  rho*dum;
  LdU(3) =  drho - dp/(a*a);
  LdU(4) = (dp + rho*a*dun )/(2*a*a);

%Wave Speeds for Roe flux part.
    ws(1) = un-a;
    ws(2) = un;
    ws(3) = un;
    ws(4) = un+a;
  abs_ws  = abs(ws);

%Harten's Entropy Fix JCP(1983), 49, pp357-393:
% only for the nonlinear fields.
  dws(1)=1/5; if ws(1) < dws(1); ws(1) = 0.5*( ws(1)*ws(1)/dws(1)+dws(1) ); end
  dws(4)=1/5; if ws(4) < dws(4); ws(4) = 0.5*( ws(4)*ws(4)/dws(4)+dws(4) ); end

%HLL wave speeds, evaluated with [nx1,ny1] (=tangent wrt n2).
   SRp = max( zero, umR + aR, um + a);
   SLm = min( zero, umL - aL, um - a);

%Modified wave speeds for the Rotated-RHLL flux: (5.12) in the original paper.
   ws = alpha2*abs_ws - ( alpha2*(SRp+SLm)*ws + two*alpha1*SRp*SLm )/ (SRp-SLm);

%Right Eigenvectors: with n2 as normal and n1 as tangent.
  mx = nx1;
  my = ny1;
  
  Rv(1,1) = 1;
  Rv(2,1) = u - a*nx2;
  Rv(3,1) = v - a*ny2;
  Rv(4,1) = H - a*un;

  Rv(1,2) = 0;
  Rv(2,2) = mx;
  Rv(3,2) = my;
  Rv(4,2) = um;

  Rv(1,3) = 1;
  Rv(2,3) = u;
  Rv(3,3) = v;
  Rv(4,3) = 0.5*(u*u+v*v);

  Rv(1,4) = 1;
  Rv(2,4) = u + a*nx2;
  Rv(3,4) = v + a*ny2;
  Rv(4,4) = H + a*un;

%Dissipation Term: Roe dissipation with the modified wave speeds.
  diss = zeros(1,4);
  for i=1:4
   for j=1:4
    diss(i) = diss(i) + ws(j)*LdU(j)*Rv(i,j);
   end
  end

%Compute the Rotated-RHLL flux.
  flux = (SRp*fL - SLm*fR)/(SRp-SLm) - 0.5*diss;

  wsn = 0.5*( abs(un) + abs(um) + a);  %Normal max wave speed times half

end