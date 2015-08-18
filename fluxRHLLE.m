function [RHLLE]= fluxRHLLE(qL,qR,gamma,normal)
    % Evaluate Rotate-HLLE fluxes

    % normal vectors
    nx = normal(1);
    ny = normal(2);
    
    % Tangent vectors
    tx = -ny;
    ty = nx;
    
    % Left state
    rL = qL(1);
    uL = qL(2)/rL;
    vL = qL(3)/rL;
    vnL = uL*nx+vL*ny;
    vtL = uL*tx+vL*ty;
    pL = (gamma-1)*( qL(4) - rL*(uL^2+vL^2)/2 );
    aL = sqrt(gamma*pL/rL);
    HL = ( qL(4) + pL ) / rL;
    
    % Right state
    rR = qR(1);
    uR = qR(2)/rR;
    vR = qR(3)/rR;
    vnR = uR*nx+vR*ny;
    vtR = uR*tx+vR*ty;
    pR = (gamma-1)*( qR(4) - rR*(uR^2+vR^2)/2 );
    aR = sqrt(gamma*pR/rR);
    HR = ( qR(4) + pR ) / rR;
    
    % Define n1 and n2, and compute alpha1 and alpha2: (4.2) in the original paper.
    % (NB: n1 and n2 may need to be frozen at some point during
    %     a steady calculation to fully make it converge. For time-accurate
    %     calculation, this is fine.)
    % NB: For a boundary face, set (nx2,ny2)=(nx,ny), (nx1,ny1)=(-ny,nx).
    
    eps = 1.0e-12; % * M_inf;
    abs_dq = sqrt( (uR-uL)^2 + (vR-vL)^2 );
    
    if ( abs_dq > eps)
        nx1 = (uR-uL)/abs_dq;
        ny1 = (vR-vL)/abs_dq;
    else
        nx1 = -ny;
        ny1 =  nx;
    end
    
    % Rey = 1000.0_p2
    % temp = ( tanh(Rey*(abs_dq-eps)) - tanh(-Rey*eps) ) &
    %       /( tanh(Rey*(   one-eps)) - tanh(-Rey*eps) )
    % nx1 = temp*(uR-uL)/(abs_dq + eps) + (one-temp)*(-ny)
    % ny1 = temp*(vR-vL)/(abs_dq + eps) + (one-temp)*( nx)
    
    alpha1 = nx * nx1 + ny * ny1;
    % To make alpha1 always positive.
    temp = sign(one,alpha1);
    nx1 = temp * nx1;
    ny1 = temp * ny1;
    alpha1 = temp * alpha1;
    
    % Take n2 as perpendicular to n1.
    nx2 = -ny1;
    ny2 =  nx1;
    alpha2 = nx * nx2 + ny * ny2;
    % To make alpha2 always positive.
    temp = sign(one,alpha2);
    nx2 = temp * nx2;
    ny2 = temp * ny2;
    alpha2 = temp * alpha2;
    
    %Now we are going to compute the Roe flux with n2 as the normal
    %and n1 as the tagent vector, with modified wave speeds (5.12)
    
    % First compute the Roe Averages
    RT = sqrt(rR/rL);
    r = RT*rL;
    u = (uL+RT*uR)/(1+RT);
    v = (vL+RT*vR)/(1+RT);
    H = ( HL+RT* HR)/(1+RT);
    a = sqrt( (gamma-1)*(H-(u^2+v^2)/2) );
    vn = u*nx2+v*ny2;
    vt = u*nx1+v*ny1;
    
    % Wave Strengths
    dr = rR - rL;     dp = pR - pL;     dvn= vnR - vnL;     dvt= vtR - vtL;
    dV = [(dp-r*a*dvn )/(2*a^2); r*dvt/a; dr-dp/(a^2); (dp+r*a*dvn)/(2*a^2)];
    
    % Wave Speed
    ws = [vn-a; vn; vn; vn+a]; abs_ws = abs(ws);
    
    % Harten's Entropy Fix JCP(1983), 49, pp357-393:
    % only for the nonlinear fields.
    dws(1)=1/5; if ws(1)<dws(1); abs_ws(1)=( abs_ws(1)*abs_ws(1)/dws(1)+dws(1) )/2; end
    dws(4)=1/5; if ws(4)<dws(4); abs_ws(4)=( abs_ws(4)*abs_ws(4)/dws(4)+dws(4) )/2; end

    % Wave speed estimates, evaluated with [nx1,ny1] (= tangent wrt n2)
    SLm = min([ vtL-aL, vt-a, 0]);
    SRp = max([ vtR+aR, vt+a, 0]);
    
    % Modifed wave speed for the Rotated-RHLL flux (5.12) in the paper.
    ws = alpha2*abs_ws - ( alpha2*(SRp+SLm)*ws + 2*alpha1*SRp*SLm )/ (SRp-SLm);

    %Right Eigenvectors: with n2 as normal and n1 as tangent.
    tx = nx1;
    ty = ny1;
    
    % Right Eigenvectors       
    Rv = [  1   ,  0   ,    1      ,  1   ;
          u-a*nx, a*tx ,    u      ,u+a*nx;
          u-a*ny, a*ty ,    u      ,u+a*ny;
          H-vn*a, vt*a ,(u^2+v^2)/2,H+vn*a];
    
    % Left and Right fluxes
    FL=[rL*vnL; rL*vnL*uL + pL*nx; rL*vnL*vL + pL*ny; rL*vnL*HL];
    FR=[rR*vnR; rR*vnR*uR + pR*nx; rR*vnR*vR + pR*ny; rR*vnR*HR];
    
    % Compute the HLL flux.
    RHLLE = ( SRp*FL - SLm*FR )/(SRp-SLm) - 0.5*Rv*(ws.*dV);
end