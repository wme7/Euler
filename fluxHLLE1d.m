function [HLLE]= fluxHLLE1d(qL,qR,gamma)
    % evaluate HLLE fluxes 1d: 
    % no normal/tangential vectors assumed here!

    % Left state
    rL = qL(1);
    uL = qL(2)/qL(1);
    pL = (gamma-1)*( qL(3) - 0.5*rL*uL*uL );
    aL = sqrt(gamma*pL/rL);
    HL = ( qL(3) + pL ) / rL;
    
    % Right state
    rR = qR(1);
    uR = qR(2)/qR(1);
    pR = (gamma-1)*( qR(3) - 0.5*rR*uR*uR );
    aR = sqrt(gamma*pR/rR);
    HR = ( qR(3) + pR ) / rR;

    % Evaluate the two wave speeds: Einfeldt.
    RT = sqrt(rR/rL); %r = RT*rL;
    u = (uL+RT*uR)/(1+RT);
    H = (HL+RT*HR)/(1+RT);
    a = sqrt( (gamma-1)*(H-0.5*u*u) );
    
    % Wave estimates for unified definition
    SLm = min([ uL-aL, u-a, 0]);
    SRp = max([ uR+aR, u+a, 0]);
    
    % Left and Right fluxes
    FL=[rL.*uL; rL.*uL.^2+pL; uL.*rL.*HL];
    FR=[rR.*uR; rR.*uR.^2+pR; uR.*rR.*HR];

    % Compute the HLL flux: (using single equation definition)
    % Right-going supersonic flow, Subsonic flow, Left-going supersonic flow
    HLLE = (SRp*FL-SLm*FR + SLm*SRp*(qR-qL))/(SRp-SLm);
end