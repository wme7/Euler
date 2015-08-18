function [mu,k] = transportCoefs(q,gamma,Pr,T_inf,Re_inf,M_inf,L_inf)
    % Load Macroscopic data
    r=q(1,:);     u=q(2,:)./r;      E=q(3,:)./r; 
    p=(gamma-1)*r.*(E-0.5*u.^2);    T=gamma*p./r;

    % Compute transport Coeficients
    C = 110.5; %[K] : for air
    mu= (1+C/T_inf)./(T+C/T_inf).*T.^(3/2) * M_inf/Re_inf;
    k = L_inf*gamma*mu/(Pr*(gamma-1)) * M_inf/Re_inf/gamma;
    
end