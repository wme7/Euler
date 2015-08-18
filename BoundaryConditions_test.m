%-------------------------------------------------------------------------
% Close with the boundary flux using the element-based formula that is
% exact for linear fluxes (See Nishikawa AIAA2010-5093 for boundary weights
% that ensure the linear exactness for 2D/3D elements).
%
%      |  Interior Domain          |
%      |        .........          |
%      |        .       .          |
%      |        .       .          |
%      o--o--o-----o---------o--o--o  <- Boundary segment
%                  n1   |   n2
%                       v
%                     n12 (unit face normal vector)
%
% NOTE: We visit each boundary face, defined by the nodes n1 and n2,
%       and compute the flux across the boundary face: left half for node1,
%       and the right half for node2. In the above figure, the dots indicate
%       the control volume around the node n1. Clearly, the flux across the
%       left half of the face contributes to the node n1. Similarly for n2.
%
%--------------------------------------------------------------------------
%  BC: Upwind flux via freestream values
%
%      NOTE: If the final solution at the boundary node is far from
%            the freestream values, then the domain is probably is not
%            large enough. 
%--------------------------------------------------------------------------
%  BC: Solid body and Supersonic outflow
%
%      NOTE: Basically, simply compute the physical flux, which
%            can be done by calling the Roe flux with wR = wL.
%            It is equivalent to the interior-extrapolation condition.
%      NOTE: Tangency condition for solid body will be applied later.
%--------------------------------------------------------------------------
%  BC: Subsonic Outflow - Fixed Back Pressure
%
%      NOTE: Fix the pressure as freestream pressure
%            on the right side of the face (outside the domain).
%            Assumption is that the outflow boundary is far from the body.
%--------------------------------------------------------------------------

for i = 1:nbound
    for  j = 1:bound(i).nbfaces
        % Get Left and Right base data:
        n1 = bound(i).bnode( j );       %Left node
        n2 = bound(i).bnode(j+1);       %Right node
        n12(1) = bound(i).bfnx(j);      %x-component of the unit face normal vector
        n12(2) = bound(i).bfny(j);      %y-component of the unit face normal vector
        mag_e12 = bound(i).bfn(j)*0.5;  %Half length of the boundary face, j.
        switch bound(i).bcname
            case 'freestream' % bnodes_numerical_flux_via_freestream
                %   1. Left node
                wL = node(n1).w;    wR = Winf;
                [num_flux,wsn] = fluxfunc2d(wL,wR,n12,fType);
                bfluxL = num_flux;
                node(n1).wsn = node(n1).wsn + wsn*mag_e12;
                
                %   2. Right node
                wL = node(n2).w;    wR = Winf;
                [num_flux,wsn] = fluxfunc2d(wL,wR,n12,fType);
                bfluxR = num_flux;
                node(n2).wsn = node(n2).wsn + wsn*mag_e12;
                
            case {'slip_wall','outflow_supersonic'} % bnodes_slip_wall
                %   1. Left node
                wL = node(n1).w;    wR = wL;
                [num_flux,wsn] = fluxfunc2d(wL,wR,n12,fType);
                bfluxL = num_flux;
                node(n1).wsn = node(n1).wsn + wsn*mag_e12;
                
                %   2. Right node
                wL = node(n2).w;    wR = wL;
                [num_flux,wsn] = fluxfunc2d(wL,wR,n12,fType);
                bfluxR = num_flux;
                node(n2).wsn = node(n2).wsn + wsn*mag_e12;
                
            case 'outflow_back_pressure' % bnodes_outflow
                %   1. Left node
                wL = node(n1).w;    wR = wL;    wR(4) = Winf(4); %Fix the pressure
                [num_flux,wsn] = fluxfunc2d(wL,wR,n12,fType);
                bfluxL = num_flux;
                node(n1).wsn = node(n1).wsn + wsn*mag_e12;
                
                %   2. Right node
                wL = node(n2).w;    wR = wL;    wR(4) = Winf(4); %Fix the pressure
                [num_flux,wsn] = fluxfunc2d(wL,wR,n12,fType);
                bfluxR = num_flux;
                node(n2).wsn = node(n2).wsn + wsn*mag_e12;
        end
        %   3. Add contributions to the two nodes (See Nishikawa AIAA2010-5093)
        switch elm( bound(i).belm(j) ).nV
            case 3	% Triangle
                node(n1).res = node(n1).res + (5*bfluxL+bfluxR)/6*mag_e12;
                node(n2).res = node(n2).res + (5*bfluxR+bfluxL)/6*mag_e12;
            case 4	% Quad
                node(n1).res = node(n1).res + bfluxL.mag_e12;
                node(n2).res = node(n2).res + bfluxR.mag_e12;
            otherwise
                error('Element is neither tria nor quad. Stop. ');
        end
    end
end