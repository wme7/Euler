function va = vanAlbada(da,db,h)
%**************************************************************************
% -- vanAlbada Slope Limiter Function--
%
% 'A comparative study of computational methods in cosmic gas dynamics', 
% Van Albada, GD, B. Van Leer and W.W. Roberts, Astronomy and Astrophysics,
% 108, p76, 1982
% -------------------------------------------------------------------------
%  Input:   da, db: two differences
% Output:   va:     limited difference
% -------------------------------------------------------------------------
%**************************************************************************
eps2=(0.3*h)^3; va=0.5*(sign(da)*sign(db)+1)*...
    ( (db^2+eps2)*da + (da^2+eps2)*db )/(da^2+db^2+2*eps2);