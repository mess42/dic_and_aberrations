% This program calculates the aberration dependence of contrast in DIC microscopy
% Copyright (C) 2015 Moritz Esslinger
% moritz.esslinger@web.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
% 





    function otf = otf_calc( cvek, xp_mesh, yp_mesh )
        % Calculates the OTF value at pupil position(s) xp_mesh, yp_mesh.
        %
        % Arguments:
        % cvek    : List of Zernike coefficents (1 x nzern 2d array of float)
        % xp_mesh : List of pupil x coordinates (1d array of float)
        % yp_mesh : List of pupil y coordinates (1d array of float)
        %           Pupil coordinates are normalized from -1 to 1.
        %           Larger values will be cut by the aperture.
        %
        %
        % Return:
        % otf : optical transfer function 
        %       at pupil position(s) xp_mesh, yp_mesh
        %       (1d array of complex)
    
        s = size(xp_mesh);
        npx = s(1);
        npy = s(2);
        nzern = length(cvek);
        
        apodisation = ones(npx,npy);
        
        % get all Zernike Fringe polynomials on the xp_mesh, yp_mesh grid
        zernike_fringes = zeros([ npx, npy, nzern ]);
        for i = 1:nzern
            zernike_fringes(:, :, i) = zernike(i, xp_mesh, yp_mesh);
        end
     
        % multiply Zernike coefficients with polynomials
        wavefront = matrix_multiplication(zernike_fringes, cvek, 3,2);     
        otf = apodisation .* exp( 2i * pi * wavefront ) ;
         
        % diffraction limit
        otf(xp_mesh.^2 + yp_mesh.^2 > 1) = 0;
    end
