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






function fringe = zernike(jzern, xp, yp)
% Returns a single Zernike Fringe polynomial on the mesh xp,yp.
% 
% Arguments:
% jzern : Zernike Fringe index (int)
% xp    : List of x coordinates (1d array of float)
% yp    : List of y coordinates (1d array of float)
%
%
% Return:
% fringe : Zernike Fringe coefficient values of term number jzern at
%          positions xp, yp
%         (1d array of float)






% fringe = RHO * sin(omega * phi)
% or
% fringe = RHO * cos(omega * phi)
%
% RHO(r) = sum_{n=1}^{N} ( R(n) * r^(n-1) )


% get Zernike Standard m,n indices
next_sq = (ceil(sqrt(jzern)))^2; % next larger square number

m_plus_n = 2*sqrt(next_sq) - 2;

m = ceil( (next_sq - jzern)/2);
n = m_plus_n - m;


% Zernike Radial Part R
R = zeros(1,n+1); 
for a = m:2:n
    k = (n-a)/2;
    R(a+1) = (-1)^k * factorial(n-k) / ...
        ( factorial(k) * factorial((n+m)/2 - k) * factorial((n-m)/2-k) );
end
omega = m;

% numerically apply on xp,yp mesh
r = sqrt( xp.^2 + yp.^2 );
phi = atan2(xp,yp);

RHO = 0 .* r;
for n = 1:length(R)
    RHO = RHO + R(n) * r.^(n-1);
end

if mod(next_sq - jzern,2) == 1 % sin
    fringe = RHO .* sin( omega * phi);
else % cos
    fringe = RHO .* cos( omega * phi);
end
