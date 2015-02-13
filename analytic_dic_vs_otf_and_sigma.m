
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



%%%%% Begin editable Section %%%%%

k_gs = linspace(0.001,2*pi ./ 0.9E-3,16); % grating wavevectors in 1/mm
%k_gs = linspace(0.001,2*pi ./ 0.9E-3,128);

NA = 0.6; % objective NA
lambda = 550E-6; % wavelength in mm
sigmas = linspace(0,1,16); % relative pupil filling, degree of coherence according to van Cittert-Zernike
%sigmas = linspace(0,1,128);

dw = 1/8 * lambda / NA; % Wollaston prism shift

t0 = 0.5 + 0.5i; % zero order grating transmission
a = 0.5 - 0.5i;  % grating modulation amplitude

% Zernike Fringe coefficients of objectives in waves
cveks = zeros(9,9); % (obj#, zernike#)
cveks(1,1) = .25;
cveks(2,2) = .25;
cveks(3,3) = .25;
cveks(4,4) = .25; % defocus
cveks(5,5) = .25; % astigmatism
cveks(6,6) = .25; % astigmatism
cveks(7,7) = .25; % coma
cveks(8,8) = .25; % coma
cveks(9,9) = .25; % spherical aberration

% Assumptions:
%  - condensor illuminates sample with a set of coherent plane waves
%  - condensor provides same amplitude for every plane wave
%  - perfect Wollaston prisms, no aberrations introduced
%  - both Wollaston prism and grating are oriented along the y-axis
%  - no apodisation
%  - phase difference between both paths extincts DIC signal of a uniform sample

%%%%% End editable Section %%%%%







k_obj_max = 2*pi*NA / lambda; % max transmitted in-plane wavevector by objective stop


% We raster grating constant, coherence, and objective aberrations
ng = length(k_gs);                  % number of gratings
nsig = length(sigmas);              % number of sigmas
nobjectives = length(cveks(:,1));   % number of objectives


% Visibilities
vis_plasdic_1 = zeros( nobjectives, ng, nsig );
vis_plasdic_2 = zeros( nobjectives, ng, nsig );
vis_dic_2     = zeros( nobjectives, ng, nsig );

tic

for isig = 1:nsig
    sigma = sigmas(isig);
    disp (['    calculating sigma = ', num2str(sigma)]);
    
    % determine required numer of illumination points 
    nillu_1d = round( 60 * sigma^2 + 20 );
    
    % grid of illumination plane wave wavevectors
    killumax = k_obj_max * sigma; % max wavevector in sigma spectrum
    [killux, killuy] = meshgrid(  linspace(-killumax,killumax,nillu_1d)  );
    clear nillu_1d;
    
    % The condensor sigma-aperture cuts too large illumination angles
    ind = (( killux.^2 + killuy.^2 ) <= killumax.^2 );
    killux = killux(ind);
    killuy = killuy(ind);

    nillu = length(killux); % number of illumination sources

    for ig = 1:ng
        k_g = k_gs(ig);   

        for iob = 1:nobjectives
            % calculate the OTF values 

            %OTF(killu)
            xp = killux * lambda / (2*pi*NA); % pupil coordinates in [-1,1]
            yp = killuy * lambda / (2*pi*NA); % pupil coordinates in [-1,1]
            OTF_killu = otf_calc(cveks(iob,:), xp, yp );
            
            % OTF(killu - k_g)
            xp = killux * lambda / (2*pi*NA);
            yp = (killuy-k_g) * lambda / (2*pi*NA);
            OTF_killu_minus_k_g = otf_calc(cveks(iob,:), xp, yp );

            % OTF(killu + k_g)
            xp = killux * lambda / (2*pi*NA);
            yp = (killuy+k_g) * lambda / (2*pi*NA);
            OTF_killu_plus_k_g = otf_calc(cveks(iob,:), xp, yp );

            
            
            % MTF = |OTF|
            MTF_killu = abs(OTF_killu);
            MTF_killu_plus_k_g = abs(OTF_killu_plus_k_g);
            MTF_killu_minus_k_g = abs(OTF_killu_minus_k_g);

            % amplitudes of PlasDIC intensities at [0,1,-1,2] * k_g + kil
            A0_PlasDIC =   4   * abs(t0)^2 * ( sin(killuy*dw) ).^2 .* ( MTF_killu ).^2 + ...
                 + 0.5 * abs(a)^2  * 4 * ( sin((k_g+killuy)*dw) ).^2 .* ( MTF_killu_plus_k_g ).^2 + ...
                 + 0.5 * abs(a)^2  * 4 * ( sin((killuy-k_g)*dw) ).^2 .* ( MTF_killu_minus_k_g ).^2 ;

            A1_PlasDIC      = -4 * t0 * conj(a) * sin(killuy*dw) .* sin((k_g-killuy)*dw) .* OTF_killu .* conj(OTF_killu_minus_k_g);
            Aminus1_PlasDIC = -4 * t0 * conj(a) * sin(killuy*dw) .* sin((k_g+killuy)*dw) .* OTF_killu .* conj(OTF_killu_plus_k_g);

            A2_PlasDIC = 2* abs(a)^2 * ( cos(2*killuy*dw) - cos(2*k_g*dw) ) .* OTF_killu_plus_k_g .* conj(OTF_killu_minus_k_g);

            % amplitudes of DIC intensities at [0,2] * k_g + kil
            A0_DIC = 0.5 * abs(a)^2 * ( sin(k_g*dw) )^2 * ( MTF_killu_plus_k_g.^2 +  MTF_killu_minus_k_g.^2 );
            A2_DIC =       abs(a)^2 * ( sin(k_g*dw) )^2 * OTF_killu_plus_k_g .* conj(OTF_killu_minus_k_g);

            %visibility of first and second grating harmonic in intensity image
            vis_plasdic_1(iob,ig,isig) = abs( sum( A1_PlasDIC- conj(Aminus1_PlasDIC)) ) / abs( sum(A0_PlasDIC) );
            vis_plasdic_2(iob,ig,isig) = abs( sum( A2_PlasDIC )               ) / abs( sum(A0_PlasDIC) );
            vis_dic_2(iob,ig,isig) = abs( sum( A2_DIC )                       ) / abs( sum(A0_DIC) );
        end
    end
end

toc


% Display results
for i = 1:nobjectives
    
    figure(i)
    imagesc(sigmas,k_gs,squeeze(vis_dic_2(i,:,:) ) ,[0,1])
    title(['DIC Visibility of second grating period, objective' num2str(i) ]);
    ylabel('grating wavevector ( 1 / mm )')
    xlabel('sigma')
    colormap('hot')
    colorbar()
    cm=get(gcf,'colormap');
    cmneu=zeros(200,3);
    cmneu(:,1)=interp1(1:length(cm),cm(:,1),linspace(1,length(cm),200));
    cmneu(:,2)=interp1(1:length(cm),cm(:,2),linspace(1,length(cm),200));
    cmneu(:,3)=interp1(1:length(cm),cm(:,3),linspace(1,length(cm),200));
    colormap(cmneu)
    
    figure(i + nobjectives + 100)  
    imagesc(sigmas,k_gs, squeeze(vis_plasdic_2(i,:,:)), [0,1])
    title(['PlasDIC Visibility of second grating period, objective' num2str(i) ]);
    ylabel('grating wavevector ( 1 / mm )')
    xlabel('sigma')
    colormap('hot')
    colorbar()
    cm=get(gcf,'colormap');
    cmneu=zeros(200,3);
    cmneu(:,1)=interp1(1:length(cm),cm(:,1),linspace(1,length(cm),200));
    cmneu(:,2)=interp1(1:length(cm),cm(:,2),linspace(1,length(cm),200));
    cmneu(:,3)=interp1(1:length(cm),cm(:,3),linspace(1,length(cm),200));
    colormap(cmneu)

    figure(i + 2* nobjectives + 1000)  
    imagesc(sigmas,k_gs, squeeze(vis_plasdic_1(i,:,:) ),[0,1] )
    title(['PlasDIC Artefact Visibility of first grating period, objective' num2str(i) ]);
    ylabel('grating wavevector ( 1 / mm )')
    xlabel('sigma')
    colormap('hot')
    colorbar()
    cm=get(gcf,'colormap');
    cmneu=zeros(200,3);
    cmneu(:,1)=interp1(1:length(cm),cm(:,1),linspace(1,length(cm),200));
    cmneu(:,2)=interp1(1:length(cm),cm(:,2),linspace(1,length(cm),200));
    cmneu(:,3)=interp1(1:length(cm),cm(:,3),linspace(1,length(cm),200));
    colormap(cmneu)

    
end
