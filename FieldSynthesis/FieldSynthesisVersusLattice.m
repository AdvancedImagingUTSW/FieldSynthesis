function [varargout] = FieldSynthesisVersusLattice(n,w,r,offset,dispRange)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%Simulation for field synthesis
%
%   compares field synthesis vs square lattice
%
%
%   Reto, May 2017
%   Mark Kittisopikul, August 2018
%
%   INPUT
%   n - Defines the size of the image and mask to be n x n
%   w - Width of the mask components
%   r - Radius of the annulus (width is centered on the annulus)
%   offset - Offset of the side components of the square lattice
%   dispRange - Set which part of mask to display in figures
%
%   OUTPUT
%   out - struct containing workspace of this function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Field Synthesis Demonstration -
% MATLAB code to demonstrate field synthesis light sheet microscopy
% Copyright (C) 2018 Reto Fioka,
%               University of Texas Southwestern Medical Center
% Copyright (C) 2018 Mark Kittisopikul,
%               Northwestern University
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.

%% Parameters
% Size of the mask
if(nargin < 1)
    n=4096;
end
% Width of the annulus
if(nargin < 2)
    w=5;
end
% Radius of the annulus
if(nargin < 3)
    r=256;
    % r = 200
end
% Offset for side slits
if(nargin < 4)
    offset = r;
    % offset=198;
end
% Display range
if(nargin < 5)
    dispRange = (-600:600)+floor(n/2)+1;
end

%% Create clean annulus
% We do not need to initialize
% since we will create the matrix with createAnnulus
% annulus = zeros(n);

% Vector for x and y, which should be symmetric
v = 1:n;
% zeroth order coefficient is at n/2+1,n/2+1 due to fftshift/ifftshift
v = v-floor(n/2)-1;

% Create an annulus of radius r with width w centered in an n x n matrix
annulus = createAnnulus(v, r, w);

% Select columns for mask
abs_v = abs(v);
% Select three sets of frequency columns
% 1) Group of columns centered on the offset to the left of width w
% 2) Group of columns in the center of width w
% 3) Group of columns centered on the offset to the right of width w
selected_columns = (abs_v < offset+w/2 & abs_v > offset-w/2) | ...
    (v < w/2 & v > -w/2);

% Remove unselected columns from mask
latticeFourierMask = annulus;
latticeFourierMask(:,~selected_columns) = false;
latticeFourierMask = double(latticeFourierMask);
% latticeFourierMask is now the Fourier mask of a square lattice

%% Field Synthesis

% The field synthesis is process is equivalent to summing over a
% 1D Fourier Transform of the mask
% 1) Shift so the 0th frequency is at 1,1
% 2) Do the 1D inverse FT
% 3) Shift so the center pixel is the center of the image
fieldSynthesisProfile = fftshift(ifft(ifftshift(latticeFourierMask)));
fieldSynthesisProfile = sum(abs(fieldSynthesisProfile).^2,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lattice simulation

% The electric field of lattice is the 2D Fourier Transform of the mask
lattice_efield=fftshift(ifft2(ifftshift(latticeFourierMask)));
% Take the square modulus to get the intensity
lattice=abs(lattice_efield).^2;
% Perform the dithering operation
latticeLineProfile=sum(lattice,2);
% Scale by n, due ifft2 normalization
latticeLineProfile=latticeLineProfile*n;

%% Plot: Compare lattice profile to field synthesis profile
figure;
% Show the convention lattice profile
subplot(3,1,1);
plot(dispRange-n/2+1,latticeLineProfile(dispRange));
xlim([min(dispRange) max(dispRange)]-n/2+1);
title('Conventional Lattice Profile');

subplot(3,1,2);
% Show the field synthesis profile
plot(dispRange-n/2+1,fieldSynthesisProfile(dispRange),'r')
xlim([min(dispRange) max(dispRange)]-n/2+1);
title('Field Synthesis Profile');

subplot(3,1,3);
% Compare the two profiles
plot(dispRange-n/2+1,latticeLineProfile(dispRange));
hold on;
plot(dispRange-n/2+1,fieldSynthesisProfile(dispRange),'r--')
xlim([min(dispRange) max(dispRange)]-n/2+1);
title('Comparison of Lattice and Field Synthesis Profiles');

%% Analysis of all interference patterns in lattice

% lattice is the intensity of the pattern as per above
% lattice=abs(B).^2;
%Fourier transform of lattice
lattice_hat=fftshift(fft2(ifftshift(lattice)));


%% Dithering lattice: lattice pattern is shifted by subpixel steps and added
% Calculate time average by dithering over the period
period = n/offset;

% To dither, we average over one period of the lattice by shifting
if(period == round(period))
    % The shifting operation can be done via a 2D convolution
    latticeDithered = conv2(lattice,ones(1,period)/period,'same');
    % % The following block of code is equivalent to the above line
    % latticeDithered = zeros(size(lattice));
    % for s=floor(-period/2):floor(period/2)-1
    %     latticeDithered = latticeDithered + circshift(lattice,s,2);
    % end
    % latticeDithered = latticeDithered / period;
else
    % Above, we assume that the period is of integer units.
    % If it were not of integer units, we can use the following code
    % Use the convolution theorem to do convolution in Fourier space
    latticeDithered = bsxfun(@times,lattice_hat,sinc(v/period));
    latticeDithered = fftshift(ifft2(ifftshift(latticeDithered)));
    % % We could also approximate the the dithering via subpixel steps
    % subpixelFactor = 1/(period-floor(period));
    % subpixelFactor = ceil( subpixelFactor );
    % subpixelFactor = min(subpixelFactor,10);
    % period = floor(period*subpixelFactor);
    % latticeDithered = conv2( interpft(lattice,n*subpixelFactor,2), ...
    %                         ones(1,period)/period,'same');
    % latticeDithered = interpft(latticeDithered,n,2);
end

%Fourier transform of dithered lattice
latticeDithered_hat=fftshift(fft2(ifftshift(latticeDithered)));

%% Plot 2x3 

h = figure;
% Make figure full screen
set(h,'Units','normalized','Position',[0 0 1 1]);

% Show the mask
subplot(2,3,1)
imshow(latticeFourierMask(dispRange,dispRange),[0 1e-6]);colormap hot
title('Electric field in pupil');

% Show the Fourier transform of the _intensity_ of the lattice
subplot(2,3,2)
imshow(abs(lattice_hat(dispRange,dispRange)),[0 1e-6]);colormap hot
title('Fourier components of lattice intensity');

% Show the Fourier transform of the dithered lattice intensity
subplot(2,3,3)
imshow(abs(latticeDithered_hat(dispRange,dispRange)),[0 1e-6]);colormap hot
title('Fourier components of dithered lattice intensity');

% Show the electric field of the lattice at the focal plane
subplot(2,3,4);
imshow(lattice_efield(dispRange,dispRange),[]);
title('Electric field of lattice at focal plane');
% Zoom in so we can see the details of the lattice
xlim([-75 75]+length(dispRange)/2+1);
ylim([-75 75]+length(dispRange)/2+1);

% Show the intensity of the lattice
subplot(2,3,5)
imshow(lattice(dispRange,dispRange),[]);
title('Intensity of lattice');
% Zoom in so we can see the details of the lattice
xlim([-75 75]+length(dispRange)/2+1);
ylim([-75 75]+length(dispRange)/2+1);

% Show the Fourier transform of the dithered lattice intensity
subplot(2,3,6)
imshow(latticeDithered(dispRange,dispRange),[]);
title('Averaged Intensity of dithered lattice');
xlim([-75 75]+length(dispRange)/2+1);
ylim([-75 75]+length(dispRange)/2+1);

%% Output
if(nargout > 0)
    % If output is requested, pack workspace into a struct
    varnames = who;
    out = struct;
    for varIdx = 1:length(varnames)
        out.(varnames{varIdx}) = eval(varnames{varIdx});
    end
    varargout{1} = out;
end