function [ annulus ] = createAnnulus( n, r, w )
%createAnnulus Create a binary annular mask
%
% INPUT (all optional)
% n - size of the annular mask as a scalar, or vector with coordinates
% r - radius of the annulus in pixels
% w - width of the annulus in pixels
%
% OUTPUT
% annulus - n x n matrix with the annulus marked with ones
%
% USAGE
% figure;
% imshow(createAnnulus(256,32,4),[]);
%
% Create Bessel beam 2D profile
% figure;
% imshow(log(abs(fftshift(ifft2(ifftshift(createAnnulus)))).^2+1),[]);
% colormap(gca,hot);
% caxis([0 6e-4]);
%
% REMARKS
% This could be streamlined using the bresenham circle algorithm

% Mark Kittisopikul, August 25th, 2018
% Lab of Robert D. Goldman;
% Northwestern University

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


if(nargin < 1)
    n = 256;
end
if(nargin < 2)
    r = 32;
end
if(nargin < 3)
    w = 4;
end

if(isscalar(n))
    v = 1:n;
    % zeroth order coefficient is at n/2+1,n/2+1 due to fftshift/ifftshift
    v = v-floor(n/2)-1;
else
    % non-scalar given. Use n as coordinates
    v = n;
end

% Calculate radial position in polar coordinate system
% Pre-bsxfun expansion code (pre 2017a):
[Y,X] = meshgrid(v,v);
Q = hypot(X,Y);

% Bsxfun expansion code (post-2017a)
% Q = hypot(v,v.');

% Create an annulus with radius r and width w
annulus = abs(Q -r) < w;

end

