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

