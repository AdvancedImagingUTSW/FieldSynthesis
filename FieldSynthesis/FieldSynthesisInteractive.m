function [ hfig ] = FieldSynthesisInteractive( mask, doshift, lineProfile )
%FieldSynthesisInteractive Create an interactive line scan demonstration of
%field synthesis
%
% INPUT
% mask - mask at the pupil, which is the Fourier transform of electrical
%        field at the focal plane. zeroth frequency should be in the
%        middle. ifftshift will be applied for calcualtions.
% doshift - if true, shift the Fourier transform of the mask so the first
%           pixel is in the center of the image rather than the upper left
% lineProfile  -  line profile for the scan in the pupil mask
%                 EITHER:
%                 1) 0 for a delta function line scan
%                 2) a positive double value indicating the sigma of the
%                      gaussianLine in pixels
%                 3) a line profile vector the same width as mask. The main
%                    peak is expected to be in the center and ifftshift
%                    will be applied
%
% OUTPUT
% hfig - handle for the display figure
%
% INTERACTIVE
% The button in the lower left plays / pauses the movie.
% The arrow buttons on the slider will move the scan by one column.
% Clicking on the trough of the slider will move the scan by five columns.
% The button in the lower right labeled R will reset the cumulative view.
%
% DISPLAY
% The display consists of 6 panels
% 1 2 3
% 4 5 6
% 1. The pupil mask, |\hat{F}|^2 in log scale
% 2. The object domain, |F|^2, scanning left to right
%    Line plot indicates beam intensity
% 3. Dithered, averaged intensity. Cumulative sum of display #2
% 4. Display of the real component of the electric field of an insteaneous
%    scan, Real{T_a}
% 5. Instaneous scan intensity, |T_a|^2
% 6. Cumulative scan intensity of display #5
%
% EXAMPLE
% FieldSynthesisInteractive; % default demonstration with cameraman
% FieldSynthesisInteractive(createAnnulus(),true); % demonstrate a Bessel beam
% Create a sinc profile to emulate a scan over a finite range
% N = 128;
% x = -ceil(N/2):floor(N/2-1)
% L_hat = fftshift(fft(ifftshift(abs(x) < 30)));
% FieldSynthesisInteractive(createAnnulus(),true,L_hat);
%
% Mark Kittisopikul , August 2018
% Goldman Lab
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
    mask = fftshift(fft2(double(imread('cameraman.tif'))));
end
if(nargin < 2)
    doshift = false;
end
if(nargin < 3)
    % delta function
    lineProfile = 0;
end

%% Setup helper functions
% Show the intensity (square modulus) using a log transform
dispLogAbs2PlusOne = @(x) mat2gray(log(abs(x).^2 + 1));
% Do a 1D smear (dithering)
smear1D = @(I,dim) repmat(mean(I,dim),circshift([1 size(I,dim)],dim,2));
% Prepare Image for 1D DFT Display
disp1DFFT = @(I,dim) dispLogAbs2PlusOne(fftshift(fft(I,[],dim),dim));
% Prepare Image for 2D DFT Display
disp2DFFT = @(I) dispLogAbs2PlusOne(fftshift(fft2(I)));

%% Setup line profile
if(isscalar(lineProfile))
    % If lineProfile is scalar, interpret it to be gaussianLineSigma
    gaussianLineSigma = lineProfile;
    assert(gaussianLineSigma >= 0,...
        'lineProfile scalar must be nonnegative');
    
    % 1D delta function
    L_hat = zeros(size(mask));
    L_hat(:,1) = 1;

    % If lineProfile is 0, then line profile is delta function
    if(gaussianLineSigma > 0)
        L_hat = circshift(L_hat,[0 ceil(5*gaussianLineSigma)]);
        L_hat = imgaussfilt(L_hat,gaussianLineSigma);
        L_hat = circshift(L_hat,[0 -ceil(5*gaussianLineSigma)]);
    end
elseif(isvector(lineProfile))
    assert(length(lineProfile) == size(mask,2), ...
        'lineProfile vector must match size(mask,2)');
    % 1-D arbitrary line profile is provided
    L_hat = lineProfile;
    L_hat = ifftshift(L_hat);
    L_hat = repmat(L_hat(:).',size(mask,1),1);
else
    assert(all(size(lineProfile) == size(mask)), ...
        'lineProfile field must match mask size');
    % 2-D arbitrary line profile is provided
    L_hat = lineProfile;
    L_hat = ifftshift(L_hat);
end

L_hat_sqmod = abs(L_hat).^2;

L = ifft2(L_hat);
L_sqmod = fftshift(abs(L).^2);

center = floor(size(mask,2)/2+1);

%% Modulate the mask so that the image at the focal plane is centered
% if(doshift)
%     shifter = zeros(size(mask));
%     shifter(ceil(size(mask,1)/2+1),ceil(size(mask,1)/2+1)) = 1;
%     mask = mask .* fft2(shifter);
% end
% mask = ifftshift(mask);

mask_unshifted = ifftshift(mask);
% startCol = find(any(mask,1),1);
% nCols = find(any(mask,1),1,'last');

startCol = 1;
nCols = size(mask,2);
F = ifft2(mask_unshifted);
F_sqmod = abs(F).^2;

if(doshift)
    F_sqmod = fftshift(F_sqmod);
end

% smeared = sum(abs(F).^2,2);
% smeared = repmat(smeared,1,size(F,2));

% Use circulation convolution due to potential x-boundary effects with small
% sigma. Linear convolution (conv2) will work in large sigma case
cconv2 = @(A,B,~) fftshift(ifft2(fft2(ifftshift(A)).*fft2(ifftshift(B))));
smeared = cconv2(F_sqmod,L_sqmod);

%% Begin display

hfig = figure;

% Mask at pupil
subplot(2,3,1);
himMask = imshowpair(dispLogAbs2PlusOne(mask),zeros(size(mask)));
title('Mask: $\log(|\hat{F}|^2+1)$','interpreter','latex');
% hline = patch([1 1],[1 size(mask,2)],0,'EdgeColor','g','EdgeAlpha',0.5);

% Focal plane
subplot(2,3,2);
hFocalPlane = imshow(F_sqmod,[]);
title('Intensity: $|F|^2$','interpreter','latex');
hold on;
hlsqmod = plot(1:size(mask,2),-mat2gray(L_sqmod(center,:))*size(mask,1)/2+size(mask,1));

% Smeared
subplot(2,3,3)
hSmeared = imshow(smeared,[]);
hold on;
plot(mat2gray(smeared(:,center))*(size(mask,2)-1)/2+1,1:size(mask,1),'m')
title('Dithered Intensity: $\sum_x |F|^2$','interpreter','latex');

% Electric field at focal plane
subplot(2,3,4);
hreal = imshow(zeros(size(mask)),[]);
hreal_title = title('Electric field: $Real\{T_a\}$','interpreter','latex');

% Instaneous intensity at line scan
subplot(2,3,5);
hsqmod = imshow(zeros(size(mask)),[]);
hold on;
hsqmod_line =  plot(zeros(1,size(mask,1)),1:size(mask,1),'m');
hsqmod_title = title('Scan Intensity: $|T_a|^2$','interpreter','latex');

% Cumulative intensity of line scans
subplot(2,3,6);
hcumulative = imshow(zeros(size(mask)),[]);
hold on;
hcumulative_line = plot(zeros(1,size(mask,1)),1:size(mask,1),'m');
title('Cum. Intensity: $\sum_a |T_a|^2$','interpreter','latex');

% Cumulative matrix
cumulativeScan = zeros(size(mask));
cumulative = zeros(size(mask));

% Play button in the lower left
hplay = uicontrol('Style','togglebutton','Units','normalized', ...
    'Position',[0 0 0.05 0.05],'String','||','Value',1, ...
    'Callback',@toggleButton);
% Reset button in the lower right
hreset = uicontrol('Style','pushbutton','Units','normalized', ...
    'Position',[0.95 0 0.05 0.05],'String','R', ...
    'Callback',@resetCumulative);
% Slide control
hslider = uicontrol('Style','slider','Units','normalized', ...
    'Position',[0.05 0 0.90 0.05],'String','Scan Position', ...
    'Min',startCol,'Max',nCols,'Value',startCol, ...
    'SliderStep',[1 5]/nCols, ...
    'Callback',@updateSlider);
% Text label for slider
uicontrol('Style','text','Units','normalized', ...
    'Position',[0 0.05 1 0.05], ...
    'String','Scan Position','HorizontalAlignment','left');

% Axis annotation
annotation('textarrow','Color','m','Position',[0.1 0.17 0 0.1],'String','z');
annotation('textarrow','Color','m','Position',[0.15 0.12 0.1 0],'String','x');

% Play on start
play();

    function play()
        % Loop from min to max values
        for aa=round(get(hslider,'Value')):nCols
            set(hslider,'Value',aa);
            updateSlider(hslider,[]);
            if(~get(hplay,'Value'))
                break;
            end
            pause(0.1);
        end
        set(hplay,'Value',0);
        toggleButton(hplay,[]);
    end
    function toggleButton(source,event)
        % Toggle the play button
        switch(get(source,'Value'))
            case 0
                set(source,'String','>');
            case 1
                set(source,'String','||');
                % If at the end, reset the line scan position on play
                if(get(hslider,'Value') == get(hslider,'Max'))
                    set(hslider,'Value',get(hslider,'Min'));
                end
                play();
        end
    end
    function updateSlider(source,event)
        % Update the slider and the corresponding images
        a = round(source.Value);
%         hline.XData = [a a];
        Ta = ifft2(ifftshift(mask.*circshift(L_hat,a-1,2)));
        Ta_sqmod = abs(Ta).^2;
        if(doshift)
            Ta = fftshift(Ta);
            Ta_sqmod = fftshift(Ta_sqmod);
        end
        F_sqmod_scanning = circshift(F_sqmod,-center+a-1,2);
        
        pupil = get(himMask,'CData');
        pupil(:,:,1) = mat2gray(circshift(L_hat_sqmod,a-1,2))*255;
        pupil(:,:,3) = pupil(:,:,1);

        set(himMask,'CData',pupil);
        
        set(hFocalPlane,'CData',F_sqmod_scanning);
        cumulativeScan = cumulativeScan + F_sqmod_scanning.*L_sqmod(center,a);
        set(hSmeared,'CData',cumulativeScan);
        
        
        set(hreal,'CData',mat2gray(real(Ta)));
        set(hsqmod,'CData',mat2gray(Ta_sqmod));
        set(hsqmod_line,'XData',mat2gray(Ta_sqmod(:,center))*(size(mask,2)-1)/2+1);
        cumulative = cumulative + Ta_sqmod;
        set(hcumulative_line,'XData', ...
            mat2gray(cumulative(:,center))*(size(mask,2)-1)/2+1);
        set(hcumulative,'CData',mat2gray(cumulative));
    end
    function resetCumulative(source,event)
        % Zero out the cumulative matrix
        cumulative = zeros(size(mask));
        set(hcumulative_line,'XData', ...
            mat2gray(cumulative(:,1))*(size(mask,2)-1)/2+1);
        set(hcumulative,'CData',mat2gray(cumulative));
        
        cumulativeScan = zeros(size(mask));
        set(hSmeared,'CData',mat2gray(cumulativeScan));
%         set(hcumulative_line,'XData', ...
%             mat2gray(cumulative(:,1))*(size(mask,2)-1)/2+1);
%         set(hcumulative,'CData',mat2gray(cumulative));

    end

end

