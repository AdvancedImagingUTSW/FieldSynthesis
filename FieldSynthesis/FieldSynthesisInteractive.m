function [ hfig ] = FieldSynthesisInteractive( mask, doshift )
%FieldSynthesisInteractive Create an interactive line scan demonstration of
%field synthesis
%
% INPUT
% mask - mask at the pupil, which is the Fourier transform of electrical
%        field at the focal plane
% doshift - if true, shift the Fourier transform of the mask so the first
%           pixel is in the center of the image rather than the upper left
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
% EXAMPLE
% FieldSynthesisInteractive; % default demonstration with cameraman
% FieldSynthesisInteractive(createAnnulus(),true); % demonstrate a Bessel beam 
%
% Mark Kittisopikul , August 2018
% Goldman Lab
% Northwestern University

if(nargin < 1)
    mask = fftshift(fft2(double(imread('cameraman.tif'))));
end
if(nargin < 2)
    doshift = false;
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

%% Modulate the mask so that the image at the focal plane is centered
if(doshift)
    shifter = zeros(size(mask));
    shifter(ceil(size(mask,1)/2+1),ceil(size(mask,1)/2+1)) = 1;
    mask = mask .* fft2(shifter);
end

mask_unshifted = ifftshift(mask);
startCol = find(any(mask,1),1);
nCols = find(any(mask,1),1,'last');
% nCols = size(mask,2);
F = ifft2(mask_unshifted);

hfig = figure;

% Mask at pupil
subplot(2,3,1);
himMask = imshow(dispLogAbs2PlusOne(mask),[]);
title('Mask: $\log(|\hat{F}|^2+1)$','interpreter','latex');
hline = patch([1 1],[1 size(mask,2)],0,'EdgeColor','g','EdgeAlpha',0.5);

% Focal plane
subplot(2,3,2);
imshow(abs(F).^2,[]);
title('Intensity: $|F|^2$','interpreter','latex');

% Smeared
subplot(2,3,3)
smeared = sum(abs(F).^2,2);
imshow(repmat(smeared,1,size(F,2)),[]);
hold on;
plot(mat2gray(smeared)*(size(mask,2)-1)/2+1,1:size(mask,1),'m')
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
title('Cum. Intensity: $\sum^a |T_a|^2$','interpreter','latex');

% 1D delta function, used for line scan display
delta = zeros(size(mask,2));
delta(:,1) = 1;

% Cumulative matrix
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
        hline.XData = [a a];
        Ta = ifft2(ifftshift(mask.*circshift(delta,a-1,2)));
        Ta_sqmod = abs(Ta).^2;
        set(hreal,'CData',mat2gray(real(Ta)));
        set(hsqmod,'CData',mat2gray(Ta_sqmod));
        set(hsqmod_line,'XData',mat2gray(Ta_sqmod(:,1))*(size(mask,2)-1)/2+1);
        cumulative = cumulative + Ta_sqmod;
        set(hcumulative_line,'XData', ...
            mat2gray(cumulative(:,1))*(size(mask,2)-1)/2+1);
        set(hcumulative,'CData',mat2gray(cumulative));
    end
    function resetCumulative(source,event)
        % Zero out the cumulative matrix
        cumulative = zeros(size(mask));
        set(hcumulative_line,'XData', ...
            mat2gray(cumulative(:,1))*(size(mask,2)-1)/2+1);
        set(hcumulative,'CData',mat2gray(cumulative));
    end

end

