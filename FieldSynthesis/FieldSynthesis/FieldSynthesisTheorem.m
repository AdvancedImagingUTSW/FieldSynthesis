function [efield,slice,smear,Q,T] = FieldSynthesisTheorem(efield)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Small program to illustrate a new Field Synthesis Theorem.
%
%   In essence it says that the projection of the absolute modulus of a
%   complex field is the same as when one takes a sliding window in the
%   Fourier domain, makes an inverse FFT of each slice, take the absolute
%   modulus of that and sum it up while moving the window through the
%   spectrum.  This has important applications for scanned light-sheets and
%   how to generate them.
%
%   Reto Fiolka, May 2017
%   Mark Kittisopikul, May 2017 - Aug 2018
%
%   INPUT
%   efield - electric field at the focal plane, may be real or complex
%            valued
%   
%   OUTPUT
%   efield - electric field at the focal plane
%   slice  - intensity of illumination pattern by field synthesis
%   smear  - intensity of illumination pattern by dithering
%   Q      - Fourier transform of individual line scan without phasing,a=10
%   T      - Fourier transform of individual line scan with phasing,a=10
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

%% Initialize input
% assume this is the e-field in real space (real or complex valued)
if(nargin < 1)
    % Use stock image as an arbitrary efield
    efield = imread('cameraman.tif');
    imaginaryAmplitude = 1;
    efield=double(efield)+1j*rand(size(efield))*imaginaryAmplitude;
end

% efield is F(x,y) in the proof and can be complex valued
% get efield properties
sz = size(efield);
N = sz(2);

% then this is the spectrum of the efield
% spectrum is F_hat(k_x,k_y) in the proof
spectrum_unshifted = fft2(efield);
spectrum=fftshift(spectrum_unshifted);
% the intensity is the modulus squared of the efield
% this is a real valued image         
I=abs(efield).^2;                         
% slice is an array that is created by superposition of
% inverse FFTs of spectral slices 
slice=zeros(sz);
% smear is created by scanning the image in real space                  
smear=zeros(sz);                      
 
%% Do line scan
for a=1:N
    % take one slice of spectrum and take inverse FFT
    T_hat=zeros(sz);
    T_hat(:,a)=spectrum(:,a);
    T=ifft2(ifftshift(T_hat));
    % superimpose intensities (modulus squared) of
    % inverse FFTs of spectral slices          
    slice=slice+abs(T).^2;           
    % smearing the image in x-direction
    % intensity is superimposed at every position
    smear=smear+(circshift(I,[0,a]));
end

% We could also create smear as follows:
% smear = repmat(sum(I,2),1,256);

%% Line profiles of the smeared intensity images
 
sliceProfile=slice(:,1);
 
smearProfile=smear(:,1)/N;
 
%alternatively, one can also just project the intensity image I in x
projectionProfile=sum(I,2)/N;

% The profile can be calculated directly by taking a 1 dimensional fourier
% transform
oneDFT = ifft(spectrum_unshifted)/N;
% This line is equivalent to the above
% oneDFT = fft(efield,[],2)/N;
oneDFT = abs(oneDFT).^2;
oneDFT = sum(oneDFT,2);

 
%% Plot Figures;
figure;
subplot(2,3,1);
imshow(I,[]);
title('Original Intensity');
subplot(2,3,2);
imshow(slice,[]);
title('Slice');
subplot(2,3,3);
imshow(smear,[]);
title('Smear (~Dither)');

 
% All profiles are identical
subplot(2,3,4:6);
plot(projectionProfile,'ko','DisplayName','Projection')
hold on;
xlim([1 256]);

plot(sliceProfile,'b+','DisplayName','Slice');

plot(smearProfile,'rx','DisplayName','Smear');

plot(oneDFT,'g.','DisplayName','1D FT');

grid on;
legend('show','Location','southwest');
title('Vertical Profiles are the Same');
xlabel('z position');

disp('Press any key');
pause;

%% Explanation of the profile of individual line scans
% T represents a selected column in the spectral field selected by the scan

efield_xft = fft(efield,[],2);
 
% for k=1:N

a = 10;
hfig = figure('units','normalized','outerposition',[0 0 1 1]);

% T is constructed similarly to above in the for loop
% The only difference is how k is indexed in the unshifted spectrum
T_hat=zeros(sz);
T_hat(:,a)=spectrum_unshifted(:,a);
T = ifft2(T_hat);

% Q differs from T because the selected column is copied into k_x = 0
Q_hat = zeros(N);
Q_hat(:,1) = spectrum_unshifted(:,a);
Q = ifft2(Q_hat);

subplot(4,3,1);
imshow(real(Q),[]);
title('Real(Q)');

subplot(4,3,2);
imshow(imag(Q),[]);
title('Imag(Q)');

subplot(4,3,3);
imshow(abs(Q).^2,[]);
title('abs(Q)^2');

subplot(4,3,4);
imshow(real(T),[]);
title('Real(T)');

subplot(4,3,5);
imshow(imag(T),[]);
title('Imag(T)');

subplot(4,3,6);
imshow(abs(T).^2,[]);
title('abs(T)^2');

subplot(4,3,7:9);
plot(real(Q(:,1))*N,'r-');
hold on;
plot(real(efield_xft(:,a)),'ro');
 
plot(imag(Q(:,1))*N,'b-');
plot(imag(efield_xft(:,a)),'bo');
 
legend({'real(Q)','Real 1D FFT of E','imag(Q)','Imag 1D FFT of E'}, ...
    'Location','southoutside','Orientation','horizontal');
title('Q is the Fourier Transform of a Line Scan in the Spectrum');
xlim([1 N]);
xlabel('z position (pixels)');

 
% pause(1);
 
subplot(4,3,10:12);
plot(real(exp(1i*2*pi*(a-1)/N.*(0:N-1))),'r-');
hold on;
plot(real(T(1,:)./Q(1,:)),'ro');
plot(imag(exp(1i*2*pi*(a-1)/N.*(0:N-1))),'b-');
plot(imag(T(1,:)./Q(1,:)),'bo');
title(['T is Q With Complex Modulation ' ...
    'Due to the Location of the Line Scan']);
xlabel('z position (pixels)');
ylabel('T/Q');
xlim([1 N]);


 
% pause(1);
 
% close(hfig);
% end;

end
