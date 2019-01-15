function [F,fieldSynthesis,dithered,Q,T] = FieldSynthesisTheorem(F,L)
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
%   F - electric field at the front focal plane, may be real or complex
%            valued
%   
%   OUTPUT
%   F - electric field at the front focal plane
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
if(nargin < 1 || isempty(F))
    % Use stock image as an arbitrary efield
    F = imread('cameraman.tif');
    imaginaryAmplitude = 1;
    F=double(F)+1j*rand(size(F))*imaginaryAmplitude;
end
if(nargin < 2 || isempty(L))
    L = ones(1,size(F,1));
end

% F(x,y) is the eletric field and can be complex valued
% get F properties. F was previously called efield
sz = size(F);
N = sz(2);
center = floor(N/2+1);

% then this is the spectrum of the efield
% spectrum is F_hat(k_x,k_y) in the proof
F_hat_unshifted = fft2(F);
F_hat=fftshift(F_hat_unshifted);
% the intensity is the modulus squared of the efield
% this is a real valued image         
F_sqmod=abs(F).^2;                         
% slice is an array that is created by superposition of
% inverse FFTs of spectral slices 
fieldSynthesis=zeros(sz);
% smear is created by scanning the image in real space                  
dithered=zeros(sz);                      

%% Calculate line scan profile at back focal plane
L_hat_unshifted = fft(ifftshift(L));
L_hat = fftshift(L_hat_unshifted);
L_hat = repmat(L_hat,size(F,1),1);

%% Calculate convolution kernel for front focal plane
L_sqmod = zeros(size(F));
L_sqmod(center,:) = abs(L).^2;
 
%% Do line scan
for a=1:N
    
    % If L were a delta function, then this would work
    % take one slice of spectrum and take inverse FFT
    % T_hat=zeros(sz);
    % T_hat(:,a)=F_hat(:,a);
    
    T_hat = F_hat.*circshift(L_hat,[0 a-center]);

    T=ifft2(ifftshift(T_hat));
    % superimpose intensities (modulus squared) of
    % inverse FFTs of spectral slices          
    fieldSynthesis=fieldSynthesis+abs(T).^2;           
    % smearing the image in x-direction
    % intensity is superimposed at every position
    dithered=dithered+(circshift(F_sqmod,[0,a-center])).*L_sqmod(center,a);
end



%% Line profiles of the smeared intensity images
 
fieldSynthesisCentralZSlice=fieldSynthesis(:,center);
 
ditheredCentralZSlice=dithered(:,center);
 
%alternatively, one can also just project the intensity image I in x
% if L_sqmod were delta(z)
% projectionProfile=sum(F_sqmod,2)/N;
% more generally, this can be done by circular convolution
smearByCircularConvolution = fftshift(ifft2(fft2(ifftshift(F_sqmod)).*fft2(ifftshift(L_sqmod))));


% The profile can be calculated directly by taking a 1 dimensional fourier
% transform
F_unshifted = ifftshift(F,2);
if(mod(N,2) == 0)
    % even
    L_flipped_unshifted = [L(1) fliplr(L(2:end))];
else
    % odd
    L_flipped_unshifted = fliplr(L);
end
L_flipped_unshifted = repmat(L_flipped_unshifted,N,1);
L_flipped_unshifted = ifftshift(L_flipped_unshifted,2);

oneDFT = fftshift(fft(F_unshifted.*L_flipped_unshifted,[],2),2);
oneDFT = abs(oneDFT).^2;
oneDFT = sum(oneDFT,2);

 
%% Plot Figures;
figure;

subplot(2,3,1);
imshow(F_sqmod,[]);
title('Original Intensity');

subplot(2,3,2);
imshow(fieldSynthesis,[]);
title('Field Synthesis');

subplot(2,3,3);
imshow(dithered,[]);
title('Smear / Dither');

 
% All profiles are identical
subplot(2,3,4:6);
plot(smearByCircularConvolution(:,center),'ko','DisplayName','Projection')
hold on;
xlim([1 N]);

plot(fieldSynthesisCentralZSlice/N,'b+','DisplayName','Field Synthesis');

plot(ditheredCentralZSlice,'rx','DisplayName','Dithered');

plot(oneDFT/N,'g.','DisplayName','1D FT');

grid on;
legend('show','Location','southwest');
title('Vertical Z Profiles are the Identical');
xlabel('z position');

% disp('Press any key');
% pause;
return;

%% Explanation of the profile of individual line scans
% T represents a selected column in the spectral field selected by the scan

% efield_xft = fft(F,[],2);
 
% for k=1:N

a = center-10;
hfig = figure('units','normalized','outerposition',[0 0 1 1]);

% T is constructed similarly to above in the for loop
% The only difference is how k is indexed in the unshifted spectrum
% T_hat=zeros(sz);
% T_hat(:,a)=F_hat_unshifted(:,a);
T_hat = F_hat_unshifted.*circshift(L_hat_unshifted,[0 a-center]);
T = ifft2(T_hat);

% Q differs from T because the selected column is copied into k_x = 0
% Q_hat = zeros(N);
% Q_hat(:,1) = F_hat_unshifted(:,a);
Q_hat = circshift(T_hat,[0 center-a]);
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
% plot(real(efield_xft(:,a)),'ro');
 
plot(imag(Q(:,1))*N,'b-');
% plot(imag(efield_xft(:,a)),'bo');
 
% legend({'real(Q)','Real 1D FFT of E','imag(Q)','Imag 1D FFT of E'}, ...
%     'Location','southoutside','Orientation','horizontal');
legend({'real(Q)','imag(Q)'}, ...
     'Location','southoutside','Orientation','horizontal');
title('Q is the Fourier Transform of a Line Scan in the Spectrum');
xlim([1 N]);
xlabel('z position (pixels)');

 
% pause(1);
 
subplot(4,3,10:12);
% plot(real(exp(1i*2*pi*(a-1)/N.*(0:N-1))),'r-');
hold on;
plot(real(T(center,:)./Q(center,:)),'ro');
% plot(imag(exp(1i*2*pi*(a-1)/N.*(0:N-1))),'b-');
plot(imag(T(center,:)./Q(center,:)),'bo');
title(['T is Q With Complex Modulation ' ...
    'Due to the Location of the Line Scan']);
xlabel('z position (pixels)');
ylabel('T/Q');
xlim([1 N]);


 
% pause(1);
 
% close(hfig);
% end;

end
