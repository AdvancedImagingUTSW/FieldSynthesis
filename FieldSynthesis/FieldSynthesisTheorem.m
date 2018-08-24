%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Small program to illustrate a new Fourier slice theorem I discovered
%   accidentally. Analytical proof is still missing.
%
%   in essence it says that the projection of the absolute modulus of a
%   complex field is the same as when one takes a sliding window in the
%   Fourier domain, makes an inverse FFT of each slice, take the absolute
%   modulus of that and sum it up while moving the window through the
%   spectrum.  This has important applications for scanned light-sheets and
%   how to generate them.
%
%   Reto, May 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
%% assume this is the e-field in real space (real or complex valued)
efield=double(imread('cameraman.tif'))+j*rand(256,256);
%% then this is the spectrum of the efield
spectrum=fftshift(fft2(efield));
%% the intensity is the modulus squared of the efield  -> a real valued image         
I=abs(efield).^2;                         
%% slice is an array that is created by superposition of inverse FFTs of spectral slices 
slice=zeros(256,256);
%% smear is created by scanning the image in real space                  
smear=zeros(256,256);                      
 

for k=1:256
    temp=zeros(256,256);
    temp(:,k)=spectrum(:,k);
    %% take one slice and take inverse FFT
    temp=ifft2(temp);
    %% superimpose intensities (modulus squared) of inverse FFTs of spectral slices          
    slice=slice+abs(temp).^2;           
    %% smearing the image in x-direction, intensity is superimposed at every position
    smear=smear+(circshift(I,[0,k]));      
end
 
figure;imshow(I,[])
figure;imshow(slice,[])
figure;imshow(smear,[])
%% profiles of the smeared intensity images
 
 
slice=slice(:,128);
slice=slice/max(slice);
 
smear=smear(:,128);
smear=smear/max(smear);
 
%alternatively, one can also just project the intensity image I in x
projection=sum(I');
projection=projection/max(projection);
 
 
%% all profiles are identical, ignoring scaling factors
figure;plot(slice,'b+');hold on;plot((smear),'r*');hold on;plot(projection,'ko')


%% The following shows modulation added to each line scan
efield_xft = fft(efield,[],2);
spectrum_unshifted = ifftshift(spectrum);
N = 256;
 
for k=1:N
 
temp=zeros(256,256);
temp(:,k)=spectrum_unshifted(:,k);
temp = ifft2(temp);
 
Q = zeros(N);
Q(:,1) = spectrum_unshifted(:,k);
Q = ifft2(Q);
 
hfig = figure;
plot(real(Q(:,1))*N,'r-');
hold on;
plot(real(efield_xft(:,k)),'ro');
 
plot(imag(Q(:,1))*N,'b-');
plot(imag(efield_xft(:,k)),'bo');
 
legend({'real(Q)','Real 1D FFT of efield','imag(Q)','Imag 1D FFT of efield'});
 
pause(1);
 
hfig2 = figure;
plot(real(exp(i*2*pi*(k-1)/N.*(0:N-1))),'r-');
hold on;
plot(real(temp(1,:)./Q(1,:)),'ro');
plot(imag(exp(i*2*pi*(k-1)/N.*(0:N-1))),'b-');
plot(imag(temp(1,:)./Q(1,:)),'bo');
 
pause(1);
 
close(hfig);
close(hfig2);
end;




?
