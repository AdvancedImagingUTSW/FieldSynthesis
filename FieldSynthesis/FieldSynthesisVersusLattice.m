%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%Simulation for field synthesis
%
%   compares field synthesis vs lattice

%
%   Reto, May 2017
%   Mark Kittisopikul, August 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;clear all

%% Parameters
% Size of the mask
n=4096;
% Width of the annulus
w=5;
% Radius of the annulus
r=258;
% Offset for side slits
offset = 256;
% offset=198;
% Display range
dispRange = (-600:600)+n/2+1;

%% Create clean annulus
% Unneeded to initialize since we will create the matrix with hypot
% A=zeros(n,n);


% Vector for x and y, which should be symmetric
v = 1:n;
% zeroth order coefficient is at n/2+1,n/2+1 due to fftshift/ifftshift
v = v-floor(n/2)-1;

annulus = createAnnulus(v, r, w);

% Select columns for mask
abs_v = abs(v);
selected_columns = (abs_v < offset & abs_v > offset-w) | (v < w/2 & v > -w/2);

% Remove unselected columns from mask
A = annulus;
A(:,~selected_columns) = false;
A = double(A);

%% Obsolete for loop version
% tic;
% for k=1:n
%     for l=1:n
%         
%         x=k-n/2-1;
%         y=l-n/2-1;
%         q=sqrt(x^2+y^2);
%         
%         if q<r+w && q>r-w && y<2 && y>-3
%             A(k,l)=1;
%           
%         end   
%         
%         if q<r+w && q>r-w && y<offset && y>offset-5
%             A(k,l)=1;
%           
%         end  
%         
%                 if q<r+w && q>r-w && y>-offset && y<-(offset-5)
%             A(k,l)=1;
%           
%         end  
%         
%     end
% end
% toc;
%% Field Synthesis

E = sum(abs(ifftshift(ifft(ifftshift(A)))).^2,2);
% E=E/max(E);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lattice simulation

% B=fftshift(fft2(A));
B=fftshift(ifft2(ifftshift(A)));
%figure;imshow(abs(B).^2,[])
latticeLineProfile=abs(B).^2;

latticeLineProfile=sum(latticeLineProfile,2);
latticeLineProfile=latticeLineProfile*n;
% latticeLineProfile=latticeLineProfile/max(latticeLineProfile);

%% compare lattice profile to field synthesis profile
figure;
subplot(3,1,1);
plot(dispRange-n/2+1,latticeLineProfile(dispRange));
xlim([min(dispRange) max(dispRange)]-n/2+1);
title('Conventional Lattice Profile');

subplot(3,1,2);
plot(dispRange-n/2+1,E(dispRange),'r')
xlim([min(dispRange) max(dispRange)]-n/2+1);
title('Field Synthesis Lattice Profile');

subplot(3,1,3);
plot(dispRange-n/2+1,latticeLineProfile(dispRange));
hold on;
plot(dispRange-n/2+1,E(dispRange),'r--')
xlim([min(dispRange) max(dispRange)]-n/2+1);
title('Comparison');

%% Analysis of all interference patterns in lattice

lattice=abs(B).^2;
test=fftshift(fft2(lattice)); %Fourier transform of lattice


%% dithering lattice: lattice pattern is shifted by subpixel steps and added
latPro=zeros(n,n);
steps=64;stepsize=0.25;
for k=1:steps
temp=lattice;



temp=fftshift(fft2(temp));
u=0:1:n-1; 
i=sqrt(-1);
wa=exp(-i*2*pi.*(u*stepsize*k)./(n));
wa=repmat(wa,[n,1]);

temp=temp.*wa;

temp=abs(ifft2(temp));

latPro=latPro+temp;
end

latPro=latPro/steps;
figure;imshow(latPro,[])
test2=fftshift(fft2(latPro)); %Fourier transform of dithered lattice

%%



h = figure;
set(h,'Position',[625 532 1165 626]);
subplot(2,3,1)
imshow(A(dispRange,dispRange),[]);colormap hot
title('Electric field in pupil');

% figure;
subplot(2,3,2)
imshow((abs(test(dispRange,dispRange))),[0,1E-6]);colormap hot
title('Fourier components of lattice intensity');

% figure;
subplot(2,3,3)
imshow((abs(test2(dispRange,dispRange))),[0,1E-6]);colormap hot
title('Fourier components of dithered lattice intensity');

% figure;
subplot(2,3,4);
imshow(B(dispRange,dispRange),[]);
% imshowpair(real(B(dispRange,dispRange)),imag(B(dispRange,dispRange)));
title('Electric field of lattice at focal plane');
% figure;
subplot(2,3,5)
imshow(lattice(dispRange,dispRange),[]);
title('Intensity of lattice');
% figure;
subplot(2,3,6)
imshow(latPro(dispRange,dispRange),[]);
title('Averaged Intensity of dithered lattice');

%% Major Line Scan Components
% figure;
% % (abs_v < offset & abs_v > offset-5) | (v < 3 & v > -3);
% idx_minus = -(offset-4:offset-1);
% idx_center = -2:2;
% idx_plus = +(offset-4:offset-1);
% accumPattern = zeros(size(A));
% 
% subplot(3,4,1);
% idx = idx_minus;
% selection = ismember(v,idx);
% lineSelection = repmat(selection,n,1);
% subPattern = fftshift(ifft2(ifftshift(A.*lineSelection)));
% accumPattern = accumPattern + abs(subPattern).^2;
% imshowpair(lineSelection(dispRange,dispRange), ...
%            A(dispRange,dispRange));
%        
% subplot(3,4,5);
% stem(dispRange,A(dispRange,selection));
% xlim([min(dispRange) max(dispRange)]);
% 
% subplot(3,4,9);
% imshow(abs(subPattern(dispRange,dispRange)).^2,[]);
%        
% subplot(3,4,2);
% idx = idx_center;
% selection = ismember(v,idx);
% lineSelection = repmat(selection,n,1);
% subPattern = fftshift(ifft2(ifftshift(A.*lineSelection)));
% accumPattern = accumPattern + abs(subPattern).^2;
% imshowpair(lineSelection(dispRange,dispRange), ...
%            A(dispRange,dispRange));
%        
% subplot(3,4,6);
% stem(dispRange,A(dispRange,selection));
% xlim([min(dispRange) max(dispRange)]);
% 
% subplot(3,4,10);
% imshow(abs(subPattern(dispRange,dispRange)).^2,[]);
%        
% subplot(3,4,3);
% idx = idx_plus;
% selection = ismember(v,idx);
% lineSelection = repmat(selection,n,1);
% subPattern = fftshift(ifft2(ifftshift(A.*lineSelection)));
% accumPattern = accumPattern + abs(subPattern).^2;
% imshowpair(lineSelection(dispRange,dispRange), ...
%            A(dispRange,dispRange));
%        
% subplot(3,4,7);
% stem(dispRange,A(dispRange,selection));
% xlim([min(dispRange) max(dispRange)]);
% 
% subplot(3,4,11);
% imshow(abs(subPattern(dispRange,dispRange)).^2,[]);
% 
% subplot(3,4,4);
% selection = ismember(v,[idx_minus,idx_plus]);
% lineSelection = repmat(selection,n,1);
% subPattern = fftshift(ifft2(ifftshift(A.*lineSelection)));
% imshowpair(lineSelection(dispRange,dispRange), ...
%            A(dispRange,dispRange));
%        
% subplot(3,4,8);
% stem(dispRange,A(dispRange,selection));
% xlim([min(dispRange) max(dispRange)]);
% 
% subplot(3,4,12);
% imshow(accumPattern,[]);