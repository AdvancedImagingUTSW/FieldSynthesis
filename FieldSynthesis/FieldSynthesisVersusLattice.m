%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%Simulation for field synthesis
%
%   compares field synthesis vs lattice
 
%
%   Reto, May 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
clear all;clear all
n=4096;
w=5;
r=200
offset=158
A=zeros(4096,4096);
% create clean annulus
tic;
for k=1:n
    for l=1:n
        
        x=k-n/2-1;
        y=l-n/2-1;
        q=sqrt(x^2+y^2);
        
        if q<r+w && q>r-w && y<2 && y>-3
            A(k,l)=1;
          
        end   
        
        if q<r+w && q>r-w && y<offset && y>offset-5
            A(k,l)=1;
          
        end  
        
                if q<r+w && q>r-w && y>-offset && y<-(offset-5)
            A(k,l)=1;
          
        end  
        
    end
end
toc; 
%% Field Synthesis
 
E = sum(abs(ifftshift(ifft(ifftshift (A)))).^2,2);
E=E/max(E);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lattice simulation
 
% B=fftshift(fft2(A));
B=fftshift(ifft2(ifftshift (A)));
%figure;imshow(abs(B).^2,[])
lattice=abs(B).^2;
 
lattice=sum(lattice');
lattice=lattice/max(lattice);
 
%% compare lattice profile to field synthesis profile
figure;plot(lattice(1500:2500));hold on;plot(E(1500:2500),'r')
 
%% Analysis of all interference patterns in lattice
 
lattice=abs(B).^2;
latPro=sum(lattice'); %projected lattice pattern
latPro=repmat(latPro',[1,4096]); %smear out lattice profile to a light-sheet: dithering
 
test=fftshift(fft2(lattice)); %Fourier transform of lattice
test2=fftshift(fft2(latPro)); %Fourier transform of dithered lattice
 
 
figure;imshow(A(1500:2500,1500:2500),[]);title('electric field in pupil')
 
figure;imshow(log(abs(test(1500:2500,1500:2500))),[]);title('Fourier components of lattice intensity')
 
figure;imshow(log(abs(test2(1500:2500,1500:2500))),[]);title('Fourier components of dithered lattice intensity')
 
figure;imshow(lattice(1500:2500,1500:2500),[])
figure;imshow(latPro(1500:2500,1500:2500),[])


