



clc;
clear;


N_1=5;% number of points in x_1-1
N_2=N_1;

%% Mesh parameters
h_1=2*pi/(N_1); % step in x_2
h_2=2*pi/(N_2); % step in x_2

%% Coordinates
x=zeros(N_2,N_1,2);
[x(:,:,1),x(:,:,2)]=meshgrid(-pi+h_1/2:h_1:pi-h_1/2,-pi+h_2/2:h_2:pi-h_2/2);

Pixels = imread('structure_3.png');
%figure;
%Pixels = im2double(Pixels);
% freqIm=fftshift(ifft2(ifftshift(Pixels)));
% Pixels=real(fftshift(fft2(ifftshift(freqIm))));
% 
% %imshow(flipud(Pixels))
% Pixels=(Pixels+flipud(Pixels));
% Pixels(abs(Pixels)<270)=0;
% Pixels(abs(Pixels)>280)=255;

%figure
%imshow(Pixels)

%freqIm=fftshift(ifft2(ifftshift(Pixels)));
%freqIm(abs(freqIm)<1e-1)=0;

%v = nonzeros(freqIm);
%figure
%imshow(freqIm)


%freqIm([1:450 550:end],:) = 0;
%freqIm(:,[1:495 505:end]) = 0;
freqIm=fftshift(ifft2(ifftshift(Pixels)));

freqIm=freqIm./5;
freqIm(500,500)=freqIm(500,500);
freqIm(500,501)=freqIm(500,501).*5;
freqIm(500,501)=freqIm(500,501).*10;


figure
newimage=real(fftshift(fft2(ifftshift(freqIm))));
%newimage(abs(newimage)<mean(newimage))=0;
newimage(newimage<0)=0;
imshow(newimage,'DisplayRange',[0,255])
save('structure11.mat','newimage');

 
 pixa=round(linspace(1,size(Pixels,2),N_1));
 piya=round(linspace(1,size(Pixels,1),N_2));
 
 
 
 