clc;clear;

N_1=5% number of points in x_1-1

%% Mesh parameters
h_1=2*pi/(N_1) ;
h_bar=2*pi/(N_1*2-1);
x(:)=-pi+h_1/2:h_1:pi-h_1/2;
x_bar(:)=-pi+h_bar/2:h_bar:pi-h_bar/2;
%m(:)=x(:)
%%  frequencies 
k=-(N_1-1)/2:(N_1-1)/2;
m=-(N_1-1)/2:(N_1-1)/2;
%  frequencies 
N_bf=N_1;
%% Frequences
F=zeros(N_bf,N_bf);
for i=1:N_bf
    for j=1:N_bf 
       F(i,j)=exp(-2*pi*1i*(k(i)*m(j)')/N_1);
    end
end
F=F/sqrt(N_bf);

%% materials

cpar=2
A=zeros(N_bf,N_bf);
 for i=1:N_1   
         A(i,i)=x(i)^2;
         a(i)=x(i)^2;
         C(i,i)=cpar;
 end


%% Frequences
In=zeros(N_bf,N_bf);
for i=1:N_bf
    for j=1:N_bf 
       In(i,j)=exp(-1i*(k(i)*x(j)'));
    end
end
%% Gradients
for i=1:N_bf
       G(i,i)=1i*(k(i));
end
%GG=G'*G
G((N_1+1)/2,(N_1+1)/2,:)=1;

PrecIvana=G*((G'*C*G)^-1)*G';
PrecIvana((N_1+1)/2,(N_1+1)/2,:)=0
%GG=GG^-1
%G*GG*G'

%% Projection Operator
for i=1:N_bf
       Gn(i,i)=1i*(k(i))./sqrt(cpar*k(i).^2);
end
Gn((N_1+1)/2,(N_1+1)/2,:)=0;

Gamma=Gn*(Gn')

Project((N_1+1)/2,(N_1+1)/2,:)=0;

%% 
 %AA=F*a'
 Fa=fftshift(fft(ifftshift(a')));
% a_hat=fft(ifftshift(a'))
 
 %LargeFa = reshape([Fa'; zeros(size(Fa'))],[],1);
 %Fa(2*numel(Fa)-1)=0
LargeFa=2.*[zeros(1,(numel(Fa)-1)/2),Fa',zeros(1,(numel(Fa)-1)/2)]

 
 La=fftshift(ifft(ifftshift(LargeFa')))
 a 
 %In*A*In'
 figure 
 plot(x,a)
 hold on
 plot(x_bar,La)
 