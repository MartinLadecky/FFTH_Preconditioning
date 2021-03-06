%function [A_0,st,t,Nbf]=Hom_solver(N)
%% HOMOGENIZATION PROBLEM
%  sol PDE -div(A(x)grad(u))=div(A(x)E))
%  then homogenized A_0E=int(A(x)(E+grad(u)))dx/(volume of omega)
%% Input
% N   [1]    -number of points in sample
%% Output
% A_0 [2,2] -Homogenized material parameters
% st  [1]   -Number of iteration steps
% t   [1]   -Time
%% BEGIN
clc;clear;

tic
E_0=[1 0;0 1];
A_0=zeros(2,2);

N_1=51%round(sqrt(N)./2).*2+1; % number of points in x_1-1
N_2=51; % number of points in x_2
% Nbf=N_1*N_2;
%% Mesh parameters
h_1=2*pi/(N_1); % step in x_2
h_2=2*pi/(N_2); % step in x_2

%% Coordinates
x=zeros(N_2,N_1,2);
[x(:,:,1),x(:,:,2)]=meshgrid(-pi+h_1/2:h_1:pi-h_1/2,-pi+h_2/2:h_2:pi-h_2/2);
%% Derivatives
G_n=G_matrix(N_1,N_2);

%% Material coeficient matrix
 Min_eig = 100000;
 Max_eig= 0;
 A=zeros(N_2,N_1,2,2);
 for i=1:N_2
     for j=1:N_1    
         A(i,j,:,:)=a_matrix(x(i,j,:)); % functional input data
         %A(i,j,:,:)=a_matrix_img(Pixels(piya(i),pixa(j)));% image input data
         pom = zeros(2);
         pom(1,1) = A(i,j,1,1); pom(1,2) = A(i,j,1,2);
         pom(2,1) = A(i,j,2,1); pom(2,2) = A(i,j,2,2);        
         p = eig(pom);
         min_eig = min(p);
         max_eig = max(p);
         if (min_eig<Min_eig) 
             Min_eig=min_eig; 
         end
         if (max_eig>Max_eig)
             Max_eig=max_eig;
         end       
     end
 end
 Max_kappa = Max_eig/Min_eig;
 [Min_eig,Max_eig,Max_kappa]
 

%% Preconditioning with constant frequencie
d=[mean(mean(A(:,:,1,1))) mean(mean(A(:,:,1,2))) ;...
   mean(mean(A(:,:,2,1))) mean(mean(A(:,:,2,2)))];

 M_f_const= d(1,1).*(G_n(:,:,1).^2)+d(2,2).*(G_n(:,:,2).^2)...% M_f by lemma3.2 (frequencies)
     +2*d(1,2).*(G_n(:,:,1).*(G_n(:,:,2)));  % predpodmienovacia matica
 M_f_const((end+1)/2,(end+1)/2)=1; % doplnenie stredneho clenu kôli inverzi "tip"


%% Preconditioning with significant frequencie

% 
A_freq11=fft2(ifftshift(A(:,:,1,1)));
A_freq22=fft2(ifftshift(A(:,:,2,2)));
A_freq12=fft2(ifftshift(A(:,:,1,2)));
    % 0 frequency
d11 = A_freq11(1,2); % nad diagonalu v matici FAF !
d22 = A_freq22(1,2);
d12 = A_freq12(1,2);

dd11 = A_freq11(1,1);
dd22 = A_freq22(1,1);
dd12 = A_freq12(1,1);

G_1 = G_n(:,:,1);
G_2 = G_n(:,:,2);
G_1p = [G_1(:,2:end),G_1(:,1)];% posunuta matica derivaci
G_2p = [G_2(:,2:end),G_2(:,1)];

GG0 = (dd11*G_1.*G_1 + dd22*G_2.*G_2 + dd12*(G_2.*G_1 + G_1.*G_2))/N_1/N_2;
GG0((end+1)/2,(end+1)/2)=-1; %Predpodmienovacia 0 f
GG0 = -GG0;

GG1 = -conj((d11*conj(G_1).*G_1p + d22*conj(G_2).*G_2p + ...
      d12*conj(G_1).*G_2p + d12*conj(G_2).*G_1p))/N_1/N_2;
GG1 = -GG1;%Predpodmienovacia 1 f

P = zeros(N_1,N_1,N_2); % sada predpodminovacich matic
for r = 1:N_2    
    m1 = GG0(r,:); % radek Predpod matice 0f
    m2 = GG1(r,:); % radek Predpod matice 1f  
    p = spdiags([m2(1:N_1-1)'],[-1],N_1,N_1); % zoradi m2 riadok do poddiagonaly, riedka matica
%     p(1,N_1) = m2(N_1);  % ZDE OPRAVA NEBO NE
    p = p + p' + diag(m1);     
    P(:,:,r) = p;% zaradi to do glob matice
end  

U1 = zeros(N_2,N_1);
U2 = U1;
for r = 1:N_2    
    m1 = GG0(r,:);
    m2 = GG1(r,:);    
    p = spdiags([m2(1:N_1-1)'],[-1],N_1,N_1);
%     p(1,N_1) = m2(N_1);  % ZDE OPRAVA NEBO NE
    p = p + p' + diag(m1); 
    pom = chol(p);  % p = pom'*pom % chol rozklad
    U1(r,:) = diag(pom); % diagonala chol
    U2(r,1:N_1-1) = diag(pom,1); % mimodiagonala chol
end

%% ZACATEK POCITANI:
steps = 100;
toler = 1e-6;

%% SOLVER without preconditionig

c_000=rand(N_2,N_1);
% c_000 = c_000-mean(c_000);
% c_000 = zeros(N_2,N_1);
c_0 = c_000;

for k=1:2 
E=E_0(:,k);
tic; 
[C,st]=CG_solver(A,G_n,c_0,E,steps,toler); % without preconditioning
T1(k) = toc;
S1(k) = st;
A_0(:,k)=Hom_parameter(C,A,G_n,E); % Compute homogenized parameter
end
A_01 = A_0
T1
S1

%% SOLVER with constant preconditionig

c_0=c_000;

for k=1:2
E=E_0(:,k);
tic;
[C,st]=CGP_solver_constant(A,G_n,c_0,E,steps,toler,GG0);% with preconditioning
T2(k) = toc;
S2(k) = st;
A_0(:,k)=Hom_parameter(C,A,G_n,E);% Compute homogenized parameter
end
A_02 = A_0
T2
S2

%% SOLVER with efficient enhanced preconditionig

c_0=c_000;

for k=1:2
E=E_0(:,k);
tic;
[C,st]=CGP_solver_1f_v3(A,G_n,c_0,E,steps,toler,GG0,GG1,U1,U2); % with better preconditioning
T4(k)=toc;
S4(k) = st;
A_0(:,k)=Hom_parameter(C,A,G_n,E);% Compute homogenized parameter
end
A_04 = A_0
T4
S4



