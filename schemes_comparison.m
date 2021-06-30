
%% BEGIN
clc;clear;

tic
E_0=[1 0;0 1];
A_0=zeros(2,2);


counter=1;
steps = 500;
for loop=[5]%(4:1:10)%[10]%(4:1:10)%[8]%(4:1:10)%%[8]%(1:1:12)%[8]%(1:1:12)%[1]%(5:1:12)%[1,5,10,11]%,20,21,22
N_1=2*(loop^2)+1%2*(loop^2)+1% number of points in x_1-1

N_2=N_1; % number of points in x_2

%% Mesh parameters
h_1=2*pi/(N_1); % step in x_2
h_2=2*pi/(N_2); % step in x_2

%% Coordinates
x=zeros(N_2,N_1,2);
[x(:,:,1),x(:,:,2)]=meshgrid(-pi+h_1/2:h_1:pi-h_1/2,-pi+h_2/2:h_2:pi-h_2/2);

%% Material coeficient matrix

 Pixels = imread('structure_4.png');

 pixa=round(linspace(1,size(Pixels,2),N_1));
 piya=round(linspace(1,size(Pixels,1),N_2));
 par=100
  
  Min_eig = 100000;
 Max_eig= 0;
 A=zeros(N_2,N_1,2,2);
  C_ref=zeros(N_2,N_1,2,2); 
  for i=1:N_2
     for j=1:N_1    
          A(i,j,:,:)=a_matrix_img_aniso(Pixels(piya(i),pixa(j)),par);
         % A(i,j,:,:)=a_matrix(x(i,j,:));
          %A(i,j,:,:)=a_matrix_img_aniso(Pixels(piya(i),pixa(j)),par);          
     end       
  end

 

%% Material ananlysis
d=[mean(mean(A(:,:,1,1))) mean(mean(A(:,:,1,2)));...
   mean(mean(A(:,:,2,1))) mean(mean(A(:,:,2,2)))];

d_inv=d^-1;
for i=1:N_2
     for j=1:N_1    
          C_ref(i,j,:,:)=d;
          C_ref_inv(i,j,:,:)=d_inv;
     end
end

%% Derivatives
G = G_clasic(N_1,N_2);
G_n=G_mean(N_1,N_2,[1 0;0 1]);
G_m=G_mean(N_1,N_2,d);

%% Preconditioning
[M_m_half_inv] = M_mean_half_inv(N_1,N_2,d); %M^-1/2

M_m= -(d(1,1).*(G(:,:,1).^2)+d(2,2).*(G(:,:,2).^2)...
                   +2*d(1,2).*(G(:,:,1).*(G(:,:,2))));  
               
M_m((end+1)/2,(end+1)/2)=1;

%% Conditions:
tau=0.25
toler = 1e-6;
c_000=zeros(N_2,N_1);


%% GBection based solvers
c_0 = c_000;

disp('Gradient-Based Conjugage Gradient solver')
for k=1:1
    E=E_0(:,k); 
    tic;
   [C_GB_CG,st_GB_CG, norm_evol_GB_CG_grad, norm_evol_GB_CG_energy, estim_GB_CG,  sol_norm_GB_CG]=solver_GP_projection_left(A,G,c_0,E,steps,toler,M_m,d,tau,G_m,C_ref_inv);
     
    T_GB_CG(k,counter) = toc;
    stepses=size(norm_evol_GB_CG_grad);

    S_GB_CG(k,counter) = stepses(2);
    A_(:,k)=Hom_parameter_grad(C_GB_CG,A,G,E) % Compute homogenized parameter
    A_GB_CG(counter)=A_(1,1);
end

disp('Modifield Gradient-Based Conjugage Gradient solver')
for k=1:1
    E=E_0(:,k); 
    tic;
    [C_GB_CG_mod,st_GB_CG_mod,norm_evol_GB_CG_mod_grad,norm_evol_GB_CG_mod_energy, estim_GB_CG_mod,  sol_norm_GB_CG_mod]=solver_GP_projection_left_modif(A,G,c_0,E,steps,toler,M_m,d,tau,G_m,C_ref_inv);

    T_GB_CG_mod(k,counter) = toc;
    stepses=size(norm_evol_GB_CG_mod_grad);

    S_GB_CG_mod(k,counter) = stepses(2);

    A_(:,k)=Hom_parameter_grad(C_GB_CG_mod,A,G,E) % Compute homogenized parameter
    A_GB_CG_mod(counter)=A_(1,1);
end

disp('Gradient-Based Preconditioned Conjugage Gradient solver')
for k=1:1
    E=E_0(:,k); 
    tic;
    [C_GB_PCG,st_GB_PCG, norm_evol_GB_PCG_grad, norm_evol_GB_PCG_energy, estim_GB_PCG,  sol_norm_GB_PCG]=solver_PCG_projection_left(A,G,c_0,E,steps,toler,M_m,d,tau);

    T_GB_PCG(k,counter) = toc;
         stepses=size(norm_evol_GB_PCG_grad);

    S_GB_PCG(k,counter) = stepses(2);

    A_(:,k)=Hom_parameter_grad(C_GB_PCG,A,G,E) % Compute homogenized parameter
    A_GB_PCG(counter)=A_(1,1);
end

%% Displacement-Based Preconditioned Conjugage Gradient solver
disp('Displacement-Based Preconditioned Conjugage Gradient solver')

c_0=c_000;
%toler = 1e-10;
for k=1:1
    E=E_0(:,k);
    tic;
    [C_DB_PCG,st_DB_PCG,norm_evol_DB_PCG_rr,norm_evol_DB_PCG_energy, norm_evol_DB_PCG_grad, estim_DB_PCG]=solver_PCG_left(A,G,c_0,E,steps,toler,M_m,tau);% with preconditioning

    T_DB_PCG(k,counter)=toc;
    S_DB_PCG(k,counter) = st_DB_PCG+1;
    A_(:,k)=Hom_parameter(C_DB_PCG,A,G,E)% Compute homogenized parameter
    A_DB_PCG(counter)=A_(1,1);
    
end

% error 
  NoP(counter)=N_1*N_2;
  counter=counter+1;
end

% % Plot estimates
 rel_estim_DB_PCG=(estim_DB_PCG./estim_DB_PCG(1));%.^(1/2);
% rel_estimg=estimg./estimg(1);
 rel_estim_GB_CG=estim_GB_CG./estim_GB_CG(1);


%% Plot residuals
 figure 
 hold on

 plot((1:S_DB_PCG(1,end)),norm_evol_DB_PCG_rr,'--')
 plot((1:S_DB_PCG(1,end)),norm_evol_DB_PCG_grad,'--')
 plot((1:S_DB_PCG(1,end)),norm_evol_DB_PCG_energy,'--')
 

  plot((1:S_GB_CG(1,end)),norm_evol_GB_CG_grad,'-.x')
  plot((1:S_GB_CG(1,end)),norm_evol_GB_CG_energy,'-.o')

    plot((1:S_GB_CG_mod(1,end)),norm_evol_GB_CG_mod_grad,'-.x')
  plot((1:S_GB_CG_mod(1,end)),norm_evol_GB_CG_mod_energy,'-.o')
  
  plot((1:S_GB_PCG(1,end)),norm_evol_GB_PCG_grad,'-.x')
  plot((1:S_GB_PCG(1,end)),norm_evol_GB_PCG_energy,'-.o')
 
set(gca, 'XScale', 'linear', 'YScale', 'log');
legend('DB || r ||','DB  ||Dr||', 'DB  ||r||_M', ...
         'GB CG ||Dr||', 'GB CG ||r||_M',...
          'GB CG mod ||Dr||', 'GB CG mod ||r||_M',...
        'GB  PCG ||Dr||', 'GB PCG ||r||_M')
title('Residuals ')
%% Plot residuals
 figure 
 hold on
 plot((1:numel(estim_DB_PCG)),abs(rel_estim_DB_PCG),'.')
  plot((1:numel(estim_DB_PCG)),abs(rel_estim_DB_PCG./(1-tau)),'.')
 plot((1:S_DB_PCG(1,end)),norm_evol_DB_PCG_rr,'--')
 plot((1:S_DB_PCG(1,end)),norm_evol_DB_PCG_grad,'--')
 plot((1:S_DB_PCG(1,end)),norm_evol_DB_PCG_energy,'--')
 
 
 plot((1:numel(estim_GB_CG)),abs(rel_estim_GB_CG),'.')
 plot((1:numel(estim_GB_CG)),abs(rel_estim_GB_CG./(1-tau)),'.')
  plot((1:S_GB_CG(1,end)),norm_evol_GB_CG_grad,'-.x')
  plot((1:S_GB_CG(1,end)),norm_evol_GB_CG_energy,'-.o')
 
set(gca, 'XScale', 'linear', 'YScale', 'log');
legend('DB estim || e_k ||_K',' DB estim || e_k ||_K','DB || r ||','DB  ||Dr||', 'DB  ||r||_M', ...
    'GB estim || e_k ||_K',' GB estim || e_k ||_K','GB  ||Dr||', 'GB ||r||_M')
title('Residuals 2')

%% Plot residuals
 figure 
 hold on
 plot((1:S_DB_PCG(1,end)),norm_evol_DB_PCG_rr,'--')
 plot((1:S_DB_PCG(1,end)),norm_evol_DB_PCG_grad,'--')
 plot((1:S_DB_PCG(1,end)),norm_evol_DB_PCG_energy,'--')

 
 %plot((1:Sg(1,end)),norm_evolg,'-o') %'DB modif',
 
  plot((1:S_GB_CG(1,end)),norm_evol_GB_CG_grad,'-.')
  plot((1:S_GB_CG(1,end)),norm_evol_GB_CG_energy,'-.')
  
%  plot((1:Spc(1,end)),norm_evolpc,'-.')%'GB modif',
  plot((1:S_GB_CG_mod(1,end)),norm_evol_GB_CG_mod_grad,'-.')

set(gca, 'XScale', 'linear', 'YScale', 'log');
legend('DB PCG || r ||','DB PCG  ||Dr||', 'DB PCG ||r||_M' ...
    ,'GB orig ||Dr|| ','GB orig ||r||_M','GB modif 2')
title('Residuals with modif CG')
%% Plot solution difference
%  figure 
%  hold on
 %plot((1:Sp(1,end)),abs(sol_normp-sol_normg),'-.*')
% plot((1:Ss(1,end)),abs(sol_normp(1,end-1)-sol_norms),'-.o')
 %
 %plot((1:Ss(1,end)),abs(sol_norms-sol_normg(1,end-1)),'-.^')
 % plot((1:Sg(1,end)),sol_normg,'-.^')
 %plot((1:Sp(1,end)),abs(norm_evolp-norm_evolg),'-.*')

% set(gca, 'XScale', 'linear', 'YScale', 'log');
% legend('sol_norm  GB- grad ','sol_norm GB - sym')
% title('solution difference')
%% Plot solution norm

 figure 
 hold on
 
%  plot((1:Sg(1,end)),abs(sol_normg(end)-sol_normg),'x')
  
 plot((1:S_GB_CG(1,end)),abs(sol_norm_GB_CG(end)-sol_norm_GB_CG),'-.*')
 
  plot((1:S_GB_CG_mod(1,end)),abs(sol_norm_GB_CG_mod(end)-sol_norm_GB_CG_mod),'-.o')
  
plot((1:S_GB_PCG(1,end)),abs(sol_norm_GB_PCG(end)-sol_norm_GB_PCG),'-.^')
 
set(gca, 'XScale', 'linear', 'YScale', 'log');
legend('GB CG orig','GB CG modif','GB PCG ')
title('solution norm')

%% Plot steps

% figure 
%  hold on
%  plot(NoP,S3(1,:),'^')
%  plot(NoP,Sp,'-.*')
%     plot(NoP,Sg,'-.^') 
%     
% legend('DB','DB modif','GB orig','GB modif')
% title(' Number of steps')
%% Plot hom_mat prop
% A_refer=real(A3);%9.403399951113226;
% figure 
%  hold on
% 
%  plot(NoP,abs(A_refer-real(A3))','^')
%   plot(NoP,abs(A_refer-real(Ag))','-.^') 
%   
%  plot(NoP,abs(A_refer-real(Ap))','-.*')
%  plot(NoP,abs(A_refer-real(Apc))','-.>') 
% 
%  set(gca, 'XScale', 'linear', 'YScale', 'log');
% 
% legend('DB','DB modif','GB orig','GB modif')
%  title(' Hom parameter error')
% figure 
%  hold on
%  plot(NoP,imag(A1)','x')
%  plot(NoP,imag(A2)','o')
%  plot(NoP,imag(A3)','^')
%  plot(NoP,imag(Ap)','-.*')
%  plot(NoP,imag(Ag)','-.^') 
% legend('in G','Symetric','Left','GB','grad norm')


close all

 %% Plot A function