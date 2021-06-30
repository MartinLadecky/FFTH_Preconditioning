
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


%% Projection based solver
c_0 = c_000;

disp('Projection based solver')
for k=1:1
    E=E_0(:,k); 
    tic;
   %[Cp,st,norm_evol_proj_grad, norm_evol_proj_energy, estimp,  sol_normp]=solver_GP_projection_left(A,G,c_0,E,steps,toler,M_m,d,tau,G_m,C_ref_inv);
   disp('Projection based solver')
    [Cp,st,norm_evol_proj_grad, norm_evol_proj_energy, estimp,  sol_normp]=solver_PCG_projection_left(A,G,c_0,E,steps,toler,M_m,d,tau);
    
    %sol_normp
    norm_evol_proj_grad
    Tp(k,counter) = toc;
         stepses=size(norm_evol_proj_grad);

    Sp(k,counter) = stepses(2);
  %  Sp(k,counter) = st+1;
    A_p(:,k)=Hom_parameter_grad(Cp,A,G,E) % Compute homogenized parameter
    Ap(counter)=A_p(1,1);
end

disp('Projection based solver: modified with C_ref new')
for k=1:1
    E=E_0(:,k); 
    tic;
    [Cpm,st,norm_evolpm, estimpm,  sol_normpm]=solver_GP_projection_left_modif(A,G,c_0,E,steps,toler,M_m,d,tau,G_m,C_ref_inv);
    %sol_normp
    norm_evolpm;
    Tpm(k,counter) = toc;
    stepses=size(norm_evolpm);

    Spm(k,counter) = stepses(2);
  %  Sp(k,counter) = st+1;
    A_p(:,k)=Hom_parameter_grad(Cpm,A,G,E) % Compute homogenized parameter
    Apm(counter)=A_p(1,1);
end

% % Projection based solver: modified with C_ref 
% c_0 = c_000;
% 
% disp('Projection based solver: modified with C_ref ')
% for k=1:1
%     E=E_0(:,k); 
%     tic;
%     [Cpc,st,norm_evolpc, estimpc,  sol_normpc]=solver_GP_projection_left_Cref(A,G,c_0,E,steps,toler,M_m,d,tau,G_m,C_ref_inv);
%     norm_evolpc;
%     Tpc(k,counter) = toc;
%      stepses=size(norm_evolpc);
% 
%     Spc(k,counter) = stepses(2);
%     A_pc(:,k)=Hom_parameter_grad(Cpc,A,G,E) % Compute homogenized parameter
%     Apc(counter)=A_pc(1,1);
% end


% % SOLVER with constant preconditionig from left grad error measure
% c_0=c_000;
% disp('SOLVER with constant preconditionig from left hand side ::: grad error measure')
% for k=1:1
%     E=E_0(:,k);
%     tic;
%     [Cg,st,norm_evolg, estimg,  sol_normg]=solver_PCG_left_grad_norm(A,G,c_0,E,steps,toler,M_m,tau);% with preconditioning
%     sol_normg
%     Tg(k,counter)=toc;
%     stepses=size(norm_evolg);
%     Sg(k,counter) =stepses(2);
%     A_g(:,k)=Hom_parameter(Cg,A,G,E)% Compute homogenized parameter
%     Ag(counter)=A_g(1,1);
%     
% end

%% SOLVER with constant preconditionig from left hand side
c_0=c_000;
disp('SOLVER with constant preconditionig from left hand side')
%toler = 1e-10;
for k=1:1
    E=E_0(:,k);
    tic;
    [C,st,norm_evol_PCG_rr,norm_evol_PCG_rz,norm_evol_PCG_DrDr, estim3 ]=solver_PCG_left(A,G,c_0,E,steps,toler,M_m,tau);% with preconditioning
    norm_evol_PCG_rr
    norm_evol_PCG_rz
    T3(k,counter)=toc;
    S3(k,counter) = st+1;
    A_3(:,k)=Hom_parameter(C,A,G,E)% Compute homogenized parameter
    A3(counter)=A_3(1,1);
    
end

% error 
  NoP(counter)=N_1*N_2;
  counter=counter+1;
end

% % Plot estimates
 rel_estim3=(estim3./estim3(1));%.^(1/2);
% rel_estimg=estimg./estimg(1);
 rel_estimp=estimp./estimp(1);
% 
%  figure 
%  hold on
% plot((1:numel(estim3)),abs(rel_estim3),'.')
%    plot((1:numel(estim3)),abs(rel_estim3./(1-tau)),'.') 
%   
%  plot((1:numel(estimg)),abs(rel_estimg),'--x')
%   plot((1:numel(estimg)),abs(rel_estimg./(1-tau)),'--x')
%   
%   plot((1:numel(estimp)),abs(rel_estimp),'--o')
%   plot((1:numel(estimp)),abs(rel_estimp./(1-tau)),'--o')
%   
% set(gca, 'XScale', 'linear', 'YScale', 'log');
% legend('rel estim3 lower','rel estim3 upper','rel estimg lower','rel estimg upper','rel estimp lower','rel estimp upper')
% title('Energetic norm estimates')

%% Plot residuals
 figure 
 hold on
 plot((1:numel(estim3)),abs(rel_estim3),'.')
  plot((1:numel(estim3)),abs(rel_estim3./(1-tau)),'.')
 plot((1:S3(1,end)),norm_evol_PCG_rr,'--')
 plot((1:S3(1,end)),norm_evol_PCG_DrDr,'--')
 plot((1:S3(1,end)),norm_evol_PCG_rz,'--')
 
 
 plot((1:numel(estimp)),abs(rel_estimp),'.')
 plot((1:numel(estimp)),abs(rel_estimp./(1-tau)),'.')
  plot((1:Sp(1,end)),norm_evol_proj_grad,'-.x')
  plot((1:Sp(1,end)),norm_evol_proj_energy,'-.o')
 
set(gca, 'XScale', 'linear', 'YScale', 'log');
legend('DB estim || e_k ||_K',' DB estim || e_k ||_K','DB || r ||','DB  ||Dr||', 'DB  ||r||_M', ...
    'GB estim || e_k ||_K',' GB estim || e_k ||_K','GB  ||Dr||', 'GB ||r||_M')
title('Residuals PCG left')

%% Plot residuals
 figure 
 hold on
 plot((1:S3(1,end)),norm_evol_PCG_rr,'--')
 plot((1:S3(1,end)),norm_evol_PCG_DrDr,'--')
 plot((1:S3(1,end)),norm_evol_PCG_rz,'--')

 
 %plot((1:Sg(1,end)),norm_evolg,'-o') %'DB modif',
 
  plot((1:Sp(1,end)),norm_evol_proj_grad,'-.')
  plot((1:Sp(1,end)),norm_evol_proj_energy,'-.')
  
%  plot((1:Spc(1,end)),norm_evolpc,'-.')%'GB modif',
  plot((1:Spm(1,end)),norm_evolpm,'-.')

set(gca, 'XScale', 'linear', 'YScale', 'log');
legend('DB PCG || r ||','DB PCG  ||Dr||', 'DB PCG ||r||_M' ...
    ,'GB orig ||Dr|| ','GB orig ||r||_M','GB modif 2')
title('Residuals')
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
% legend('sol_norm  proj- grad ','sol_norm proj - sym')
% title('solution difference')
%% Plot solution norm

 figure 
 hold on
 
  %plot((1:Sg(1,end)),abs(sol_normg(end)-sol_normg),'x')
  
 plot((1:Sp(1,end)),abs(sol_normp(end)-sol_normp),'-.*')
  plot((1:Spm(1,end)),abs(sol_normpm(end)-sol_normpm),'-.o')

 
set(gca, 'XScale', 'linear', 'YScale', 'log');
legend('DB modif','GB orig','GB modif')
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
% legend('in G','Symetric','Left','proj','grad norm')


close all

 %% Plot A function