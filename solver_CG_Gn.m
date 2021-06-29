function [c_1,st,norm_evol] = solver_CG_Gn(A,G_n,c_0,E,steps,toler,M_m)
%% input
% A   [N_bf2,N_bf1,2,2]- matrix of material parameters in every point of grid
% G_n [N_bf2,N_bf1,2]   -matrix of coeficients of 1st derivative
% c_0 [N_bf2,N_bf1]     -initial solution
% E   [2,1]             -vector [0;1] or [1;0]
% toler [1]             -relative tolerance
% steps [1]             -max number of steps
%% Output
% c_0 [N_bf,1]          -vector of solution
% st  [1]               -number of stepts
%% 
% toler           % relative tolerance
% steps           % -max number of steps

c_1=c_0./M_m;
M_0 = LHS_freq(A,c_1,G_n);
b_0 = RHS_freq(A,E,G_n);

r_0 = b_0-M_0; % 

%nr0 = norm(r_0.*(M_m.^-2).*r_0,'fro');
nr0 = sqrt(sum(sum(abs(r_0.*(M_m.^-2).*r_0))));
% z_0 = r_0;%./M_f;
% 
% grad_z_0=G_plain.*(z_0);
% nr0 =sqrt(scalar_product_grad(grad_z_0,grad_z_0))
norm_evol(1)=nr0/nr0;


p_0 = r_0;
for st = 1:steps
    
    M_1 = LHS_freq(A,p_0,G_n);
    alfa_0 =sum(sum((r_0.')'.*r_0))/sum(sum((p_0.')'.*M_1));
    x_1 = c_1 + alfa_0*p_0;
   
    r_1 = r_0 - alfa_0*M_1;
    
%     grad_Mr_1=G_plain.*(r_1);
%     norm_evol(st)=sqrt(scalar_product_grad(grad_Mr_1,grad_Mr_1))/nr0;
    norm_evol(st+1)=sqrt(sum(sum(abs(r_1.*(M_m.^-2).*r_1))))/nr0;  %sum(sum(r_1.*(M_m.^-2).*r_1))/nr0;
                    
    if ( norm_evol(st+1)<toler)
         c_1 = x_1; 
         break; 
    end
    
    beta_0 =sum(sum((r_1.')'.*r_1))/sum(sum((r_0.')'.*r_0));
    p_1 = r_1 + beta_0*p_0;
    %%
    p_0 = p_1;
    r_0 = r_1;
    c_1 = x_1 ;
end
c_1 =c_1.*M_m;
%x_0=fftshift(ifft2(ifftshift(c_0)));
end