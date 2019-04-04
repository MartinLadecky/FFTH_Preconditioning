function [c_1,st,norm_evol] = CGP_solver_left(A,G,c_0,E,steps,toler,M_f)
%% input
% A   [N_bf2,N_bf1,2,2] -matrix of material parameters in every point of grid
% G_n [N_bf2,N_bf1,2]   -matrix of coeficients of 1st derivative
% c_0 [N_bf2,N_bf1]     -initial solution
% E   [2,1]             -vector [0;1] or [1;0]
% toler [1]             -relative tolerance
% steps [1]             -max number of steps
%% Output
% c_0 [N_bf2,N_bf1]     -vector of solution
% st  [1]               -number of stepts
%% 
% toler           % relative toleranceF
% steps           % -max number of steps
norm_evol=0;
M_0 = LHS_freq(A,c_0,G); % System matrix * initial solution(x_0=c_0)
b_0 = RHS_freq(A,E,G);     % Right hand side vector

r_0 = b_0-M_0; % x_0=0
 nr0 = norm(r_0.*r_0,'fro')
 
z_0 = r_0./M_f; % solve lin system rM_0=M_f^(-1)*r_0: rM_0 is idagonal matrix

p_0 = z_0;

for st = 1:steps
    Ap_0 = LHS_freq(A,p_0,G);
    
    
    alfa_0 = real(sum(sum((z_0.')'.*r_0))/sum(sum((p_0.')'.*Ap_0)));
    
  %  M_1 = M_1./M_f;
%     M_1 = LHS_freq(A,p_0,G_n) ;
%     M_1 = M_1./M_f;
%     alfa_01= real(sum(sum((r_0.')'.*r_0))/sum(sum((p_0.')'.*M_1)));
%    
    c_1 = c_0 + alfa_0.*p_0;
    r_1 = r_0-alfa_0*Ap_0;
    
        if (norm(r_1.*r_1,'fro')/nr0<toler)
            % c_1 = c_0; 
            break; 
        end    
    z_1=r_1./M_f;
    
    beta_1 = real(sum(sum((z_1.')'.*r_1))/sum(sum((z_0.')'.*r_0)));
    p_1 = z_1 + beta_1*p_0;    
    %% 
    p_0 = p_1;
    r_0 = r_1;
    z_0 = z_1; 
    c_0 = c_1;
   % norm_evol(st)=norm(c_1,'fro');
end
%c_1=(M_f.^(-1/2)).*c_1

%x_0=fftM_f.^(-1/2)shift(ifft2(ifftshift(c_0)));
end