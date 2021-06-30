function [c_1,st,norm_evol_rr, norm_evol_rz, estim,  norm_sol,e_norm_error] = solver_GP_projection_left(A,G,c_0,E,steps,toler,M_f,C_ref,tau,G_m,C_ref_inv)
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
    %norm_evol=0;
    
    grad_c_0=G_m.*c_0; % Gradient [N_bf2,N_bf1,2]
    c_0=grad_c_0;
    %grad_c_0=grad_c_0-mean(mean(grad_c_0));
     % Gradient in real space
    norm_sol(1)=sqrt(scalar_product_grad(c_0,c_0));
    
    M_0 = Projection(A,grad_c_0,G,M_f); % Projection*Material*grad
    b_0 = Projection_rhs(A,E,G,M_f) ;% Right hand side vector

    r_0 = b_0-M_0; % x_0=0 

    
    %nr0=sqrt(scalar_product_grad(r_0,r_0));

%     Fr_0=fftshift(ifft2(ifftshift(r_0)));
% 
%     CFr_0=cat(3,C_ref_inv(:,:,1,1).*Fr_0(:,:,1)+C_ref_inv(:,:,1,2).*Fr_0(:,:,2),...
%             C_ref_inv(:,:,2,1).*Fr_0(:,:,1)+C_ref_inv(:,:,2,2).*Fr_0(:,:,2));
%     z_0=fftshift(fft2(ifftshift(CFr_0)));    
        
    %z_0 = r_0./C_ref_full; %gradient of z_0, precon is incorp. in the projection
     z_0=r_0;
    nz0=sqrt(scalar_product_grad(z_0,z_0));
    norm_evol_rr(1)=nz0/nz0;
    
    nz0r0 =sqrt(scalar_product_grad_energy_ref(r_0,r_0,C_ref));
    norm_evol_rz(1)=nz0r0/nz0r0;
    
    
    p_0 = z_0;
    
    k=1;
    d=0;
    C_precise=load('experiment_data/sol_10_10.mat');
    estim=0;
    for st = 1:steps
        Ap_0 = Projection(A,p_0,G,M_f);%LHS_freq(A,p_0,G); 
        %alfa_0 = real(sum(sum((z_0.')'.*r_0))/sum(sum((p_0.')'.*Ap_0)));
%         disp('Mean grad_FAp_0')
%         mean(mean(fftshift(ifft2(ifftshift(Ap_0)))))
        
        
        z_0r_0=scalar_product_grad(z_0,r_0);
        p_0Ap_0=scalar_product_grad(p_0,Ap_0);
        alfa_0 = z_0r_0/p_0Ap_0;
        
        c_1 = c_0 + alfa_0.*p_0; % [N_bf2,N_bf1] +  scalar*[N_bf2,N_bf1,2]
%         c_1((end+1)/2,(end+1)/2,:)=0;
        norm_sol(st+1)=sqrt(scalar_product_grad_energy(c_1,c_1,A));
       % e_norm_error(st)=sqrt(scalar_product_grad_energy(c_1-G.*C_precise.C,c_1-G.*C_precise.C,A));
         %norm_sol(st)=sqrt(scalar_product_grad_energy(c_1,c_1,A));
         %disp('Mean c1')
         %mean(mean(fftshift(ifft2(ifftshift(c_1)))))
%         imag(sum(sum(fftshift(ifft2(ifftshift(c_1))))))
% %         
        r_1 = r_0 - alfa_0*Ap_0 ;
        
        z_1=r_1;
%         Fr_1=fftshift(ifft2(ifftshift(r_1)));
% 
%         CFr_1=cat(3,C_ref_inv(:,:,1,1).*Fr_1(:,:,1)+C_ref_inv(:,:,1,2).*Fr_1(:,:,2),...
%             C_ref_inv(:,:,2,1).*Fr_1(:,:,1)+C_ref_inv(:,:,2,2).*Fr_1(:,:,2));
%         z_1=fftshift(fft2(ifftshift(CFr_1)));
        
        nz1=sqrt(scalar_product_grad(z_1,z_1));
        norm_evol_rr(st+1)=nz1/nz0;
        
        
        nz1r1=sqrt(scalar_product_grad_energy_ref(r_1,r_1,C_ref));
        norm_evol_rz(st+1)=nz1r1/nz0r0;
        
            if ( norm_evol_rr(st+1)<toler)
                % c_1 = c_0; 
                break; 
            end  
        
    
        z_1r_1=scalar_product_grad(z_1,r_1);
        beta_1 =  z_1r_1/z_0r_0;
        p_1 = z_1 + beta_1*p_0;
        
        %error estimates
        Delta(st)=real(alfa_0*z_1r_1);
        curve(st)=0;
        curve=curve+Delta(st);
        if st >1
            S=findS(curve,Delta,k);
            num = S*Delta(st);
            den = sum(Delta(k:st-1));
            
            while (d>= 0) && (num/den<= tau)
                delay(k)=d;
                estim(k)=den;
                k=k+1;
                d=d-1;
                den=sum(Delta(k:st-1));
            end
            d=d+1;
        end
        %% 
        p_0 = p_1;
        r_0 = r_1;
        z_0 = z_1; 
        c_0 = c_1;
    end
end


function [S]=findS(curve,Delta,k)
    ind=find((curve(k)./curve) <= 1e-4,1,'last');
    if isempty(ind)
       ind = 1 ;
    end
    S = max(curve(ind:end-1)./Delta(ind:end-1));
end