function [c_1,st,norm_evol_rr, norm_evol_rz, estim,  norm_sol,e_norm_error] = solver_PCG_projection_left(A,G,c_0,E,steps,toler,M_f,C_ref,tau)

    grad_c_0=G.*c_0; % Gradient [N_bf2,N_bf1,2]
    c_0=grad_c_0;
    %grad_c_0=grad_c_0-mean(mean(grad_c_0));
     % Gradient in real space
    norm_sol(1)=sqrt(scalar_product_grad(c_0,c_0));
    
    
    M_0 =Material_data_mul(A,grad_c_0)
    Fb_0=cat(3,A(:,:,1,1).*E(1)+A(:,:,1,2).*E(2),...
              A(:,:,2,1).*E(1)+A(:,:,2,2).*E(2));
          
    b_0=fftshift(fft2(ifftshift(Fb_0)));
    
    
    r_0 = b_0-M_0;
    
    Dr_0=G.*r_0;
    Dr_0=Dr_0(:,:,1)+Dr_0(:,:,2)
        
%     M_0 = Projection_plain(grad_c_0,G,M_f); % Projection*Material*grad
%     b_0 = Projection_plain(E,G,M_f); % Right hand side vector
  % MM=M_f.^-1;
      
    MDr=Dr_0./M_f
  %  MDr=MM.*Dr_0
    nz0r0_test = sqrt(scalar_product(Dr_0,MDr))
    
    DMDr=G.*MDr
    nz0r0_test = sqrt(scalar_product_grad(r_0,DMDr))
    
    rDMDr=(r_0(:,:,1).')'.*DMDr(:,:,1)+(r_0(:,:,2).')'.*DMDr(:,:,2)
    sqrt(sum(sum(rDMDr)))
    
    z_0=Projection_plain(r_0,G,M_f);
    
    nz0=sqrt(scalar_product_grad(r_0,z_0))
    norm_evol_rr(1)=nz0/nz0
    
    nz0r0 =sqrt(scalar_product_grad_energy_ref(z_0,r_0,C_ref));
    norm_evol_rz(1)=nz0r0/nz0r0;
    
    
    p_0 = z_0;
    
    k=1;
    d=0;
    C_precise=load('experiment_data/sol_10_10.mat');
    estim=0;
    for st = 1:steps
        Ap_0 = Material_data_mul(A,p_0);%Projection(A,p_0,G,M_f);%LHS_freq(A,p_0,G); 

        z_0r_0=scalar_product_grad(r_0,z_0);
        p_0Ap_0=scalar_product_grad(p_0,Ap_0);
        alfa_0 = z_0r_0/p_0Ap_0;
        
        c_1 = c_0 + alfa_0.*p_0; % [N_bf2,N_bf1] +  scalar*[N_bf2,N_bf1,2]

        norm_sol(st+1)=sqrt(scalar_product_grad_energy(c_1,c_1,A));

        r_1 = r_0 - alfa_0*Ap_0 ;
        
        z_1=Projection_plain(r_1,G,M_f);

        nz1=sqrt(scalar_product_grad(r_1,z_1));
        norm_evol_rr(st+1)=nz1/nz0;
        
        
        nz1r1=sqrt(scalar_product_grad_energy_ref(z_1,r_1,C_ref));
        norm_evol_rz(st+1)=nz1r1/nz0r0;
        
            if ( norm_evol_rr(st+1)<toler)
                % c_1 = c_0; 
                break; 
            end  
        
    
        z_1r_1=scalar_product_grad(r_1,z_1)
        beta_1 =  z_1r_1/z_0r_0
        p_1 = z_1 + beta_1*p_0;
        
%         %error estimates
%         Delta(st)=real(alfa_0*z_1r_1);
%         curve(st)=0;
%         curve=curve+Delta(st);
%         if st >1
%             S=findS(curve,Delta,k);
%             num = S*Delta(st);
%             den = sum(Delta(k:st-1));
%           
%             while (d>= 0) && (num/den<= tau)
%                
%                 delay(k)=d;
%                 estim(k)=den;
%                 k=k+1;
%                 d=d-1;
%                 den=sum(Delta(k:st-1));
%              
%             end
%             d=d+1;
%         end
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