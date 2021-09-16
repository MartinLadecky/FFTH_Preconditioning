function [c_1,st,norm_evol_rr, norm_evol_rz,   norm_sol, e_norm_error]...
    = solver_GB_PCG_right(A,G,c_0,E,M,C_ref,steps,toler)

    grad_c_0=G.*c_0; % 
    c_0=grad_c_0;

    norm_sol(1)=sqrt(scalar_product_grad_energy(c_0,c_0,A));
    
    
    M_0 =apply_material_data(A,grad_c_0);
    b_0= apply_material_data(A,E);
    
    
    r_0 = b_0-M_0;

    
    z_0=Projection_plain_Cref(r_0,G,M,C_ref);

   
    Dr_0=conj(G).*r_0;
    Dr_0=Dr_0(:,:,1)+Dr_0(:,:,2) ;
    nr0=sqrt(scalar_product(Dr_0,Dr_0));
    
    %nr0=sqrt(scalar_product_grad(r_0,r_0));
    norm_evol_rr(1)=nr0/nr0;
    
    nz0r0 =sqrt(scalar_product_grad(r_0,z_0));
    norm_evol_rz(1)=nz0r0/nz0r0;
    
    
    p_0 = z_0;
    
    for st = 1:steps
        Ap_0 = apply_material_data(A,p_0);
        
        z_0r_0=scalar_product_grad(r_0,z_0);
        p_0Ap_0=scalar_product_grad(p_0,Ap_0);
        alfa_0 = z_0r_0/p_0Ap_0;
        
        c_1 = c_0 + alfa_0.*p_0;
        
        norm_sol(st+1)=sqrt(scalar_product_grad_energy(c_1,c_1,A));

        r_1 = r_0 - alfa_0*Ap_0;
        
        z_1=Projection_plain_Cref(r_1,G,M,C_ref);

        Dr_1=conj(G).*r_1;
        Dr_1=Dr_1(:,:,1)+Dr_1(:,:,2) ;
        nr1=sqrt(scalar_product(Dr_1,Dr_1));
         
        %nr1=sqrt(scalar_product_grad(r_1,r_1));
        norm_evol_rr(st+1)=nr1/nr0;   
        
        nz1r1=sqrt(scalar_product_grad(r_1,z_1));
        norm_evol_rz(st+1)=nz1r1/nz0r0;
        
            if ( norm_evol_rr(st+1)<toler)
                break; 
            end           
    
        z_1r_1=scalar_product_grad(r_1,z_1);
        beta_1 =  z_1r_1/z_0r_0;
        p_1 = z_1 + beta_1*p_0;
        
        %% reiniticialisation
        p_0 = p_1;
        r_0 = r_1;
        z_0 = z_1; 
        c_0 = c_1;
        
    end
end
