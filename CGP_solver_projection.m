function [grad_c_1,st,norm_evol] = CGP_solver_projection(A,G,c_0,E,steps,toler,M_f,C_ref)

    grad_c_0=G.*c_0; % Gradient
    FGFx=fftshift(ifft2(ifftshift(grad_c_0))); % Gradient in real space
    
    
    M_0 = Projection(A,FGFx,G,C_ref); % Projection*Material*grad
    
    b_0 = Projection_rhs(A,E,G,C_ref);  % Right hand side vector

    r_0 = b_0-M_0; % x_0=0
    nr0_1 =norm(r_0(:,:,1),'fro');
    nr0_2=norm(r_0(:,:,2),'fro');
    nr0=sqrt(nr0_1^2+nr0_2^2);
    
    z_0 = r_0;

    p_0 = z_0;

    for st = 1:steps
        Ap_0 = Projection(A,p_0,G,C_ref); 
        z_0r_0=scalar_product_grad(z_0,r_0);
        p_0Ap_0=scalar_product_grad(p_0,Ap_0);
        alfa_0 = z_0r_0/p_0Ap_0;

        grad_c_1 = grad_c_0 + alfa_0.*p_0;
        r_1 = r_0-alfa_0*Ap_0;

        nr1_1 =norm(r_1(:,:,1),'fro');
        nr1_2=norm(r_1(:,:,2),'fro');
        nr1=sqrt(nr1_1^2+nr1_2^2);
        norm_evol(st)=nr1/nr0;
            if ( norm_evol(st)<toler)
                % c_1 = c_0; 
                break; 
            end    
        z_1=r_1;

        beta_1 = real(sum(sum(sum((z_1).*r_1)))/sum(sum(sum((z_0).*r_0))));
        p_1 = z_1 + beta_1*p_0;    
        %% 
        p_0 = p_1;
        r_0 = r_1;
        z_0 = z_1; 
        grad_c_0 = grad_c_1;
    end
end

