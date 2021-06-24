function [U1,U2] = precon_decomposition(A,G,N_1,N_2)
% This function create a decomposition of the
% preconditioner with one significant frequence
% into two matrices U1 and U2

%% BEGIN

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

G_1 = G(:,:,1);
G_2 = G(:,:,2);
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
    p =p + p' + diag(m1);     
    P(:,:,r) = p;% zaradi to do glob matice
end  
U1 = zeros(N_2,N_1);

U2 = U1;
for r = 1:N_2    
    m1 = GG0(r,:);
    m2 = GG1(r,:);    
    p = spdiags([m2(1:N_1-1)'],[-1],N_1,N_1);
%     p(1,N_1) = m2(N_1);  % ZDE OPRAVA NEBO NE
    p =   p + p' + diag(m1);

    pom = chol(p);  % p = pom'*pom % chol rozklad
    U1(r,:) = diag(pom); % diagonala chol
    U2(r,1:N_1-1) = diag(pom,1); % mimodiagonala chol
end

end