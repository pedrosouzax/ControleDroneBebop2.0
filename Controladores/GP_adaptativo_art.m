function [mu,cov, Zt,zt_1,yt] = GP_adaptativo_art(sigma_quad,l,zt_1,Zt,yt,i,omega_quad,etol, p_max)

for j=1:4
    for k=1:length(Zt(j,:))
        %kt(k,k) = kernel(sigma_quad, zt_1(j), Zt(:,k),l, length(zt_1(j)));
        Ktt(:,:) = kernel(sigma_quad,Zt(j,k),Zt,l, length(Zt(j,:)));
        %kt_s(k,k) = kernel(sigma_quad, zt_1(j), zt_1(j),l, length(zt_1(j)));
    end
    
        kt(:,:) = kernel(sigma_quad, zt_1(j), Zt,l, length(Zt(j,:)));
        %Ktt(k,k) = kernel(sigma_quad,Zt(j,k),Zt(:,k),l, length(Zt(j,k)));
        kt_s(:,:) = kernel(sigma_quad, zt_1(j), zt_1(j),l, length(zt_1(j)));

    Ct = diag(Ktt) + omega_quad*eye(i);
    Ct  = chol(Ct);
    Bt = Ct*yt;
    mu(j) = Bt'*kt;
    cov = kt_s -kt'*Ct*kt;
    
end


    
%[Zt] = KLITest(Zt, Ktt, kt, kt_s, etol, p_max);

% Ktt = chol(Ktt);
% alpha = Ktt*kt;
% gama = kt_s - kt'*alpha;
 
%    
%    mu(j,:) = Bt'*kt;
% %Bt = Ct*yt' % yt = valor conhecido da função 
% 
% cov = kt_s -kt'*Ct*kt;

%Zt(:,i) = zt_1; 
%yt(:,i) = mu;

mu;
end

