% var_bd.m
% Computes expected posterior covariance determinant using Eq. (10) from
% Barendregt et al., 2022.
function V = var_bd(T,mu,lambda,N_trunc,pn,Xp,Tp)
C = zeros(2,2);
h = [mu;lambda]; ind = [1 2];
[M,L] = meshgrid(mu,lambda); H = cat(3,M,L);
% Calculate transition probabilities conditioned on previous measurement:
p = NaN(length(mu),length(lambda),Xp+N_trunc);
for i = 1:length(mu)
    for j = 1:length(lambda)
        for n = 0:(Xp+N_trunc)
            p(i,j,n+1)= bd_trans_prob(mu(i),lambda(j),n,Xp,T,Tp);
        end
    end
end

for i = 1:2
    for j = 1:2
        for n = 1:(Xp+N_trunc+1)
            p_k = (p(:,:,n))';
            if i==j
                % Compute integrands for diagonal entries:
                Int = H(:,:,i)-trapz(h(i,:),H(:,:,i).*trapz(h(ind(ind~=i),:),p_k.*pn.^2,i),ind(ind~=i))/...
                    trapz(h(1,:),trapz(h(2,:),p_k.*pn.^2,1),2);
                % Compute diagonal entries of covariance matrix:
                C(i,j) = C(i,j)+trapz(h(1,:),trapz(h(2,:),Int.^2.*p_k.*pn.^2,1),2);
            else
                % Compute integrands for off-diagonal entries:
                Int_i = H(:,:,i)-trapz(h(i,:),H(:,:,i).*trapz(h(j,:),p_k.*pn.^2,i),j)/...
                    trapz(h(1,:),trapz(h(2,:),p_k.*pn.^2,1),2);
                Int_j = H(:,:,j)-trapz(h(j,:),H(:,:,j).*trapz(h(i,:),p_k.*pn.^2,j),i)/...
                    trapz(h(1,:),trapz(h(2,:),p_k.*pn.^2,1),2);
                % Compute off-diagonal entries of covariance matrix:
                C(i,j) = C(i,j)+trapz(h(1,:),trapz(h(2,:),Int_i.*Int_j.*p_k.*pn.^2,1),2);
            end
        end
    end
end

V = det(C);
end