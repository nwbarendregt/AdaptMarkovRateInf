% var_ring.m
% Calculates expected posterior covariance for a ring chain of state using
% Eq. (10) from Barndregt et al., 2022.
function V = var_ring(T,hp,hm,A,m,pn,Xp,Tp)
C = zeros(2,2);
h = [hp;hm]; ind = [1 2];
[Hp,Hm] = meshgrid(hp,hm); H = cat(3,Hp,Hm);
% Calculate state transition probabilities:
p = NaN(length(hp),length(hm),m);
for i = 1:length(hp)
    for j = 1:length(hm)
        for n = 1:m
            p_ij = ring_trans_prob(A(hp(i),hm(j)),T-Tp,Xp);
            p(i,j,n) = p_ij(n);
        end
    end
end

for i = 1:2
    for j = 1:2
        for n = 1:m
            p_k = (p(:,:,n))';
            if i==j
                % Compute integrands for diagonal entries of covariance
                % matrix:
                Int = H(:,:,i)-trapz(h(i,:),H(:,:,i).*trapz(h(ind(ind~=i),:),p_k.*pn.^2,i),ind(ind~=i))/...
                    trapz(h(1,:),trapz(h(2,:),p_k.*pn.^2,1),2);
                % Compute diagonal entries of covariance matrix:
                C(i,j) = C(i,j)+trapz(h(1,:),trapz(h(2,:),Int.^2.*p_k.*pn.^2,1),2);
            else
                % Compute integrands for off-diagonal entries of covariance
                % matrix:
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