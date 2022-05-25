% var_bin_struct.m
% Computes expected posterior covariance for arbitrary Markov chain with 
% binary transition rates using Eq.(10) from Barendregt et al., 2022. 
function V = var_bin_struct(T,A,pn,Xp,Tp,states,n_nodes,max_edges)

C  = zeros(max_edges); h_ind = 1:max_edges;
% Calculate conditional transition probabilities:
for k = 1:n_nodes
    p_k = NaN(length(states),1);
    for j = 1:length(states)
        p_k(j) = bin_struct_trans_prob(A(states(j,:)),T-Tp,k,Xp);
    end
    p_k = reshape(p_k,2*ones(1,max_edges));

    for j = 1:length(states(1,:))
        if j ~= 2
            % Calculate integrand of covariance (taking transpose as necessary):
            Int = [0 1]-sum([0 1].*squeeze(sum(p_k.*pn.^2,h_ind(h_ind ~= j)))')/...
                sum(p_k.*pn.^2,'all');
            % Calculate entries of covariance matrix (taking transpose as necessary): 
            C(j,j) = C(j,j)+sum(Int.^2.*squeeze(sum(p_k.*pn.^2,h_ind(h_ind ~= j)))');
        else
            % Calculate integrands of covariance:
            Int = [0 1]-sum([0 1].*squeeze(sum(p_k.*pn.^2,h_ind(h_ind ~= j))))/...
                sum(p_k.*pn.^2,'all');
            % Calculate entries of covariance matrix:
            C(j,j) = C(j,j)+sum(Int.^2.*squeeze(sum(p_k.*pn.^2,h_ind(h_ind ~= j))));
        end
    end
end
V = det(C);
end