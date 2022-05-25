% bin_struct_connect.m
% Constructs infinitesimal generator matrix for an arbitrary finite-state
% Markov chain with binary transition rates:
function A = bin_struct_connect(h,n_nodes)
A = NaN(n_nodes);
for i = 1:n_nodes
    % Set off-diagonal entries based on connectivity vector:
    A(i,1:n_nodes ~= i) = h(1+(n_nodes-1)*(i-1):1+(n_nodes-1)*(i-1)+(n_nodes-2));
    A(i,i) = -sum(A(i,:),'omitnan');
end
end