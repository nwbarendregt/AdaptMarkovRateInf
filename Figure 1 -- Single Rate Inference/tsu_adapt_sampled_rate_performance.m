clear
load('gamma_samples.mat');

prior = @(x) beta^alpha*x.^(alpha-1).*exp(-beta*x)/gamma(alpha);
Var_tol = [0.1 0.05 0.01 0.005 0.001];
h0_mesh = linspace(0,10);
pn_0 = prior(h0_mesh);

N_samples = NaN(length(h0),length(Var_tol));
MSE = NaN(length(h0),length(Var_tol));

poolObj = parpool(20);
parfor i = 1:length(h0)
    N_samples_i = NaN(1,length(Var_tol));
    MSE_i = NaN(1,length(Var_tol));
    for j = 1:length(Var_tol)
        pn = pn_0;
        E_h0 = trapz(h0_mesh,h0_mesh.*pn);
        V = trapz(h0_mesh,(h0_mesh-E_h0).^2.*pn);
        T = 0;
        while V >= Var_tol(j)
            ub = 10/E_h0;
            Tn = fminbnd(@(T)var_tsu(T,h0_mesh,pn),0,ub);
            T = [T T(end)+Tn];
            if Tn < -log(rand)/h0(i)
                pn = exp(-h0_mesh*Tn).*pn;
                pn = pn/trapz(h0_mesh,pn);
            else
                pn = (1-exp(-h0_mesh*Tn)).*pn;
                pn = pn/trapz(h0_mesh,pn);
            end
            E_h0 = trapz(h0_mesh,h0_mesh.*pn);
            V = trapz(h0_mesh,(h0_mesh-E_h0).^2.*pn);
        end
        N_samples_i(j) = length(T);
        MSE_i(j) = trapz(h0_mesh,(h0_mesh-h0(i)).^2.*pn);
    end
    N_samples(i,:) = N_samples_i;
    MSE(i,:) = MSE_i;
    disp(i)
end
delete(poolObj);

N_samples_adapt = N_samples;
MSE_adapt = MSE;

save('tsu_adapt_sampled_rate_performance_data.mat');

figure
loglog(Var_tol,mean(N_samples_adapt,1),'k-o','linewidth',5)
xlabel('Variance Tolerance'); ylabel('Average Number of Samples');
figure
loglog(Var_tol,mean(MSE_adapt,1),'k-o','linewidth',5)
xlabel('Variance Tolerance'); ylabel('Average MSE');