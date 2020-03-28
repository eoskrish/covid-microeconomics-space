function dydt = ODEsolver(t,y,N,n_age_cat,n_eco_cat)
        
    % extract rate parameters     
    R0 = 2.7;
%     zeta = ones(1,n_age_cat);
    zeta = [0.05,   0.05,   0.10,    0.15,    0.20,    0.25,    0.30,    0.40,    0.50];
    zeta = (1-zeta).';
    
    tl = 5;
    ti = 3;
    th = 4;
    tc = 14;
    
%     all porbabilitues equal
%     m = 0.75*ones(n_age_cat,n_eco_cat);
%     c = 0.24*ones(n_age_cat,n_eco_cat);
%     f = 0.36*ones(n_age_cat,n_eco_cat);
    
    % demographic dependent proababilites
    confirmed = [5 5 10 15 20 25 30 40 50] *1/100;
    m = [1 3 3 3 6 10 25 35 50] * 1/100;
    c = [5 10 10 15 20 25 25 35 45 55] * 1/100;
    f = [30 30 30 30 30 40 40 50 50] * 1/100;
    m = ((1-confirmed.*m)).';
%     m = (1-m).';
    c = (c).';
    f = (f).';    
    
    tmax = -70;
    minM = 0.8;
    tau = 10;
    epsilon = 0;
    %with mitigation
%     beta = R0 * zeta * (minM + (1-minM)*exp(-t/tau)) * (1 + epsilon * cos(2*pi*(t-tmax))) / ti;
    % without mitigation
    beta = R0 * zeta * (1 + epsilon * cos(2*pi*(t-tmax)/ti));
    
    sp = reshape(y,7,n_age_cat,n_eco_cat);
    
    S = zeros(n_age_cat,n_eco_cat);
    E = zeros(n_age_cat,n_eco_cat);
    I = zeros(n_age_cat,n_eco_cat);
    R = zeros(n_age_cat,n_eco_cat);
    C = zeros(n_age_cat,n_eco_cat);
    H = zeros(n_age_cat,n_eco_cat);
    D = zeros(n_age_cat,n_eco_cat);
    
    S(:,:) = sp(1,:,:);
    E(:,:) = sp(2,:,:);
    I(:,:) = sp(3,:,:);
    R(:,:) = sp(4,:,:);
    H(:,:) = sp(5,:,:);
    C(:,:) = sp(6,:,:);
    D(:,:) = sp(7,:,:);
    
    % ODEs
    dSdt = zeros(n_age_cat,n_eco_cat);
    dEdt = zeros(n_age_cat,n_eco_cat);
    dIdt = zeros(n_age_cat,n_eco_cat);
    dRdt = zeros(n_age_cat,n_eco_cat);
    dHdt = zeros(n_age_cat,n_eco_cat);
    dCdt = zeros(n_age_cat,n_eco_cat);
    dDdt = zeros(n_age_cat,n_eco_cat);
    
    totalI = sum(sum(I));
    
%     N = sum(sum(y));
    import_rate = 0;
    
    for j = 1 : n_eco_cat
        for i = 1 : n_age_cat
            dSdt(i,j) = - beta(i) / N * S(i,j) * totalI ;
            dEdt(i,j) = beta(i) / N * S(i,j) * totalI - E(i,j) / tl + import_rate;
            dIdt(i,j) = E(i,j)/tl - I(i,j)/ti ;
            dRdt(i,j) = (1-c(i,j)) * H(i,j)/th + m(i,j) * I(i,j)/ti ;
            dHdt(i,j) = (1-m(i,j))*I(i,j)/ti + (1-f(i,j)) * C(i,j)/tc - H(i,j)/th ;
            dCdt(i,j) = c(i,j)*H(i,j)/th - C(i,j)/tc ;
            dDdt(i,j) = f(i,j) * C(i,j)/tc ;
        end
    end
    
    dydt = [];
    for j = 1 : n_eco_cat
        for i = 1 : n_age_cat
            tmp = [dSdt(i,j) dEdt(i,j) dIdt(i,j) dRdt(i,j) dHdt(i,j) dCdt(i,j) dDdt(i,j)];
            dydt = [dydt tmp];
        end
    end
    dydt = dydt.';
    
end

