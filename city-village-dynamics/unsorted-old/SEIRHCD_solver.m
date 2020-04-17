function dydt = SEIRHCD_solver(t,y,n_age_cat,n_eco_cat,import_rate)
        
    % extract rate parameters     
    R0 = 2.7;
    
    %if degree of isolation is taken to 1 for all age categories
    zeta = ones(1,n_age_cat);
     
    % if degree of isolation are different for different age categories
    % zeta = [0.05, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50];
    % zeta = (1-zeta).';
    
    % time spent in different states (jump rates from a state)
    ti = 3;     % time spent in infectious state
    tp = 15;    % 1/rate of S ->  D due to impoverished conditions
    
    % demographic dependent proababilites
    conf = [5 5 10 15 20 25 30 40 50] *1/100;    % don't fully understand why this is necessary
    m = [1 3 3 3 6 10 25 35 50] * 1/100;         % values taken from COVID scenarios
    c = [5 10 10 15 20 25 35 45 55] * 1/100;     % values taken from COVID scenarios
    f = [30 30 30 30 30 40 40 50 50] * 1/100;    % values taken from COVID scenarios
    
    m = (1 - conf.*m);        % still don't understand this but 
                              % borrowed directly from neher lab's code
                              % possibly bcoz only the ones that are
                              % confirmed as Covid+ve  are admitted to
                              % a hospital
    c = (c);               % converting to column vector
    f = (f).';             % converting to column vector
    
    f = repmat(f,1,n_eco_cat);
    
    p = [1 1 0];
    p = repmat(p,n_age_cat,1);
    
    % seasonal forcing parameters
    tmax = -70;     % when was the peak (days), adjust according 
                    % to start date
    epsilon = 0;    % forcing amplitude the 
    
    % beta with mitigation - exponential 
    % minM = 0.4;     % by what factor does mitigation finally reduces 
    % tau = 5;        % time period of relaxation    
    % beta = R0 * zeta * (minM + (1-minM)*exp(-t/tau)) * (1 + epsilon * cos(2*pi*(t-tmax))) / ti;
    
    % beta with mitigation - linear 
%     tau = 30;
%     minM = 0.4;
%     if t < tau
%         beta = R0 * zeta * (1 - (1-minM)*t/tau) * (1 + epsilon * cos(2*pi*(t-tmax))) / ti;
%     else 
%         beta = R0 * zeta * minM * (1 + epsilon * cos(2*pi*(t-tmax))) / ti;
%     end
%     
    contact_map = R0*ones(n_age_cat,n_eco_cat,n_age_cat,n_eco_cat);
    
    % without mitigation
    % beta = R0 * zeta * (1 + epsilon * cos(2*pi*(t-tmax)))/ti;       %why are they dividing R0 by ti  
    
    sp = reshape(y,4,n_age_cat,n_eco_cat);
    
    S = zeros(n_age_cat,n_eco_cat);     % susceptible
    E = zeros(n_age_cat,n_eco_cat);     % exposed
    I = zeros(n_age_cat,n_eco_cat);     % infectious
    R = zeros(n_age_cat,n_eco_cat);     % recovered
    H = zeros(n_age_cat,n_eco_cat);     % hospitalised
    C = zeros(n_age_cat,n_eco_cat);     % critical
    D = zeros(n_age_cat,n_eco_cat);     % casualties
    
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
    N = sum(y);
    
    import_S = zeros(n_age_cat,n_eco_cat);
    import_E = zeros(n_age_cat,n_eco_cat);
    import_I = zeros(n_age_cat,n_eco_cat);
    import_R = zeros(n_age_cat,n_eco_cat);
    
    import_S(2:4,2) = 0.8*import_rate*[0.50; 0.30; 0.20];
    import_E(2:4,2) = 0.1*import_rate*[0.50; 0.30; 0.20];
    import_I(2:4,2) = 0.1*import_rate*[0.50; 0.30; 0.20];
    
    beta = zeros(n_age_cat,n_eco_cat);
    
    for j = 1 : n_eco_cat
        for i = 1 : n_age_cat
            beta(:,:) = contact_map(i,j,:,:);
            dSdt(i,j) = - beta(i) / N * S(i,j) * totalI + import_S(i,j);
            dEdt(i,j) = beta(i) / N * S(i,j) * totalI - E(i,j) / tl + import_E(i,j);
            dIdt(i,j) = E(i,j)/tl - I(i,j)/ti + import_I(i,j);
            dRdt(i,j) = (1-c(i,j)) * H(i,j)/th + m(i,j) * I(i,j)/ti + import_R(i,j);
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
    
    check = y + dydt*0.1;
    dydt(check<0) = 0;
    
end

