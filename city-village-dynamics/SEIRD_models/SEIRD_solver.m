function [dydt,fluxSD,fluxID] = SEIRD_solver(economic_module,dt,t,y,CM,zeta,n,f,imports)
      
    n_age_cat = 3;
    n_eco_cat = 3;
    
    % Epidemeologic parameters    
%     R0 = 2.7;
    
    % time spent in different states (jump rates from a state)
    te = 5;     % time spent in exposed state (delay in becoming infectious)
    ti = 42;     % time spent in infectious state
    th = 4;     % time spent in hospital before recovery or transfer to ICU
    tc = 14;    % time spent in ICU 

    prob_factor = zeta * 0.019;               % this is the beta in Ronojoy's work
    gamma = 0;
    % redefining ti
%     ti = 0.8*ti + 0.15*(ti+th) + 0.05*(ti+th+tc);
    
    rate_SD = (2.5 / 365 *10^6) / (0.8 * 1.4 * 10^9 );
    tp = 1/rate_SD;    % 1/rate of S ->  D due to impoverished conditions

    % effect of mitigation on economic S -> D rate 
    alpha = (1/zeta)^n;
    if economic_module == 0
        alpha = 0;  % if you want to turn off effect on economic situation
    end
    % who's affected by poverty? all poor people equally
    p = alpha * [1/6 1/6 0; 1/6 1/6 0; 1/6 1/6 0];
     
    sp = reshape(y,5,n_age_cat,n_eco_cat);
    
    S = zeros(n_age_cat,n_eco_cat);     % susceptible
    E = zeros(n_age_cat,n_eco_cat);     % exposed
    I = zeros(n_age_cat,n_eco_cat);     % infectious
    R = zeros(n_age_cat,n_eco_cat);     % recovered
    D = zeros(n_age_cat,n_eco_cat);     % casualties
    
    S(:,:) = sp(1,:,:);
    E(:,:) = sp(2,:,:);
    I(:,:) = sp(3,:,:);
    R(:,:) = sp(4,:,:);
    D(:,:) = sp(5,:,:);
   
    N_ae = zeros(3,3);
    N_ae(:,:) = sum(sp,1);
%     size(N_ae)
    
    % ODEs
    dSdt = zeros(n_age_cat,n_eco_cat);
    dEdt = zeros(n_age_cat,n_eco_cat);
    dIdt = zeros(n_age_cat,n_eco_cat);
    dRdt = zeros(n_age_cat,n_eco_cat);
    dDdt = zeros(n_age_cat,n_eco_cat);
    
    fluxSD = zeros(n_age_cat,n_eco_cat);
    fluxID = zeros(n_age_cat,n_eco_cat);
    
    totalE = sum(sum(E));
    totalI = sum(sum(I));
    N = sum(y);
    
    no_of_ICUs = 7000;
    if totalI > no_of_ICUs
        f = 1.25*f;
    end
    
    import_S = zeros(n_age_cat,n_eco_cat);
    import_E = zeros(n_age_cat,n_eco_cat);
    import_I = zeros(n_age_cat,n_eco_cat);
    import_R = zeros(n_age_cat,n_eco_cat);
    
    import_S(:,:) = imports(:,:,1);
    import_E(:,:) = imports(:,:,2);
    import_I(:,:) = imports(:,:,3);
    import_R(:,:) = imports(:,:,4);

    beta = zeros(n_age_cat,n_eco_cat);
    
    for j = 1 : n_eco_cat
        for i = 1 : n_age_cat
            beta(:,:) = prob_factor * CM(i,j,:,:);
            fluxSD(i,j) = p(i,j)* S(i,j)/tp;
%             fluxSE(i,j) = S(i,j) * ( sum(sum(beta.*E./ N_ae)) + sum(sum(gamma*beta.*I./ N_ae)) );
            fluxID(i,j) = f(i,j) * I(i,j)/ti ;
            dSdt(i,j) = - S(i,j) * ( sum(sum(beta.*E./ N_ae)) + sum(sum(gamma*beta.*I./ N_ae)) ) - p(i,j) * S(i,j)/tp + import_S(i,j);
            dEdt(i,j) = S(i,j) * ( sum(sum(beta.*E./ N_ae)) + sum(sum(gamma*beta.*I./ N_ae)) ) - E(i,j)/te + import_E(i,j);
            dIdt(i,j) = E(i,j)/te - I(i,j)/ti + import_I(i,j);
            dRdt(i,j) = ( 1 - f(i,j) ) * I(i,j)/ti + import_R(i,j);
            dDdt(i,j) = f(i,j) * I(i,j)/ti + p(i,j)* S(i,j)/tp ;
        end
    end
    
    dydt = [];  
    for j = 1 : n_eco_cat
        for i = 1 : n_age_cat
            tmp = [dSdt(i,j) dEdt(i,j) dIdt(i,j) dRdt(i,j) dDdt(i,j)];
            dydt = [dydt tmp];
        end
    end
    
    check = y + dydt*dt;
    dydt(check<0) = 0;
end

