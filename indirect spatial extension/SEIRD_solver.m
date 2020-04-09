function [dydt,fluxSD,fluxSE] = SEIRD_solver(t,y,contact_map,zeta,n,imports)
      
    n_age_cat = 3;
    n_eco_cat = 3;
    
    % Epidemeologic parameters
    
    R0 = 2.7;
    
    % time spent in different states (jump rates from a state)
    te = 5;     % time spent in exposed state (delay in becoming infectious)
    ti = 3;     % time spent in infectious state
    th = 4;     % time spent in hospital before recovery or transfer to ICU
    tc = 14;    % time spent in ICU 
    
    % redefining ti
    ti = 0.8*ti + 0.15*(ti+th) + 0.05*(ti+th+tc);
    
    tp = 10*365;    % 1/rate of S ->  D due to impoverished conditions
    
    % demographic dependent proababilites
    conf = [5 5 10 15 20 25 30 40 50] *1/100;    % don't fully understand why this is necessary
    m = [1 3 3 3 6 10 25 35 50] * 1/100;         % values taken from COVID scenarios
    c = [5 10 10 15 20 25 35 45 55] * 1/100;     % values taken from COVID scenarios
    f = [30 30 30 30 30 40 40 50 50] * 1/100;    % values taken from COVID scenarios
    
    m = (1 - conf.*m);          % still don't understand this but 
                                % borrowed directly from neher lab's code
                                % possibly bcoz only the ones that are
                                % confirmed as Covid+ve  are admitted to
                                % a hospital
                                
   
    
    c = (m.*c);                 % converting to column vector
    f = (c.*f).';               % converting to column vector
    
    %reconstructing fatality for the purposes of 3 age categories, 3 economic categories
    age_dist = [17 18.3 17.4 15.6 12.3 9.3 6.3 2.8 0]/100;
    age_dist(end) = 1 - sum(age_dist);
    age_dist = age_dist.';
    f(1) = (f(1)*age_dist(1) + f(2)*age_dist(2))/(age_dist(1) + age_dist(2)); 
    f(2) = sum(f(2:6).*age_dist(2:6))/sum(age_dist(2:6));
    f(3) = sum(f(7:9).*age_dist(7:9))/sum(age_dist(7:9));
    f(4:end) = [];
    f = repmat(f,1,n_eco_cat);
    
    
    
    % effect of mitigation on economic S -> D rate 
%     alpha = (1/zeta)^n;
    alpha = 0;  % if you want to turn off effect on economic situation

    % who's affected by poverty
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
    fluxSE = zeros(n_age_cat,n_eco_cat);
    
    totale = sum(sum(E));
    totalI = sum(sum(I));
    N = sum(y);
    
    import_S = zeros(n_age_cat,n_eco_cat);
    import_E = zeros(n_age_cat,n_eco_cat);
    import_I = zeros(n_age_cat,n_eco_cat);
    import_R = zeros(n_age_cat,n_eco_cat);
    
    import_S(:,:) = imports(:,:,1);
    import_E(:,:) = imports(:,:,2);
    import_I(:,:) = imports(:,:,3);
    import_R(:,:) = imports(:,:,4);

    beta = zeros(n_age_cat,n_eco_cat);
    prob_factor = zeta * 0.019;               % this is the beta in Ronojoy's work
    
    for j = 1 : n_eco_cat
        for i = 1 : n_age_cat
            beta(:,:) = prob_factor * contact_map(i,j,:,:);
            fluxSD(i,j) = p(i,j) * S(i,j)/tp;
            fluxSE(i,j) = S(i,j) * ( sum(sum(beta.*E./ N_ae)) + sum(sum(beta.*I./ N_ae)) );
            dSdt(i,j) = - S(i,j) * ( sum(sum(beta.*E./ N_ae)) + sum(sum(beta.*I./ N_ae)) ) - p(i,j) * S(i,j)/tp + import_S(i,j);
            dEdt(i,j) = S(i,j) * ( sum(sum(beta.*E./ N_ae)) + sum(sum(beta.*I./ N_ae)) )  - E(i,j)/te + import_E(i,j);
            dIdt(i,j) = E(i,j)/te - I(i,j)/ti + import_I(i,j);
            dRdt(i,j) = ( 1 - f(i,j) ) * I(i,j)/ti + import_R(i,j);
            dDdt(i,j) = f(i,j) * I(i,j)/ti + p(i,j) * S(i,j)/tp ;
        end
    end
    
    dydt = [];  
    for j = 1 : n_eco_cat
        for i = 1 : n_age_cat
            tmp = [dSdt(i,j) dEdt(i,j) dIdt(i,j) dRdt(i,j) dDdt(i,j)];
            dydt = [dydt tmp];
        end
    end
    
    check = y + dydt*0.1;
    dydt(check<0) = 0;
end

