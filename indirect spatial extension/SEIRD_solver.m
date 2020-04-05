function dydt = SEIRD_solver(t,y,n_age_cat,n_eco_cat,imports)
        
    % extract rate parameters     
    R0 = 2.7;
    
    %if degree of isolation is taken to 1 for all age categories
%     zeta = ones(1,n_age_cat);
    zeta = 0.75;	
    
    n = 1;
    alpha = (1/zeta)^n;
%     alpha = 0;  % if you want to turn off effect on economic situation
    
    % if degree of isolation are different for different age categories
    % zeta = [0.05, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50];
    % zeta = (1-zeta).';
    
    % time spent in different states (jump rates from a state)
    te = 5;         % time spent in infectious state
    ti = 14;        % time spent in infectious state
    tp = 10*365;    % 1/rate of S ->  D due to impoverished conditions
    
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
    c = (m.*c);                   % converting to column vector
    f = (c.*f).';                 % converting to column vector
    
    %reconstructing fatality for the purposes of 3 age categories, 3
    %economic categories
    age_dist = [17 18.3 17.4 15.6 12.3 9.3 6.3 2.8 0]/100;
    age_dist(end) = 1 - sum(age_dist);
    age_dist = age_dist.';
    f(1) = (f(1)*age_dist(1) + f(2)*age_dist(2))/(age_dist(1) + age_dist(2)); 
    f(2) = sum(f(2:6).*age_dist(2:6))/sum(age_dist(2:6));
    f(3) = sum(f(7:9).*age_dist(7:9))/sum(age_dist(7:9));
    f(4:end) = [];
    f = repmat(f,1,n_eco_cat);
    
    % who's affected by poverty
    p = alpha * [1/6 1/6 0; 1/6 1/6 0; 1/6 1/6 0];
    
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
    contact_map = R0 * zeta * ones(n_age_cat,n_eco_cat,n_age_cat,n_eco_cat) / ti;
    
    % without mitigation
    % beta = R0 * zeta * (1 + epsilon * cos(2*pi*(t-tmax)))/ti;       %why are they dividing R0 by ti  
    
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
    
    % ODEs
    dSdt = zeros(n_age_cat,n_eco_cat);
    dEdt = zeros(n_age_cat,n_eco_cat);
    dIdt = zeros(n_age_cat,n_eco_cat);
    dRdt = zeros(n_age_cat,n_eco_cat);
    dDdt = zeros(n_age_cat,n_eco_cat);
    
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
    
% import_S(1:2,2) = *import_rate*[0.50; 0.30; 0.20];
% import_I(1:2,2) = 0.9*import_rate*[0.50; 0.30; 0.20];
% import_R(1:2,2) = 0.1*import_rate*[0.50; 0.30; 0.20];

    beta = zeros(n_age_cat,n_eco_cat);
    
    for j = 1 : n_eco_cat
        for i = 1 : n_age_cat
            beta(:,:) = contact_map(i,j,:,:);
            dSdt(i,j) = - 1/N * S(i,j) * ( sum(sum(beta.*E)) + sum(sum(beta.*I)) ) - p(i,j) * S(i,j)/tp + import_S(i,j);
            dEdt(i,j) = 1/N * S(i,j) * ( sum(sum(beta.*E)) + sum(sum(beta.*I)) )  - E(i,j)/te + import_E(i,j);
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
%     dydt = dydt.';
    
    check = y + dydt*0.1;
    dydt(check<0) = 0;
    
end

