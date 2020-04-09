clear all
close all

data = load('IndiaData.csv');
data(1:42) = [];
length(data)
tdata = [0:1:length(data)-1];

%% Simulation parameters

n_days = 360;         % observation period in number of days

N = 1.4*10^9;            % city population
InitNumCases = data(1);      % initial number of cases in the city

Nv = 5000;              % village population
num_villages = 10;      % number of villages

T_lockdown = 21;         % time point of lockdown initiation
zeta_old = 1;           % without mitigation
zeta_new = 0.65;         % by how much is R0 reduced with mitigation

% total people migrated from city to villages in response to migration 
% = import_rate * tp_migration
import_rate = 0;     % rate of migration from the city after lockdown (#people per day)
tp_migration = 14;      % time period of migration from city to villages 

% new rate_SD = rate_SD / (zeta)^n
n = 0;                  % effect if mitigation on death rates due to economic reasons

% contact map
dummy_age_dist = [0.17/2 0.17/2 0.183/2 0.183/2 0.174/2 0.174/2 0.156/2 0.156/2 0.123/2 0.123/2 0.093/2 0.093/2 0.063/2 0.063/2 0.028/2 0.024];

age_CM = zeros(3,3);
base_CM = load('India-contact-map.csv');                % work by Kiesha Prem, Alex R Cook, Mark Jit
base_CM = base_CM .* repmat(dummy_age_dist.',1,16);
age_CM(1,1) = sum(sum(base_CM(1:3,1:3)));
age_CM(1,2) = sum(sum(base_CM(1:3,4:12)));
age_CM(1,3) = sum(sum(base_CM(1:3,13:16)));
age_CM(2,1) = sum(sum(base_CM(4:12,1:3)));
age_CM(2,2) = sum(sum(base_CM(4:12,4:12)));
age_CM(2,3) = sum(sum(base_CM(4:12,13:16)));
age_CM(3,1) = sum(sum(base_CM(13:16,1:3)));
age_CM(3,2) = sum(sum(base_CM(13:16,4:12)));
age_CM(3,3) = sum(sum(base_CM(13:16,13:16)));

% definition used by Neher Lab
%     contact_map = R0 * zeta * ones(n_age_cat,n_age_cat,n_eco_cat,n_eco_cat) / ti;

% definition used by Ronojoy Adhikari
% 1. within an economic group, contact map is as given by age_CM
% 2. rich-poor interactions are limited to young
contact_map = zeros(3,3,3,3);
contact_map(:,1,:,1) = age_CM;
contact_map(:,1,:,2) = age_CM;               % poor - poor
contact_map(2,1,2,3) = 0.1*age_CM(2,2);      % young,rich - young,poor
contact_map(:,2,:,1) = age_CM;
contact_map(:,2,:,2) = age_CM;
contact_map(2,2,2,3) = 0.1*age_CM(2,2);
contact_map(2,3,2,1) = 0.1*age_CM(2,2);
contact_map(2,3,2,2) = 0.1*age_CM(2,2);
contact_map(:,3,:,3) = age_CM;

%% City and Village population 

% age distribution is the same in the city and villages (taken from cumulative age dist of India)
age_dist = [0.26 0.64 0.1];
age_dist = age_dist.';

n_age_cat = length(age_dist);   % 0-14, 15-59,  60+ 
n_eco_cat = 3;                  % immobile poor, mobile poor, rich

% declare S, E, I, R, D variables in the city and villages
% City:
S = zeros(n_age_cat, n_eco_cat);
E = zeros(n_age_cat, n_eco_cat);
I = zeros(n_age_cat, n_eco_cat);
R = zeros(n_age_cat, n_eco_cat);
D = zeros(n_age_cat, n_eco_cat);
% Village:
Sv = zeros(n_age_cat, n_eco_cat);
Ev = zeros(n_age_cat, n_eco_cat);
Iv = zeros(n_age_cat, n_eco_cat);
Rv = zeros(n_age_cat, n_eco_cat);
Dv = zeros(n_age_cat, n_eco_cat);

% Initialization: 
% 1. the infections are only in the city 
% 2. only the young (15-59) rich were infected 
% 3. S in villages consisted of 80% immobile poor, 20% rich
% 4. S in city consisted of 35% immobile poor, 35% mobile poor, 30% rich

I(:,1) = 0 * transpose([0 InitNumCases 0]);
I(:,2) = 0 * transpose([0 InitNumCases 0]);
I(:,3) = 1 * transpose([0 InitNumCases 0]);

Sv(:,1) = 0.8 * age_dist * Nv;
Sv(:,2) = 0 * age_dist * Nv;
Sv(:,3) = 0.2 * age_dist * Nv;

S(:,1) = 0.4 * age_dist * (N-InitNumCases);
S(:,2) = 0.4 * age_dist * (N-InitNumCases);
S(:,3) = 0.2 * age_dist * (N-InitNumCases);

% putting all into a vector
y0 = [];
for j = 1 : n_eco_cat
    for i = 1 : n_age_cat
        tmp = [S(i,j) E(i,j) I(i,j) R(i,j) D(i,j)];
        y0 = [y0 tmp];
    end
end
y0 = y0.';

yv0 = [];
for j = 1 : n_eco_cat
    for i = 1 : n_age_cat
        tmp = [Sv(i,j) Ev(i,j) Iv(i,j) Rv(i,j) Dv(i,j)];
        yv0 = [yv0 tmp];
    end
end
yv0 = yv0.';

dt = 0.1;               % time step of integration
T = [0:dt:n_days];      % time vector

% initializing vector for city
Y = zeros(length(T),5*n_age_cat*n_eco_cat);
Y(1,:) = y0.';

% initializing vector for village
Yv = zeros(length(T),5*n_age_cat*n_eco_cat);
Yv(1,:) = yv0.';

% tracking fluxes of S -> D and S -> E
city_flux_SD = zeros(length(T),3,3);
city_flux_SE = zeros(length(T),3,3);
vill_flux_SD = zeros(length(T),3,3);
vill_flux_SE = zeros(length(T),3,3);

for t_ind = 1 : length(T)-1
    
    if T(t_ind) < T_lockdown
        zeta = zeta_old;
    else 
        zeta = zeta_new;
    end
    
    if ( T(t_ind) >= T_lockdown && T(t_ind) < T_lockdown+tp_migration+1 && zeta < zeta_old )
        sp = reshape(Y(t_ind,:),5,n_age_cat,n_eco_cat);
        S(:,:) = sp(1,:,:);
        E(:,:) = sp(2,:,:);
        I(:,:) = sp(3,:,:);
        R(:,:) = sp(4,:,:);
        D(:,:) = sp(5,:,:);
    
        % distributing migrating people of people into categories: 
        % 1. only the mobile poor  from the cities migrate to villages 
        % 2. each category 
        total_mobile = S(1:2,2) + E(1:2,2) + I(1:2,2) + R(1:2,2);
        imports = zeros(n_age_cat,n_eco_cat,4);
        for e_ind = 1 : 2
            if total_mobile(e_ind) ~= 0
                imports(e_ind,2,1) = import_rate * S(e_ind,2) ./ total_mobile(e_ind) .* age_dist(e_ind)/sum(age_dist(1:2));
                imports(e_ind,2,2) = import_rate * E(e_ind,2) ./ total_mobile(e_ind) .* age_dist(e_ind)/sum(age_dist(1:2));
                imports(e_ind,2,3) = import_rate * I(e_ind,2) ./ total_mobile(e_ind) .* age_dist(e_ind)/sum(age_dist(1:2));
                imports(e_ind,2,4) = import_rate * R(e_ind,2) ./ total_mobile(e_ind) .* age_dist(e_ind)/sum(age_dist(1:2));
            end
        end

    else
        imports = zeros(n_age_cat,n_eco_cat,4);
    end
    
    [dYdt,city_flux_SD(t_ind,:,:),city_flux_SE(t_ind,:,:)] = SEIRD_solver(T(t_ind),Y(t_ind,:),contact_map,zeta,n,-imports);
    Y(t_ind+1,:) = Y(t_ind,:) + dYdt*dt;
    
    [dYvdt,vill_flux_SD(t_ind,:,:),vill_flux_SE(t_ind,:,:)] = SEIRD_solver(T(t_ind),Yv(t_ind,:),contact_map,zeta,n,imports/num_villages);
    Yv(t_ind+1,:) = Yv(t_ind,:) + dYvdt*dt;
    
end

%% extracting information: City

S = zeros(length(T),n_age_cat,n_eco_cat);
E = zeros(length(T),n_age_cat,n_eco_cat);
I = zeros(length(T),n_age_cat,n_eco_cat);
R = zeros(length(T),n_age_cat,n_eco_cat);
D = zeros(length(T),n_age_cat,n_eco_cat);

Z = zeros(1,5*n_age_cat*n_eco_cat);
for i = 1 : length(T)
    
    Z(:) = Y(i,:);
    sp = reshape(Z,5,n_age_cat,n_eco_cat);

    S(i,:,:) = sp(1,:,:);
    E(i,:,:) = sp(2,:,:);
    I(i,:,:) = sp(3,:,:);
    R(i,:,:) = sp(4,:,:);
    D(i,:,:) = sp(5,:,:);
    
end

%% summing up

totalS_a = sum(S,3);
totalS_e = sum(S,2);
totalS = sum(totalS_a,2);

totalE_a = sum(E,3);
totalE_e = sum(E,2);
totalE = sum(totalE_a,2);

totalI_a = sum(I,3);
totalI_e = sum(I,2);
totalI = sum(totalI_a,2);

totalR_a = sum(R,3);
totalR_e = sum(R,2);
totalR = sum(totalR_a,2);

totalD_a = sum(D,3);
totalD_e = sum(D,2);
totalD = sum(totalD_a,2);

disp('city scene')
text = sprintf('max infectious: %.f', max(totalI));
disp(text)
text = sprintf('susceptible left %.f', totalS(end));
disp(text)
text = sprintf('exposed left %.f', totalE(end));
disp(text)
text = sprintf('infectious left %.f', totalI(end));
disp(text)
text = sprintf('total receoverd %.f', totalR(end));
disp(text)
text = sprintf('cumulative deaths %.f', totalD(end));
disp(text)

finalN = totalS + totalE + totalI + totalR + totalD;

%% plot each of the 9 (a,e) categories separately
% total over a and e
% figure, plot(T,totalS, T,totalE, T,totalI, T,totalR, T,totalD,'LineWidth',3)
% figure, plot(T,totalE+totalI, T,totalR, T,totalD,'LineWidth',3)

figure, plot(T,totalI, 'LineWidth',3)
hold on
plot(tdata,data,'o')
xlabel('time (days)')
ylabel('number')
% legend('S','E','I','R','D')
% legend('E+I','R','D')
legend('I')
title('City - total')
% saveas(gcf,'seird-city.png')

%% individual (a,e)
figure
count = 0;
for ind_a = 1 : 3
    for ind_e = 1 : 3
        count = count + 1;
        subplot(3,3,count)
        plot(T,S(:,ind_a,ind_e),T,E(:,ind_a,ind_e),T,I(:,ind_a,ind_e),T,R(:,ind_a,ind_e),T,D(:,ind_a,ind_e))
        text = catlabel(ind_a,ind_e);
        title(text)
    end
end
suptitle('City')

% total over a
figure
count = 0;
for ind_e = 1 : 3
    count = count + 1;
    subplot(3,1,count)
    plot(T,totalS_e(:,1,ind_e),T,totalE_e(:,1,ind_e),T,totalI_e(:,1,ind_e),T,totalR_e(:,1,ind_e),T,totalD_e(:,1,ind_e),'Linewidth',2)
    text = catlabel(0,ind_e);
    title(text)
end
suptitle('City - summed over age categories')
% saveas(gcf,'seird-city.png')

% total over e
figure
count = 0;
for ind_a = 1 : 3
    count = count + 1;
    subplot(3,1,count)
    plot(T,totalS_a(:,ind_a,1),T,totalE_a(:,ind_a,1),T,totalI_a(:,ind_a,1),T,totalR_a(:,ind_a,1),T,totalD_a(:,ind_a,1),'Linewidth',2)
    text = catlabel(ind_a,0);
    title(text)
end
suptitle('City - summed over econominc categories')
% saveas(gcf,'seird-city.png')

%% extracting information: Village

Z = zeros(1,5*n_age_cat*n_eco_cat);
for i = 1 : length(T)
    
    Z(:) = Yv(i,:);
    sp = reshape(Z,5,n_age_cat,n_eco_cat);

    S(i,:,:) = sp(1,:,:);
    E(i,:,:) = sp(2,:,:);
    I(i,:,:) = sp(3,:,:);
    R(i,:,:) = sp(4,:,:);
    D(i,:,:) = sp(5,:,:);
    
end

%% summing up

totalS_a = sum(S,3);
totalS_e = sum(S,2);
totalS = sum(totalS_a,2);

totalE_a = sum(E,3);
totalE_e = sum(E,2);
totalE = sum(totalE_a,2);

totalI_a = sum(I,3);
totalI_e = sum(I,2);
totalI = sum(totalI_a,2);

totalR_a = sum(R,3);
totalR_e = sum(R,2);
totalR = sum(totalR_a,2);

totalD_a = sum(D,3);
totalD_e = sum(D,2);
totalD = sum(totalD_a,2);

disp('city scene')
text = sprintf('max infectious: %.f', max(totalI));
disp(text)
text = sprintf('susceptible left %.f', totalS(end));
disp(text)
text = sprintf('exposed left %.f', totalE(end));
disp(text)
text = sprintf('infectious left %.f', totalI(end));
disp(text)
text = sprintf('total receoverd %.f', totalR(end));
disp(text)
text = sprintf('cumulative deaths %.f', totalD(end));
disp(text)

finalN = totalS + totalE + totalI + totalR + totalD;

%% plot each of the 9 (a,e) categories separately

% total over a and e
figure, plot(T,totalS, T,totalE, T,totalI, T,totalR, T,totalD,'LineWidth',3)
xlabel('time (days)')
ylabel('number')
legend('S','E','I','R','D')
title('Village - total')
% saveas(gcf,'seird-city.png')

% individual (a,e)
figure
count = 0;
for ind_a = 1 : 3
    for ind_e = 1 : 3
        count = count + 1;
        subplot(3,3,count)
        plot(T,S(:,ind_a,ind_e),T,E(:,ind_a,ind_e),T,I(:,ind_a,ind_e),T,R(:,ind_a,ind_e),T,D(:,ind_a,ind_e),'Linewidth',2)
        text = catlabel(ind_a,ind_e);
        title(text)
    end
end
suptitle('Village')

% total over a
figure
count = 0;
for ind_e = 1 : 3
    count = count + 1;
    subplot(3,1,count)
    plot(T,totalS_e(:,1,ind_e),T,totalE_e(:,1,ind_e),T,totalI_e(:,1,ind_e),T,totalR_e(:,1,ind_e),T,totalD_e(:,1,ind_e),'Linewidth',2)
    text = catlabel(0,ind_e);
    title(text)
end
suptitle('Village - summed over age categories')
% saveas(gcf,'seird-city.png')

% total over e
figure
count = 0;
for ind_a = 1 : 3
    count = count + 1;
    subplot(3,1,count)
    plot(T,totalS_a(:,ind_a,1),T,totalE_a(:,ind_a,1),T,totalI_a(:,ind_a,1),T,totalR_a(:,ind_a,1),T,totalD_a(:,ind_a,1),'Linewidth',2)
    text = catlabel(ind_a,0);
    title(text)
end
suptitle('Village - summed over econominc categories')
% saveas(gcf,'seird-city.png')