%% Clear all variables and close all figures

clear all
close all
clc

%% Define Tasks 

% compare to real data or not ?
% 0 for NO
% 1 for YES
extract_data = 0;

% consider villages around?
% 0 for NO
% 1 for YES
toggle_village = 1;

% include economic module or not?
% 0 for NO
% 1 for YES
economic_module = 1;

% plots for individual cases while scanning n and zeta?
case_plots = 1;

% set range of n and zeta to be covered
n_range = linspace(0,5,6);
zeta_range = flip(linspace(0.75,1,11));

% extrct final state defined as rate of dI/dt < 0.001;
% extract time of peak infections
% extract final time 

%% extract data from JHU CSSE

if extract_data == 1
   
    data = load('IndiaData.csv');
    data(1:42) = [];
    length(data)
    tdata = [0:1:length(data)-1];

    recov_data = load('IndiaRecovered.csv');
    recov_data(1:42) = [];

    death_data = load('IndiaDeaths.csv');
    death_data(1:42) = [];
    
end

%% Simulation parameters

n_days = 2*365;         % observation period in number of days

N = 12*10^6;            % city population
InitNumCases = 5;       % initial number of cases in the city

Nv = 1500;              % average village population
num_villages = 5000;    % number of villages: there are about 25k villages in Karnataka and 
                        % let's say Bangalore takes in from a fifth of
                        % these villages, therefore 5000 villages

T_lockdown = 21;        % time point of lockdown initiation
duration_lockdown = 30;
zeta_old = 1;           % without mitigation


% total people migrated from city to villages in response to migration 
% = import_rate * tp_migration
tp_migration = 60;                      % time duration of migration period sfrom city to villages 
import_rate = 0.2*N/tp_migration;       % rate of migration from the city after lockdown (#people per day)

if extract_data == 1
    N = 1.4*10^9;               % city population
    InitNumCases = data(1);     % initial number of cases in the city
    import_rate = 0;
end

if toggle_village == 0 
    num_villages = 0;
    Nv = 0;
    import_rate = 0;
end

%% Extracting contact map from Keisha et al

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
contact_map(2,1,2,3) = 0.5*age_CM(2,2);      % young,rich - young,poor
contact_map(:,2,:,1) = age_CM;
contact_map(:,2,:,2) = age_CM;
contact_map(2,2,2,3) = 0.5*age_CM(2,2);
contact_map(2,3,2,1) = 0.5*age_CM(2,2);
contact_map(2,3,2,2) = 0.5*age_CM(2,2);
contact_map(:,3,:,3) = age_CM;

age_CM_rural = zeros(3,3);
work_CM = load('India-Work-CM.csv');
base_CM_rural = base_CM - work_CM;
base_CM_rural = base_CM_rural .* repmat(dummy_age_dist.',1,16);
age_CM_rural(1,1) = sum(sum(base_CM_rural(1:3,1:3)));
age_CM_rural(1,2) = sum(sum(base_CM_rural(1:3,4:12)));
age_CM_rural(1,3) = sum(sum(base_CM_rural(1:3,13:16)));
age_CM_rural(2,1) = sum(sum(base_CM_rural(4:12,1:3)));
age_CM_rural(2,2) = sum(sum(base_CM_rural(4:12,4:12)));
age_CM_rural(2,3) = sum(sum(base_CM_rural(4:12,13:16)));
age_CM_rural(3,1) = sum(sum(base_CM_rural(13:16,1:3)));
age_CM_rural(3,2) = sum(sum(base_CM_rural(13:16,4:12)));
age_CM_rural(3,3) = sum(sum(base_CM_rural(13:16,13:16)));

% definition used by Neher Lab
%     contact_map = R0 * zeta * ones(n_age_cat,n_age_cat,n_eco_cat,n_eco_cat) / ti;

% definition used by Ronojoy Adhikari
% 1. within an economic group, contact map is as given by age_CM
% 2. rich-poor interactions are limited to young
contact_map_rural = zeros(3,3,3,3);
contact_map_rural(:,1,:,1) = age_CM_rural;
contact_map_rural(:,1,:,2) = age_CM_rural;               % poor - poor
contact_map_rural(2,1,2,3) = 0.5*age_CM_rural(2,2);      % young,rich - young,poor
contact_map_rural(:,2,:,1) = age_CM_rural;
contact_map_rural(:,2,:,2) = age_CM_rural;
contact_map_rural(2,2,2,3) = 0.5*age_CM_rural(2,2);
contact_map_rural(2,3,2,1) = 0.5*age_CM_rural(2,2);
contact_map_rural(2,3,2,2) = 0.5*age_CM_rural(2,2);
contact_map_rural(:,3,:,3) = age_CM_rural;

% age distribution is the same in the city and villages (taken from cumulative age dist of India)
age_dist = [0.26 0.64 0.1];
age_dist = age_dist.';

n_age_cat = length(age_dist);   % 0-14, 15-59,  60+ 
n_eco_cat = 3;                  % immobile poor, mobile poor, rich


%% Epidemiological parameters

% data from Neherlab.org/covid19-scenarios
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
age_dist_9 = [17 18.3 17.4 15.6 12.3 9.3 6.3 2.8 0]/100;
age_dist_9(end) = 1 - sum(age_dist_9);
age_dist_9 = age_dist_9.';
f(1) = (f(1)*age_dist_9(1) + f(2)*age_dist_9(2))/(age_dist_9(1) + age_dist_9(2)); 
f(2) = sum(f(2:6).*age_dist_9(2:6))/sum(age_dist_9(2:6));
f(3) = sum(f(7:9).*age_dist_9(7:9))/sum(age_dist_9(7:9));
f(4:end) = [];
f = repmat(f,1,n_eco_cat);

% adjustment to fit India data
f = 4*f;

%%

finalS = zeros(length(n_range),length(zeta_range));
finalE = zeros(length(n_range),length(zeta_range));
finalI = zeros(length(n_range),length(zeta_range));
finalR = zeros(length(n_range),length(zeta_range));
finalD = zeros(length(n_range),length(zeta_range));

deaths_SD = zeros(length(n_range),length(zeta_range));
deaths_SE = zeros(length(n_range),length(zeta_range));

city_maxI = zeros(length(n_range),length(zeta_range));
city_tmaxI = zeros(length(n_range),length(zeta_range));
vill_maxI = zeros(length(n_range),length(zeta_range));
vill_tmaxI = zeros(length(n_range),length(zeta_range));

for n_ind = 1 : length(n_range)
    
folder = [pwd, '\', 'n=',num2str(n_range(n_ind))];
    
for zeta_ind = 1 : length(zeta_range)

subfolder = [folder,'\','zeta=',num2str(zeta_range(zeta_ind))];

if ~exist(subfolder, 'dir')
   mkdir(subfolder)
end
    
zeta_new = zeta_range(zeta_ind);           % by how much is contact or R0 reduced with mitigation

% new rate_SD = rate_SD / (zeta)^n
n = n_range(n_ind);                         % effect if mitigation on death rates due to economic reasons

text = sprintf('n = %f, zeta = %f', n, zeta_new);
disp(text)
    
% Initializing City and Village populations

dt = 0.1;               % time step of integration
T = [0:dt:n_days];       % time vector

% declare S, E, I, R, D variables in the city and villages
% City:
S = zeros(n_age_cat, n_eco_cat);
E = zeros(n_age_cat, n_eco_cat);
I = zeros(n_age_cat, n_eco_cat);
R = zeros(n_age_cat, n_eco_cat);
D = zeros(n_age_cat, n_eco_cat);

% Initialization: 
% 1. the infections are only in the city 
% 2. only the young (15-59) rich were infected 
% 3. S in villages consisted of 80% immobile poor, 20% rich
% 4. S in city consisted of 35% immobile poor, 35% mobile poor, 30% rich

I(:,1) = 0 * transpose([0 InitNumCases 0]);
I(:,2) = 0 * transpose([0 InitNumCases 0]);
I(:,3) = 1 * transpose([0 InitNumCases 0]);

E(:,3) = 0.4 * transpose([0 InitNumCases 0]);

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

% initializing vector for city
Y = zeros(length(T),5*n_age_cat*n_eco_cat);
Y(1,:) = y0.';

% tracking fluxes of S -> D and S -> E
city_flux_SD = zeros(length(T),3,3);
city_flux_ID = zeros(length(T),3,3);

if toggle_village == 1

    % Village:
    Sv = zeros(n_age_cat, n_eco_cat);
    Ev = zeros(n_age_cat, n_eco_cat);
    Iv = zeros(n_age_cat, n_eco_cat);
    Rv = zeros(n_age_cat, n_eco_cat);
    Dv = zeros(n_age_cat, n_eco_cat);

    Sv(:,1) = 0.75 * age_dist * Nv;
    Sv(:,2) = 0.05 * age_dist * Nv;
    Sv(:,3) = 0.2 * age_dist * Nv;

    yv0 = [];
    for j = 1 : n_eco_cat
        for i = 1 : n_age_cat
            tmp = [Sv(i,j) Ev(i,j) Iv(i,j) Rv(i,j) Dv(i,j)];
            yv0 = [yv0 tmp];
        end
    end
    yv0 = yv0.';

    % initializing vector for village
    Yv = zeros(length(T),5*n_age_cat*n_eco_cat);
    Yv(1,:) = yv0.';

    vill_flux_SD = zeros(length(T),3,3);
    vill_flux_ID = zeros(length(T),3,3);

end

% Run the simulation

for t_ind = 1 : length(T)-1
    
    if T(t_ind) < T_lockdown
        zeta = zeta_old;
    else 
        zeta = zeta_new;
    end
    
    % if you want to implement mitigation as in R Adhikari's work
    % u = -tanh((T(t_ind)-T_lockdown)/1) + tanh((T(t_ind)-51)/1);
    
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
    
    [dYdt,city_flux_SD(t_ind,:,:),city_flux_ID(t_ind,:,:)] = SEIRD_solver(economic_module,dt,T(t_ind),Y(t_ind,:),contact_map,zeta,n,f,-imports);
    Y(t_ind+1,:) = Y(t_ind,:) + dYdt*dt;
    
    if toggle_village == 1
    [dYvdt,vill_flux_SD(t_ind,:,:),vill_flux_ID(t_ind,:,:)] = SEIRD_solver(economic_module,dt,T(t_ind),Yv(t_ind,:),contact_map_rural,zeta,n,f,imports/num_villages);
    Yv(t_ind+1,:) = Yv(t_ind,:) + dYvdt*dt;
    end
    
%     if sum(sum(city_flux_ID(t_ind,:,:))) < 0.01 && sum(sum(vill_flux_ID(t_ind,:,:))) < 0.01
%         break
%     end
    
end

% extracting information: City
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

% summing up

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

% disp('city scene')
% text = sprintf('max infectious: %.f', max(totalI));
% disp(text)
% text = sprintf('susceptible left %.f', totalS(end));
% disp(text)
% text = sprintf('exposed left %.f', totalE(end));
% disp(text)
% text = sprintf('infectious left %.f', totalI(end));
% disp(text)
% text = sprintf('total receoverd %.f', totalR(end));
% disp(text)
% text = sprintf('cumulative deaths %.f', totalD(end));
% disp(text)

% test with real data to adjust parameters

if extract_data ==1 && economic_module == 0

    figure, plot(T,totalI,'b', T,totalR,'r', T,totalD,'g', T,totalE,'k', 'LineWidth',3)
    hold on
    plot(tdata,data,'ob', tdata,recov_data,'or', tdata,death_data,'og')
    ylim = get(gca,'YLim');
    rectangle('Position',[0 ylim(1) T_lockdown ylim(2)-ylim(1)],'EdgeColor',[0.5 0.5 0.5 0.2],'FaceColor',[0.5 0.5 0.5 0.2])
    xlabel('time (days)')
    ylabel('number')
    legend('I','R','D','E')
    title('Number of cases in the country (Mar 4 - Apr 15, 2020)')

end

% individual case plots 

if case_plots == 1 
    
    % fluxes 
    figure, plot(T,sum(sum(city_flux_SD(:,:,:),3),2), T,sum(sum(city_flux_ID(:,:,:),3),2), 'LineWidth',3)
    xlabel('flux_{SD} (per day)')
    ylabel('flux_{SE} (per day)')
    legend('flux_SD','flux_ID')
    figname = 'City_fluxes.png';
    fullFileName = fullfile(subfolder, figname);
    saveas(gcf, fullFileName, 'png')
    close
    
    % plot each of the 9 (a,e) categories separately

    % total over a and e
    figure, plot(T,totalS, T,totalE, T,totalI, T,totalR, T,totalD,'LineWidth',3)
    xlabel('time (days)')
    ylabel('number')
    legend('S','E','I','R','D')
    title('City - total')
    figname = 'City_total.png';
    fullFileName = fullfile(subfolder, figname);
    saveas(gcf, fullFileName, 'png')
    close

    % individual (a,e)
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
    figname = 'City.png';
    fullFileName = fullfile(subfolder, figname);
    saveas(gcf, fullFileName, 'png')
    close

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
    figname = 'City_econ.png';
    fullFileName = fullfile(subfolder, figname);
    saveas(gcf, fullFileName, 'png')
    close

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
    figname = 'City_age.png';
    fullFileName = fullfile(subfolder, figname);
    saveas(gcf, fullFileName, 'png')
    close

end 

% extracting information: Village

if toggle_village == 1
    
    Sv = zeros(length(T),n_age_cat,n_eco_cat);
    Ev = zeros(length(T),n_age_cat,n_eco_cat);
    Iv = zeros(length(T),n_age_cat,n_eco_cat);
    Rv = zeros(length(T),n_age_cat,n_eco_cat);
    Dv = zeros(length(T),n_age_cat,n_eco_cat);

    Zv = zeros(1,5*n_age_cat*n_eco_cat);
    for i = 1 : length(T)

        Z(:) = Yv(i,:);
        sp = reshape(Z,5,n_age_cat,n_eco_cat);

        Sv(i,:,:) = sp(1,:,:);
        Ev(i,:,:) = sp(2,:,:);
        Iv(i,:,:) = sp(3,:,:);
        Rv(i,:,:) = sp(4,:,:);
        Dv(i,:,:) = sp(5,:,:);

    end

    % summing up

    totalSv_a = sum(Sv,3);
    totalSv_e = sum(Sv,2);
    totalSv = sum(totalSv_a,2);

    totalEv_a = sum(Ev,3);
    totalEv_e = sum(Ev,2);
    totalEv = sum(totalEv_a,2);

    totalIv_a = sum(Iv,3);
    totalIv_e = sum(Iv,2);
    totalIv = sum(totalIv_a,2);

    totalRv_a = sum(Rv,3);
    totalRv_e = sum(Rv,2);
    totalRv = sum(totalRv_a,2);

    totalDv_a = sum(Dv,3);
    totalDv_e = sum(Dv,2);
    totalDv = sum(totalDv_a,2);

%     disp('city scene')
%     text = sprintf('max infectious: %.f', max(totalI));
%     disp(text)
%     text = sprintf('susceptible left %.f', totalS(end));
%     disp(text)
%     text = sprintf('exposed left %.f', totalE(end));
%     disp(text)
%     text = sprintf('infectious left %.f', totalI(end));
%     disp(text)
%     text = sprintf('total receoverd %.f', totalR(end));
%     disp(text)
%     text = sprintf('cumulative deaths %.f', totalD(end));
%     disp(text)
% 
%     finalN = totalS + totalE + totalI + totalR + totalD;

    % plot each of the 9 (a,e) categories separately

    if case_plots == 1

        % total over a and e
        figure, plot(T,totalSv, T,totalEv, T,totalIv, T,totalRv, T,totalDv,'LineWidth',3)
        xlabel('time (days)')
        ylabel('number')
        legend('S','E','I','R','D')
        title('Village - total')
        figname = 'Village_total.png';
        fullFileName = fullfile(subfolder, figname);
        saveas(gcf, fullFileName, 'png')
        close

        % individual (a,e)
        figure
        count = 0;
        for ind_a = 1 : 3
            for ind_e = 1 : 3
                count = count + 1;
                subplot(3,3,count)
                plot(T,Sv(:,ind_a,ind_e),T,Ev(:,ind_a,ind_e),T,Iv(:,ind_a,ind_e),T,Rv(:,ind_a,ind_e),T,Dv(:,ind_a,ind_e),'Linewidth',2)
                text = catlabel(ind_a,ind_e);
                title(text)
            end
        end
        suptitle('Village')
        figname = 'Village.png';
        fullFileName = fullfile(subfolder, figname);
        saveas(gcf, fullFileName, 'png')
        close

        % total over a
        figure
        count = 0;
        for ind_e = 1 : 3
            count = count + 1;
            subplot(3,1,count)
            plot(T,totalSv_e(:,1,ind_e),T,totalEv_e(:,1,ind_e),T,totalIv_e(:,1,ind_e),T,totalRv_e(:,1,ind_e),T,totalDv_e(:,1,ind_e),'Linewidth',2)
            text = catlabel(0,ind_e);
            title(text)
        end
        suptitle('Village - summed over age categories')
        figname = 'Village_econ.png';
        fullFileName = fullfile(subfolder, figname);
        saveas(gcf, fullFileName, 'png')
        close

        % total over e
        figure
        count = 0;
        for ind_a = 1 : 3
            count = count + 1;
            subplot(3,1,count)
            plot(T,totalSv_a(:,ind_a,1),T,totalEv_a(:,ind_a,1),T,totalIv_a(:,ind_a,1),T,totalRv_a(:,ind_a,1),T,totalDv_a(:,ind_a,1),'Linewidth',2)
            text = catlabel(ind_a,0);
            title(text)
        end
        suptitle('Village - summed over econominc categories')
        figname = 'Village_age.png';
        fullFileName = fullfile(subfolder, figname);
        saveas(gcf, fullFileName, 'png')
        close

    end
    
end
% plot final state as a function of mitigation strategy

finalNv = totalSv + totalEv + totalIv + totalRv + totalDv;

finalS(n_ind, zeta_ind) = finalS(n_ind, zeta_ind) + num_villages * totalSv(end);
finalE(n_ind, zeta_ind) = finalE(n_ind, zeta_ind) + num_villages * totalEv(end);
finalI(n_ind, zeta_ind) = finalI(n_ind, zeta_ind) + num_villages * totalIv(end);
finalR(n_ind, zeta_ind) = finalR(n_ind, zeta_ind) + num_villages * totalRv(end); 
finalD(n_ind, zeta_ind) = finalD(n_ind, zeta_ind) + num_villages * totalDv(end);

deaths_SD(n_ind,zeta_ind) = sum(sum(sum(city_flux_SD*dt))) + sum(sum(sum(vill_flux_SD*dt)));
deaths_SE(n_ind,zeta_ind) = sum(sum(sum(city_flux_ID*dt))) + sum(sum(sum(vill_flux_ID*dt)));

[val, pos] = max(totalI);

city_maxI(n_ind,zeta_ind) = val;
city_tmaxI(n_ind,zeta_ind) = T(pos);

[val, pos] = max(totalIv);

vill_maxI(n_ind,zeta_ind) = val;
vill_tmaxI(n_ind,zeta_ind) = T(pos);

end

end

%% final state results

f1 = figure();
f2 = figure();
f3 = figure();
f4 = figure();
f5 = figure();
f6 = figure();
f7 = figure();

for n_ind = 1 : length(n_range)

    txt = ['n = ',num2str(n_range(n_ind))];
    
    figure(f1)
    plot(zeta_range, deaths_SD(n_ind,:),'DisplayName',txt)
    xlabel('\zeta')
    ylabel('deaths_{SD}')
    hold on

    figure(f2)
    plot(zeta_range, deaths_SE(n_ind,:),'DisplayName',txt)
    xlabel('\zeta')
    ylabel('deaths_{SE}')
    hold on

    figure(f3)
    plot(zeta_range, finalD(n_ind,:),'DisplayName',txt)
    xlabel('\zeta')
    ylabel('total deaths')
    hold on
    
    figure(f4)
    plot(zeta_range, city_maxI(n_ind,:),'DisplayName',txt)
    xlabel('\zeta')
    ylabel('max infectons (city)')
    hold on

    figure(f5)
    plot(zeta_range, city_tmaxI(n_ind,:),'DisplayName',txt)
    xlabel('\zeta')
    ylabel('time of peak infections (city)')
    hold on
    
    figure(f6)
    plot(zeta_range, vill_maxI(n_ind,:),'DisplayName',txt)
    xlabel('\zeta')
    ylabel('max infections (village)')
    hold on
    
    figure(f7)
    plot(zeta_range, vill_tmaxI(n_ind,:),'DisplayName',txt)
    xlabel('\zeta')
    ylabel('time of peak infections (village)')
    hold on
    
end 

figure(f1)
hold off 
legend show

figure(f2)
hold off 
legend show

figure(f3)
hold off 
legend show

figure(f4)
hold off 
legend show

figure(f5)
hold off 
legend show

figure(f6)
hold off 
legend show

figure(f7)
hold off 
legend show
%% Plot the age_CM

figure()
flipped_age_CM = flip(age_CM,1);
imagesc(flipped_age_CM)
yticks([1 2 3])
yticklabels({'3','2','1'})
xticks([1 2 3])
title('Contact map C^a_{a^\prime} within the same economic category e^\prime = e')
xlabel('age category, a')
ylabel('age category, a^\prime')