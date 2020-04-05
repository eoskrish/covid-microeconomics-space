clear all
close all
% clc

n_days = 12*30;
import_rate = 100;

N = 12*10^6;
numCases = 10;

Nv = 5000;
num_villages = 1;

age_dist = [0.26 0.64 0.1];
age_dist = age_dist.';

n_age_cat = length(age_dist);   % 0-14, 15-59,  60+ 
n_eco_cat = 3;                  % immobile poor, mobile poor, rich

% declare S, I, R, D variables
% Pop = zeros(n_age_cat,n_eco_cat,4);
S = zeros(n_age_cat, n_eco_cat);
E = zeros(n_age_cat, n_eco_cat);
I = zeros(n_age_cat, n_eco_cat);
R = zeros(n_age_cat, n_eco_cat);
D = zeros(n_age_cat, n_eco_cat);

Sv = zeros(n_age_cat, n_eco_cat);
Ev = zeros(n_age_cat, n_eco_cat);
Iv = zeros(n_age_cat, n_eco_cat);
Rv = zeros(n_age_cat, n_eco_cat);
Dv = zeros(n_age_cat, n_eco_cat);

% initialization
% Pop(:,1,2) = 0 * transpose([0 numCases 0]);
% Pop(:,2,2) = 0 * transpose([0 numCases 0]);
% Pop(:,3,2) = 1 * transpose([0 numCases 0]);
I(:,1) = 0 * transpose([0 numCases 0]);
I(:,2) = 0 * transpose([0 numCases 0]);
I(:,3) = 1 * transpose([0 numCases 0]);

Sv(:,1) = 0.8 * age_dist * Nv;
Sv(:,2) = 0 * age_dist * Nv;
Sv(:,3) = 0.2 * age_dist * Nv;

S(:,1) = 0.35 * age_dist * (N-numCases);
S(:,2) = 0.35 * age_dist * (N-numCases);
S(:,3) = 0.30 * age_dist * (N-numCases);

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


dt = 0.1;
T = [0:dt:n_days];

Y = zeros(length(T),5*n_age_cat*n_eco_cat);
Y(1,:) = y0.';

Yv = zeros(length(T),5*n_age_cat*n_eco_cat);
Yv(1,:) = yv0.';


for t_ind = 1 : length(T)-1
    
    if T(t_ind) < 14
    sp = reshape(Y(t_ind,:),5,n_age_cat,n_eco_cat);
    S(:,:) = sp(1,:,:);
    E(:,:) = sp(2,:,:);
    I(:,:) = sp(3,:,:);
    R(:,:) = sp(4,:,:);
    D(:,:) = sp(5,:,:);
    
    % distributing movement of people into categories
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
    
    dYdt = SEIRD_solver(T(t_ind),Y(t_ind,:),n_age_cat,n_eco_cat,-imports);
    Y(t_ind+1,:) = Y(t_ind,:) + dYdt*dt;
    
    dYvdt = SEIRD_solver(T(t_ind),Yv(t_ind,:),n_age_cat,n_eco_cat,imports/num_villages);
    Yv(t_ind+1,:) = Yv(t_ind,:) + dYvdt*dt;
    
end
    
% options1 = odeset('MaxStep',0.1);
% options2 = odeset(options1,'NonNegative',1);
% [T,Y] = ode45(@(t,y)SIR_solver(t,y,n_age_cat,n_eco_cat,imports),[0:n_days],y0,options1);

S = zeros(length(T),n_age_cat,n_eco_cat);
E = zeros(length(T),n_age_cat,n_eco_cat);
S = zeros(length(T),n_age_cat,n_eco_cat);
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

totalS = sum(S,3);
totalS = sum(totalS,2);
totalE = sum(E,3);
totalE = sum(totalE,2);
totalI = sum(I,3);
totalI = sum(totalI,2);
totalR = sum(R,3);
totalR = sum(totalR,2);
totalD = sum(D,3);
totalD = sum(totalD,2);

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

figure, plot(T,totalS, T,totalE, T,totalI, T,totalR, T,totalD,'LineWidth',3)
xlabel('time (days)')
ylabel('number')
legend('S','E','I','R','D')
title('city scene')
% saveas(gcf,'seird-city.png')

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

totalS = sum(S,3);
totalS = sum(totalS,2);
totalE = sum(E,3);
totalE = sum(totalE,2);
totalI = sum(I,3);
totalI = sum(totalI,2);
totalR = sum(R,3);
totalR = sum(totalR,2);
totalD = sum(D,3);
totalD = sum(totalD,2);

disp('village scene')
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

figure, plot(T,totalS, T,totalE, T,totalI, T,totalR, T,totalD,'LineWidth',3)
xlabel('time (days)')
ylabel('number')
legend('S','E','I','R','D')
title('village scene')
% saveas(gcf,'seird-village.png')

% figure, loglog(T,totalS, T,totalI, T,totalR, T,totalD)
% xlabel('time (days)')
% ylabel('number')
% legend('S','I','R','D')
% saveas(gcf,'log-log-plot.png')
