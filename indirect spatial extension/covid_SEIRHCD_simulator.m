clear all 
close all

% declare S, E, I, R, H, C, D variables

% N = 66165886;
N = 100000;
numCases = 10;

age_dist = [17 18.3 17.4 15.6 12.3 9.3 6.3 2.8 0]/100;
age_dist(end) = 1 - sum(age_dist);
age_dist = age_dist.';

n_age_cat = length(age_dist); % 0-9, 1-19,..., 80+ 
n_eco_cat = 3; % immobile poor, moile poor, rich

S = zeros(n_age_cat, n_eco_cat);
I = zeros(n_age_cat, n_eco_cat);
R = zeros(n_age_cat, n_eco_cat);
D = zeros(n_age_cat, n_eco_cat);

I(:,1) = 0 * transpose([0 0 0 0 numCases 0 0 0 0]);
I(:,2) = 0 * transpose([0 0 0 0 numCases 0 0 0 0]);
I(:,3) = 1 * transpose([0 0 0 0 numCases 0 0 0 0]);

S(:,1) = 0.3 * age_dist * (N-numCases);
S(:,2) = 0.3 * age_dist * (N-numCases);
S(:,3) = 0.4 * age_dist * (N-numCases);

y0 = [];
for j = 1 : n_eco_cat
    for i = 1 : n_age_cat
        tmp = [S(i,j) I(i,j) R(i,j) D(i,j)];
        y0 = [y0 tmp];
    end
end
y0 = y0.';  

import_rate = 5;

options1 = odeset('MaxStep',0.1);
options2 = odeset(options1,'NonNegative',1);
[T,Y] = ode45(@(t,y)SIR_solver(t,y,n_age_cat,n_eco_cat,-import_rate),[0:50],y0,options1);

S = zeros(length(T),n_age_cat,n_eco_cat);
I = zeros(length(T),n_age_cat,n_eco_cat);
R = zeros(length(T),n_age_cat,n_eco_cat);
D = zeros(length(T),n_age_cat,n_eco_cat);

Z = zeros(1,4*n_age_cat*n_eco_cat);
for i = 1 : length(T)
    
    Z(:) = Y(i,:);
    sp = reshape(Z,4,n_age_cat,n_eco_cat);

    S(i,:,:) = sp(1,:,:);
    I(i,:,:) = sp(2,:,:);
    R(i,:,:) = sp(3,:,:);
    D(i,:,:) = sp(4,:,:);
    
end

totalS = sum(S,3);
totalS = sum(totalS,2);
totalI = sum(I,3);
totalI = sum(totalI,2);
totalR = sum(R,3);
totalR = sum(totalR,2);
totalD = sum(D,3);
totalD = sum(totalD,2);

sprintf('%.f, %.f, %.f, %.f', totalS(end), totalI(end), totalR(end), totalD(end))
finalN = totalS + totalI + totalR + totalD;

figure, plot(T,totalS, T,totalI, T,totalR, T,totalD)
xlabel('time (days)')
ylabel('number')
legend('S','I','R','D')
saveas(gcf,'real-time-plot.png')

% figure, loglog(T,totalS, T,totalI, T,totalR, T,totalD)
% xlabel('time (days)')
% ylabel('number')
% legend('S','I','R','D')
% saveas(gcf,'log-log-plot.png')    