function covid_SIRD_simulator(flag,init_pop,init_cases,import_rate,n_days)

N = init_pop;
numCases = init_cases;

age_dist = [0.26 0.64 0.1];
age_dist = age_dist.';

n_age_cat = length(age_dist);   % 0-14, 15-59,  60+ 
n_eco_cat = 3;                  % immobile poor, mobile poor, rich

% declare S, I, R, D variables
Pop = zeros(n_age_cat,n_eco_cat,4);
% S = zeros(n_age_cat, n_eco_cat);
% I = zeros(n_age_cat, n_eco_cat);
% R = zeros(n_age_cat, n_eco_cat);
% D = zeros(n_age_cat, n_eco_cat);

% initialization
I(:,1) = 0 * transpose([0 numCases 0]);
I(:,2) = 0 * transpose([0 numCases 0]);
I(:,3) = 1 * transpose([0 numCases 0]);

if flag == 0
    % for village
    S(:,1) = 1 * age_dist * (N-numCases);
    S(:,2) = 0 * age_dist * (N-numCases);
    S(:,3) = 0 * age_dist * (N-numCases);
else
    % for city
    S(:,1) = 0.35 * age_dist * (N-numCases);
    S(:,2) = 0.35 * age_dist * (N-numCases);
    S(:,3) = 0.30 * age_dist * (N-numCases);
end

% y0 = [];
% for j = 1 : n_eco_cat
%     for i = 1 : n_age_cat
%         tmp = [S(i,j) I(i,j) R(i,j) D(i,j)];
%         y0 = [y0 tmp];
%     end
% end
% y0 = y0.';

% distributing movement of people into categories
total_mobile = S(1:2,2) + I(1:2,2) + R(1:2,2);
imports = zeros(n_age_cat,n_eco_cat,3);
imports(1:2,2,1) = import_rate * S(1:2,2) ./ total_mobile(1:2) .* age_dist(1:2)/sum(age_dist(1:2));
imports(1:2,2,2) = import_rate * I(1:2,2) ./ total_mobile(1:2) .* age_dist(1:2)/sum(age_dist(1:2));
imports(1:2,2,3) = import_rate * R(1:2,2) ./ total_mobile(1:2) .* age_dist(1:2)/sum(age_dist(1:2));

options1 = odeset('MaxStep',0.1);
% options2 = odeset(options1,'NonNegative',1);
[T,Y] = ode45(@(t,y)SIR_solver(t,y,n_age_cat,n_eco_cat,imports),[0:n_days],y0,options1);

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

text1 = sprintf('max infectious: %.f', max(totalI));
disp(text1)
text2 = sprintf('susceptible left %.f', totalS(end));
disp(text2)
text3 = sprintf('infectious left %.f', totalI(end));
disp(text3)
text4 = sprintf('total receoverd %.f', totalR(end));
disp(text4)
text5 = sprintf('cumulative deaths %.f', totalD(end));
disp(text5)

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

end