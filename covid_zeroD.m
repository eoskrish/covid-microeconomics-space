% declare S, E, I, R, H, C, D variables
N = 66165886;

n_age_cat = 9;
n_eco_cat = 1;

S = zeros(n_age_cat, n_eco_cat);
E = zeros(n_age_cat, n_eco_cat);
I = zeros(n_age_cat, n_eco_cat);
R = zeros(n_age_cat, n_eco_cat);
H = zeros(n_age_cat, n_eco_cat);
C = zeros(n_age_cat, n_eco_cat);
D = zeros(n_age_cat, n_eco_cat);

I = [0 0 0 0 3 0 0 0 0];
E = [0 0 0 0 7 0 0 0 0];
S = [17 18.3 17.4 15.6 12.3 9.3 6.3 2.8 0];
S(end) = 100 - sum(S);
S = S /100 * (N-sum(E)-sum(I));

S = S.';
I = I.';
E = E.';

y0 = [];
for j = 1 : n_eco_cat
    for i = 1 : n_age_cat
        tmp = [S(i,j) E(i,j) I(i,j) R(i,j) H(i,j) C(i,j) D(i,j)];
        y0 = [y0 tmp];
    end
end
y0 = y0.';  

[T,Y] = ode45(@(t,y)ODEsolver(t,y,N,n_age_cat,n_eco_cat),[0:182],y0);

S = zeros(length(T),n_age_cat,n_eco_cat);
E = zeros(length(T),n_age_cat,n_eco_cat);
I = zeros(length(T),n_age_cat,n_eco_cat);
R = zeros(length(T),n_age_cat,n_eco_cat);
C = zeros(length(T),n_age_cat,n_eco_cat);
H = zeros(length(T),n_age_cat,n_eco_cat);
D = zeros(length(T),n_age_cat,n_eco_cat);

Z = zeros(1,7*n_age_cat*n_eco_cat);
for i = 1 : length(T)
    
    Z(:) = Y(i,:);
    sp = reshape(Z,7,n_age_cat,n_eco_cat);

    S(i,:,:) = sp(1,:,:);
    E(i,:,:) = sp(2,:,:);
    I(i,:,:) = sp(3,:,:);
    R(i,:,:) = sp(4,:,:);
    H(i,:,:) = sp(5,:,:);
    C(i,:,:) = sp(6,:,:);
    D(i,:,:) = sp(7,:,:);
    
end

totalS = 0;
totalE = 0;
totalI = 0;
totalR = 0;
totalH = 0;
totalC = 0;
totalD = 0;
for i = 1 : 9
    totalS = totalS + Y(:,(i-1)*7+1);
    totalE = totalE + Y(:,(i-1)*7+2);
    totalI = totalI + Y(:,(i-1)*7+3);
    totalR = totalR + Y(:,(i-1)*7+4);
    totalH = totalH + Y(:,(i-1)*7+5);
    totalC = totalC + Y(:,(i-1)*7+6);
    totalD = totalD + Y(:,i*7);
end
loglog(T,totalS, T,totalE, T,totalI, T,totalR, T,totalH, T,totalC, T,totalD)