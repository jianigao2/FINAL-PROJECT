% PROGRAM NAME: ps4huggett.m
clear all;
close all;
clc;

% PARAMETERS
beta = .9932; %discount factor 
sigma = 1.5; % coefficient of risk aversion
b = 0.5; % replacement ratio (unemployment benefits)

w1=0.4539; %rural weight
w2=0.5461; %urban weight

y_s = [1, b]; % endowment in employment states for rural
y_s1 = [2, b]; % endowment in employment states for urban

PI = [.97 .03; .5 .5]; % transition matrix

d=0.3;    %Housing consumption percentage
eta=0.25; %transaction cost

% ASSET VECTOR
a_lo = -2; %lower bound of grid points
a_hi = 5;%upper bound of grid points
a_num = 10;
a = linspace(a_lo, a_hi, a_num); % asset (row) vector

% Housing Vector
h_lo=0.1;
h_hi=5;
h_num=10;
h=linspace(h_lo,h_hi,h_num);

% INITIAL GUESS FOR q
q_min = 0.93;
q_max = 1.2;


% ITERATE OVER ASSET PRICES
aggsav = 1 ;
a=repmat(a,[a_num 1 1]);
a=a(:);
while (abs(aggsav) >= 0.01) ;
    q_guess = (q_min + q_max) / 2;
    % CURRENT RETURN (UTILITY) FUNCTION
    con = bsxfun(@minus, a', q_guess * a);
    ho=repmat(h',[1 h_num 1]);
    hous=repmat(ho,h_num);
    con=con-hous-eta.*(abs(hous'-hous));
    cons = bsxfun(@plus, con, permute(y_s, [1 3 2])); %rural consumption
    cons1 = bsxfun(@plus, con, permute(y_s1, [1 3 2])); %urban consumption
    hous=repmat(hous,[1 1 2]);
    ret=((cons.^(1-d)).*(hous.^d)).^(1-sigma)./(1-sigma); % current period utility/rural
    ret1=((cons1.^(1-d)).*(hous.^d)).^(1-sigma)./(1-sigma); % current period utility/urban
    ret(cons<0)=-Inf;
    ret1(cons1<0)=-Inf;
    % INITIAL VALUE FUNCTION GUESS
    v_guess = zeros(2,a_num*h_num);
    v_guess1 = zeros(2,a_num*h_num);
    
    % VALUE FUNCTION ITERATION
    v_tol = 1;
    while v_tol >.000001;
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
        v_mat=ret+beta*repmat(permute((PI*v_guess),[3 2 1]), [a_num 1 1]);
        v_mat1=ret+beta*repmat(permute((PI*v_guess1),[3 2 1]), [a_num 1 1]);
        
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
       [vfn, pol_indx] = max(v_mat, [], 2); %max for each row/rural
       [vfn1, pol_indx1] = max(v_mat1, [], 2); %max for each row/rural
       
       v_tol = max(max(abs(permute(vfn, [3 1 2])-v_guess)));
       v_tol1 = max(max(abs(permute(vfn1, [3 1 2])-v_guess1)));
 
       v_tol=max(v_tol,v_tol1);
       
       v_guess = permute(vfn, [3 1 2]);
       v_guess1 = permute(vfn1, [3 1 2]);
    end;
    
    % KEEP DECSISION RULE
    ab=a';
    pol_indx=permute(pol_indx, [3 1 2]);
    pol_fn = ab(pol_indx);
    
    pol_indx1=permute(pol_indx1, [3 1 2]);
    pol_fn1 = ab(pol_indx1);
    
    % SET UP INITITAL DISTRIBUTION
    Mu=ones (2,a_num*h_num); % Mu is the distribution for a(a,h,y)
    Mu=Mu/sum(Mu(:));
    
    Mu1=ones (2,a_num*h_num); % Mu is the distribution for a(a,h,y)
    Mu1=Mu1/sum(Mu1(:));
    
    mu_tol=1;
  while mu_tol>1e-07 
      % ITERATE OVER DISTRIBUTIONS
      MuNew = zeros(size(Mu));
      MuNew1= zeros(size(Mu1));
     [emp_ind, a_ind, mass] = find(Mu); % find non-zero indices
     [emp_ind1, a_ind1, mass1] = find(Mu1); % find non-zero indices
     
    for ii = 1:length(emp_ind)
        apr_ind = pol_indx(emp_ind(ii),a_ind(ii)); % which a prime does the policy fn prescribe?
        MuNew(:, apr_ind) = MuNew(:, apr_ind) + ... % which mass of households goes to which exogenous state?
            (PI(emp_ind(ii), :)*mass(ii))';    
    end
    
    for jj = 1:length(emp_ind1)
        apr_ind1 = pol_indx1(emp_ind1(jj),a_ind1(jj)); % which a prime does the policy fn prescribe?
        MuNew1(:, apr_ind1) = MuNew1(:, apr_ind1) + ... % which mass of households goes to which exogenous state?
            (PI(emp_ind1(jj), :)*mass1(jj))';    
    end
    mu_tol= max(max(abs(Mu-MuNew)));
    mu_tol1= max(max(abs(Mu1-MuNew1)));
    Mu=MuNew;
    Mu1=MuNew1;
  end
   %Market clears
   aggsav = w1*sum( pol_fn(:) .* Mu(:) )+w2*sum( pol_fn1(:) .* Mu1(:) ); % Aggregate future assets;
   if aggsav>0;
       q_min=q_guess;
   else q_max=q_guess;
   end
end
%%
%Q2
figure(1)
subplot(1,2,1)
plot(a,v_guess(1,:),'blue')
hold on
plot(a,v_guess(2,:),'red')
legend('Employed','Unemployed','location','southeast')
title(['Value Function'])
hold off
subplot(1,2,2)
plot(a,pol_fn(1,:),'blue',a,pol_fn(2,:),'red')
legend('Employed','Unemployed','location','southeast')
title(['Policy Function'])


%Q3. gini coefficient and lorenz curve
pop=reshape(mu',[2*a_num,1]);
wealth=reshape([a+y_s(1);a+y_s(2)]',[2*a_num,1]);
earning=reshape([repmat(y_s(1),[1 a_num]);repmat(y_s(2),[1 a_num])]',[2*a_num,1]);

%%%%%%% Distribution Plot%%%%%
figure(2)
bar(wealth(1:a_num),mu(1,:),'blue')
hold on
bar(wealth(a_num+1:2*a_num),mu(2,:),'red')
xlim([-2 2.5]);
legend('Employed','Unemployed','location','northwest')
title('Distribution')
hold off
%%%%%%% Gini coefficient and lorenz curve%%%%
figure(4)
suptitle('Lorenz Curve II' )
subplot(1,2,1)
area(WEALTH(:,2),WEALTH(:,3),'FaceColor',[0.5,0.5,1.0])
hold on
plot([0,1],[0,1],'--k')
axis square
title(['Wealth, Gini=',num2str(gini_wealth2)])
hold off
subplot(1,2,2)
area(EARNING(:,2),EARNING(:,3),'FaceColor',[0.5,0.5,1.0])
hold on
plot([0,1],[0,1],'--k')
axis square
title(['Earning, Gini=',num2str(gini_earning2)])
hold off

PI_stat=PI^a_num; %long-run probability
c=PI_stat(1,:)*y_s';  
W_FB=((c^(1-sigma))/(1-sigma))/(1-beta);%first best welfare
lambda=((v_guess.^(-1)).*W_FB).^(1/(1-sigma))-1;%consumption equivalent
Hh_B=sum(sum((lambda>0).*mu));%franction W_FB is better than v(s,a)
WG=sum(sum(lambda.*mu)); % welfare gain

figure(5)
plot(a,lambda(1,:),'blue',a,lambda(2,:),'r')
legend('Employed','Unemployed','location','northeast')
suptitle(['Consumption Equivalent'])
title( ['Welfare Gain=',num2str(WG)])
