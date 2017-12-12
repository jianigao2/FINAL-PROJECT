% PROGRAM NAME: ps4huggett.m
clear all;
close all;
clc;

% PARAMETERS
beta = .9932; %discount factor 
sigma = 1.5; % coefficient of risk aversion
b = 0.5; % replacement ratio (unemployment benefits)
y_s = [1, b]; % endowment in employment states
PI = [.97 .03; .5 .5]; % transition matrix


% ASSET VECTOR
a_lo = -2; %lower bound of grid points
a_hi = 5;%upper bound of grid points
num_a = 1000;

a = linspace(a_lo, a_hi, num_a); % asset (row) vector

% INITIAL GUESS FOR q
q_min = 0.98;
q_max = 1.2;


% ITERATE OVER ASSET PRICES
aggsav = 1 ;
while (abs(aggsav) >= 0.01) ;
    q_guess = (q_min + q_max) / 2;
    % CURRENT RETURN (UTILITY) FUNCTION
    cons = bsxfun(@minus, a', q_guess * a);
    cons = bsxfun(@plus, cons, permute(y_s, [1 3 2]));
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
    ret(cons<0)=-Inf;
    % INITIAL VALUE FUNCTION GUESS
    v_guess = zeros(2, num_a);
    
    % VALUE FUNCTION ITERATION
    v_tol = 1;
    while v_tol >.000001;
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
        value_mat=ret+beta*repmat(permute((PI*v_guess),[3 2 1]), [num_a 1 1]);

        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
       [vfn, pol_indx] = max(value_mat, [], 2); %max for each row
        
       v_tol = max(abs(permute(vfn, [3 1 2])-v_guess));
       v_tol = max(v_tol(:));
       
       v_guess = permute(vfn, [3 1 2]);
    end;
    
    % KEEP DECSISION RULE
    pol_indx=permute(pol_indx, [3 1 2]);
    pol_fn = a(pol_indx);
    
    % SET UP INITITAL DISTRIBUTION
    mu=zeros(2,num_a);
    mu(:)=1/(2*num_a);
    
    dis=1;
  while dis>0.0000001 
      % ITERATE OVER DISTRIBUTIONS
      MuNew = zeros(size(mu));
     [emp_ind, a_ind, mass] = find(mu); % find non-zero indices
    
    for ii = 1:length(emp_ind)
        apr_ind = pol_indx(emp_ind(ii),a_ind(ii)); % which a prime does the policy fn prescribe?
        
        MuNew(:, apr_ind) = MuNew(:, apr_ind) + ... % which mass of households goes to which exogenous state?
            (PI(emp_ind(ii), :)*mass(ii))';
        
    end
    dis = max(max(abs(mu-MuNew)));
    mu=MuNew;
  end
   %Market clears
   aggsav=mu(1,:)*pol_fn(1,:)'+mu(2,:)*pol_fn(2,:)';
   if aggsav>0;
       q_min=q_guess;
   else q_max=q_guess;
   end
end

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
pop=reshape(mu',[2*num_a,1]);
wealth=reshape([a+y_s(1);a+y_s(2)]',[2*num_a,1]);
earning=reshape([repmat(y_s(1),[1 num_a]);repmat(y_s(2),[1 num_a])]',[2*num_a,1]);

%%%%%%% Distribution Plot%%%%%
figure(2)
bar(wealth(1:num_a),mu(1,:),'blue')
hold on
bar(wealth(num_a+1:2*num_a),mu(2,:),'red')
xlim([-2 2.5]);
legend('Employed','Unemployed','location','northwest')
title('Distribution')
hold off
%%%%%%% Gini coefficient and lorenz curve%%%%
% Type I. without non-negative wealth
wealth1=wealth;
wealth1(wealth1<0)=0;
figure(3)
suptitle('Lorenz Curve I' )
subplot(1,2,1)
gini_wealth1=gini(pop, wealth1,true);% download gini.m online
title(['Wealth, gini=',num2str(gini_wealth1)])
subplot(1,2,2)
gini_earning1=gini(pop, earning,true);
title(['Wealth, gini=',num2str(gini_earning1)])
% Type II. with negative wealth
WEALTH=sortrows([wealth,pop,pop.*wealth]);
WEALTH=cumsum(WEALTH);
pw=WEALTH(:,2);
pw=pw(end);
WEALTH(:,2)=WEALTH(:,2)/pw;
w=WEALTH(:,3);
w=w(end);
WEALTH(:,3)=WEALTH(:,3)/w;
gini_wealth2 = 1 - sum((WEALTH(1:end-1,3)+WEALTH(2:end,3)) .* diff(WEALTH(:,2)));

EARNING=sortrows([earning,pop,pop.*earning]);
EARNING=cumsum(EARNING);
pe=EARNING(:,2);
pe=pe(end);
EARNING(:,2)=EARNING(:,2)/pe;
e=EARNING(:,3);
e=e(end);
EARNING(:,3)=EARNING(:,3)/e;
gini_earning2 = 1 - sum((EARNING(1:end-1,3)+EARNING(2:end,3)) .* diff(EARNING(:,2)));

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

PI_stat=PI^num_a; %long-run probability
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
