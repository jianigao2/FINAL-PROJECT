% PROGRAM NAME: ps4huggett.m
clear all;
close all;
clc;

% PARAMETERS
beta = .95; %discount factor 
sigma = 1.5; % coefficient of risk aversion
b = 0.5; % replacement ratio (unemployment benefits)

w1=0.4539; %rural weight
w2=1-w1; %urban weight

y_s = [1, b]; % endowment in employment states for rural
y_s1 = [1.5, b]; % endowment in employment states for urban

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
q_max = 1;


% ITERATE OVER ASSET PRICES
aggsav = 1 ;
a=repmat(a,[a_num 1 1]);
a=a(:);
while abs(aggsav)>= 0.01
    q_guess = (q_min + q_max) / 2;
    % CURRENT RETURN (UTILITY) FUNCTION
    con = bsxfun(@minus, a, q_guess * a');
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
    while v_tol >=1e-03;
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
        v_mat=ret+beta*repmat(permute((PI*v_guess),[3 2 1]), [a_num*h_num 1 1]);
        v_mat1=ret1+beta*repmat(permute((PI*v_guess1),[3 2 1]), [a_num*h_num 1 1]);
        
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
       [vfn, pol_indx] = max(v_mat, [], 2); %max for each row/rural
       [vfn1, pol_indx1] = max(v_mat1, [], 2); %max for each row/rural
       
        vfn=permute(vfn,[3 1 2]);
        vfn1=permute(vfn1,[3 1 2]);
        
        v_tol=max(abs(v_guess(:)-vfn(:)));
        v_tol1=max(abs(v_guess1(:)-vfn1(:)));
        v_tol=max(v_tol,v_guess);
        
        v_guess=vfn;
        v_guess1=vfn1;
    end;
    
    % KEEP DECSISION RULE
    pol_indx=permute(pol_indx, [3 1 2]);
    pol_indx1=permute(pol_indx1, [3 1 2]);
    
    ab=a';
    pol_fn = ab(pol_indx);
    pol_fn1 = ab(pol_indx1);
    
    % SET UP INITITAL DISTRIBUTION
    Mu=ones(2,a_num*h_num); % Mu is the distribution for a(a,h,y)
    Mu=Mu/sum(Mu(:));
    
    Mu1=ones (2,a_num*h_num); % Mu is the distribution for a(a,h,y)
    Mu1=Mu1/sum(Mu1(:));
    
    mu_tol=1;
  while mu_tol>1e-05 
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
    mu_tol = max(abs(MuNew(:) - Mu(:)));
    mu_tol1 = max(abs(MuNew1(:) - Mu1(:)));
    
    mu_tol=max(mu_tol,mu_tol1);
    
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

