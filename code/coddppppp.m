%%%%% value function iteration, two-dimensional grid over k,n
clear;
close all;

%%%%%parameters%%%%%
beta=0.95;
sigma=1.5;
d=0.3;
eta=0.25;
b=0.5;
y=[1,b];
pi=[0.97,0.03;0.5,0.5];
pi_inv=pi^1000;
pi_inv=pi_inv(1,:);

%%%%%set up the grids%%%
h_min=y(2)/5;
h_max=5;
h_num=11;
h=linspace(h_min,h_max,h_num);

a_min=-2;
a_max=3;
a_num=11;
a=linspace(a_min,a_max,a_num);

%%%initial guess for q
q_min=0.93;
q_max=1;
q_guess=(q_min+q_max)/2;

 a=repmat(a,[a_num 1 1]);
 a=a(:);

%%%%iterate over assets prices
aggsav=1;
while abs(aggsav)>=0.01
    
    cons=bsxfun(@minus,a,q_guess*a');
    hou=repmat(h',[1 h_num 1]);
    hous=repmat(hou,h_num);
    cons=cons-hous-eta.*(abs(hous'-hous));
    cons=bsxfun(@plus,cons,permute(y,[1 3 2]));
    hous=repmat(hous,[1 1 2]);
    ret=((cons.^(1-d)).*(hous.^d)).^(1-sigma)./(1-sigma);
    ret(cons<0)=-Inf;
    
    v_guess=zeros(2,a_num*h_num);
    
    %%%value function iteration
    v_tol=1;
    while v_tol>=0.001
        v_mat=ret+beta*repmat(permute(pi*v_guess,[3 2 1]),[a_num*h_num 1 1]);
        [vfn,pol_ind]=max(v_mat,[],2);
        vfn=permute(vfn,[3 1 2]);
        v_tol=abs(max(v_guess(:)-vfn(:)));
        v_guess=vfn;
    end
    pol_ind=permute(pol_ind,[3,1,2]);
    aa=a';
    pol_fn=aa(pol_ind);
    
   
    
    %%%%distribution
    Mu=ones (2,a_num*h_num);
    Mu=Mu/sum(Mu(:));
    
    mu_tol = 1;
   
   while mu_tol > 1e-8
   [emp_ind, a_ind,mass] = find(Mu); % find non-zero indices
    
   MuNew = zeros(size(Mu));
      for ii = 1:length(emp_ind)
        apr_ind = pol_ind(emp_ind(ii), a_ind(ii)); 
        MuNew(:, apr_ind) = MuNew(:, apr_ind) + ...
            (pi(emp_ind(ii), :) * mass(ii))';
      end
    

    mu_tol = max(abs(MuNew(:) - Mu(:)));
    
    Mu = MuNew ;
   end
    aggsav = sum( pol_fn(:) .* Mu(:) ); % Aggregate future assets

    if aggsav > 0 ;
    q_min = q_guess ;
    end ;
    if aggsav < 0;
    q_max = q_guess ;
    end ;
    q_guess = (q_max + q_min)/2 ;
    
    
end



    



