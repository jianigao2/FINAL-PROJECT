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