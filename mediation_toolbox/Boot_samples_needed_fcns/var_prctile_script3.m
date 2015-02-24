%% plot B needed againt target alpha inflation for various n

alpha = .05; n = 20; targetu = [.1:.1:1./alpha];

tor_fig;
legstr = {'10' '20' '50' '100' '500' '1000'};
for n = [10 20 50 100 500 1000]

    [B,min_n,Blowerlimit] = var2b(v_from_targetu(targetu,alpha),n,alpha);
    
    plot(targetu,B,'Color',rand(1,3),'LineWidth',2);
    
end
legend(legstr);
xlabel('Target alpha inflation'); ylabel('B needed');


%% plot B needed againt target alpha inflation for various n

alpha = .05; n = 30; targetu = 20;
n = [10 20 30 50 100];

tor_fig;
legstr = {'.05' '.01' '.005' '.001'};
for alpha = [.05 .01 .005 .001]

    [B,min_n,Blowerlimit] = var2b(v_from_targetu(targetu,alpha),n,alpha);
    
    plot(n,B,'Color',rand(1,3),'LineWidth',2);
    
end
legend(legstr);
xlabel('N in sample'); ylabel('B needed');
title(['Target alpha inflate = ' num2str(targetu)]);

%% plot B needed against acceptable 95% upper limit on alpha for various n

alpha = .05; n = 30; abs_perr = .2;    % absolute error in p-value
n = [10:10:100];

tor_fig;
legstr = {'.05' '.01' '.005' '.001'};
for alpha = [.05 .01 .005 .001]

    [B,min_n,Blowerlimit] = var2b(v_from_alphaaccept(abs_perr,alpha),n,alpha);
    
    plot(n,B,'Color',rand(1,3),'LineWidth',2);
    
end
legend(legstr);
xlabel('N in sample'); ylabel('B needed');

title(['Target alpha inflate = ' num2str(targetu)]);

%% plot B against relative contribution of B to error
tor_fig;
B = 200:10:10000;
for n = [10 20 30 50 100]
    [relB,propB,alphainflate] = boot_rel_contrib_B(.05,n,B);
    plot(B,propB,'Color',rand(1,3),'LineWidth',2);
end

legstr = {'10' '20' '30' '50' '100'};
legend(legstr);
xlabel('Bootstrap samples'); ylabel('Prop. contrib. of bootstrap to prctile variance');

%% plot B against relative contribution of B to alpha error
tor_fig;
B = 200:10:10000;
for n = [10 20 30 50 100]
    [relB,propB,alphainflate] = boot_rel_contrib_B(.05,n,B);
    plot(B,alphainflate,'Color',rand(1,3),'LineWidth',2);
end

legstr = {'10' '20' '30' '50' '100'};
legend(legstr);
xlabel('Bootstrap samples'); ylabel('Prop. contrib. of bootstrap to prctile variance');

%% plot B needed so that bootstrap contribution of error to p-value is
%utarget of the alpha value 95% of the time

alpha = [.001:.001:.05];
tor_fig;
legstr = {'.2' '.15' '.10' '.05'};
for utarget = [.2 .15 .1 .05]

    B = Bneeded(alpha,utarget);
    
    plot(alpha,B,'Color',rand(1,3),'LineWidth',2);
    
end
legend(legstr);
xlabel('Alpha value'); ylabel('B needed');

%% plot B needed so that bootstrap contribution of error to p-value is
%utarget of the alpha value on average
colors = {[1 0 0] [1 .5 0] [0 0 1] [0 .6 .1]};
alpha = [.001:.001:.05];
tor_fig;
legstr = {'.2' '.15' '.10' '.05'};
i = 1;
for utarget = [.2 .15 .1 .05]

    [tmp,B] = Bneeded(alpha,utarget);
    
    plot(alpha,B,'Color',colors{i},'LineWidth',2);
    i = i+1;
end
legend(legstr);
xlabel('Alpha value'); ylabel('B needed');

