% Boldrin-Montes model with family altruism (Dec 1, 2019)
% Written by Michael Hatcher (m.c.hatcher@soton.ac.uk)

clc
clear
close all

A = 10; 
B = 2.5;
betta = 0.30;
deltta = 2/3; deltta1 = 2/3;
mu = 1;
zetta = 0.60; zetta1 = 0.90;
n = 0.5;

%Credit constraint
lambda = 0.1;
lambda1 = 0.1;

lambda_crit = (1-zetta)*(1-deltta)/(1-deltta*(1-zetta));
zetta_crit = deltta*(1-deltta)*(1-lambda)/(1-deltta^2*(1-lambda));
lambda_crit1 = (1-zetta)*(1-deltta1)/(1-deltta1*(1-zetta));
zetta_crit1 = deltta1*(1-deltta1)*(1-lambda1)/(1-deltta1^2*(1-lambda1));

k_init = 1;
h_init = 1;

gama_star = betta*deltta*(zetta-lambda)/(1-deltta*(1-lambda));
gama_star1 = betta*deltta1*(zetta1-lambda1)/(1-deltta1*(1-lambda1));
ngama = 2000;
gama_max = 5;

for j = 1:ngama
    
    gama = gama_max*j/ngama;
    gama_vec(j) = gama;
    
    %Complete markets

coef_k(j) = A*(1-deltta)*(deltta*betta*(1-zetta) + gama)/( (1-deltta*(1-zetta))*(1+betta+gama) );

coef_h(j) = ( A*B^(1/zetta)*zetta*deltta*(deltta*betta*(1-zetta) + gama)/( (1-deltta*(1-zetta))*(1+betta+gama)*(1+n) ) )^zetta;

xstar(j) = (coef_k(j)/ ((1+n)*coef_h(j)) )^(1/(deltta+zetta*(1-deltta)));

Gstar(j) = (coef_k(j)/(1+n))*xstar(j)^(-deltta);

Gstar_check(j) = (coef_k(j)/(1+n))^( zetta*(1-deltta)/(deltta+zetta*(1-deltta)) )*coef_h(j)^(deltta/(deltta+zetta*(1-deltta)));

Gstar_ann(j) = 100*(Gstar(j)^(1/30) - 1);

Ratio(j) = (1-deltta)*A*xstar(j)^(-deltta)/(Gstar(j)*(1+n));  

R(j) = Ratio(j)*(1+n)*Gstar(j);

    %Incomplete markets
    
coef_k_ic(j)  = A*betta*deltta*(1-lambda)/( (1+betta+gama)*(1+deltta*lambda/(1-deltta)) );

coef_h_ic(j)  = ( ( A*B^(1/zetta)*(1-lambda)*deltta*( gama + (betta+gama)*deltta*lambda/(1-deltta) )  ) / ( (1+betta+gama)*(1+deltta*lambda/(1-deltta))*(1+n) ) )^zetta;

xstar_ic(j)  = (coef_k_ic(j) / ((1+n)*coef_h_ic(j) ) )^(1/(deltta+zetta*(1-deltta)));

Gstar_ic(j)  = (coef_k_ic(j) /(1+n))*xstar_ic(j) ^(-deltta);

Gstar_check_ic(j)  = (coef_k_ic(j) /(1+n))^( zetta*(1-deltta)/(deltta+zetta*(1-deltta)) )*coef_h_ic(j)^(deltta/(deltta+zetta*(1-deltta)));

Gstar_ann_ic(j) = 100*(Gstar_ic(j)^(1/30) - 1);

Ratio_ic(j)  = (1-deltta)*A*xstar_ic(j)^(-deltta)/(Gstar_ic(j) *(1+n));  

Multiplier(j)  = Gstar_ic(j)*(1+n)*(  ( zetta*(1+betta+gama)*(1+deltta*lambda/(1-deltta)) )/( (1-lambda)*(gama + (gama+betta)*lambda*deltta/(1-deltta)) ) - (1-deltta)*A/coef_k_ic(j) );

Unconstrained(j) = 0;

Gstar_ic2(j)  = (coef_k_ic(j) /(1+n))*xstar_ic(j) ^(-deltta);

Gstar_check_ic2(j)  = (coef_k_ic(j) /(1+n))^( zetta*(1-deltta)/(deltta+zetta*(1-deltta)) )*coef_h_ic(j)^(deltta/(deltta+zetta*(1-deltta)));

Gstar_ann_ic2(j) = 100*(Gstar_ic(j)^(1/30) - 1);

Ratio_ic2(j)  = (1-deltta)*A*xstar_ic(j)^(-deltta)/(Gstar_ic(j) *(1+n));

R_ic(j) = Ratio_ic2(j)*(1+n)*Gstar_ic2(j);


if Multiplier(j)  <= 0
    Gstar_ic2(j)  = Gstar(j);
    Gstar_check_ic2(j)  = Gstar_check(j);
    Gstar_ann_ic2(j) = Gstar_ann(j);
    Unconstrained(j)  = 1;
    Ratio_ic2(j) = Ratio(j);
    R_ic2(j) = R(j);
end


%With different parameter values

    %Complete markets

coef_k1(j) = A*(1-deltta1)*(deltta1*betta*(1-zetta1) + gama)/( (1-deltta1*(1-zetta1))*(1+betta+gama) );

coef_h1(j) = ( A*B^(1/zetta)*zetta1*deltta1*(deltta1*betta*(1-zetta1) + gama)/( (1-deltta1*(1-zetta1))*(1+betta+gama)*(1+n) ) )^zetta1;

xstar1(j) = (coef_k1(j)/ ((1+n)*coef_h1(j)) )^(1/(deltta1+zetta1*(1-deltta1)));

Gstar1(j) = (coef_k1(j)/(1+n))*xstar1(j)^(-deltta1);

Gstar_check1(j) = (coef_k1(j)/(1+n))^( zetta1*(1-deltta1)/(deltta1+zetta1*(1-deltta1)) )*coef_h1(j)^(deltta1/(deltta1+zetta1*(1-deltta1)));

Gstar_ann1(j) = 100*(Gstar1(j)^(1/30) - 1);

Ratio1(j) = (1-deltta1)*A*xstar1(j)^(-deltta1)/(Gstar1(j)*(1+n));  

    %Incomplete markets
    
coef_k_ic1(j)  = A*betta*deltta1*(1-lambda1)/( (1+betta+gama)*(1+deltta1*lambda1/(1-deltta1)) );

coef_h_ic1(j)  = ( ( A*B^(1/zetta1)*(1-lambda1)*deltta1*( gama + (betta+gama)*deltta1*lambda1/(1-deltta1) )  ) / ( (1+betta+gama)*(1+deltta1*lambda1/(1-deltta1))*(1+n) ) )^zetta1;

xstar_ic1(j)  = (coef_k_ic1(j) / ((1+n)*coef_h_ic1(j) ) )^(1/(deltta1+zetta1*(1-deltta1)));

Gstar_ic1(j)  = (coef_k_ic1(j) /(1+n))*xstar_ic1(j)^(-deltta1);

Gstar_check_ic1(j)  = (coef_k_ic1(j) /(1+n))^( zetta1*(1-deltta1)/(deltta1+zetta1*(1-deltta1)) )*coef_h_ic1(j)^(deltta1/(deltta1+zetta1*(1-deltta1)));

Gstar_ann_ic1(j) = 100*(Gstar_ic1(j)^(1/30) - 1);

Ratio_ic1(j)  = (1-deltta1)*A*xstar_ic1(j)^(-deltta1)/(Gstar_ic1(j)*(1+n));  

Multiplier1(j)  = Gstar_ic1(j)*(1+n)*(  ( zetta1*(1+betta+gama)*(1+deltta1*lambda1/(1-deltta1)) )/( (1-lambda1)*(gama + (gama+betta)*lambda1*deltta1/(1-deltta1)) ) - (1-deltta1)*A/coef_k_ic1(j) );

Unconstrained1(j) = 0;

Gstar_ic21(j)  = (coef_k_ic1(j)/(1+n))*xstar_ic1(j) ^(-deltta1);

Gstar_check_ic21(j)  = (coef_k_ic1(j) /(1+n))^( zetta1*(1-deltta1)/(deltta1+zetta1*(1-deltta1)) )*coef_h_ic1(j)^(deltta1/(deltta1+zetta1*(1-deltta1)));

Gstar_ann_ic21(j) = 100*(Gstar_ic1(j)^(1/30) - 1);

Ratio_ic21(j)  = (1-deltta1)*A*xstar_ic1(j)^(-deltta1)/(Gstar_ic1(j)*(1+n));


if Multiplier1(j)  <= 0
    Gstar_ic21(j)  = Gstar1(j);
    Gstar_check_ic21(j)  = Gstar_check1(j);
    Gstar_ann_ic21(j) = Gstar_ann1(j);
    Unconstrained1(j)  = 1;
    Ratio_ic21(j) = Ratio1(j);
    
end


end

figure(1)
title('Growth rates versus altruism: complete mkts vs incomplete mkts (red)')
hold on,
plot(gama_vec, Gstar), plot(gama_vec, Gstar_ic, 'r')
xlabel('Altruism parameter (gamma)')
ylabel('Growth rate')
hold off

figure(2)
title('Growth rates versus altruism: complete mkts vs incomplete mkts (red)')
hold on,
plot(gama_vec, Gstar), plot(gama_vec, Gstar_ic2, '--r') 
xlabel('Altruism parameter (gamma)')
ylabel('Growth rate')
hold off

%%Check
%figure(2)
%hold on,
%plot(gama_vec, Gstar_check), plot(gama_vec, Gstar_check_ic, 'r')
%hold off

figure(2)
title('Growth rates versus altruism: complete mkts vs incomplete mkts (red)')
hold on,
subplot(2,2,1)
title('Growth rates versus altruism: complete mkts vs incomplete mkts (red)')
hold on,
plot(gama_vec(1:end), Gstar(1:end)), hold on, 
plot(gama_vec(1:end), Gstar_ic2(1:end), '--r')
%axes('position', [0.25,0.25,0.25,0.25])
%hold on, box on
your_index = gama_vec > 0 & gama_vec <= 0.7;
hold on, 
%plot(gama_vec(your_index), Gstar_ann(your_index)), plot(gama_vec(your_index), Gstar_ann_ic(your_index), 'r')
axis tight
subplot(2,2,2)
title('Growth rates versus altruism: complete mkts vs incomplete mkts (red)')
hold on,
plot(gama_vec(1:end), Ratio(1:end)), hold on, 
plot(gama_vec(1:end), Ratio_ic2(1:end), '--r')
%axes('position', [0.15,0.15,0.15,0.15])
%hold on, box on
your_index2 = gama_vec > 0 & gama_vec <= 2;
hold on, 
%plot(gama_vec(your_index2), Ratio(your_index2)), plot(gama_vec(your_index2), Ratio_ic(your_index2), 'r')
axis tight
hold off
hold on,
subplot(2,2,3)
title('Growth rates versus altruism: complete mkts vs incomplete mkts (red)')
hold on,
plot(gama_vec(1:end), Gstar1(1:end)), hold on, 
plot(gama_vec(1:end), Gstar_ic21(1:end), '--r')
%axes('position', [0.25,0.25,0.25,0.25])
%hold on, box on
your_index = gama_vec > 0 & gama_vec <= 0.7;
%hold on, 
%plot(gama_vec(your_index), Gstar_ann1(your_index)), plot(gama_vec(your_index), Gstar_ann_ic21(your_index), 'r')
axis tight
subplot(2,2,4)
title('Growth rates versus altruism: complete mkts vs incomplete mkts (red)')
hold on,
plot(gama_vec(1:end), Ratio1(1:end)), hold on, 
plot(gama_vec(1:end), Ratio_ic21(1:end), '--r')
%axes('position', [0.15,0.15,0.15,0.15])
%hold on, box on
your_index2 = gama_vec > 0 & gama_vec <= 2;
%hold on, 
%plot(gama_vec(your_index2), Ratio1(your_index2)), plot(gama_vec(your_index2), Ratio_ic21(your_index2), 'r')
axis tight
hold off

figure(3)
title('Growth rates versus altruism: complete mkts vs incomplete mkts (red)')
hold on,
plot(gama_vec, Gstar), plot(gama_vec, Gstar_ic, 'r'), plot(gama_vec, Gstar_ic2, '--g') 
xlabel('Altruism parameter (gamma)')
ylabel('Growth rate')
hold off

figure(4)
title('Dynamic efficiency check: complete mkts vs incomplete mkts (red)')
hold on,
plot(gama_vec, Ratio), plot(gama_vec, Ratio_ic2, '--r')
xlabel('Altruism parameter (gamma)')
ylabel('Ratio of MPK to (1+g)(1+n)')
hold off

figure(5)
title('Growth rates versus altruism: complete mkts vs incomplete mkts (red)')
hold on,
subplot(1,2,1)
title('Growth rates versus altruism: complete mkts vs incomplete mkts (red)')
hold on,
plot(gama_vec(1:end), Gstar(1:end)), hold on, 
plot(gama_vec(1:end), Gstar_ic2(1:end), '--r')
%axes('position', [0.25,0.25,0.25,0.25])
%hold on, box on
your_index = gama_vec > 0 & gama_vec <= 0.7;
hold on, 
%plot(gama_vec(your_index), Gstar_ann(your_index)), plot(gama_vec(your_index), Gstar_ann_ic(your_index), 'r')
axis tight
subplot(1,2,2)
title('Growth rates versus altruism: complete mkts vs incomplete mkts (red)')
hold on,
plot(gama_vec(1:end), Gstar1(1:end)), hold on, 
plot(gama_vec(1:end), Gstar_ic21(1:end), '--r')
%axes('position', [0.25,0.25,0.25,0.25])
%hold on, box on
your_index = gama_vec > 0 & gama_vec <= 0.7;
%hold on, 
%plot(gama_vec(your_index), Gstar_ann1(your_index)), plot(gama_vec(your_index), Gstar_ann_ic21(your_index), 'r')
axis tight
hold off
