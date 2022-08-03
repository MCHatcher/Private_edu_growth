% Boldrin-Montes model with family altruism (Dec 1, 2019)
% Written by Michael Hatcher (m.c.hatcher@soton.ac.uk)

clc
clear
close all

A = 10; B = 2.5;
A1 = 10; B1 = 2.5;
betta = 0.30;
deltta = 2/3;
gama = 0.15;
mu = 1; zetta = 0.60;
n = 0.5;

k_init = 1;
h_init = 1;

ngama = 80000;
gama_max = 10;

for j = 1:ngama
    
    gama = gama_max*j/ngama;
    gama_vec(j) = gama;
    
    %Complete markets

coef_k(j) = A*(1-deltta)*(deltta*betta*(1-zetta) + gama)/( (1-deltta*(1-zetta))*(1+betta+gama) );

coef_h(j) = ( A*B^(1/zetta)*zetta*deltta*(deltta*betta*(1-zetta) + gama)/( (1-deltta*(1-zetta))*(1+betta+gama)*(1+n) ) )^zetta;

xstar(j) = (coef_k(j)/ ((1+n)*coef_h(j)) )^(1/(deltta+zetta*(1-deltta)));

Gstar(j) = (coef_k(j)/(1+n))*xstar(j)^(-deltta);

Gstar_check(j) = (coef_k(j)/(1+n))^( zetta*(1-deltta)/(deltta+zetta*(1-deltta)) )*coef_h(j)^(deltta/(deltta+zetta*(1-deltta)));

Ratio(j) = (1-deltta)*A*xstar(j)^(-deltta)/(Gstar(j)*(1+n));

Ratio_check(j) = (1-deltta)/coef_k(j)/A;

Gstar_ann(j) = 100*(Gstar(j)^(1/30) - 1);

Line(j) = Ratio(j);

if Ratio(j) >= 1 && gama_vec(j) < deltta*zetta*betta 
    %Line(j) = Ratio(j);
    Line(j) = NaN;
end

if Ratio(j) < 1 
    Line(j) = NaN;
end

if gama_vec(j) >= deltta*zetta*betta  && Ratio(j) > gama_vec(j)/(deltta*zetta*betta) || Ratio(j) <1
    Line(j) = NaN;
end

Line1(j) = NaN;

if Ratio(j) >= 1 && gama_vec(j) < deltta*zetta*betta 
    %Line1(j) = NaN;
    Line1(j) = Ratio(j);
end

if Ratio(j) < 1 
    Line1(j) = Ratio(j);
end

if gama_vec(j) >= deltta*zetta*betta  && Ratio(j) > gama_vec(j)/(deltta*zetta*betta)  || Ratio(j) <1
    Line1(j) = Ratio(j);
end

    %Incomplete markets
    
coef_k_ic(j) = A1*betta*deltta/( (1+betta+gama) );

coef_h_ic(j) = ( A1*B1^(1/zetta)*gama*deltta/ ( (1+betta+gama)*(1+n) ) )^zetta;

xstar_ic(j) = (coef_k_ic(j)/ ((1+n)*coef_h_ic(j)) )^(1/(deltta+zetta*(1-deltta)));

Gstar_ic(j) = (coef_k_ic(j)/(1+n))*xstar_ic(j)^(-deltta);

Gstar_check_ic(j) = (coef_k_ic(j)/(1+n))^( zetta*(1-deltta)/(deltta+zetta*(1-deltta)) )*coef_h_ic(j)^(deltta/(deltta+zetta*(1-deltta))); 

Ratio_ic(j) = (1-deltta)*A1*xstar_ic(j)^(-deltta)/(Gstar_ic(j)*(1+n));

Ratio_ic_check(j) = (1-deltta)/coef_k_ic(j)/A1;

Gstar_ann_ic(j) = 100*(Gstar_ic(j)^(1/30) - 1);

gama1 = ( betta*deltta - (1-deltta)*(1+betta) )/(1-deltta);
gama2 = zetta*(1-deltta*(1+betta))/(1-zetta*(1-deltta));

Line2(j) = Ratio_ic(j);

if gama_vec(j) < max(gama1,gama2) 
    %Line2(j) = Ratio_ic(j);
    Line2(j) = NaN;
end

if Ratio_ic(j) < 1 
    Line2(j) = NaN;
end

if gama_vec(j) >= max(gama1,gama2)  
    Line2(j) = Ratio_ic(j);
end

Line3(j) = NaN;

if gama_vec(j) < max(gama1,gama2) 
    Line3(j) = Ratio_ic(j);
end

if Ratio_ic(j) < 1 
    Line3(j) = Ratio_ic(j);
end

if gama_vec(j) >= max(gama1,gama2)   
    Line3(j) = NaN;
end

end

figure(1)
title('Growth rate along BGP')
hold on,
plot(gama_vec, Gstar), hold on, 
plot(gama_vec, Gstar_ic, 'r')
axes('position', [0.25,0.25,0.25,0.25])
box on
your_index = gama_vec > 0 & gama_vec <= 0.7;
hold on, 
plot(gama_vec(your_index), Gstar(your_index)), plot(gama_vec(your_index), Gstar_ic(your_index), 'r')
axis tight

figure(2)
%title('Growth rate along BGP')
hold on,
subplot(1,2,1)
title('Growth rate along BGP')
hold on,
plot(gama_vec, Gstar), hold on, 
plot(gama_vec, Gstar_ic, 'r')
axes('position', [0.25,0.25,0.25,0.25])
hold on, box on
your_index = gama_vec >= 0 & gama_vec <= 0.6;
hold on, 
plot(gama_vec(your_index), Gstar(your_index)), plot(gama_vec(your_index), Gstar_ic(your_index), 'r')
axis tight
subplot(1,2,2)
title('Dynamic efficiency of BGP: solid line')
hold on,
plot(gama_vec, Line), hold on, 
plot(gama_vec, Line2, 'r'), hold on,
plot(gama_vec, Line1, '--b'), hold on, 
plot(gama_vec, Line3, '--r')
axes('position', [0.15,0.15,0.15,0.15])
hold on, box on
your_index2 = gama_vec >= 0 & gama_vec <= 0.6;
hold on, 
plot(gama_vec(your_index2), Line(your_index2)), plot(gama_vec(your_index2), Line2(your_index2), 'r'), hold on,
plot(gama_vec(your_index2), Line1(your_index2), '--b'), plot(gama_vec(your_index2), Line3(your_index2), '--r')
axis tight
hold off

figure(3)
title('Growth rate along BGP')
hold on,
plot(gama_vec, Gstar), hold on, 
plot(gama_vec, Gstar_ic, 'r')
axes('position', [0.25,0.25,0.25,0.25])
box on
your_index = gama_vec > 0 & gama_vec <= 0.7;
hold on, 
plot(gama_vec(your_index), Gstar(your_index)), plot(gama_vec(your_index), Gstar_ic(your_index), 'r')
axis tight

%%Check
%figure(2)
%hold on,
%plot(gama_vec, Gstar_check), plot(gama_vec, Gstar_check_ic, 'r')
%hold off

figure(4)
title('Dynamic efficiency of BGP: solid line')
hold on,
plot(gama_vec, Ratio), plot(gama_vec, Ratio_ic, 'r')
xlabel('Altruism parameter (gamma)')
ylabel('Ratio of MPK to (1+g)(1+n)')
hold off

figure(5)
plot(gama_vec, Line), hold on, plot(gama_vec, Line1, '--r'), hold on, plot(gama_vec, Line2, 'g'), hold on, plot(gama_vec, Line3, '--y'), 
