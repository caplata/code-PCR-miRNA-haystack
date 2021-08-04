%% Experimental & Computational Conditions
clear all
close all

cNa=0.15; %NaCl concentration
Lo=1000;%length of the target
num_ol=11000; %number of target chains
N_tar=300; %number of target sites (inside a chain) 
L=8; %length of  seed
T_body=37; % °C

R=1.987; %gas costant in cal K^-1 mol^-1
Tk = 273.15; % kelvin-celsius

Tmin=0; %min temperature
Tmax=100; %max temperature
NT = 200; %temperature steps


%% Fig. 7: L primer dependence
L_max=20;
L_prim=[1:L_max];

for i=1:L_max
    [z_perf_L(i), z_med_L(i)]=weights(L_prim(i), cNa, T_body);
end

p_perf_L=N_tar*z_perf_L./(N_tar*z_perf_L+Lo*num_ol*z_med_L); 

p_perf_L_noT=z_perf_L./(z_perf_L+Lo*num_ol*z_med_L);  
p_perf_L_100=100*z_perf_L./(100*z_perf_L+Lo*num_ol*z_med_L); 

figure (7)
plot(L_prim(:), p_perf_L(:),'-o', 'LineWidth',3);
hold on
plot(L_prim(:), p_perf_L_100(:),'-o', 'LineWidth',3);
plot(L_prim(:), p_perf_L_noT(:),'-d', 'LineWidth',3);


text(L-3,0.64, sprintf('L = %i ',L), 'FontSize',15);
val=L*ones(6,1);
valn=(0:5)/5;
plot(val(:), valn(:),':','Color', [0.2, 0.2, 0.2], 'LineWidth',1.5, 'Handlevisibility', 'off');

grid on
xlabel('L')
ylabel('\phi_0')
xlim([0  L_max])
title('miRNA Pairing Statistics')  
legend('N_{tar}=300','N_{tar}=100','N_{tar}=1','Location','southeast')
set(gca, 'FontSize', 15)


%% Fig. 8: dependence on the system degeneracy and T

g_body=Lo*num_ol/N_tar;

g_body_var=round((logspace(0, 6.8)));

[z_perf_Lo, z_med_Lo]=weights(L, cNa, T_body);

p_perf_Lo_n=z_perf_Lo./(z_perf_Lo+g_body_var.*z_med_Lo); 

figure (8)
plot(g_body_var(:), p_perf_Lo_n(:), 'LineWidth',3);
hold on


text(g_body*1.5,0.75, sprintf('L_{eff} = %0.1e ',g_body), 'FontSize',15);
val=g_body*ones(2,1);
valn=[10^-4 , 1];
plot(val(:), valn(:),':','Color', [0.2, 0.2, 0.2], 'LineWidth',1.5);

grid on
xlabel('L_{eff} = L_0 / N_{tar}')
ylabel('\phi_0')
xlim([g_body_var(1) g_body_var(end)])
title('miRNA Pairing Statistics')  
set(gca, 'Xscale', 'Log','FontSize', 15)

T_var=(linspace(Tmin, Tmax, NT));

for i=1:NT
    [z_perf(i), z_med(i)]=weights(L, cNa, T_var(i));
end

p_perf=N_tar*z_perf./(N_tar*z_perf+Lo*num_ol*z_med); 
p_med=Lo*num_ol*z_med./(N_tar*z_perf+Lo*num_ol*z_med); 


axes('Position',[.25 .23 .35 .30])
box on
plot(T_var(:), p_perf(:), 'LineWidth',3);
hold on
grid on
text(T_body-5, 0.54, sprintf('T = %i °C',T_body), 'FontSize',15);
val=T_body*ones(6,1);
valn=(0:5)/5;
plot(val(:), valn(:),':','Color', [0.2, 0.2, 0.2], 'LineWidth',1.5);

set(gca, 'FontSize', 15)
xlabel('T (°C)','FontSize', 15)
ylabel('\phi_0')
xlim([30  42])
ylim([0.52  0.61])

%% functions: DH, DS and weights
function [z_perf, z_med]=weights(L, cNa, T)

    R=1.987; %gas constant in cal K^-1 mol^-1
    Tk = 273.15; % kelvin-celsius

    z_perf = exp(-(DH_f(L,0,0,0)./(T+Tk)-DS_f(L,0,0,0,cNa))/R);
    z_rand = 0;

    for ae1 = 0:(L-1)
        for ae2 = 0:(L-1-ae1)
            for ai = 0:max(0, L-2-ae1-ae2)
                z_rand=z_rand+3^(ae1+ae2+ai)*nchoosek(max(0,L-2-ae1-ae2),ai)...
                    *exp(-(DH_f(L, ae1, ae2, ai)./(T+Tk)-DS_f(L,ae1,ae2,ai,cNa))/R);
            end
        end
    end
    
    for ae1=0:(L-1)
        z_rand=z_rand+3^(L-1)...
            *exp(-(DH_f(L,ae1,L-1-ae1,0)./(T+Tk)-DS_f(L,ae1,L-1-ae1,0,cNa))/R);
    end    
    
    z_med=z_rand/4^L;

end



%total enthaply and entropy: ATCG averaged
function DH = DH_f(L, ae1, ae2, ai)

%enthalpic contributions in kcal mol^-1 :
DH0 = 7.33;
DHq = -10.783125; %average
DHs = -2.534;%average
DHi = 0.154166667;

%DHpen=2.2; %AT penalty

DH=1000*(DH0+(L-1-ae1-ae2-ci_f(ai))*DHq+cs_f(ae1, ae2)*DHs+ci_f(ai)*DHi);
 
end

function DS = DS_f(L, ae1, ae2, ai, cNa)

%entropic contributions in cal K^-1 mol^-1 :
DS0 = 9;
DSq = -27.8875; %media
DSi = -0.841666667;
DSs=-6.749755621;

DS=DS0+(L-1-ae1-ae2-ci_f(ai))*DSq+cs_f(ae1, ae2)*DSs+ci_f(ai)*DSi;

f_CG=0.5;
DS=DS+DH_f(L, ae1, ae2, ai)*((4.29*f_CG-3.95)*1e-5*log(cNa) +9.4*1e-6*(log(cNa))^2 ); %Owczarczy
end


%number of shift/external penalities (0, 2, 3, 4)
function cs = cs_f(ae1, ae2)

cs=4;
              
end

%number of internal penalities 
function ci = ci_f(ai)

%ci=ai;% consecutive approx (0, 2, 3, 4, 5,...)

%if ai>0
%    ci=ci+1;
%end

ci=2*ai;% non-consecutive approx (0, 2, 4, 6,...)
 
end

