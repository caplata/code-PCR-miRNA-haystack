%% Experimental & Computational Conditions
clear all
close all

cNa=0.055; %NaCl concentration
Lo=2*(3*10^9); %Genome length, x2
L=20; %Primer Length
T_pcr=55; % °C

R=1.987; %gas constant in cal K^-1 mol^-1
Tk = 273.15; %kelvin-celsius

Tmin=0; %min temperature
Tmax=100; %max temperature
NT = 200; %temperature steps

%% Pairing statistics: all contributions (including wrong 3')

%% Fig. 3: L primer dependence
L_prim=[1:40];

for i=1:40
    [z_perf_L(i), z_med_L(i), z_e1_L(i), z_no_e1_L(i)]=weights(L_prim(i), cNa, T_pcr);
end

p_perf_L=z_perf_L./(z_perf_L+Lo*z_med_L); 
p_med_L=Lo*z_med_L./(z_perf_L+Lo*z_med_L); 
p_e1_L=Lo*z_e1_L./(z_perf_L+Lo*z_med_L); 
p_no_e1_L=Lo*z_no_e1_L./(z_perf_L+Lo*z_med_L); 

figure (3)
plot(L_prim(:), p_perf_L(:),'-o', 'LineWidth',3);
hold on
plot(L_prim(:), p_med_L(:), '-.o','LineWidth',3);
plot(L_prim(:), p_e1_L(:),'-d','LineWidth',2);
plot(L_prim(:), p_no_e1_L(:),'-d','LineWidth',2);

text(L-6,0.3, sprintf('L = %i ',L), 'FontSize',15);
val=L*ones(6,1);
valn=(0:5)/5;
plot(val(:), valn(:),':','Color', [0.2, 0.2, 0.2], 'LineWidth',1.5);

grid on
xlabel('L')
ylabel('\phi')
xlim([0  40])
title('PCR Pairing Statistics')  
legend('\phi_0','1-\phi_0', '\phi_{(\alpha_{e1}>0)}','\phi_{(\alpha_{e1}=0)}','Location','northeast')
set(gca, 'FontSize', 15)


%% distrib CG (solo sulle quartine): solo 3' giusto (t.c. la polimerasi produca)

%% Fig. 4: L primer dependence
c_p=[0 0.4470 0.7410];
c_err=[0.8500 0.3250 0.0980];
L_prim=[1:40];

for i=1:40
    [z_perf_L_CG(i), z_no_e1_L_CG(i)]=weights_CG_quart(L_prim(i), cNa, T_pcr, L_prim(i)/2+2);
    [z_perf_L_00(i), z_no_e1_L_00(i)]=weights_CG_quart(L_prim(i), cNa, T_pcr, L_prim(i)/2);
    [z_perf_L_AT(i), z_no_e1_L_AT(i)]=weights_CG_quart(L_prim(i), cNa, T_pcr, L_prim(i)/2-2);
end

p_perf_L_00=z_perf_L_00./(z_perf_L_00+Lo*z_no_e1_L_00); 
p_no_e1_L_00=Lo*z_no_e1_L_00./(z_perf_L_00+Lo*z_no_e1_L_00); 
p_perf_L_CG=z_perf_L_CG./(z_perf_L_CG+Lo*z_no_e1_L_CG); 
p_no_e1_L_CG=Lo*z_no_e1_L_CG./(z_perf_L_CG+Lo*z_no_e1_L_CG); 
p_perf_L_AT=z_perf_L_AT./(z_perf_L_AT+Lo*z_no_e1_L_AT); 
p_no_e1_L_AT=Lo*z_no_e1_L_AT./(z_perf_L_AT+Lo*z_no_e1_L_AT); 

figure (4)
plot(L_prim(:), p_perf_L_00(:), '-o', 'LineWidth',3, 'Color', c_p);
hold on
plot(L_prim(:), p_perf_L_CG(:), '--o', 'LineWidth',2,'Color', c_p);
plot(L_prim(:), p_perf_L_AT(:), ':o', 'LineWidth',2,'Color', c_p);

text(L+2.5,0.25, sprintf('L = %i ',L), 'FontSize',15);
val=L*ones(6,1);
valn=(0:5)/5;
plot(val(:), valn(:),':','Color', [0.2, 0.2, 0.2], 'LineWidth',1.5);

grid on
xlabel('L')
ylabel('\phi_{0}'' ')
xlim([5  30])
title('PCR Pairing Statistics')  
legend('n_{CG}=L/2','n_{CG}=L/2+2','n_{CG}=L/2-2','Location','northwest')
set(gca, 'FontSize', 15)


%% Fig. 5: T dependence

T_var=(linspace(Tmin, Tmax, NT));

for i=1:NT
    [z_perf_CG(i), z_no_e1_CG(i)]=weights_CG_quart(L, cNa, T_var(i), L/2+2);
    [z_perf_00(i), z_no_e1_00(i)]=weights_CG_quart(L, cNa, T_var(i), L/2);
    [z_perf_AT(i), z_no_e1_AT(i)]=weights_CG_quart(L, cNa, T_var(i), L/2-2);
end

p_perf_00=z_perf_00./(z_perf_00+Lo*z_no_e1_00); 
p_no_e1_00=Lo*z_no_e1_00./(z_perf_00+Lo*z_no_e1_00); 

p_perf_CG=z_perf_CG./(z_perf_CG+Lo*z_no_e1_CG); 
p_no_e1_CG=Lo*z_no_e1_CG./(z_perf_CG+Lo*z_no_e1_CG); 

p_perf_AT=z_perf_AT./(z_perf_AT+Lo*z_no_e1_AT); 
p_no_e1_AT=Lo*z_no_e1_AT./(z_perf_AT+Lo*z_no_e1_AT); 


figure (5)
plot(T_var(:), p_perf_00(:), 'LineWidth',3, 'Color', c_p);
hold on
plot(T_var(:), p_perf_CG(:), '--', 'LineWidth',2,'Color', c_p);
plot(T_var(:), p_perf_AT(:), ':', 'LineWidth',2,'Color', c_p);

text(40,0.35, sprintf('T = %i °C',T_pcr), 'FontSize',15);
val=T_pcr*ones(6,1);
valn=(0:5)/5;
plot(val(:), valn(:),':','Color', [0.2, 0.2, 0.2], 'LineWidth',1.5);

grid on
xlabel('T (°C)')
ylabel('\phi_{0}'' ')
xlim([20  80])
title('PCR Pairing Statistics')  
legend('n_{CG}=L/2','n_{CG}=L/2+2','n_{CG}=L/2-2','Location','southwest')
set(gca, 'FontSize', 15)



%% Fig. 6: L0 oligomer dependence

Lo_var=round((logspace(7, 11)));
Lo_var(35)=Lo;

[z_perf_Lo_CG, z_no_e1_Lo_CG]=weights_CG_quart(L, cNa, T_pcr, L/2+2);
[z_perf_Lo_00, z_no_e1_Lo_00]=weights_CG_quart(L, cNa, T_pcr, L/2);
[z_perf_Lo_AT, z_no_e1_Lo_AT]=weights_CG_quart(L, cNa, T_pcr, L/2-2);

p_perf_Lo_CG=z_perf_Lo_CG./(z_perf_Lo_CG+Lo_var.*z_no_e1_Lo_CG); 
p_no_e1_Lo_CG=Lo_var.*z_no_e1_Lo_CG./(z_perf_Lo_CG+Lo_var.*z_no_e1_Lo_CG);
p_perf_Lo_00=z_perf_Lo_00./(z_perf_Lo_00+Lo_var.*z_no_e1_Lo_00); 
p_no_e1_Lo_00=Lo_var.*z_no_e1_Lo_00./(z_perf_Lo_00+Lo_var.*z_no_e1_Lo_00); 
p_perf_Lo_AT=z_perf_Lo_AT./(z_perf_Lo_AT+Lo_var.*z_no_e1_Lo_AT); 
p_no_e1_Lo_AT=Lo_var.*z_no_e1_Lo_AT./(z_perf_Lo_AT+Lo_var.*z_no_e1_Lo_AT); 

figure (6)
plot(Lo_var(:), p_perf_Lo_00(:), '-', 'LineWidth',3, 'Color', c_p);
hold on
plot(Lo_var(:), p_perf_Lo_CG(:), '--', 'LineWidth',2,'Color', c_p);
plot(Lo_var(:), p_perf_Lo_AT(:), ':', 'LineWidth',2,'Color', c_p);

text(Lo/7,0.11, 'L_0 = 6 10^9 ', 'FontSize',15);
val=Lo*ones(6,1);
valn=(0:5)/5;
plot(val(:), valn(:),':','Color', [0.2, 0.2, 0.2], 'LineWidth',1.5);

grid on
xlabel('L_0')
ylabel('\phi_{0}'' ')
xlim([Lo_var(1) Lo_var(end)])
title('PCR Pairing Statistics')  
legend('n_{CG}=L/2','n_{CG}=L/2+2','n_{CG}=L/2-2','Location','southwest')
set(gca, 'Xscale', 'Log','FontSize', 15)%


%% SM Fig. 1: bounded primer per genome
massa_gen=10*10^-9; %10 ng of genome
vol=25*10^-6; %reaction volume  25 microL
c_primer=400*10^-9; % primer concentration 400 nM
MW_ssDNA=3*10^9*330; %molecular weight genomic ssDNA  in g/mol
c_gen=0.8e-16;

L_melt=[20, 21, 22];
T_melt=linspace(40,60 ,100);

z_prim_per=zeros(size(L_melt,2),size(T_melt,2));
z_prim_a=zeros(size(L_melt,2),size(T_melt,2));

for i=1:100
    for j=1:3
        [z_perf, z_prim_a(j,i), z_e1, z_no_e1]=weights(L_melt(j), cNa,T_melt(i));
    end
end


figure (1)
hold on
grid on
frac=0;
for i=1:size(L_melt,2)
    z_prim_per(i,:)= c_gen*(Lo*z_prim_a(i,:)'+exp(-(DH_f(L_melt(i),0,0,0)./(T_melt(:)+Tk)-DS_f(L_melt(i),0,0,0,cNa))/R));
    rho_gen(i,:)=2./(1+(1-2*frac)*z_prim_per(i,:)+sqrt(1+2*z_prim_per(i,:)+((1-2*frac)*z_prim_per(i,:)).^2));
    plot(T_melt(:), c_primer/c_gen*(1-rho_gen(i,:)),  'LineWidth',2);
end

xlabel('T(°C)')
ylabel('n_b ')
xlim([T_melt(1) T_melt(end)])
legend('L=20','L=21','L=22','Location','northwest')
set(gca,'FontSize', 15)%



T_melt=linspace(53,57 ,100);

z_prim_per=zeros(size(L_melt,2),size(T_melt,2));
z_prim_a=zeros(size(L_melt,2),size(T_melt,2));

for i=1:100
    for j=1:3
        [z_perf, z_prim_a(j,i), z_e1, z_no_e1]=weights(L_melt(j), cNa,T_melt(i));
    end
end



axes('Position',[.25 .23 .35 .30])
box on
grid on
frac=0;
for i=1:size(L_melt,2)
    z_prim_per(i,:)= c_gen*(Lo*z_prim_a(i,:)'+exp(-(DH_f(L_melt(i),0,0,0)./(T_melt(:)+Tk)-DS_f(L_melt(i),0,0,0,cNa))/R));
    rho_gen(i,:)=2./(1+(1-2*frac)*z_prim_per(i,:)+sqrt(1+2*z_prim_per(i,:)+((1-2*frac)*z_prim_per(i,:)).^2));
    plot(T_melt(:), c_primer/c_gen*(1-rho_gen(i,:)),  'LineWidth',2);
hold on
end

xlabel('T(°C)')
ylabel('n_b ')
xlim([T_melt(1) T_melt(end)])
legend('L=20','L=21','L=22','Location','northwest')
set(gca,'FontSize', 15)%


%% functions: DH, DS and weights


function [z_perf, z_med, z_e1, z_no_e1]=weights(L, cNa, T)

R=1.987; %gas constant in cal K^-1 mol^-1
Tk = 273.15; % kelvin-celsius

    z_perf = exp(-(DH_f(L,0,0,0)./(T+Tk)-DS_f(L,0,0,0,cNa))/R);
    z_rand = 0;
    z_e1 = 0;
    z_e1_plus1 = 0;
    z_no_e1 = 0;

    for ae1 = 0:(L-1)
        for ae2 = 0:(L-1-ae1)
            for ai = 0:max(0, L-2-ae1-ae2)
                z_rand=z_rand+3^(ae1+ae2+ai)*nchoosek(max(0,L-2-ae1-ae2),ai)...
                    *exp(-(DH_f(L, ae1, ae2, ai)./(T+Tk)-DS_f(L,ae1,ae2,ai,cNa))/R);
                if ae1>0
                    z_e1=z_e1+3^(ae1+ae2+ai)*nchoosek(max(0,L-2-ae1-ae2),ai)...
                    *exp(-(DH_f(L, ae1, ae2, ai)./(T+Tk)-DS_f(L,ae1,ae2,ai,cNa))/R);
                    if ae1>1
                        z_e1_plus1=z_e1_plus1+3^(ae1+ae2+ai)*nchoosek(max(0,L-2-ae1-ae2),ai)...
                            *exp(-(DH_f(L, ae1, ae2, ai)./(T+Tk)-DS_f(L,ae1,ae2,ai,cNa))/R);                   
                    end
                else
                    z_no_e1=z_no_e1+3^(ae1+ae2+ai)*nchoosek(max(0,L-2-ae1-ae2),ai)...
                    *exp(-(DH_f(L, ae1, ae2, ai)./(T+Tk)-DS_f(L,ae1,ae2,ai,cNa))/R);
                end
            end
        end
    end
    
    z_med=z_rand/4^L;
    z_no_e1=z_no_e1/4^L;
    z_e1=z_e1/4^L;

    
end

function [z_perf, z_no_e1]=weights_CG(L, cNa, T, n_CG)

R=1.987; %gas constant in cal K^-1 mol^-1
Tk = 273.15; % kelvin-celsius

    z_perf = exp(-(DH_f_cg(L,0,0,0, n_CG)./(T+Tk)-DS_f_cg(L,0,0,0,cNa, n_CG))/R);
    z_no_e1 = 0;

    ae1 = 0;
    for ae2 = 0:(L-1-ae1)
        for ai = 0:max(0, L-2-ae1-ae2)
            z_no_e1=z_no_e1+3^(ae1+ae2+ai)*nchoosek(max(0,L-2-ae1-ae2),ai)...
                *exp(-(DH_f_cg(L, ae1, ae2, ai, n_CG)./(T+Tk)-DS_f_cg(L,ae1,ae2,ai,cNa, n_CG))/R);
        end
    end


    z_no_e1=z_no_e1/4^L;

end

function [z_perf, z_no_e1]=weights_CG_quart(L, cNa, T, n_CG)

R=1.987; %gas costant in cal K^-1 mol^-1
Tk = 273.15; %kelvin-celsius

    z_perf = exp(-(DH_f_cg_quart(L,0,0,0, n_CG)./(T+Tk)-DS_f_cg_quart(L,0,0,0,cNa, n_CG))/R);
    z_no_e1 = 0;

    ae1 = 0;
    for ae2 = 0:(L-1-ae1)
        for ai = 0:max(0, L-2-ae1-ae2)
            z_no_e1=z_no_e1+3^(ae1+ae2+ai)*nchoosek(max(0,L-2-ae1-ae2),ai)...
                *exp(-(DH_f_cg_quart(L, ae1, ae2, ai, n_CG)./(T+Tk)-DS_f_cg_quart(L,ae1,ae2,ai,cNa, n_CG))/R);
        end
    end

    z_no_e1=z_no_e1/4^L;

end


%total enthaply and entropy: ATCG averaged
function DH = DH_f(L, ae1, ae2, ai)

%enthalpic contributions in kcal mol^-1 :
DH0 = 0.2;
DHq = -8.2375; %average
DHs = -2.534;%average

DHi = 0.154166667;

DHpen=2.2; %AT penalty

DH=1000*(DH0+0.5*2*DHpen+(L-1-ae1-ae2-ci_f(ai))*DHq+cs_f(ae1, ae2)*DHs+ci_f(ai)*DHi);
 
end

function DS = DS_f(L, ae1, ae2, ai, cNa)

%entropic contributions in cal K^-1 mol^-1 :
DS0 = -5.7;
DSq = -22.01875; %media
DSs = -6.749755621;
DSi = -0.84166667;

DSpen = 6.9; %penalty AT end

DS=DS0+0.5*2*DSpen+(L-1-ae1-ae2-ci_f(ai))*DSq+cs_f(ae1, ae2)*DSs+ci_f(ai)*DSi;
% DS=DS  + (L-abs(as)-ae1-ae2-ci_f(ai))*0.5*0.368*log(cNa);% SantaLucia
f_CG=0.5;
DS=DS+DH_f(L, ae1, ae2, ai)*((4.29*f_CG-3.95)*1e-5*log(cNa) +9.4*1e-6*(log(cNa))^2 ); %Owczarczy
end

%total enthaply and entropy: f_CG detailed

function DH_cg = DH_f_cg(L, ae1, ae2, ai, n_CG)

%enthalpic contributions in kcal mol^-1 :
DH0 = 0.2;
%DHq = -8.2375; %average
DHq_AT=-7.4; %just AT
DHq_CG=-9.1; %just CG
DHq_mean = (DHq_AT*(L-n_CG)+DHq_CG*n_CG)/L; 

%DHs = -2.534; %average
DHs_CG =-3.9875; %just AT
DHs_AT =-1.0813; %just CG
DHs_mean = (DHs_AT*(L-n_CG)+DHs_CG*n_CG)/L; 

DHi = 0.154166667;

DHpen = 2.2; %penalty AT end
DHpen_mean=DHpen*(L-n_CG)/L; 

DH_cg=1000*(DH0+DHpen_mean+(L-1--ae1-ae2-ci_f(ai))*DHq_mean+cs_f(ae1, ae2)*DHs_mean+ci_f(ai)*DHi);
 
end

function DS_cg = DS_f_cg(L, ae1, ae2, ai, cNa, n_CG)

%entropic contributions in cal K^-1 mol^-1 :
DS0 = -5.7;
%DSq = -22.01875; %average
DSq_AT=-21.075; %just AT
DSq_CG=-22.85;  %just CG
DSq_mean = (DSq_AT*(L-n_CG)+DSq_CG*n_CG)/L; 

%DSs = -6.749755621; %average
DSs_CG =-11.1498; %  just CG (from: DGs_CG =-0.5294)
DSs_AT =-2.7648; %   just AT (from:  DGs_AT = -0.2238)
DSs_mean=(DSs_AT*(L-n_CG)+DSs_CG*n_CG)/L; 

DSi = -0.84166667;

DSpen = 6.9; %penalty AT end
DSpen_mean=DSpen*(L-n_CG)/L; 

DS_cg=DS0+DSpen_mean+(L-1-ae1-ae2-ci_f(ai))*DSq_mean+cs_f(ae1, ae2)*DSs_mean+ci_f(ai)*DSi;
% DS_cg= DS_cg + (L-abs(as)-ae1-ae2-ci_f(ai))*0.368*log(cNa); % SantaLucia
DS_cg= DS_cg +DH_f_cg(L,ae1, ae2, ai, n_CG)*((4.29*n_CG/L-3.95)*1e-5*log(cNa) +9.4*1e-6*(log(cNa))^2 ); %Owczarczy
end


%total enthaply and entropy: f_CG detailed

function DH_cg = DH_f_cg_quart(L, ae1, ae2, ai, n_CG)

%enthalpic contributions in kcal mol^-1 :
DH0 = 0.2;
%DHq = -8.2375; %average
DHq_AT=-7.4; %just AT
DHq_CG=-9.1; %just CG
DHq_mean = (DHq_AT*(L-n_CG)+DHq_CG*n_CG)/L; 

DHs = -2.534;%average

DHi = 0.154166667;

DHpen = 2.2; %penalty AT end
DHpen_mean=DHpen*(L-n_CG)/L; 

DH_cg=1000*(DH0+DHpen_mean+(L-1--ae1-ae2-ci_f(ai))*DHq_mean+cs_f(ae1, ae2)*DHs+ci_f(ai)*DHi);
 
end

function DS_cg = DS_f_cg_quart(L, ae1, ae2, ai, cNa, n_CG)

%entropic contributions in cal K^-1 mol^-1 :
DS0 = -5.7;
%DSq = -22.01875; %average
DSq_AT=-21.075; %just AT
DSq_CG=-22.85;  %just CG
DSq_mean = (DSq_AT*(L-n_CG)+DSq_CG*n_CG)/L; 

DSs = -6.749755621; %average
DSi = -0.84166667;

DSpen = 6.9; %penalty AT end
DSpen_mean=DSpen*(L-n_CG)/L; 

DS_cg=DS0+DSpen_mean+(L-1-ae1-ae2-ci_f(ai))*DSq_mean+cs_f(ae1, ae2)*DSs+ci_f(ai)*DSi;
% DS_cg= DS_cg + (L-abs(as)-ae1-ae2-ci_f(ai))*0.368*log(cNa); % SantaLucia
DS_cg= DS_cg +DH_f_cg(L,ae1, ae2, ai, n_CG)*((4.29*n_CG/L-3.95)*1e-5*log(cNa) +9.4*1e-6*(log(cNa))^2 ); %Owczarczy
end


%number of shift/external penalities (0, 2, 3, 4)
function cs = cs_f(ae1, ae2)

cs=2;

if ae1>0
    cs=cs+1;
end

if ae2>0
    cs=cs+1;
end
   
              
end

%number of internal penalities 
function ci = ci_f(ai)

%ci=ai;% consecutive approx (0, 2, 3, 4, 5,...)

%if ai>0
%    ci=ci+1;
%end

ci=2*ai;% non-consecutivi pprox (0, 2, 4, 6,...)
 
end
