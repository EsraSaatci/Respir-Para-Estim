clc
% close all
clear all

% channels: ai21,ai29,ai20,ai28,ai19,ai27,ai18 -> lung sounds 1:7 (LS7 is
% tracheal sound)
% channels: ai17, ai26 -> temperature and humidity -> 8,9
% channels: ai0,ai1,ai2 -> flow, pressure, volume -> 10,11,12

fs=8000; % Sampling frequency

FEV1=[4.98;4.98;4.28;4.28;4.28;3.92;3.92;3.26;3.29;4.45;4.45;3.33;...
    3.33;3.51;3.51;4.82;4.82;4.58;4.58]; 
FVC=[5.73;5.73;5.69;5.69;5.69;4.24;4.24;4.55;4.55;5.11;5.11;4.43;...
    4.43;3.67;3.67;5.05;5.05;5.47;5.47];
PEF=[10.18;10.18;12.35;12.35;12.35;7.76;7.76;6.87;6.87;9.91;9.91;8.84;...
    8.84;7.02;7.02;10.58;10.58;7.50;7.50];

RV=[1.96;1.96;1.23;1.23;1.23;1.34;1.34;1.87;1.87;1.47;1.47;1.27;...
    1.27;1.14;1.14;1.45;1.45;0.57;0.57];
TLC=[7.51;7.51;6.57;6.57;6.57;5.34;5.34;6.05;6.05;6.20;6.20;5.29;...
    5.29;4.56;4.56;6.02;6.02;5.81;5.81];
sRaw=[0.75;0.75;0.88;0.88;0.88;0.81;0.81;1.52;1.52;0.59;0.59;1.32;...
    1.32;0.89;0.89;0.25;0.25;0.69;0.69];
Raw=[0.17;0.17;0.23;0.23;0.23;0.25;0.25;0.40;0.40;0.12;0.12;0.35;...
    0.35;0.39;0.39;0.13;0.13;0.23;0.23];

R5=[0.25;0.25;0.53;0.53;0.53;0.40;0.40;0.66;0.66;0.35;0.35;0.63;...
    0.63;0.43;0.43;0.28;0.28;0.37;0.37];
X5=[-0.09;-0.09;-0.22;-0.22;-0.22;-0.09;-0.09;-0.32;-0.32;-0.06;-0.06;-0.26;...
    -0.26;-0.09;-0.09;-0.07;-0.07;-0.14;-0.14];
R35=[0.21;0.21;0.40;0.40;0.40;0.39;0.39;0.34;0.34;0.43;0.43;0.36;...
    0.36;0.37;0.37;0.24;0.24;0.35;0.35];
X35=[0.10;0.10;0.08;0.08;0.08;0.16;0.16;0.00;0.00;0.27;0.27;0.06;...
    0.06;0.10;0.10;0.25;0.25;0.16;0.16];
AX=[0.32;0.32;1.45;1.45;1.45;0.37;0.37;4.32;4.32;0.09;0.09;2.47;...
    2.47;0.56;0.56;0.19;0.19;0.46;0.46];

PARA=[FEV1 FVC PEF/10 RV TLC sRaw 10*Raw 10*R5 -10*X5 10*R35 10*X35 AX];
titlenames={'FEV1';'FVC';'PEF';'RV';'TLC';'sRaw';'Raw';'R5';'-X5';'R35';'X35';'AX'};
PARAread=PARA;

A=load('E5.mat');
AA=A.AA_R_stop;

g=[];Rcbox=[];
dic=0;
for di=1:23 %23 recordings 
    dic=dic+1; % use dic for the indexing only 

Rc=AA{di}; Rc=Rc(:);
Rcbox=[Rcbox;Rc];
Rc_fea1_(dic)=mean(abs(Rc));  %feature1 not used
Rc_fea2_(dic)=var(Rc); %feature2
Rc_fea3_(dic)=skewness(Rc); %feature3
g=[g;dic*ones(length(Rc), 1)];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig11 Features

Rc_fea2_quite=[Rc_fea2_(1),Rc_fea2_(3),Rc_fea2_(6:2:end)];
Rc_fea2_force=[Rc_fea2_(2),Rc_fea2_(4),Rc_fea2_(5:2:end)];
figure, bar([1;3;[6:2:23]'],Rc_fea2_quite,0.2,'k'),hold on
bar([2;4;[5:2:23]'],Rc_fea2_force,0.4,'c')
set(gca,'FontSize',24)
xlabel('Recordings','FontSize',24)
ylabel('$G_{RC}$','Interpreter','latex','FontSize',24)

Rc_fea3_quite=[Rc_fea3_(1),Rc_fea3_(3),Rc_fea3_(6:2:end)];
Rc_fea3_force=[Rc_fea3_(2),Rc_fea3_(4),Rc_fea3_(5:2:end)];
figure, bar([1;3;[6:2:23]'],Rc_fea3_quite,0.2,'k'),hold on
bar([2;4;[5:2:23]'],Rc_fea3_force,0.4,'c')
set(gca,'FontSize',24)
xlabel('Recordings','FontSize',24)
ylabel('$S_{RC}$','Interpreter','latex','FontSize',24)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%exclude recordings 12,13,16 and 17

index=[1:11,14:15,18:23];
Rc_fea2=Rc_fea2_(index);
Rc_fea3=Rc_fea3_(index);
N=length(Rc_fea2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correlation coefficient part

for ii=1:12
 [Cor_Rc_fea2_ P2]=corrcoef(Rc_fea2,PARAread(:,ii));
 [Cor_Rc_fea3_ P3]=corrcoef(Rc_fea3,PARAread(:,ii));
 Cor_Rc_fea2(ii)=Cor_Rc_fea2_(1,2);
 Cor_Rc_fea3(ii)=Cor_Rc_fea3_(1,2);
 Corr_Rc_fea2(ii)=P2(1,2);
 Corr_Rc_fea3(ii)=P3(1,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig12 Correlation coefficients

figure, stem(Cor_Rc_fea2,'LineWidth',2),grid on
set(gca,'xtick',[0:12],'xticklabel',['0   ';'FEV1';'FVC ';'PEF ';'RV  ';'TLC ';'sRaw';'Raw ';'R5  ';'X5  ';'R35 ';'X35 ';'AX  '])
ylabel('${\hat \rho _{G_{RC},RP}}$','Interpreter','latex','FontSize',24)
figure, stem(Cor_Rc_fea3,'LineWidth',2),grid on
set(gca,'xtick',[0:12],'xticklabel',['0   ';'FEV1';'FVC ';'PEF ';'RV  ';'TLC ';'sRaw';'Raw ';'R5  ';'X5  ';'R35 ';'X35 ';'AX  '])
ylabel('${\hat \rho _{S_{RC},RP}}$','Interpreter','latex','FontSize',24)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Blant-Altman part

Bland_Rc_fea2=Rc_fea2'; % feature power
% Bland_Rc_fea2=Rc_fea3'; % feature skewness

for ii=1:12
BAdif=Bland_Rc_fea2-PARAread(:,ii); 
BAdifmean=mean(BAdif); BAdifmeanv=BAdifmean*ones(N,1);
BAdifstd=std(BAdif); BAdifstdv1=BAdifstd*ones(N,1); 
BAdifstdv2=BAdifmeanv-1.96*BAdifstdv1; BAdifstdv3=BAdifmeanv+1.96*BAdifstdv1;


BAmean=(Bland_Rc_fea2+PARAread(:,ii))/2;

BAmean_quite=[BAmean(1);BAmean(3);BAmean(6);BAmean(8);BAmean(10);BAmean(12);BAmean(14);BAmean(16);BAmean(18);];
BAmean_force=[BAmean(2);BAmean(4);BAmean(5);BAmean(7);BAmean(9);BAmean(11);BAmean(13);BAmean(15);BAmean(17);BAmean(19);];

BAdif_quite=[BAdif(1);BAdif(3);BAdif(6);BAdif(8);BAdif(10);BAdif(12);BAdif(14);BAdif(16);BAdif(18);];
BAdif_force=[BAdif(2);BAdif(4);BAdif(5);BAdif(7);BAdif(9);BAdif(11);BAdif(13);BAdif(15);BAdif(17);BAdif(19);];


BAxa=linspace(min(BAmean),max(BAmean),N);
figure, scatter(BAmean_quite,BAdif_quite,'k','filled','LineWidth',4),hold on  
scatter(BAmean_force,BAdif_force,'c','filled','LineWidth',4)
plot(BAxa,BAdifmeanv,'g','LineWidth',4),plot(BAxa,BAdifstdv2,'m','LineWidth',4),plot(BAxa,BAdifstdv3,'b','LineWidth',4), grid on, 
set(gca,'FontSize',24)
xlabel('$M_j$','Interpreter','latex','FontSize',24)
ylabel('$DF_j$','Interpreter','latex','FontSize',24)
title(titlenames(ii),'FontSize',24)
axis([0.5,6.5,-6.5,8.0]) % feature power
% axis([0.2,4.8,-7,3.8]) % feature skewness


[LRRR S]=polyfit(BAmean,BAdif,1); saveR(ii)=S.normr;
FFTRR_est=polyval(LRRR,BAmean);
plot(BAmean,FFTRR_est,'--r','LineWidth',4)

end



