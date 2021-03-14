clc
% close all
clear all

% channels: ai21,ai29,ai20,ai28,ai19,ai27,ai18 -> lung sounds 1:7 (LS7 is
% tracheal sound)
% channels: ai17, ai26 -> temperature and humidity -> 8,9
% channels: ai0,ai1,ai2 -> flow, pressure, volume -> 10,11,12

fs=8000; % Sampling frequency

%% Which data
upcross1 = load('upcross.mat');
downcross1 = load('downcross.mat');
namedata1='C:\MATLAB\Olcumler\esra0\';
% namedata1='C:\Users\esra.saatci\Documents\Olcumler\';
% namedata1='C:\Users\Esra\Documents\Olcumler\';
names_data={'aziz1.mat','aziz2.mat','dogancan1.mat','Dogancan2.mat',...
    'Dogancan3.mat','duygu2.mat','duygu1.mat','feray2.mat','feray1.mat',...
    'ferhan1.mat','ferhan2.mat','ghazal3.mat','ghazal2.mat',...
    'kenan1.mat','kenan2.mat','melek2.mat','melek1.mat',...
    'noushin1.mat','noushin2.mat','onur2.mat','onur1.mat','umut2.mat',...
    'umut1.mat'};
names_channels={'ai21','ai29','ai20','ai28','ai19','ai27','ai18',...
    'ai17','ai26','ai0','ai1','ai2'};

w_length_all = [256 512 1024 2048];

wic=0;
for wi=1:4 %4 window length
    wic=wic+1;
    w_length=w_length_all(wi);
    
    dic=0;
for di=1:23 %23 recordings 
    dic=dic+1; % use dic for the indexing only 
      
%determination of the correct periods
upcross2=upcross1.s_upcrosss{di};
downcross2=downcross1.s_downcrosss{di};

namedata2=names_data{di}; %path name for the recording di
namedata=strcat(namedata1,namedata2); 
iter1=insexpsep(namedata2,upcross2,downcross2); iter2=iter1(1); %first 3 period
datae=load(namedata); %data for recording di 
dataee=fields(datae);
data=datae.(dataee{1}).ai0; dataf=data(upcross2(iter2):upcross2(iter2+4)); %flow signal
data=datae.(dataee{1}).ai1; datap=data(upcross2(iter2):upcross2(iter2+4)); %pressure signal
timee=datae.(dataee{1}).Time;    % time axis

N=length(data);

%downsample pressure and flow
datap_down=downsample(datap,w_length); dataf_down=downsample(dataf,w_length); 
datap_down=datap_down-mean(datap_down);dataf_down=dataf_down-mean(dataf_down);
datap_down=datap_down/max(datap_down); dataf_down=dataf_down/max(dataf_down);

cic=0;
for ci=1:7  %run only for RSS
    cic=cic+1;
    cna=names_channels{ci}; %Channel name for the channel ci
    data=datae.(dataee{1}).(cna); data=data(upcross2(iter2):upcross2(iter2+4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mouth Pressure Estimator

A=load('E1.mat');
AA=A.HFD_framestop;
AAA=AA{wic,dic,cic};

HFD_frames = AAA'-mean(AAA); % zero mean needed for RLS
N_frames=length(HFD_frames);
xf1_af=HFD_frames; xf1_afn=xf1_af/max(xf1_af); %reference signal
x1_af=datap_down; x1_afn=x1_af/max(x1_af);     %desired signal

M1_RLS=20;                                     %Filter order 
lamda=0.99; lamda1=1/lamda;
P_RLS=10^10*eye(M1_RLS);
wa_RLS=0.1*ones(M1_RLS,1);
wt_RLS=zeros(M1_RLS,1);
u_RLS=[]; PM_RLS=[];

for i_RLS=0:N_frames-M1_RLS
    xf1u=xf1_afn(M1_RLS+i_RLS:-1:1+i_RLS);
    u_RLS=[u_RLS xf1u];
end
for i_RLS=1:N_frames-M1_RLS
    k_RLS=(P_RLS*lamda1*u_RLS(:,i_RLS))/(1+(lamda1*u_RLS(:,i_RLS)'*P_RLS*u_RLS(:,i_RLS)));
%     y_RLS(i_RLS)=wa_RLS'*u_RLS(:,i_RLS);
    y_RLS_R=filter(wa_RLS',1,u_RLS(:,i_RLS));
    y_RLS(i_RLS)=y_RLS_R(1);
    e_RLS(i_RLS)=x1_afn(i_RLS)-y_RLS(i_RLS);
    wa_RLS=wa_RLS+k_RLS*e_RLS(i_RLS);
    P_RLS=lamda1*P_RLS-lamda1*k_RLS*u_RLS(:,i_RLS)'*P_RLS;
    wt_RLS=[wt_RLS wa_RLS];
    PM_RLS(:,:,i_RLS)=P_RLS;
end

AA_RLS=filter(wa_RLS,1,xf1_af);
AA_RLS=AA_RLS/max(AA_RLS); 
AA_RLS_stop{wic,dic,cic}=AA_RLS;

MSE1(cic)=sum((AA_RLS(2:end)-x1_afn(2:end)).^2)/length(x1_afn(2:end));
end
MSE2(dic,:)=MSE1;
end
MSE3(wic,:,:)=MSE2;
end


save('E3.mat','AA_RLS_stop') 
% save('E4.mat','MSE3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig8 Estimated mouth pressure and measured mouth pressure

wic=2; %run only wi=2
dic=1; %run di=1:23
cic=1; %run only ci=1

figure,
[PAX,PH1,PH2]=plotyy(1:N_frames,AA_RLS,1:N_frames,datap_down); grid on 
xlabel('Data Points,i','FontSize',24)
ylabel(PAX(1),'${{{\hat P}_{ao,i}}}$','FontSize',24,'Interpreter','latex')
ylabel(PAX(2),'${{P_{ao,i}}}$','FontSize',24,'Interpreter','latex')
set(PAX(1),'YLim',[-1 1],'ytick',[-1:0.5:1],'XLim',[0 N_frames],'xtick',[0:50:N_frames],'XTickLabel',[0:50:N_frames])
set(PAX(2),'YLim',[-1 1],'ytick',[-1:0.5:1],'YTickLabel',[-1 -0.5 0 0.5 1],'XLim',[0 N_frames],'xtick',[0:50:N_frames])
set(gca,'FontSize',24)
set(PAX,'FontSize',24),set(PH1,'LineWidth',2),set(PH2,'LineWidth',2) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig3 The number of Min MSE

MSE4=shiftdim(MSE3,1);
for dic=1:23
HDD1=squeeze(MSE4(dic,:,:));
[v,loc] = min(HDD1(:));
[ii,jj] = ind2sub(size(HDD1),loc);
iii(dic)=ii; % for each subject min hata at iii channel
jjj(dic)=jj; % for each subject min hata at jjj window
end

N1 = histcounts(iii,7); %number of min error et each channel
N2 = histcounts(jjj,4); %number of min error et each window

figure, histogram(iii),grid on 
set(gca,'xtick',[1:7],'xticklabel',['LS1';'LS2';'LS3';'LS4';'LS5';'LS6';'LS7'])
xlabel('Respiratory Sound Signals','FontSize',24)
ylabel('The number of MSE','FontSize',24)
set(gca,'FontSize',24)
figure, histogram(jjj),grid on 
set(gca,'xtick',[1:4],'xticklabel',['W1';'W2';'W3';'W4'])
xlabel('Windows','FontSize',24)
ylabel('The number of MSE','FontSize',24)
set(gca,'FontSize',24)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig9 Bar graphs of MSE values

wic=2; %run only wi=2
dic=1; %run di=1:23
cic=1; %run only ci=1

MSE5=squeeze(MSE3(wic,:,cic)');
MSE5_quite=[MSE5(1);MSE5(3);MSE5(6:2:end)];
MSE5_force=[MSE5(2);MSE5(4);MSE5(5:2:end)];
figure, bar([1;3;[6:2:23]'],MSE5_quite,0.2,'k'),hold on
bar([2;4;[5:2:23]'],MSE5_force,0.4,'c')
set(gca,'FontSize',24)
xlabel('Recordings','FontSize',24)
ylabel('MSE','FontSize',24)