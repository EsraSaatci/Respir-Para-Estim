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
    disp(wi)
    
    dic=0;
for di=1:23 %23 recordings 
    dic=dic+1; % use dic for the indexing only 
    disp(di)

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
for ci=1:12  %12 channels
    cic=cic+1;
    cna=names_channels{ci}; %Channel name for the channel ci
    data=datae.(dataee{1}).(cna); data=data(upcross2(iter2):upcross2(iter2+4));
    disp(ci)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic statistics of data in ci'th channel
Meanaa=mean(data); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hurst Exponent Estimator

data1=data-Meanaa; data1=data1/max(data1);

Ndata1(dic,cic)=length(data1); N_HDD=length(data1);

if ci==8||ci==9||ci==10||ci==11||ci==12 %see Ref15,17
    data1_dif=diff(data1); % change to fractional Gauss see Ref15,17
    HFD_frames=HEE(data1_dif,w_length); %windowing
    HFD_frames=HFD_frames+1; %correction see Ref15,17
%     HFD_frames_nw=HEE(data1,N_HDD); %no windowing
%     HFD_frames_nw=HFD_frames_nw+1;
%     disp(ci)
else
    HFD_frames=HEE(data1,w_length); %windowing
%     HFD_frames_nw=HEE(data1,N_HDD); %no windowing

N_frames=length(HFD_frames);
end

HFD_framestops{wic,dic,cic}=HFD_frames;
% HFD_framestop_nw{wic,dic,cic}=HFD_frames_nw;%-mean(HFD_frames_nw);
% N_framestop{wic,dic,cic}=N_frames;

end
end
end
save('E1.mat','HFD_framestops') 
% save('E2.mat','HFD_framestop_nw')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig6 Hurst exponent time series and airflow signal

wic=2; %run only wi=2
dic=1; %run di=1:23
cic=1; %run only ci=1

figure,
[PAX,PH1,PH2]=plotyy(1:N_frames,HFD_frames,1:N_frames,dataf_down); grid on 
xlabel('Data Points,i','FontSize',24)
ylabel(PAX(1),'$H_{i}$','FontSize',24,'Interpreter','latex')
ylabel(PAX(2),'$Q_{i}$','FontSize',24,'Interpreter','latex')
set(PAX(1),'YLim',[0 2],'ytick',[0:1:2],'XLim',[0 N_frames],'xtick',[0:50:N_frames],'XTickLabel',[0:50:N_frames])
set(PAX(2),'YLim',[-1.5 1.5],'ytick',[-1.5 0 1.5],'YTickLabel',[-1.5 0 1.5],'XLim',[0 N_frames],'xtick',[0:50:N_frames])
set(gca,'FontSize',24)
set(PAX,'FontSize',24),set(PH1,'LineWidth',2),set(PH2,'LineWidth',2) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig7 Box plot of Hurst exponents

HFD_frames_save=[];g=[];

wic=2; %run only wi=2
cic=1; %run only ci=1

for dic=1:23
HFD_frames1=HFD_framestop{wic,dic,cic};
HFD_frames1=HFD_frames1(:);
HFD_frames_save=[HFD_frames_save;HFD_frames1];
g = [g;dic*ones(length(HFD_frames1), 1)];
end
figure,
boxplot(HFD_frames_save, g),grid on 
ylim([0 2])
xlabel('Recordings','FontSize',24)
ylabel('Hurst Exponents','FontSize',24)
set(gca,'FontSize',24)
