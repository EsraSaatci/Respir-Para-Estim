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

A=load('E3.mat');
AA=A.AA_RLS_stop;

w_length_all = [256 512 1024 2048];

wic=0;
for wi=2:2 %only for w_length=512
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
for ci=1:1  %only for RSS=1
    cic=cic+1;
    cna=names_channels{ci}; %Channel name for the channel ci
    data=datae.(dataee{1}).(cna); data=data(upcross2(iter2):upcross2(iter2+4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Lung Pressure Estimator

AAA=AA{wi,di,ci};

disp(dic)
u=AAA;u=u-mean(u);u=u/max(u);
y=dataf_down;y=y-mean(y);y=y/max(y);

N_kalman=length(y);
deltat=1/fs;
teta_e=2+randn(4,1);

Q = 1*[1 0;0 1]; % initial process noise covariance 
R = 0.1; % measurement error covariance
xhat=[0.1;0.1];
P = 10*eye(2);
poshata=[];
poshatb=[];
poshatc=[];
poshatd=[];
xEst=[];
Pdiciii=[0 0];

for iii=1:10

A_e=[-teta_e(1)*teta_e(3) -teta_e(3)*teta_e(2); -teta_e(4)*teta_e(2) -teta_e(4)*teta_e(2)];
A_d_e=eye(2)+deltat*A_e; 
B_e=[teta_e(3)*teta_e(2);teta_e(4)*teta_e(2)];
B_d_e=deltat*B_e;
C_e=[-teta_e(2) -teta_e(2)];
D_e=teta_e(2);
z_e=y-D_e.*u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kalman Part

F=A_d_e;
H=C_e;
saveP1=[]; saveP2=[];

for k = 1 : N_kalman
    Q=(1/k)*Q;
    % Kalman Filter Calculations   
    s = H * P * H' + R; 
    % Kalman Gain matrix
    K = P * H' * inv(s);
    % State estimate
    xhat_new = F * xhat + B_d_e * u(k);
    xhat = xhat_new + K * (z_e(k) - (H * xhat_new));
    % error covariance matrix
    P_new = F * P * F' + Q;
    P = P_new - K * H * P_new; 
    % savings
    saveP1(k)=P(1,1);
    saveP2(k)=P(2,2);
    xEst(:,k)=xhat;
end
Pdiciii=[Pdiciii;saveP1(k),saveP2(k)];
savexEstall{dic,iii}=xEst;
Q = (1/iii)*[1 0;0 1]; % new process noise covariance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optimization Part

for i=1:10

    teta0=[abs(1*randn);abs(1*randn);abs(1*randn);abs(1*randn)];
    lb=[0.01;0.01;0.01;0.01];
    ub=[100;100;100;100];
    objective = @(teta) objective_mat(xEst,u,y,teta);
    constraints = @(teta) constraints_mat(xEst,u,y,teta);
    options = optimoptions(@fmincon,'Display','off','Algorithm','interior-point',...
      'StepTolerance',1e-18,'ConstraintTolerance',1e-6,'OptimalityTolerance',1e-6);
    options.InitBarrierParam=1;
    [xxx,fval,exitflag,output,lambda,grad,hessian] = fmincon(objective,teta0,[],[],[],[],[],[],constraints,options);
    xout(:,i)=xxx;
    xxoutput{i}=output.message;
end
% meanxxx=mean(xout,2); 
meanxxx(:,iii)=median(xout,2);
savexoutall{1,iii}=xout;
teta_e=median(xout,2);
end

end

mxEst1=xEst(1,:)-mean(xEst(1,:)); mxEst1=mxEst1/max(mxEst1);
mxEst2=xEst(2,:)-mean(xEst(2,:)); mxEst2=mxEst2/max(mxEst2);

Rc(dic)=1/teta_e(2);
E2(dic)=teta_e(3);
Ec(dic)=teta_e(4);
R1(dic)=1/(teta_e(1)-teta_e(2));

Csize=size(mxEst1); Csize1=ceil(2*Csize/3);
A_R_stop=(u(Csize1:end)-(mxEst1(Csize1:end)+mxEst2(Csize1:end))')./max(y(Csize1:end)); %equation 5
AA_R_stop{dic}=A_R_stop;

end
end

save('E5.mat','AA_R_stop')
% save('E6.mat','AA_R_stop')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fig10 Estimated internal pressures

wic=2; %run only wi=2
dic=1; %run di=1:23
cic=1; %run only ci=1

figure,
[PAX,PH1,PH2]=plotyy(1:N_kalman,mxEst1,1:N_kalman,mxEst2); grid on 
xlabel('Data Points,i','FontSize',24)
ylabel(PAX(1),'${{P_{2,i}}}$','FontSize',24,'Interpreter','latex')
ylabel(PAX(2),'${{P_{1,i}}}$','FontSize',24,'Interpreter','latex')
set(PAX(1),'YLim',[-2.5 2.5],'ytick',[-2:1:2],'XLim',[0 N_kalman],'xtick',[0:50:N_kalman],'XTickLabel',[0:50:N_kalman])
set(PAX(2),'YLim',[-2.5 2.5],'ytick',[-2:1:2],'YTickLabel',[-2:1:2],'XLim',[0 N_kalman],'xtick',[0:50:N_kalman])
set(gca,'FontSize',24)
set(PAX,'FontSize',24),set(PH1,'LineWidth',2),set(PH2,'LineWidth',2) 


