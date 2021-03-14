function iter1=insexpsep(namedata2,upcross,downcross)

% determine the correct cycles
if strcmp(namedata2, 'aziz2.mat')
    iter1=[1:5,10:length(upcross)];
elseif strcmp(namedata2, 'dogancan1.mat')
    iter1=4:length(upcross);
elseif strcmp(namedata2, 'Dogancan2.mat')
    iter1=[1:15,17:length(upcross)];
elseif strcmp(namedata2, 'feray1.mat')
    iter1=10:length(upcross);
elseif strcmp(namedata2, 'ferhan1.mat')
    iter1=[1,3:13,16:17,20,24:27,29:length(upcross)];
elseif strcmp(namedata2, 'ferhan2.mat')
    iter1=[1:32,36:length(upcross)];
elseif strcmp(namedata2, 'ghazal1.mat')
    iter1=0; disp('dont use ghazal1')
elseif strcmp(namedata2, 'ghazal2.mat')
    iter1=[2:4,9:length(upcross)];
elseif strcmp(namedata2, 'ghazal3.mat')
    iter1=[1:4,8:12,14:15,17:length(upcross)];
elseif strcmp(namedata2, 'ghazal4.mat')
    iter1=0; disp('dont use ghazal4')
elseif strcmp(namedata2, 'kenan1.mat')
    iter1=[1:2,6:length(upcross)];
elseif strcmp(namedata2, 'melek1.mat')
    iter1=[4:5,7,9:length(upcross)];
elseif strcmp(namedata2, 'melek2.mat')
    iter1=[9:10,15:length(upcross)];
elseif strcmp(namedata2, 'noushin1.mat')
    iter1=[1:15,17:19,21:length(upcross)];
elseif strcmp(namedata2, 'onur2.mat')
    iter1=4:length(upcross);
elseif strcmp(namedata2, 'umut2.mat')
    iter1=[1:9,14:length(upcross)];
elseif strcmp(namedata2, 'umut1.mat')
    iter1=[1:3,5,9:12,14,17:length(upcross)];
else
iter1=1:length(upcross);
end

