
function XXX=HEE(Data_all,w_length)

% Hurst exponent (H)        Conversion                  Adjustment of Hq and Ht
% <0.2                      signal=cumsum(X-mean(X))    ?1
% 0.2–0.8                   Noconversion                0
% 0.8–1.2                   Noconversion                0
% 1.2–1.8                   signal=diff(X)              +1
% >1.8                      signal=diff(diff(X))        +2

N_data=length(Data_all);

y = buffer(Data_all,w_length);
[A numFrames]=size(y);

for i=1:numFrames
    Data=y(:,i);
    X=cumsum(Data-mean(Data));
    % X=Data;%-mean(Data);
    X=transpose(X);
    N=length(X);

scmin=16;   % min scale
scres=30;   % number of scales
scmax=floor(N/4);  %scmin+scmin*scres; % max scale
exponents=linspace(log2(scmin),log2(scmax),scres);
scale=round(2.^exponents);
% scale=[16,32,64,128,256,512,1024];  % multiple segment sizes (i.e., scales)
m=1;
for ns=1:length(scale)
    segments(ns)=floor(length(X)/scale(ns));
    for v=1:segments(ns)
        Idx_start=((v-1)*scale(ns))+1;
        Idx_stop=v*scale(ns);
        Index{v,ns}=Idx_start:Idx_stop;
        X_Idx=X(Index{v,ns});
        C=polyfit(Index{v,ns},X(Index{v,ns}),m);
        fit{v,ns}=polyval(C,Index{v,ns});
        RMS{ns}(v)=sqrt(mean((X_Idx-fit{v,ns}).^2));    % local fluctuations
        RMS_scale2{ns}(Index{v,ns})=RMS{ns}(v).*ones(size(Index{v,ns}));
    end
    F(ns)=sqrt(mean(RMS{ns}.^2));   % the scaling function 
end

C=polyfit(log2(scale),log2(F),1);
H=C(1);
RegLine=polyval(C,log2(scale));
YMatrix1=[log2(F); RegLine];
X1=log2(scale);
XXX(i)=H;
end
