function [f, feq] = constraints_mat (x,u,y,teta)
% teta=teta_correct;
fs = 8000; deltat=1/fs;
A=[-teta(1)*teta(3) -teta(3)*teta(2); -teta(4)*teta(2) -teta(4)*teta(2)];
A_d=eye(2)+deltat*A; 
B=[teta(3)*teta(2);teta(4)*teta(2)];
B_d=deltat*B;
C=[-teta(2) -teta(2)];
D=teta(2);

N=length(u);
f1=[];

for i=1:N-1
%     f1=[f1;x(1,i+1)-A_d(1,1)*x(1,i)-A_d(1,2)*x(2,i)-B_d(1)*u(i)];
    f1=f1+(x(1,i+1)-A_d(1,1)*x(1,i)-A_d(1,2)*x(2,i)-B_d(1)*u(i))^2;
end
for i=1:N-1
%     f1=[f1;x(2,i+1)-A_d(2,1)*x(1,i)-A_d(2,2)*x(2,i)-B_d(2)*u(i)];
    f1=f1+(x(2,i+1)-A_d(2,1)*x(1,i)-A_d(2,2)*x(2,i)-B_d(2)*u(i))^2;
end
% for i=1:N-1
% %     f1=[f1;y(i)-C(1)*x(1,i)-C(2)*x(2,i)-D(1)*u(i)-teta(1)*0.001-teta(3)*0.001-teta(4)*0.001];
% f1=[f1;y(i)-C(1)*x(1,i)-C(2)*x(2,i)-D(1)*u(i)];
% end

feq=f1; f=[];