load FD.dat;
load adjoint.dat;
load err.dat;
hold on;
filename='0-4';
% errorpath='./error/'
% gradientpath='gradient/'
fig=plot(FD);
fig=plot(adjoint);
legend('Finite difference gradient','Adjoint gradient');
saveas(fig,filename,'jpg');
figure;
Y=zeros(size(adjoint));
filename1=strcat('e',filename);
fig_n=errorbar(Y,err);
saveas(fig_n,filename1,'jpg');
close all;