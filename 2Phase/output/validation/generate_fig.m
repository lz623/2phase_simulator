filename='0-4'
% filename=to_string(filename);
fig=openfig(filename);
legend('Finite difference gradient','Adjoint gradient');
saveas(fig,filename,'jpg');