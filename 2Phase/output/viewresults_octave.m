clear all 
close all
NX = 50;
NY = 50;
Np=40;
load results.out
%load results_with_frac.out
%results=results_with_frac;
% 
sw = reshape( results(1:NX*NY,2),NX,NY);
po = reshape( results(1:NX*NY,1),NX,NY);
sw=sw'
po=po'
clear results

fig=imagesc(po);
title('Po');
shading interp,colorbar,colormap jet;
saveas(fig,'visulize_output/Po_profile','jpg');
figure
fig=imagesc(sw);
title('Sw');
shading interp,colorbar,colormap jet;
h=colorbar;
set(h, 'ylim', [0.2 1])
saveas(fig,'visulize_output/Sw_profile','jpg');
load well_oil_rate.dat;
load well_water_rate.dat;

oil=well_oil_rate;
water=well_water_rate;
wor=water./oil;

for i=1:Np
    fig=subplot(5,8,i),[hAx,hLine1,hLine2]=plotyy(oil(:,1),oil(:,i+1),oil(:,1),wor(:,i+1));
    str=['Oil production rate and WOR for well' num2str(i)];
    title(str);
    xlabel('Time (day)');
    ylabel(hAx(1),'Oil production rate(std/d)') % left y-axis
    ylabel(hAx(2),'Water oil ratio') % right y-axis
end

saveas(fig,'visulize_output/qo_wor','fig');
