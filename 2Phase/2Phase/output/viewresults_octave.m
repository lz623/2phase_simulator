clear all 
close all
NX = 220;
NY = 60;

load results.out
sw = reshape( results(:,2),NX,NY);
po = reshape( results(:,1),NX,NY);
clear results

figure
subplot(1,2,1)
imagesc(sw)
title('Sw')
subplot(1,2,2)
imagesc(po)
title('Po')
