nm=zeros(10000,1);
for i=1:200
   close all;
   filename=['m' int2str(i-1) '.dat'];
   m=importdata(filename);
   m(m<0,:)=0.1;
   nm(1:2500)=m(1:2500);
   x=reshape(m(1:2500),50,50);
   
   imagesc(exp(x));
   colormap jet;
   colorbar;
   shading interp;

   nm(2501:5000)=m(1:2500);
   nm(5001:10000)=m;
   output=fopen(filename,'w');
   fprintf(output,'%d\n',nm);
   fclose(output);
   i
end