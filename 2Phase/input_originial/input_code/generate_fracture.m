function [ a,b ] = generate_fracture(NX,NY,Lx,Ly,frac )
mu=frac(2,1);
phi=frac(2,2);
n_e=frac(1,1);
n_d=frac(1,2);
DX=Lx/NX;
DY=Ly/NY;

n=n_e+int8(n_d*randn(1));

for nm=1:frac(7,1)
    
close all;
%%stochastic input fracture geometry
theta_e=frac(3,1);
theta_d=frac(3,2);

theta=theta_e+theta_d*randn(n,1);
length=mu+phi*randn(n,1);
%%fracture property
b=zeros(n,3);
b(:,1)=frac(4,1)+frac(4,2)*randn(n,1);
b(:,2)=frac(5,1)+frac(5,2)*randn(n,1);
b(:,3)=frac(6,1)+frac(6,2)*randn(n,1);
a=rand(n,2);
x=Lx*0.1+a(:,1)*Lx*0.8;
y=a(:,2)*Ly;



for i=1:n
   a(i,1)=x(i)-length(i)*cos(theta(i));
   a(i,2)=y(i)-length(i)*sin(theta(i));
   a(i,3)=x(i)+length(i)*cos(theta(i));
   a(i,4)=y(i)+length(i)*sin(theta(i));  
   if a(i,1)<0
       a(i,1)=0.1;
   end
   if a(i,1)>Lx
       a(i,1)=Lx-0.1;
   end
   if a(i,3)<0
       a(i,3)=0.1;
   end
   if a(i,3)>Lx
       a(i,3)=Lx-0.1;
   end
      if a(i,2)<0
       a(i,2)=0.1;
   end
   if a(i,2)>Ly
       a(i,2)=Ly-0.1;
   end
      if a(i,4)<0
       a(i,4)=0.1;
   end
   if a(i,4)>Ly
       a(i,4)=Ly-0.1;
   end
end
for i=1:n
    fig=line([a(i,1),a(i,3)],[a(i,2),a(i,4)],'LineWidth',2);
    hold on;
end
axis([0,Lx,0,Ly]);
pattern_name=['../f_pattern/frac_pattern', num2str(nm)];
saveas(fig,pattern_name,'jpg');
fracsystem= ['../f_model/fracturesystem',num2str(nm),'.dat'];
fracsystem=fopen(fracsystem,'w');
fprintf(fracsystem,'%d\n',n);
fprintf(fracsystem,'%f \t%f \t%f \t%f \t\n',a')
fprintf(fracsystem,'%f \t%f \t%f \t\n',b');
fclose(fracsystem);

end
close all;
end

