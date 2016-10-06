function [ ] = generate_well(well,schedule,nw)

wells= ['../horizontal_well.dat'];
wells=fopen(wells,'w');
fprintf(wells,'%d\n',nw);
for i=1:nw
    fprintf(wells,'%d\n',well(i,1));
    fprintf(wells,'%d\n',well(i,2));
    fprintf(wells,'%d\n',well(i,3));
    fprintf(wells,'%f \t%f \t%f \t%f \t\n',well(i,4:7)');
    fprintf(wells,'\n\n');
end
fclose(wells);
pro_t=['../schedule.dat'];
pro_t=fopen(pro_t,'w');
fprintf(pro_t,'%d\n',schedule(:,1));
fclose(pro_t);

schedules= ['../hwellschedule.dat'];
schedules=fopen(schedules,'w');
nc=size(schedule(:,i+1));
for i=1:nw
    %schedule(:,i+1)=schedule(:,i+1)+0.1*schedule(:,i+1).*randn(nc);
    fprintf(schedules,'%d\n',schedule(:,i+1));
    fprintf(schedules,'\n\n');
end
fclose(schedules);
end

