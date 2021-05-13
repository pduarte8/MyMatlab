Prefix = 'ocean_his_';
for i =1:121
   if i < 10
    fnumber = strcat('000',num2str(i));
   end
   if i >=10 & i <100
    fnumber = strcat('00',num2str(i));
   end
   if i >= 100
    fnumber = strcat('0',num2str(i));
   end
   fname = strcat(Prefix,fnumber,'.nc');
   ocean_time = ncread(fname,'ocean_time');
   [m,n] = size(ocean_time);
   [Y1 M1 D1 h1 m1 s1] = unixsecs2date(ocean_time(1));
   strY1 = num2str(Y1);
   strM1 = num2str(M1);
   if M1 < 10
     strM1 = strcat('0',num2str(M1));
   end
   strD1 = num2str(D1);
   if D1 < 10
     strD1 = strcat('0',num2str(D1));
   end
   strh1 = num2str(h1);
   if h1 < 10
     strh1 = strcat('0',num2str(h1));
   end
   t1 = strcat(strY1,strM1,strD1,strh1);
   [Y2 M2 D2 h2 m2 s2] = unixsecs2date(ocean_time(m));
   strY2 = num2str(Y2);
   strM2 = num2str(M2);
   if M2 < 10
     strM2 = strcat('0',num2str(M2));
   end
   strD2 = num2str(D2);
   if D2 < 10
     strD2 = strcat('0',num2str(D2));
   end
   strh2 = num2str(h2);
   if h2 < 10
     strh2 = strcat('0',num2str(h2));
   end
   t2 = strcat(strY2,strM2,strD2,strh2);
   NewFname = strcat('ocean_his.nc_',t1,'-',t2);
   %movefile(fname,NewFname);
   copyfile(fname,NewFname);
end