%TVF(time,amplitude,period,group_vel,dist)
%time and amplitude are seismograms, period and group_vel is a 
%group velocity function in terms of period, and the distance between stations
%This will interpolate the smoothly varying group_vel function 
%and create a time varied filtered record that should optimize that
%signal.

function [anew,f]=tvf(t,aa,T,gvw,dist)
dt=t(2)-t(1);
f=zeros(length(aa),1);
dimx=length(aa);
dimx2=ceil(dimx/2);
f(1:dimx2)=[0:(dimx2-1)]/dimx/dt;
if mod(length(f),2)==0
f(dimx2+1:end)=flipud(f(2:dimx2-1));
else
f(dimx2+1:end)=flipud(f(2:dimx2));
end
T0=1./f;
anew=zeros(size(aa));
iii=find(T0(1:dimx2)>=(T(1)) & T0(1:dimx2)<=T(end));
Tin=T0(iii);
gvww=interp1(T,gvw,T0(iii),'pchip');
gradgvw=gradient(gvww,.5);%
%luv=find(T0(iii)>=6 & T0(iii)<8);gradgvw(luv)=300/80;
%luv=find(T0(iii)>=8 & T0(iii)<12);gradgvw(luv)=300/80;%these are used for the middle periods to expand the window
%plot(T0(iii),3.5+80*abs(gradgvw))
b=fft(aa);
for i=1:length(Tin)
%This is to make the window...
Ta=round(dist/gvww(i)-Tin(i)*(3.5+80*abs(gradgvw(i))));
Tb=round(dist/gvww(i)+Tin(i)*(3.5+80*abs(gradgvw(i))));
Tm=round(dist/gvww(i));
kk=find(t>=Ta & t<=Tb);
wind=zeros(size(aa));
wind(kk)=cos(pi*(t(kk)-Tm)/(Tb-Ta));


%This is to make the fft of the freqency of interest.
c=zeros(size(b));
jj=find(T0==Tin(i));
c(jj)=b(jj);
anew=anew+wind.*ifft(c);
end
