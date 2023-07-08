clc
clear
close all

%%
wname = 'db2';
%%
%read the seismogram
[t,tr] = readsac('2017S11.12H1810100112.181000.mseed_CSN1__SH_Z_SAC');
tr = detrend(tr,'linear');
tr = tr - mean(tr);
tr = tr + 0.0*rms(tr)*randn(size(tr));
% figure
 % plot(t,tr)
%%
%wavelet transform using wname==db4
W_tr = modwt(tr',wname);

%%
%envelope estimation using hilbert transform
env = zeros(size(W_tr,1)-1,size(W_tr,2));
a=size(W_tr);

for i=1:size(W_tr,1)-1
    env(i,:) = envelope(W_tr(i,:));
end
% figure
% imagesc(env)
%%
%plot the envelop sum
% figure
DF=sum(abs(env));
% plot(t,DF)
%%
%define the window for samples ER1
l = 100;
DF_new=[zeros(1,l),DF,zeros(1,l)];
ER1 = zeros(size(DF));
aa=.1:.1:1;
% mm=length(aa);
% for j=1:mm
%     nn(j)=aa(j);
  
for i = (l+1):(length(DF)+l)
    ER1(i-l)=sum(DF_new(i:i+l))/sum(DF_new(i-l:i));
end
ER1(1:l)=ER1(l+1)*ones(1,l);
ER1(end-l+1:end)=ER1(end-l-2)*ones(1,l);

w_peaks = diff(ER1);
w_peaks = w_peaks.*(w_peaks>0);

[peak,tt]=findpeaks(w_peaks,'MinPeakHeight',.04,'MinPeakDistance',l);%,'Threshold',.1);

 
%%
%calculating ER2
ER2=ER1.*abs(DF);
% figure;plot(t,ER2,'k')
figure
% subplot(3,2,1:2)
figure;plot(t,tr,'k')
if (wname(1:3)=='db4')
    title('TU.KSN.U.SAC ')
elseif(wname(1:3)=='db2')
    title('2017B11.12H1810100112.181000.mseed_KBAM__SH_Z_SACC')
elseif (wname(1:3) =='db3')
    title('Dubichease3')
elseif (wname(1:3)=='db1')
    title ('Dubichease 1')
elseif (wname(1:3)=='db8')
    title ('dubicheas 8')
elseif (wname=='sym2')
    title ('symlet 2 "TU.KSN.U.SAC"')
    elseif (wname =='sym3')
    title ('symlet 3')
     elseif (wname =='sym8')
    title ('symlet 8')
end
% zoom on;
% %pause();
% zoom off;
%%
hold on
  plot([t(tt),t(tt)],[max(abs(tr)),-max(abs(tr))],':r','linewidth',1.5)
 
 tt1= input('tr ra vared konid:');
 hold on
 plot([tt1,tt1],[max(abs(tr)),-max(abs(tr))],'k--','linewidth',1.5)
% subplot(3,2,3)
figure;imagesc(env)
colormap('gray')
set(gca,'YDir','normal')
% subplot(3,2,4)
figure;plot(t,DF,'r')
% subplot(3,2,5:6)
figure;plot(w_peaks,'b')
% hold on
% plot(diff(ER1),'r')
title(sprintf('MinPeakHeight=%0.1f',.4))
% subplot(3,2,6)
% plot(t,ER2)
hold on;
plot(tt,peak,'*k')
saveas(gcf,sprintf('MinPeakHeight_%0.1f_%s.jpg',.4,wname),'jpeg')