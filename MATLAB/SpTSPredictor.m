% function [xt,xl,y1,y2,A1,A2,rmse1,rmse2]=SpTSPredictor(file,ssp,sp,pp,L,tol,delta)
% Examples:
% L=1100;
% [xt,xl,y1,y2,A1,A2,rmse1,rmse2]=SpTSPredictor('../DataSets/TemperatureData.csv',1/8412,.5,1,L,3e-2,3e-2);
% figure(1)
% subplot(211),plot(1:length(xl),xl,'k',1:length(y1),y1,'r'),axis tight
% subplot(212),stem(A1),axis tight;
% figure(2)
% subplot(211),plot(1:length(xl),xl,'k',1:length(y2),y2,'r'),axis tight
% subplot(212),stem(A2),axis tight
function [xt,xl,y1,y2,A1,A2,rmse1,rmse2]=SpTSPredictor(file,ssp,sp,pp,L,tol,delta)
x=importdata(file);
sl=length(x);
ssp=floor(sl*ssp);
x=x(1:ssp:sl);
sl=length(x);
sp=floor(sp*sl);
xt=x(1:sp);
pp=max(min(floor(pp*sp),sl-sp),L);
xl=x((sp+1):end);
L0=LagEstimate(xt,1e-1);
L=max(L0(1),L);
H=hankel(xt(1:L),xt(L:end));
H0=H(:,1:(end-1));
H1=H(end,2:end);
disp('=============================================');
disp('Computing predictive model with SDSI:')
tic,A1=fliplr(SpSolver(H0.',H1.',L,tol,delta).');toc
disp('=============================================');
disp('=============================================');
disp('Computing predictive model with SINDy:')
tic,A2=fliplr(SINDy(H0.',H1.',500,delta).');toc
disp('=============================================');
pp=max(pp,length(xl));
y1=DLTIPredictor(A1,xt((end-1099):end),pp);
y2=DLTIPredictor(A2,xt((end-1099):end),pp);
rmse1=sqrt(sum((y1-xl.').^2)/length(xl));
rmse2=sqrt(sum((y2-xl.').^2)/length(xl));
