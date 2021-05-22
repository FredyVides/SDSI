function [d,rmse,A,Ap]=Example1
t=(0:256)*1/32;
y=importdata('../DataSets/Example1Signal.csv');
yp=importdata('../DataSets/Example1NoisySignal.csv');
[d,rk,q,H0,H1]=deg(25,1,yp(1:70),1e-2);
A=fliplr(SpSolver(H0.',H1(end,:).',d,4.5e-2).');
s=LTIPredictor(A,H0(:,1),length(y)-(d));
rmse=sqrt((sum((y((d+1):end)-s((d+1):end)).^2)/length(s((d+1):end))));
figure(1),
subplot(411),plot(t,y,'k',t,s,'r-.','linewidth',1.5)
legend('Original signal','Identified signal')
L=[1 7 12];
for k=1:3
Hp=hankel(yp(1:(d-L(k))),yp((d-L(k)):end));
Hp0=Hp(:,1:(end-1));
Hp1=Hp(:,2:end);
Ap{k}=fliplr(SpSolver(Hp0.',Hp1(end,:).',d-L(k),4.5e-2).');
sp=LTIPredictor(Ap{k},Hp0(:,1),length(y)-(d-L(k)));
subplot(4,1,k+1),plot(t,y,'k',t,sp,'r-.','linewidth',1.5)
legend('Original signal','Identified signal')
rmse=[rmse sqrt((sum((y((d-L(k)+1):end)-sp((d-L(k)+1):end)).^2)/length(sp((d-L(k)+1):end))))];
end
figure(2),
subplot(411),stem(A,'k'),axis tight;
for k=1:3
subplot(4,1,k+1),stem(Ap{k},'k'),axis tight;
end
