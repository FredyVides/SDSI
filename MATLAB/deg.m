function [d,rk0,q,H0,H1]=deg(L,N,data,tol)
% Example:
% data=importdata('AlmostPeriodicSignal.csv');
% sl=floor(.16*size(data,2));
% sample=data(2,1:sl);
% Lags=LagEstimate(sample,1e-1); 
% [d,rk,q,H0,H1]=deg(Lags(1),1,sample,1e-4);
% SpA=fliplr(SpSolver(H0.',H1(end,:).',d,1e-3).');
% A=fliplr(H1(end,:)/H0);
% ysp=ARForcast(SpA,H0(:,1),205-d);
% y=ARForcast(A,H0(:,1),205-d);
% subplot(411);plot(data(1,1:sl),sample,'b.-',data(1,:),data(2,:),'r',data(1,:),y,'k-.'); 
% legend('training data','signal','identified signal'),title('Standard system identification')
% subplot(412);stem(A),axis tight,title('Identified dynamics matrix')
% subplot(413);plot(data(1,1:sl),sample,'b.-',data(1,:),data(2,:),'r',data(1,:),ysp,'k-.'); 
% legend('training data','signal','identified signal'),title('Sparse system identification')
% subplot(414);stem(SpA),axis tight,title('Identified dynamics matrix')
% data=importdata('PeriodicSignal.csv');
% sl=floor(.14*size(data,2));
% sample=data(2,1:sl);
% Lags=LagEstimate(sample,1e-1); 
% [d,rk,q,H0,H1]=deg(Lags(1),1,sample,1e-6);
% SpA=fliplr(SpSolver(H0.',H1(end,:).',d,1e-3).');
% A=fliplr(H1(end,:)/H0);
% ysp=ARForcast(SpA,H0(:,1),205-d);
% y=ARForcast(A,H0(:,1),205-d);
% subplot(411);plot(data(1,1:sl),sample,'b.-',data(1,:),data(2,:),'r',data(1,:),y,'k-.'); 
% legend('training data','signal','identified signal'),title('Standard system identification')
% subplot(412);stem(A),axis tight,title('Identified dynamics matrix')
% subplot(413);plot(data(1,1:sl),sample,'b.-',data(1,:),data(2,:),'r',data(1,:),ysp,'k-.'); 
% legend('training data','signal','identified signal'),title('Sparse system identification')
% subplot(414);stem(SpA),axis tight,title('Identified dynamics matrix')
% data=importdata('NLOscillatorSignal.csv');
% t=data(1,:);
% Lt=length(t);
% sl=floor(Lt*.2);
% Lt=length(t);
% s=t(1):(t(Lt)-t(1))/(Lt-1):t(Lt);
% x1=interp1(t,data(2,:),s,'spline');
% x2=interp1(t,data(3,:),s,'spline');
% sample=x1(:,1:sl);
% Lags=LagEstimate(sample,1e-1);            
% [d,rk,q,H0,H1]=deg(Lags(1),1,sample,5e-8);
% A=fliplr(H1(end,:)/H0);
% SpA=fliplr(SpSolver(H0.',H1(end,:).',d,1e-4).');
% y=ARForcast(A,H0(:,1),Lt-d);
% ysp=ARForcast(SpA,H0(:,1),Lt-d);
% subplot(411);plot(s(1:sl),sample,'b.-',s,x1,'r',s,y,'k-.'); 
% legend('training data','signal','identified signal'),title('Standard system identification')
% subplot(412);stem(A),axis tight,title('Identified dynamics matrix')
% subplot(413);plot(s(1:sl),sample,'b.-',s,x1,'r',s,ysp,'k-.'); 
% legend('training data','signal','identified signal'),title('Sparse system identification')
% subplot(414);stem(SpA),axis tight,title('Identified dynamics matrix')
% figure
% sample=x2(:,1:sl);
% Lags=LagEstimate(sample,1e-1);            
% [d,rk,q,H0,H1]=deg(Lags(1),1,sample,5e-8);
% A=fliplr(H1(end,:)/H0);
% SpA=fliplr(SpSolver(H0.',H1(end,:).',d,1e-4).');
% y=ARForcast(A,H0(:,1),Lt-d);
% ysp=ARForcast(SpA,H0(:,1),Lt-d);
% subplot(411);plot(s(1:sl),sample,'b.-',s,x2,'r',s,y,'k-.'); 
% legend('training data','signal','identified signal'),title('Standard system identification')
% subplot(412);stem(A),axis tight,title('Autoregressive coefficients')
% subplot(413);plot(s(1:sl),sample,'b.-',s,x2,'r',s,ysp,'k-.'); 
% legend('training data','signal','identified signal'),title('Sparse system identification')
% subplot(414);stem(SpA),axis tight,title('Autoregressive coefficients')
S=size(data,2);
L=min(S-1,L);
sample=SPSampler(data,N,(S-1)/S);
[q,rk0]=rk(data,tol);
if rk0<min(size(data))
	q=q(:,1:rk0);
	d=1;
	[~,H0,H1]=HankelSMData(data(:,1:rk0),d,N);
	return;
end
for k=1:L 
	[H0,~,~]=HankelSMData(sample,k,N);
	[H1,~,~]=HankelSMData(data,k+1,N);
	s=svd(H0.');
	rk0=sum(s>tol);
	s=svd(H1.');
	rk1=sum(s>tol);
	if rk0==rk1 | k==L | rk0<min(size(H0))
		d=k;
		[~,H0,H1]=HankelSMData(data,d,N);
		[q,~,~]=svd(H0.',0);
		q=q(:,1:rk0);
		return;
	end
end
end
