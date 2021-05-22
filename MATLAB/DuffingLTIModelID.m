function [A0,x0,A,x1,G1,G2,rmse1,rmse2]=DuffingLTIModelID(L,t,data,sp,pp,tol)
	% [t,x]=DuffingNetwork(2.2,-36,0,.2,[8 7 4 15 14 9],[0 20]);
	% [A0,x0,A1,x1,G1,G2,rmse1,rmse2]=DuffingLTIModelID(80,t,x,.3,.5,1e-6);
[st,sx]=DataSpliner(t,data(:,1:3));
sample=SPSampler(sx(:,1:3).',1,sp);
Lags=LagEstimate(sample(1,:),1e-1);
E=eye(3);
g1=E([2 3 1],:);
g2=E([1 3 2],:);
spsample=[sample g1*sample g1^2*sample  g2*sample g2*g1*sample g2*g1^2*sample];
[d,rk,q,H0,H1]=deg(Lags(1),6,spsample,tol);
Lag=max(Lags(1),d);
L=min(L,Lag);
[~,H0,H1]=HankelSMData(spsample,L,6);
A0=kron(diag(ones(1,L-1),1),E);
A=A0;
E=eye(L);
G1=kron(E,g1);
G2=kron(E,g2);
disp('Running SDSI:')
tic,
A0((3*L-2):3*L,:)=SpSolver(H0.',H1((3*L-2):3*L,:).',3*L,tol).';
A0=(A0+G1'*A0*G1+(G1^2)'*A0*(G1^2))/3;
A0=(A0+G2'*A0*G2)/2;
toc
disp('Running SINDy:')
tic,
A((3*L-2):3*L,:)=SINDy(H0.',H1((3*L-2):3*L,:).',500,tol).';%H1((3*L-2):3*L,:)/H0;
A=(A+G1'*A*G1+(G1^2)'*A*(G1^2))/3;
A=(A+G2'*A*G2)/2;
toc
St=floor(sp*size(sx,1));
Pt=floor(pp*St);
x0=sx(St:(St+L-1),:).';
x0=x0(:);
x1=x0;
for k=1:Pt, x0 = [x0 A0*x0(:,k)];x1 = [x1 A*x1(:,k)];end
figure(1),
subplot(311),plot(st(1:(St+Pt+1)),sx(1:(St+Pt+1),1),'k',st(St:(St+Pt)),x0(1,:),'r-.',st(1:St),sx(1:St,1),'b','linewidth',1.5);
legend('Reference signal','Identified signal','Training data')
subplot(312),plot(st(1:(St+Pt+1)),sx(1:(St+Pt+1),2),'k',st(St:(St+Pt)),x0(2,:),'r-.',st(1:St),sx(1:St,2),'b','linewidth',1.5);
legend('Reference signal','Identified signal','Training data')
subplot(313),plot(st(1:(St+Pt+1)),sx(1:(St+Pt+1),3),'k',st(St:(St+Pt)),x0(3,:),'r-.',st(1:St),sx(1:St,3),'b','linewidth',1.5);
legend('Reference signal','Identified signal','Training data')
rmse1(1)=sqrt(sum((x0(1,:)-sx(St:(St+Pt),1).').^2)/(St+Pt+1));
rmse1(2)=sqrt(sum((x0(2,:)-sx(St:(St+Pt),2).').^2)/(St+Pt+1));
rmse1(3)=sqrt(sum((x0(3,:)-sx(St:(St+Pt),3).').^2)/(St+Pt+1));
rmse2(1)=sqrt(sum((x1(1,:)-sx(St:(St+Pt),1).').^2)/(St+Pt+1));
rmse2(2)=sqrt(sum((x1(2,:)-sx(St:(St+Pt),2).').^2)/(St+Pt+1));
rmse2(3)=sqrt(sum((x1(3,:)-sx(St:(St+Pt),3).').^2)/(St+Pt+1));
figure(2)
subplot(311),plot(st(1:(St+Pt+1)),sx(1:(St+Pt+1),1),'k',st(St:(St+Pt)),x1(1,:),'r-.',st(1:St),sx(1:St,1),'b','linewidth',1.5);
legend('Reference signal','Identified signal','Training data')
subplot(312),plot(st(1:(St+Pt+1)),sx(1:(St+Pt+1),2),'k',st(St:(St+Pt)),x1(2,:),'r-.',st(1:St),sx(1:St,2),'b','linewidth',1.5);
legend('Reference signal','Identified signal','Training data')
subplot(313),plot(st(1:(St+Pt+1)),sx(1:(St+Pt+1),3),'k',st(St:(St+Pt)),x1(3,:),'r-.',st(1:St),sx(1:St,3),'b','linewidth',1.5);
legend('Reference signal','Identified signal','Training data')
figure(3),
subplot(121),spy(A0,'k');
subplot(122),spy(A,'k');
end
