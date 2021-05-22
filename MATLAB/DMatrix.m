function D=DMatrix(x)
E=@(j,k)(j==k);
order=length(x)-1;
for j=1:(order+1)
dp(j,:)=polyder(polyfit (x,E(j,1:(order+1)),order));
end
for j=1:(order+1)
for k=1:(order+1)
D(j,k)=polyval(dp(k,:),x(j));
end
end
end
