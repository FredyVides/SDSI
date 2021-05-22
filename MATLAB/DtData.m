function dtdata=DtData(data,dt,order)
if mod(order,2)~=0,error('order must be even');end
	Dt=DMatrix((0:order)-floor(order/2));
	m=floor(order/2);
	M=size(data,2);
	for j=1:m
		dtdata(:,j)=(1/dt)*(data(:,1:(1+order))*Dt(j,:).');
		dtdata(:,M-j+1)=(1/dt)*(data(:,(M-order):M)*Dt(order-j+2,:).');
	end
	for j=(m+1):(M-m)
		dtdata(:,j)=(1/dt)*(data(:,(j-m):(j+m))*Dt(m+1,:).');
	end
end
