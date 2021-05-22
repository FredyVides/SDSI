function [t,x]=DuffingNetwork(a,b,c,n,x0,T)
	f=@(t,y)[y(4);y(5);y(6);
	        -c*y(4)-y(1).*(b+(a^2)*y(1))+n*((y(1)-y(2))+(y(1)-y(3)));
		-c*y(5)-y(2).*(b+(a^2)*y(2))+n*((y(2)-y(1))+(y(2)-y(3)));
		-c*y(6)-y(3).*(b+(a^2)*y(3))+n*((y(3)-y(1))+(y(3)-y(2)))];
	opt=odeset ("RelTol", eps);
	[t,x] = ode45 (f, [T(1), T(2)], x0,opt);
end
