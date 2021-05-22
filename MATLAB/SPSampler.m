function sample=SPSampler(data,L,sp)
	S=size(data,2);
	s=floor(S/L);
	J=1:s:(S-s+1);
	s=min(s,floor(sp*s));
	S=[];
	for j=1:L
		S=[S J(j):(J(j)+s-1)];
	end
	sample=data(:,S);
end
