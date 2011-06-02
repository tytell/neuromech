function t = timetpaps(x,y,u,v,p, n,reps)

good = find(isfinite(x) & isfinite(y) & isfinite(u) & isfinite(v));
N = length(good);

for i = 1:length(n),
    for j = 1:reps,
		k = round(rand(n(i),1)*(N-1)) + 1;
	
		fprintf('%d elements: ',n(i));
		
		tic;
		sp = tpapslong([x(k)+u(k)/2 y(k)+v(k)/2]',[u(k) v(k)]',p);
		t(i,j) = toc;
		
		fprintf('%f seconds\n',t(i,j));
	end;
	fprintf('* %d elements: mean %f +- %f sec\n', n(i), mean(t(i,:)), std(t(i,:))/sqrt(reps));
end;
