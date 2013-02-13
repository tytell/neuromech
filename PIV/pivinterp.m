function varargout = pivinterp(varargin)
% function [u2,v2,DU,sp] = pivinterp(x,y,u,v,data, p,sigmin,snrmin,magstd, comp)
% or pivinterp(in,[out,] p,sigmin,snrmin,magstd, comp)

if (iscellstr(varargin{1})),
	in = varargin{1};
	if (iscellstr(varargin{2})),
		out = varargin{2};
		i = 3;
	else
		out = 0;
		i = 2;
	end;
	p = varargin{i};
	sigmin = varargin{i+1};
	snrmin = varargin{i+2};
	magstd = varargin{i+3};
	if (nargin == i+4)
		comp = varargin{i+4};
	else
		comp = [];
	end;
else
	x = varargin{1};
	y = varargin{2};
	u = varargin{3};
	v = varargin{4};
	data = varargin{5};
	
	p = varargin{6};
	sigmin = varargin{7};
	snrmin = varargin{8};
	magstd = varargin{9};
	if (nargin == 10)
		comp = varargin{10};
	else
		comp = [];
	end;
	
	in = 0;
	out = 0;
end;

N = size(in,1);
elapsed = 0;
timeperfile = [];

for i = 1:N,
	tic;

	if (iscellstr(in)),
		load(in{i},'x','y','u','v','data');
	end;
	
	us = repmat(NaN,size(x));
	vs = repmat(NaN,size(x));
	
	DU.dudx = repmat(NaN,size(x));
	DU.dudy = repmat(NaN,size(x));
	DU.dvdx = repmat(NaN,size(x));
	DU.dvdy = repmat(NaN,size(x));
	
	mag = sqrt(u.^2 + v.^2);
	
	k = find(((data.Error == 0) | (data.Error == 16)) & (data.Signal > sigmin) & ...
                (data.SNR > snrmin) & (mag < magstd*nanstd(mag(:))) & ...
                isfinite(x) & isfinite(y) & isfinite(u) & isfinite(v));

    if (isempty(timeperfile)),
		switch lower(comp),
		case 'scombrid',
			A = 5.403e-8;
			B = 2.6468;
		case 'trout',
			A = 1.4884e-7;
			B = 2.5479;
		otherwise
			A = NaN;
			B = NaN;
		end;
		
		if (~isnan(A)),
			dur = N * A * length(k)^B;
			fprintf('Estimated duration: %f seconds\n',dur);
			
			timeperfile = dur/N;
		end;
	end;
	
	smoothspline = tpapslong([x(k)+u(k) y(k)+v(k)]', [u(k) v(k)]', p);
	
	q = find(isfinite(x) & isfinite(y));
	uv = fnval(smoothspline,[x(q) y(q)]');
	dx = fnval(fnder(smoothspline,[1 0]),[x(q) y(q)]');
	dy = fnval(fnder(smoothspline,[0 1]),[x(q) y(q)]');
	
	us(q) = uv(1,:);
	vs(q) = uv(2,:);
	DU.dudx(q) = dx(1,:);
	DU.dvdx(q) = dx(2,:);
	DU.dudy(q) = dy(1,:);
	DU.dvdy(q) = dy(2,:);

	if (iscellstr(in)),
		if (iscellstr(out)),
			fn = out{i};
		else
			fn = in{i};
		end;
		save(fn,'us','vs','DU','smoothspline','-append');
	end;
	
	dur = toc;
	elapsed = elapsed + dur;
	timeperfile = elapsed/i;
	remain = (N-i)*timeperfile;
	
	if (remain ~= 0),
		fprintf(1, 'Elapsed: %.1f sec, Remaining: %.1f sec.\n',elapsed,remain);
	end;
end;

fprintf(1, 'Total time: %.1f sec',elapsed);

if (nargout > 0),
	varargout{1} = us;
	varargout{2} = vs;
	varargout{3} = DU;
	varargout{4} = smoothspline;
end;
