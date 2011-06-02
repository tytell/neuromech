function [circ,stats] = circulation(x,y,u,v, varargin)
% function [circ,stats] = circulation(x,y,u,v, varargin)

if (ischar(varargin{1})),
    method = varargin{1};
    varargin = varargin(2:end);
else
    method = 'circle';
end;

dx = diff(x,[],2);
dy = diff(y,[],1);
dr = nanmean(abs([dx(:); dy(:)]));

switch lower(method),
  case 'manual',
    xpt = varargin{1};
    ypt = varargin{2};
    
    if ((ndims(xpt) ~= ndims(ypt)) || ...
        ~all(size(xpt) == size(ypt))),
        error('circulation contours must have the same size');
    end;

    if (size(xpt,1) == 1),
        xpt = xpt';
        ypt = ypt';
    end;
    
    %run through all the contours
    N = size(xpt,2);
    circ = NaN(1,N);
    stats = struct('xpt',{},'ypt',{},'s',{},'upar',{});
    for i = 1:N,
        good = isfinite(xpt(:,i));
        
        %get the arc length
        s = [0; cumsum(sqrt(diff(xpt(good,i)).^2 + diff(ypt(good,i)).^2))];
        len = s(end);
        
        %and the value that's closest to and smaller than dr that will complete the
        %contour exactly
        dr1 = len / ceil(len/dr);
        
        %interpolate along the contour for even spacing in arc length
        s1 = (0:dr1:len);
        xpt1 = interp1(s,xpt(good,i), s1);
        ypt1 = interp1(s,ypt(good,i), s1);
        
        %interpolate the velocities
        ur = interp2(x,y,u, xpt1,ypt1, 'linear');
        vr = interp2(x,y,v, xpt1,ypt1, 'linear');
        
        %tangent vector (magnitude = 1 by definition)
        tanx = deriv(s1, xpt1);
        tany = deriv(s1, ypt1);
        
        %velocity tangent to the contour
        upar = tanx.*ur + tany.*vr;
        
        %integral
        circ(i) = trapz(s1,upar);
        
        %save some extra details
        stats(i).xpt = xpt1;
        stats(i).ypt = ypt1;
        stats(i).ur = ur;
        stats(i).vr = vr;
        stats(i).s = s1;
        stats(i).upar = upar;
    end;
    
 case 'circle',
  ctrx = varargin{1};
  ctry = varargin{2};
  R = varargin{3};
  if (size(R,2) == 1),
      R = R(:,[1 1]);
  else
     %sort the radii so that we have major (longest) axis first
      R = -sort(-R,2);
  end;

  if (length(varargin) == 4),
      angle = varargin{4};
  else
      angle = zeros(size(ctrx));
  end;
  
  for i = 1:length(ctrx),
      %put the radii together in reverse so that we're certain
      %the maximum value really is the known maximum radius.  The
      %smallest value may well be somewhat greater than dr/2, but that's
      %OK 
      r = (max(R(i,:)):-dr:dr/2)';
      r = r(end:-1:1);
      r(:,2) = r(:,1)*R(i,2)/R(i,1);

      for j = 1:size(r,1),
         %figure out the number of points to step around the 
         %circumference using about the same spacing as the 
         %points themselves have - so we don't oversample
          npt = ceil(2*pi*r(j,1)/dr);

         %first center it on zero
          theta = linspace(0,2*pi,npt+1);
          xpt = r(j,1).*cos(theta);
          ypt = r(j,2).*sin(theta);

         %then rotate the oval, if necessary
          if (angle(i) ~= 0),
              xpt2 = cos(angle(i))*xpt - sin(angle(i))*ypt;
              ypt2 = sin(angle(i))*xpt + cos(angle(i))*ypt;
              xpt = xpt2;
              ypt = ypt2;
          end;
      
         %then center on ctrx, ctry
          xpt = xpt + ctrx(i);
          ypt = ypt + ctry(i);

         %interpolate the points
          ur = interp2(x,y,u, xpt,ypt, '*cubic');
          vr = interp2(x,y,v, xpt,ypt, '*cubic');
      
         %vector tangent to the oval contour (not unit length)
          dxdtheta = -r(j,1)*sin(theta)*cos(angle(i)) - ...
              r(j,2)*cos(theta)*sin(angle(i));
          dydtheta = -r(j,1)*sin(theta)*sin(angle(i)) + ...
              r(j,2)*cos(theta)*cos(angle(i));
          mag = sqrt(dxdtheta.^2 + dydtheta.^2);

         %use it to get contour length
          s = cumtrapz(theta, sqrt(dxdtheta.^2 + dydtheta.^2));

         %get the tangential velocity along the contour
          upar = dxdtheta./mag.*ur + dydtheta./mag.*vr;
          
         %integrate along the contour length
          rcirc1(j) = trapz(s, upar);
      end;

      [q,ind] = max(abs(rcirc1));

      stats(i).CircSteps = rcirc1;
      stats(i).rSteps = r;
      stats(i).rAtMaxCirc = r(ind);
      stats(i).MaxCirc = rcirc1(ind);
      
      circ(i) = rcirc1(ind);
  end;

 otherwise,
  error('Unsupported circulation method %s.',method);
end;




	