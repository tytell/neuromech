function varargout = hdf5err(fcn,varargin)

try
    fcnout = cell(1,nargout-1);
    [fcnout{:}] = feval(fcn,varargin{:});
        
    varargout(1:nargout-1) = fcnout(1:nargout-1);
    varargout{nargout} = [];
catch except
    switch except.identifier,
        case {'MATLAB:TooManyOutputs'},
            except = MException('hdf5err:wrongNumberOfOutputs',...
                'Wrong number of outputs for hdf5err (one more than the hdf5 function)');
            throw(except);
            
        otherwise,
            match = regexp(except.identifier, 'H5ML_hdf5', 'once');
            if (~isempty(match)),
                varargout = cell(1,nargout);
                varargout{nargout} = except;
            else
                throwAsCaller(except);
            end;
    end;
end;

        
    