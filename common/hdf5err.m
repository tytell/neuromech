function varargout = hdf5err(fcn,varargin)
% function [out...] = hdf5err(fcn,params...)
% Catches errors in the HDF5 functions and returns an error as a return
% parameter, rather than throwing an error.  Makes HDF5 error handling a
% little simpler.
%
% Example:
%   [datasetid1,err] = ...
%      hdf5err(@H5D.open,fileid,gpdata.channels(ch).hwsdata);
%   Calls H5D.open(fileid,gpdata.channels(ch).hwsdata)
%   If it can't find the dataset, returns an error number in err, otherwise
%   err is empty.
%
% Mercurial revision hash: $Revision$ $Date$
% Copyright (c) 2010, Eric Tytell <tytell at jhu dot edu>


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
            match = regexp(except.identifier, '(H5ML_hdf5|hdf5lib)', 'once');
            if (~isempty(match)),
                varargout = cell(1,nargout);
                varargout{nargout} = except;
            else
                throwAsCaller(except);
            end;
    end;
end;

        
    