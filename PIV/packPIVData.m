function data = packPIVData(data)

err = cat(3,data.Error);
iscorr = cat(3,data.IsCorrected);

data = rmfield(data,{'Error','IsCorrected'});

for i = 1:size(err,3),
	data(i).Error = sparse(err(:,:,i));
	data(i).IsCorrected = sparse(iscorr(:,:,i));
end;
