function tdata=convert_slicextime_to_timeseries(slicedata,slice_acq_order,keep_slices)

zdim = size(slicedata,1);
tdim = size(slicedata,2);
if (zdim ~= length(slice_acq_order))
  disp('Error; slice data does not match slice acq order vector');
  return
end

if ~exist('keep_slices')
  keep_slices = ones(1,zdim);
end

tdata = zeros(1,zdim*tdim);
tx = 1:zdim*tdim;
for z = 1:zdim
  if keep_slices(z)
    tdata(z:zdim:zdim*tdim) = slicedata(slice_acq_order(z),:);
  else
    tx(z:zdim:zdim*tdim)=0;
  end
end