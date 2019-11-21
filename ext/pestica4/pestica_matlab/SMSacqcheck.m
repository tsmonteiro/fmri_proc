function [MBacc zmbdim uniq_slice_timing_ms uniq_acq_order ] = SMSacqcheck(TRms, zdim, slice_timing_ms)


timegap_slice = TRms/zdim; % second unit
uniq_slice_timing_ms = unique(slice_timing_ms,'stable');
[uniq_sorted_slice_timing_ms uniq_acq_order] = sort(uniq_slice_timing_ms);
timegap_acq_ms = mean(diff(uniq_sorted_slice_timing_ms));

if length(slice_timing_ms) > length(uniq_slice_timing_ms) 
  MBacc = round(timegap_acq_ms/timegap_slice);
  disp(['Note: SMS acquisition is applied (MB acc. fac = ' num2str(MBacc) ').'])
else
  MBacc =1;
end
zmbdim = zdim/MBacc;
