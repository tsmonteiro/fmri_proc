function [peak_indices,avg_period,new_phasedata,midpoint_phasedata,midpoints]=get_phase_period(phasedata,type,threshold)

% gets the phase peaks from retroicor's algorithm and corrects them slightly
%  1. does not contain algorithm that was attempting to account for envelope modulation of resp signal
%     by having a non-smoothly increasing phase (it is now smoothed)
%  2. beginning and end phase chunks are replaced with the average period of signal
%  3. outliers in period that are near > 2 or < 0.5 times the average are replaced with 
%     either an extra peak or a removed peak to get rid of missed or extra detected resp or cardiac.
%   caveat: signal can't be more than 10 percent bad for this to work (this number isn't quantified)
%   major caveat: cannot handle properly if more than one extra or more than one missed peaks next to each other.
%
%  INPUT:  phasedata is the phase output by retroicor, type is either 'card' or 'resp'
%  OUTPUT: peak index, average period, the corrected phase waveform and the envelope if 'resp' selected
%          envelope of modulation is same timepoints as peak_indices
%          Note, the maximum phase data is no longer 2pi or pi, because then 0/2pi are counted twice
%          I fixed it so it linearly increases no matter what

if (exist('threshold','var')==0)
  threshold=0.62;
  disp(sprintf('setting threshold to %f',threshold));
end
if (exist('type','var')==0)
  disp('you must specify second param as type of phasedata, "card" or "resp" (single quotes)');
  return
end
offset=-1;
if (type=='card')
  disp('setting phase to run from 0 to 2pi (rather than -pi to pi)');
  offset=1;
end
if (type=='resp')
  disp('setting phase to run from -pi to pi');
  offset=0;
end
if (offset==-1)
  disp('did not choose "card" or "resp" as second input parameter');
  return;
end

% COUNT THE NOMINAL NUMBER OF PEAKS
tdim=length(phasedata);
ncycles=0;
for i=1:tdim-1
  if (phasedata(i)-phasedata(i+1)>3.)
    ncycles=ncycles+1;
  end
end

% FIND THE INDICES OF THOSE PEAKS
peak_indices=0;
cyclelen=0;
for pass=1:2
  if (pass==2)
    averagecycle=round(mean(cyclelen));
  end
  for i=1:ncycles
    if (i==1)
      beginning=0;
      ending=find_period_start(1,phasedata);
      cyclelen(i)=ending-beginning+1;
    else
      beginning=ending;
      ending=find_period_start(i,phasedata);
      cyclelen(i)=ending-beginning+1;
    end
    if (pass==2)
      peak_indices(i)=ending;
    end
  end
end

% GO THROUGH PEAKS AND FIND OUTLIERS THAT ARE MORE THAN 3 STD AWAY FROM MEAN AND REPLACE WITH APPROPRIATE
avg_period=round(mean(cyclelen));
cyclelen(1)=avg_period;
cyclelen(end)=avg_period;
%lessthan=length(find(cyclelen<mean(cyclelen)-threshold*std(cyclelen)));
lessthan=length(find(cyclelen<threshold*mean(cyclelen)));
if (lessthan>0)
  %bad_indices=find(cyclelen<mean(cyclelen)-threshold*std(cyclelen));
  bad_indices=find(cyclelen<threshold*mean(cyclelen));
  for i=1:lessthan
    % if less than 3 sigma, then period is spurious signal, remove this peak and save in a new matrix
    temp=peak_indices([1:bad_indices(i)-1]);
    temp=[temp peak_indices([bad_indices(i)+1:end])];
    peak_indices=temp;
    temp=cyclelen([1:bad_indices(i)-1]);
    temp=[temp cyclelen([bad_indices(i)+1:end])];
    cyclelen=temp;
    bad_indices=bad_indices-1;
  end
end
%morethan=length(find(cyclelen>mean(cyclelen)+threshold*std(cyclelen)));
morethan=length(find(cyclelen>(1/threshold)*mean(cyclelen)));
if (morethan>0)
  %bad_indices=find(cyclelen>mean(cyclelen)+threshold*std(cyclelen));
  bad_indices=find(cyclelen>(1/threshold)*mean(cyclelen));
  for i=1:morethan
    % if more than 3 sigma, then period is due to missed signal, put in new peak halfway between
    peak_indices(end+1)=round((peak_indices(bad_indices(i))-peak_indices(bad_indices(i)-1))/2.)+peak_indices(bad_indices(i)-1);
    temp=cyclelen([1:bad_indices(i)-1]);
    temp(bad_indices(i))=round((peak_indices(bad_indices(i))-peak_indices(bad_indices(i)-1))/2.);
    temp(bad_indices(i)+1)=round((peak_indices(bad_indices(i))-peak_indices(bad_indices(i)-1))/2.);
    cyclelen=[temp cyclelen([bad_indices(i)+1:end])];
  end
  peak_indices=sort(peak_indices);
end
avg_period=round(mean(cyclelen));

% CONVERT PEAKS INDICES INTO PHASE DATA FOR FIRST PEAK
new_phasedata=zeros(length(phasedata),1);
% use avg_period to estimate initial and end slope
for j=1:peak_indices(1)
  if (peak_indices(1)==1)
    new_phasedata(j)=pi;
  else
    % if peak_indices(1)>avg_period, we should assume we missed the first one, put it in here either at peak_indices(1)-avg_period,
    % or halfway between peak_indices(1)-avg_period and one (if peak_indices(1) is more than double avg_period)
    if (peak_indices(1)>avg_period)
      if (peak_indices(1)>2*avg_period)
        new_peak_index=round((peak_indices(1)-avg_period)/2);
      else
        new_peak_index=peak_indices(1)-avg_period;
      end
      len=peak_indices(1)-new_peak_index;
      if (j>new_peak_index)
        new_phasedata(j)=((j-new_peak_index-1)*(2*pi/len))-pi;
        %new_phasedata(j)=((j-new_peak_index-(peak_indices(i)+1))*(2*pi/len))-pi;
      else
        new_phasedata(j)=j*(2*pi/new_peak_index)-pi;
      end
    else
      new_phasedata(j)=((j+abs(avg_period-peak_indices(1))-1)*(2*pi/avg_period))-pi;
    end
  end
end
% CONVERT PEAKS INDICES INTO PHASE DATA FOR MID PEAKS
for i=1:length(peak_indices)-1
  len=peak_indices(i+1)-peak_indices(i);
  for j=peak_indices(i)+1:peak_indices(i+1)
     new_phasedata(j)=((j-(peak_indices(i)+1))*(2*pi/len))-pi;
  end
end
% CONVERT PEAKS INDICES INTO PHASE DATA FOR END PEAK
%len=abs(avg_period-(length(phasedata)-peak_indices(end)));
len=avg_period;
for j=peak_indices(end)+1:length(phasedata)
  if (peak_indices(end)+1>=length(phasedata))
    new_phasedata(j)=-pi;
  else
    new_phasedata(j)=((j-(peak_indices(end)+1))*(2*pi/len))-pi;
  end
end

% check pi bound
largerthanpi =find(abs(new_phasedata) > pi);
if length(largerthanpi)
  for n = 1:length(largerthanpi)
    phinx= largerthanpi(n);
    while new_phasedata(phinx) > pi
      new_phasedata(phinx) = new_phasedata(phinx) -2*pi;
    end
    while new_phasedata(phinx) < -pi
      new_phasedata(phinx) = new_phasedata(phinx) +2*pi;
    end
  end
end

% OFFSET DATA FOR CARDIAC DATA EXPECTED BY RETROICOR
if (offset==1)
  new_phasedata=new_phasedata+pi;
end
avg_period=round(mean(cyclelen));

% GET MIDPOINTS BETWEEN PEAKS DATA FOR MIDPOINT-ORIENTED PHASE DATA
midpoints(1)=avg_period;
for x=1:length(peak_indices)-1
  times(x)=peak_indices(x+1)-peak_indices(x);
  midpoints(x)=round(times(x)/2.0)+peak_indices(x);
end
% if first midpoint is greater than avg_period, push in extra midpoint
if (midpoints(1)>avg_period)
  midpoints=[midpoints(1)-avg_period midpoints];
end
if (midpoints(1)>avg_period)
  midpoints=[midpoints(1)-avg_period midpoints];
end
% if last midpoints is more than avg_period away from end, put extra midpoints
if ((length(new_phasedata)-midpoints(end))>avg_period)
  midpoints(end+1)=midpoints(end)+avg_period;
end

% CONVERT PEAKS INDICES INTO PHASE DATA FOR FIRST PEAK
midpoint_phasedata=zeros(length(phasedata),1);
% use avg_period to estimate initial and end slope
for j=1:midpoints(1)
  if (midpoints(1)==1)
    midpoint_phasedata(j)=pi;
  else
    % if midpoints(1)>avg_period, we should assume we missed the first one, put it in here either at midpoints(1)-avg_period,
    % or halfway between midpoints(1)-avg_period and one (if midpoints(1) is more than double avg_period)
    if (midpoints(1)>avg_period)
      if (midpoints(1)>2*avg_period)
        new_peak_index=round((midpoints(1)-avg_period)/2);
      else
        new_peak_index=midpoints(1)-avg_period;
      end
      len=midpoints(1)-new_peak_index;
      if (j>new_peak_index)
        midpoint_phasedata(j)=((j-new_peak_index-1)*(2*pi/len))-pi;
        %midpoint_phasedata(j)=((j-new_peak_index-(midpoints(i)+1))*(2*pi/len))-pi;
      else
        midpoint_phasedata(j)=j*(2*pi/new_peak_index)-pi;
      end
    else
      midpoint_phasedata(j)=((j+abs(avg_period-midpoints(1))-1)*(2*pi/avg_period))-pi;
    end
  end
end
% CONVERT PEAKS INDICES INTO PHASE DATA FOR MID PEAKS
for i=1:length(midpoints)-1
  len=midpoints(i+1)-midpoints(i);
  for j=midpoints(i)+1:midpoints(i+1)
    midpoint_phasedata(j)=((j-(midpoints(i)+1))*(2*pi/len))-pi;
  end
end
% CONVERT PEAKS INDICES INTO PHASE DATA FOR END PEAK
%len=abs(avg_period-(length(phasedata)-midpoints(end)));
len=avg_period;
for j=midpoints(end)+1:length(phasedata)
  if (midpoints(end)+1>=length(phasedata))
    midpoint_phasedata(j)=-pi;
  else
    midpoint_phasedata(j)=((j-(midpoints(end)+1))*(2*pi/len))-pi;
  end
end
% OFFSET DATA FOR CARDIAC DATA EXPECTED BY RETROICOR
if (offset==1)
  midpoint_phasedata=midpoint_phasedata+pi;
end
avg_period=round(mean(cyclelen));

return;

function a=find_period_start(cycle,inputphase)
  thiscycle=0;
  for q=1:length(inputphase)-1
    if (inputphase(q)-inputphase(q+1)>3.)
      thiscycle=thiscycle+1;
      if (cycle==thiscycle)
        a=q;
	return;
      end
    end
  end
  disp('did not find this cycle');
  a=0;
  return

