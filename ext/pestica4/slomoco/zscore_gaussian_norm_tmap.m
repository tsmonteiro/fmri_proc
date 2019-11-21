function [zstudt,mu,sigma]=zscore_gaussian_norm_tmap(studt,bins,plott)
%chi2=sum((gdist(f1:f2)-n(f1:f2)).^2)/sum(n(f1:f2).^2);
% routine to calculate mean and standard deviation of distribution of student's t
% from mean-1/2FWHM to mean+1/2FWHM
% Author: M. Lowe
% Date: 19-DEC-2005
% modified E. Beall Jan 2007
% modification: fits using maximum likelihood instead of normfit
% this is not the fischer's z-transform, this is instead an empirical conversion
% of a connectivity distribution to a normalized gaussian. This assumes
% that the bulk of the distribution is null, which may or may not be appropriate
% for your situation
dim=size(studt);
studt=reshape(studt,[prod(dim) 1]);
zstudt=zeros(prod(dim),1);

temp=studt(find(studt(:)~=0));

mm=min(temp);
ma=max(temp);

if (exist('bins')==0)
  if (dim(1)==128)
    bins=250;
  elseif (dim(1)==64)
    bins=62;
  else
    disp('voxel size unexplored so far...');
  end
end

i=mm:(ma-mm)/bins:ma;
sthist=hist(temp,i);

% have histogram, now find mean-1/2FWHM, mean+1/2FWHM

mast=max(sthist);

jj=find(sthist>0.5*mast);

%jj(1) is lowest index, jj(size(jj,2)) is highest index

sortstudt=sort(temp);  % sort studt to find fit range
stind1=find(sortstudt>i(jj(1))); % stind1(1) is lowest studt in FWHM
stind1=stind1(1);
stind2=find(sortstudt<i(jj(size(jj,2)))); % stind2(size(stind2,1)) is high in FWHM
stind2=stind2(size(stind2,1));
% fit over FWHM only
[mu,sigma]=normfit(sortstudt(stind1:stind2));
% setup range of data to fit over (mu+-3sigma)
rdata=sortstudt(find(sortstudt>mu-3*sigma & sortstudt<mu+3*sigma));
xaxis=linspace(min(rdata),max(rdata),bins);
yaxis=hist(rdata,xaxis);
% find FWHM to fit over
a=xaxis(find(xaxis>sortstudt(stind1),1));
b=xaxis(find(xaxis<=sortstudt(stind1))); b=[a b(end)];
f1=b(find(abs(b-sortstudt(stind1))==min(abs(b-sortstudt(stind1)))));
f1=find(xaxis==f1);
a=xaxis(find(xaxis>sortstudt(stind2),1));
b=xaxis(find(xaxis<=sortstudt(stind2))); b=[a b(end)];
f2=b(find(abs(b-sortstudt(stind2))==min(abs(b-sortstudt(stind2)))));
f2=find(xaxis==f2);
% need to truncate CDF on both sides, assuming same region
pdf_truncnorm=@(x,mu,sigma) normpdf(x,mu,sigma) ./ (normcdf(xaxis(f2),mu,sigma) - normcdf(xaxis(f1),mu,sigma));
[paramEsts,paramCIs] = mle(sortstudt(stind1:stind2), 'pdf',pdf_truncnorm, 'start',[mu,sigma], 'lower',[-Inf 0]);

mu=paramEsts(1);
sigma=paramEsts(2);

[n,xaxis]=hist(temp,bins);
gdist=(sum(n)/sum(normpdf(xaxis,mu,sigma)))*normpdf(xaxis,mu,sigma);
% find height at FWHM of studt and gdist (max/2), normalize to FWHM
a=xaxis(find(xaxis>sortstudt(stind1),1));
b=xaxis(find(xaxis<=sortstudt(stind1))); b=[a b(end)];
f1=b(find(abs(b-sortstudt(stind1))==min(abs(b-sortstudt(stind1)))));
f1=find(xaxis==f1);
a=xaxis(find(xaxis>sortstudt(stind2),1));
b=xaxis(find(xaxis<=sortstudt(stind2))); b=[a b(end)];
f2=b(find(abs(b-sortstudt(stind2))==min(abs(b-sortstudt(stind2)))));
f2=find(xaxis==f2);
%normalization=((n(f1)+n(f2))/2)/(max(gdist)/2);
normalization=max(n)/max(gdist);
gdist=gdist*normalization;
nr=hist(sortstudt(stind1:stind2),xaxis);
chi2=sum(((n(f1:f2)-gdist(f1:f2)).^2)/var(n(f1:f2)));
chi2=chi2pdf(chi2,length(f1:f2));
if (exist('plott')==0)
  plott=0;
end
if (plott==1)
  figure;
  plot(xaxis,n);
  hold on
  plot(xaxis,gdist,'r');
  plot(xaxis,nr,'k');
  legend(sprintf('mu: %.02f, sigma: %.02f, chi2: %.02f',mu,sigma,chi2));
end

% take input studt "im" and convert tscores into corr coeffs and back again, applying correction
for x=1:length(studt)
  if (studt(x)==0)
    zstudt(x)=0;
    continue;
  end
  zstudt(x)=(studt(x)-mu)/sigma;
end
if (length(dim)==2)
  zstudt=reshape(zstudt,[dim(1) dim(2)]);
elseif (length(dim)==3)
  zstudt=reshape(zstudt,[dim(1) dim(2) dim(3)]);
elseif (length(dim)==4)
  zstudt=reshape(zstudt,[dim(1) dim(2) dim(3) dim(4)]);
else
  zstudt=reshape(zstudt,[dim(1) dim(2) dim(3) dim(4) dim(5)]);
end

