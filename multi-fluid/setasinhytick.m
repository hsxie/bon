% 21-09-22 13:57
% to check

function setasinhytick(fac,ymin,ymax)
% ymin<0, ymax>0
ymax=max(ymax,1/fac);
ymin=min(ymin,-1/fac);

% close all; clear; clc;
% figure;
% n2=-10:0.01:10; nx=length(n2); x=(1:nx)/nx;
% fac=10; ymin=-100; ymax=100;
% plot(x,asinh(fac*n2),'linewidth',2);

pmax=round(log10(ymax));
pmin=round(log10(1/fac));
mmax=round(log10(-ymin));
mmin=round(log10(1/fac));

dy=2;
y1=mmax:-dy:mmin; y2=pmin:dy:pmax;
yval=[asinh(-fac*10.^y1),0,asinh(fac*10.^y2)];

nl1=length(y1);  nl2=length(y2);

for jl=1:nl1
    ystr{jl}=['-10^{',num2str(y1(jl)),'}'];
end
ystr{nl1+1}='0';
for jl=1:nl2
    ystr{nl1+1+jl}=['10^{',num2str(y2(jl)),'}'];
end

set(gca,'ytick',yval);
set(gca,'yticklabel',ystr);
ylim([min(yval),max(yval)]);

end