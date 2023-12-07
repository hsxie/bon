% 2021-08-25 22:25 Hua-sheng XIE, ENN, huashengxie@gmail.com
% Name BO-n (bon.m), a branch of BO/BO-Ray/BO-n code series.
% Solve multi-fluid plasma model with thermal effect to obtain k_per(omega)
% Which could contain more physics than in the cold plasma model.
% 21-08-26 07:47 not ok yet
% 07:58 use eig(matA/matB) instead of eig(matA,matB) seem ok.
% But still not agree with cold plasma model for T->0
% 20:50 fix a bug on sum_{s,sp}, still not ok
% 21-08-27 16:09 fix a bug of initialize matrix matA, etc.
% seem ok now! agree with previous 'coldwave_N2_scan_tnpar.m' with small T.
% 17:16 add scan density
% 21-09-22 13:39 rewrite the plot and scan (f, npara, ns0, R)
% 15:49 add the cold plasma option

close all; clear; clc;

icase=1;
icold=0; % icold=0, warm; 
         % icold=1, cold plasma
         
if(icase==1) % scan f, 21-09-22 13:39
    xx=10.^(6:0.01:9); % give f
    B0=6e0; % background B field, Tesla
    n00=2e20; % m^-3
%     f=28e9; % Hz
    jcase=2;
    if(jcase==1) % fixed npara
        npar=0.5;
    elseif(jcase==2) % fix kz
        kz=10;
    end
elseif(icase==2) % scan density, 21-09-23 17:08
    jcase=3;
    if(jcase==1)
        xx=10.^(15:0.01:21); % give density
        B0=1.5e0; % background B field, Tesla
        npar=1.5; % m^-3  % 1.4445
        f=0.5e9; % Hz
    elseif(jcase==2)
        xx=10.^(18:0.001:20); % give density
        B0=1.5e0; % background B field, Tesla
%         npar=0.77; % m^-3  % 0.77
%         f=1*28e9; % Hz
        npar=0.7; % m^-3  % 0.65
        f=2*28e9; % Hz
    elseif(jcase==3) % 21-09-25 08:23 for ICRF
        xx=10.^(16:0.001:20); % give density
        B0=3e0; % background B field, Tesla
        npar=5.0; % m^-3  % 0.65
        f=50e6; % Hz
    end
elseif(icase==3) % scan npar, 21-09-24 08:19
    xx=10.^(-2:0.01:3); % give n_para
    B0=1.5e0; % background B field, Tesla
    n00=50e18; % m^-3
    f=0.5e9; % Hz
elseif(icase==4) % scan R, 21-09-27 08:13    
    
    f=60e6; % Hz
    npar=5.0; % give n_para
    npar2=npar^2;
    
    % 2. set parameters, modify here for your own case
    R0=2.0;
    xx=(1.0:0.001:3.0);
    R=xx;
    
    B00=3e0; % background B field, Tesla
    n000=5e19;
    
    B0x=B00*R0./R;
    % den0x=n00*exp(-(R-R0).^2/(0.5^2));
    deltanR=0.5;
    den0x=n000*exp(-sqrt((R-R0).^4)/(deltanR^2));
    jcase=1;
    
%     xx=10.^(-2:0.01:3); % give R
%     B0=1.5e0; % background B field, Tesla
%     n00=50e18; % m^-3
%     f=0.5e9; % Hz
end
nx=length(xx);

% 1. default SI unit
c2=(2.99792458E8)^2; % speed of ligth c^2
epsilon0=8.854187817E-12;
mu0=1/(c2*epsilon0);
kB=1.38064852e-23;
qe=1.60217662e-19; % electron charge, coulombs
mp=1.6726219e-27; % proton mass, kg
me=9.1094e-31; % electron mass, kg
c=sqrt(c2);

% 2. set parameters, modify here for your own case
ms=[me,2*mp]; % speceis mass
qs=[-1,1]*qe; % speceis charge
ns_unit=[1,1]; % *den0, species density, m^-3
Ts_unit=[1,1]; % eV
if(icase==2  && jcase==3)
ms=[me,1*mp,2*mp]; % speceis mass
qs=[-1,1,1]*qe; % speceis charge
ns_unit=[1,0.1,0.9]; % *den0, species density, m^-3
Ts_unit=[1,1,1]; % eV
end
if(icase==4  && jcase==2)
ms=[me,1*mp,2*mp]; % speceis mass
qs=[-1,1,1]*qe; % speceis charge
ns_unit=[1,0.1,0.9]; % *den0, species density, m^-3
Ts_unit=[1,1,1]; % eV
end
S=length(ms);

T00=1e3; % eV

if(icold==0)
    NN=2*S+4;
elseif(icold==1) % cold plasma
    NN=2;
end
% matA=zeros(NN,NN); % 21-08-27 16:04 wrong!
nper2x=zeros(NN,nx);
kx1x=zeros(NN,nx);

for jx=1:nx
    
if(icase==1)
    f=xx(jx);
    if(jcase==2)
        w=f*2*pi; % omega
        npar=kz*c/w;
    end
elseif(icase==2)
    n00=xx(jx);
elseif(icase==3)
    npar=xx(jx);
elseif(icase==4)
    n00=den0x(jx);
    B0=B0x(jx);
end

w=f*2*pi; % omega

kz=npar*w/c;

ns0=ns_unit*n00;
Ts0=Ts_unit*T00*1.6022e-19/kB; % Ts, eV -> K
cs=sqrt(2*kB*Ts0./ms); % thermal velocity, note we use gamma_s=2
cs2=cs.^2;

wcs=qs*B0./ms; % the gyro frequency
wps=sqrt(ns0.*qs.^2./ms/epsilon0); % plasma frequency
wcs2=wcs.^2;
wps2=wps.^2;
wp2=sum(wps2);

if(icold==0)
matA=zeros(NN,NN);

% set matrix elements
for s=1:S
    
    % (1). dns
    indi=2*(s-1)+1;
    indj=2*(s-1)+2;
%     matA(indi,indj)=ns0(s)*(w-wcs2(s)/w)/cs2(s); % (dns,dvsx)
    matA(indi,indj)=matA(indi,indj)+ns0(s)*(w-wcs2(s)/w)/cs2(s); % (dns,dvsx),21-08-27 08:02
    for sp=1:S
        indj=2*(sp-1)+2;
%         matA(indi,indj)=matA(indi,indj)-ns0(s)/w*wps2(sp)/cs2(s); % (dns,dvspx)
        matA(indi,indj)=matA(indi,indj)-... % fix a bug, 21-08-26 20:32
            ns0(s)/w*qs(s)/ms(s)*qs(sp)*ns0(sp)/epsilon0/cs2(s); % (dns,dvspx)
    end
    indj=2*S+1;
    matA(indi,indj)=ns0(s)/w*(qs(s)*wcs(s)/ms(s))/cs2(s); % (dns,dEy)
    indj=2*S+3;
    matA(indi,indj)=-1i*ns0(s)/w*(qs(s)*kz*c2/ms(s))/cs2(s); % (dns,dBy)
       
    % (2). dvsx
    indi=2*(s-1)+2;
    indj=2*(s-1)+1;
    matA(indi,indj)=(w-kz^2*cs2(s)/w)/ns0(s); % (dvsx,dns)
    indj=2*S+2;
    matA(indi,indj)=-1i*qs(s)/ms(s)*kz*ns0(s)/w/ns0(s); % (dvsx,dEz)
    
end

% (3). dEy
indi=2*S+1;
indj=2*S+4;
matA(indi,indj)=w; % (dEy,dBz)

% (4). dEz
indi=2*S+2;
indj=2*S+3;
matA(indi,indj)=kz^2*c2/w-w; % (dEz,dBy)
for sp=1:S
    indj=2*(sp-1)+2;
    matA(indi,indj)=-1i*kz/(w*epsilon0)*ns0(sp)*qs(sp); % (dEz,dvspx)
end

% (5). dBy
indi=2*S+3;
indj=2*S+2;
matA(indi,indj)=(wp2/w-w)/c2; % (dBy,dEz)
for sp=1:S
    indj=2*(sp-1)+1;
    matA(indi,indj)=(-1i*kz/(w*epsilon0)*cs2(sp)*qs(sp))/c2; % (dBy,dnsp)
end

% (6). dBz
indi=2*S+4;
indj=2*S+1;
matA(indi,indj)=(w-kz^2*c2/w-wp2/w)/c2; % (dBz,dEy)
for sp=1:S
    indj=2*(sp-1)+2;
    matA(indi,indj)=B0/w*wps2(sp)/c2; % (dBz,dvspx)
end

kxx=eig(matA);
% d0=vpa(eig(matA),32); kxx=double(d0);

[~,ind]=sort(real(kxx.^2),'descend');
kxx=kxx(ind);

nper2=(kxx*c/w).^2;

N2=nper2+npar^2;

elseif(icold==1) % cold plasma, 21-09-22 15:56
    T00=0; Ts0=0;
    
    dielS=1-sum(wps2./(w^2-wcs2));
    dielD=sum(wcs.*wps2./(w*(w^2-wcs2)));
    dielP=1-sum(wps2./w^2);
    
    npar2=npar^2;
    
    % the polynomial coefficients
    polyc4=dielS;
    polyc2=-((dielS+dielP)*(dielS-npar2)-dielD^2);
    polyc0=dielP*((dielS-npar2)^2-dielD^2);
    
    disppolynomial=[polyc4, polyc2, polyc0];
    %         n2temp=sort(real(roots(disppolynomial)),'descend');
    %         n2temp=roots(disppolynomial); [~,ind]=sort(abs(n2temp),'descend');
    n2temp=roots(disppolynomial); [~,ind]=sort(real(n2temp),'descend');
    nper2=n2temp(ind);
    
%     nper2=[(-polyc2+sqrt(polyc2^2-4*polyc4*polyc0))/(2*polyc4);
%         (-polyc2-sqrt(polyc2^2-4*polyc4*polyc0))/(2*polyc4)];
    
%     kxx=sqrt(nper2); % wrong, 21-09-27 00:43
    kxx=sqrt(nper2*w^2/c2);
end

nper2x(:,jx)=nper2;
kx1x(:,jx)=kxx;

end
%% plot
close all;
h=figure('unit','normalized','Position',[0.01 0.05 0.6 0.7],...
    'DefaultAxesLineWidth',2,...
    'DefaultAxesFontSize',15);
%   'DefaultAxesFontWeight','bold',...

if(icase==4) % scan R, 21-09-27 08:13
    ax1=axes('position',[0.1,0.1,0.55,0.8]);
    
else
    ax1=axes('position',[0.1,0.1,0.8,0.8]);
end

jplt=2;
if(jplt==1)
    for jn=1:NN
        semilogx(xx,real(nper2x(jn,:)),'.','linewidth',3,'MarkerSize',12); hold on;
    end
    
    patch([min(xx),max(xx),max(xx),min(xx)],[0,0,-1e10,-1e10],'m');
    alpha(0.1);
    ylim([-1e2,1e2]);
    grid on;grid minor;
    box on;
    
    Clog=-2; %
    symlog(gca,'y',Clog);
elseif(jplt==2)
    
    fac=10; ylmax=1e10; ylmin=-1e10;
    semilogx(xx,0.*xx,'k-','linewidth',2); hold on;
    for jn=1:NN
        yy=squeeze(nper2x(jn,:))+0*npar^2;tmp=imag(yy);
%         semilogx(xx,asinh(fac*real(yy)),...
%             'b.','linewidth',3,'MarkerSize',12); hold on;
        
%         scatter(xx,asinh(fac*real(yy)),25,...
%             log10(abs(tmp)./abs(yy)),'filled'); hold on;
        
patch([xx,xx(end)],[asinh(fac*real(yy)),NaN],[abs(imag(yy)./abs(yy)),NaN],...
    'EdgeColor','interp','MarkerFaceColor',...
    'flat','FaceColor','none','LineWidth',2,'LineStyle','-'); hold on;

%         semilogx(xx,asinh(fac*real(sqrt(nper2x(jn,:).*xx.^2/c2))),...
%             '.','linewidth',3,'MarkerSize',12); hold on;
    end
    
    setasinhytick(fac,ylmin,ylmax);
    grid on;
    patch([min(xx),max(xx),max(xx),min(xx)],...
        [0,0,asinh(ylmin*fac),asinh(ylmin*fac)],'m');
    alpha(0.1);
%     ylim([-1e5,1e5]);
end

if(icase==1)
title(['B_0=',num2str(B0),'T, S=',num2str(S),', n_{00}=',num2str(n00,3),...
    'm^{-3}, T_{00}=',num2str(T00,3),'eV, f=',num2str(f,3),...
    'Hz,',10,'q_s/e=[',num2str(qs/qe,3),'], m_s/m_e=[',num2str(ms/me,4),'], ',...
    'n_s/n_{00}=[',num2str(ns_unit,3),'], T_s/T_{00}=[',num2str(Ts_unit,3),']']);
xlabel('f=\omega/2\pi [Hz]');
elseif(icase==2)
title(['B_0=',num2str(B0),'T, S=',num2str(S),', n_{||}=',num2str(npar,3),...
    ', T_{00}=',num2str(T00,3),'eV, f=',num2str(f,3),...
    'Hz,',10,'q_s/e=[',num2str(qs/qe,3),'], m_s/m_e=[',num2str(ms/me,4),'], ',...
    'n_s/n_{00}=[',num2str(ns_unit,3),'], T_s/T_{00}=[',num2str(Ts_unit,3),']']);
xlabel('n_{00}[m^{-3}]');
elseif(icase==3)
title(['B_0=',num2str(B0),'T, S=',num2str(S),', n_{00}=',num2str(n00,3),...
    'm^{-3}, T_{00}=',num2str(T00,3),'eV, f=',num2str(f,3),...
    'Hz,',10,'q_s/e=[',num2str(qs/qe,3),'], m_s/m_e=[',num2str(ms/me,4),'], ',...
    'n_s/n_{00}=[',num2str(ns_unit,3),'], T_s/T_{00}=[',num2str(Ts_unit,3),']']);
xlabel('n_{||}=k_{||}c/\omega');
elseif(icase==4)
    
    if(jcase==1)
        wchx=abs(qs(2)*B0x./ms(2))*2; % wc_H
        wcdx=abs(qs(2)*B0x./ms(2)); % wc_D
    else
        wchx=abs(qs(2)*B0x./ms(2)); % wc_H
        wcdx=abs(qs(3)*B0x./ms(3)); % wc_D
    end
    
    xwch=interp1(wchx,xx,w,'linear','extrap');
    xwcd=interp1(wcdx,xx,w,'linear','extrap');
%     wpex=sqrt(ns_unit(1)*den0x.*qs(1)^2/ms(1)/epsilon0); % wpe
%     wuhx=sqrt(wcex.^2+wpex.^2);
%     idwch=find(abs(wchx-w)<=3e-1*w);
%     idwcd=find(abs(wcdx-w)<=3e-1*w);
    
    plim=1;
    if(plim==1)
        
%         for jn=1:NN
%             semilogy(xx,real(nper2x(jn,:)),'.','linewidth',3,'MarkerSize',12); hold on;
%         end
        
        plot([xwch,xwch],[-1e3,1e3],'k--','linewidth',1); hold on;
        plot([xwcd,xwcd],[-1e3,1e3],'m--','linewidth',1); hold on;
        hold on;
        
        text(xwch,2e-1,'f_{c,H}','fontsize',14,'color','k','HorizontalAlignment','right');
        text(xwcd,2e2,'f_{c,D}','fontsize',14,'color','m','HorizontalAlignment','left');
        
%         ylim([1e-1,1e7]);xlim([0.2,1.5]);
    elseif(plim==2)
        
%         for jn=1:NN
%             plot(xx,real(nper2x(jn,:)),'.','linewidth',3,'MarkerSize',12); hold on;
%         end
%         plot([xx(idwch(end)),xx(idwch(end))],[-1e10,1e10],'k--','linewidth',2); hold on;
%         plot([xx(idwcd(end)),xx(idwcd(end))],[-1e10,1e10],'m--','linewidth',2); hold on;
        hold on;
        
        patch([min(xx),max(xx),max(xx),min(xx)],[0,0,-1e10,-1e10],'m');
        alpha(0.1);
        ylim([-1e8,1e8]);
        grid on;grid minor; box on;
        
        xlim([min(xx),max(xx)]);
        
        
        Clog=-2; %
        symlog(gca,'y',Clog);
        
%         text(xx(idwch(end)),-2,'f_{c,H}','fontsize',14,'color','k','HorizontalAlignment','right');
%         text(xx(idwcd(end)),2,'f_{c,D}','fontsize',14,'color','m','HorizontalAlignment','left');
    end
    title(['S=',num2str(S),', n_{||}=',num2str(npar),...
        ', T_{00}=',num2str(T00,3),'eV, f=',num2str(f,3),'Hz, n_{000}=',...
        num2str(n000,3),'m^{-3}, q_s/e=[',num2str(qs/qe,3),'],',10,...
        'n_s/n_0=[',num2str(ns_unit,3),'], m_s/m_e=[',num2str(ms/me,4),']']);
    ylabel('n_\perp^2=k_\perp^2c^2/\omega^2');
    
    xlabel('R[m]');
    grid on; box on;%grid minor;
    xlim([min(xx),max(xx)]);
    %
    
    ax2=axes('position',[0.75,0.1,0.2,0.35]);
    for js=1:S
        %     plot(R,den0x*ns_unit(js),'linewidth',2); hold on;
        semilogy(R,den0x*ns_unit(js),'linewidth',2); hold on;
    end
    xlabel('R [m]');ylabel('n_s [m^{-3}]');
    grid on;
    
%     if(plim==1)
%         xlim([0.2,1.5]);
%     end
    
    ax3=axes('position',[0.75,0.55,0.2,0.35]);
    plot(R,B0x/B00,'linewidth',2); hold on;
    plot(R,den0x/n000,'linewidth',2); hold on;
    xlabel('R [m]');ylabel('B_0 [T]');
    legend('B_0/B_{00}','n/n_{00}');legend('boxoff');
    
%     if(plim==1)
%         xlim([0.2,1.5]);
%     end
    
    title(['B_{00}=',num2str(B00),'T, R_0=',num2str(R0),...
        'm, \Delta{}R=',num2str(deltanR),'m']);
% title(['B_0=',num2str(B0),'T, S=',num2str(S),', n_{00}=',num2str(n00,3),...
%     'm^{-3}, T_{00}=',num2str(T00,3),'eV, f=',num2str(f,3),...
%     'Hz,',10,'q_s/e=[',num2str(qs/qe,3),'], m_s/m_e=[',num2str(ms/me,4),'], ',...
%     'n_s/n_{00}=[',num2str(ns_unit,3),'], T_s/T_{00}=[',num2str(Ts_unit,3),']']);
% xlabel('R [m]');
end
ylabel('n_\perp^2=k_\perp^2c^2/\omega^2'); 


%%
save(['bon_scan_icase=',num2str(icase),',S=',num2str(S),...
    ',B0=',num2str(B0,3),',npar=',num2str(npar,3),...
    ',n00=',num2str(n00,3),',T00=',num2str(T00,3),',f=',num2str(f,3),...
    ',icold=',num2str(icold),',jplt=',num2str(jplt),'.mat'],...
    'S','ms','qs','B0','ns0','Ts0','icold','f','xx','wcs',...
    'nx','wps','wps2','nper2x','kx1x','kz','npar','NN');


