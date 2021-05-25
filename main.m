% mms.db_init('local_file_db','/data/mms'); %Accessing database

% starttime_brstData = [2017 01 15 06 42 45];
% stoptime_brstData =  [2017 01 15 06 43 15];
% startepochtt_brstData = irf_time(starttime_brstData,'vector>epochtt');        % Starting date
% stopepochtt_brstData = irf_time(stoptime_brstData,'vector>epochtt');          % Stop Date
% tint_brstData  = irf.tint(startepochtt_brstData,stopepochtt_brstData);

% clc, clear all
% plot(fdist_e.data(1,:,1,1));
% print -dpng bilder/delme.png


ggrList = [5, 10, 20];
load('shockList.mat')
list = zeros(length(Sdwb.Data(1,:,1,1)),5);

kB = 1.38*1e-23;
mu0 = 4*pi*1e-7;
mp = 1.67e-27;
me = 9.11e-31;
qe = 1.602e-19;
DegRad = pi/180;

%%
for p = 1:length(Sdwb.Data(1,:,1,1))



format long;
startTime = Sdwb.Data(3,p,1,1);
stopTime = Sdwb.Data(4,p,1,1);
if startTime == 0
    continue
else
% % tint = Sdwb.Data(2,p,1,1);
% tint = stopTime-startTime;
% 
startTime_brst = irf_time(startTime,'epoch>utc');
stopTime_brst = irf_time(stopTime,'epoch>utc');


% time_int = irf_time(tint,'tint>utc')
% 
% tint_brstData = irf.tint(time_int(1),time_int(2));

% starttime_brstData = [2017 01 15 06 42 45];
% stoptime_brstData =  [2017 01 15 06 43 15];
% startepochtt_brstData = irf_time(starttime_brstData,'vector>epochtt');        % Starting date
% stopepochtt_brstData = irf_time(stoptime_brstData,'vector>epochtt');          % Stop Date
tint_brstData1  = irf.tint(startTime_brst,stopTime_brst);

mms.db_init('local_file_db','/data/mms');
ne1 = mms.db_get_ts('mms1_fpi_brst_l2_des-moms','mms1_des_numberdensity_brst',tint_brstData1);
[~,index] = max(ne1.data);
disp(index);
t = ne1.time(index)+[-30, 30];
if ne1.time(index)-ne1.time(1) < 30 || ne1.time(end)-ne1.time(index) < 30
    tint_brstData = tint_brstData1;
else
    tint_brstData = t;
end

fdist_e = mms.db_get_ts('mms1_fpi_brst_l2_des-dist','mms1_des_dist_brst',tint_brstData);
B1gse = mms.db_get_ts('mms1_fgm_brst_l2','mms1_fgm_b_gse_brst_l2',tint_brstData); %B-field in GSE-coordinates
E1gse = mms.db_get_ts('mms1_edp_brst_l2_dce','mms1_edp_dce_gse_brst_l2',tint_brstData); %E-field in GSE-coordinates
EnergyBin_e = mms.db_get_ts('mms1_fpi_brst_l2_des-dist','mms1_des_energy_brst',tint_brstData); %electron energy
EnergyBinDelta_e = mms.db_get_ts('mms1_fpi_brst_l2_des-dist','mms1_des_energy_delta_brst',tint_brstData);
Eespectomni1 =  mms.db_get_ts('mms1_fpi_brst_l2_des-moms','mms1_des_energyspectr_omni_brst',tint_brstData);
pos_gsm = mms.get_data('R_gsm',tint_brstData); %


MaPBzGSM = irf_get_data_omni(tint_brstData+[-600 600],'Ma,P,BzGSM,ByGSM,Bx', 'omni_min'); % 1st column time (seconds since 1970), 2nd column velocity
[aa nvec_gsm] = model.magnetopause_normal(pos_gsm.gsmR1(1,:)./(6372), MaPBzGSM(1,4), MaPBzGSM(1,3), 'bs', MaPBzGSM(1,2));
Bvec = [MaPBzGSM(1,6) MaPBzGSM(1,5) MaPBzGSM(1,4)];

ang = irf_ang(Bvec, nvec_gsm, 1);
if ang > 90
    ang = 180-ang;
end

list(p,1) = p;
list(p,3) = ang;
disp(ang);




B1abs = irf_abs(B1gse);
B1gseall = TSeries(B1gse.time,[B1abs.data B1gse.data(:,1) B1gse.data(:,2) B1gse.data(:,3)]); %All magnetic field values
E1gseall = TSeries(E1gse.time,[E1gse.data(:,1) E1gse.data(:,2) E1gse.data(:,3)]); %All electric field values

He = 0;
Dphi = ( 11.25*DegRad )*ones(1,32); 
theta = DegRad*linspace(11.25/2,180-11.25/2,16);
dTHETA = sin( theta ).*( 11.25*DegRad*ones(1,16) );
rotEdE_e =  2*( (EnergyBin_e.data(:,:)).^(0.5) ).*EnergyBinDelta_e.data(:,:); 
velVolElem_e = (10^6)*( (2*(qe/me)^3)^(0.5) )*rotEdE_e; %Phase space volume element in (cm/s)^3 

o=2;
disp(p);

int = 0;
e0 = 0;
e02 = 0;
for i = 1:32 % phi
    for j = 1:16 % theta
        for x=5:19 % energi
            ss = (x.*fdist_e.data(:,x,i,j)).*velVolElem_e(:,x)*dTHETA(j)*Dphi(i); % x*f(x)*dtheta*dphi v?ntev?rde
            ss2 = (EnergyBin_e.data(:,x)).*(fdist_e.data(:,x,i,j)).*velVolElem_e(:,x)*dTHETA(j)*Dphi(i); % x*f(x)*dtheta*dphi v?ntev?rde
            e0 = e0 + ss;
            e02 = e02 + ss2;
            tt = fdist_e.data(:,x,i,j).*velVolElem_e(:,x)*dTHETA(j)*Dphi(i); % vanlig integral
            int = int + tt;
        end
    end
end

expValReal = e02./int;
expVal = e0./int; % normaliserar v?ntev?rdet

evTest1 = expValReal(end);
evTest2 = expValReal(1);
if evTest1 < evTest2
    evNoll = evTest1; % sparar f?rsta v?ntev?rdet
else
end
    evNoll = evTest2;
list(p,4) = evNoll;

disp('vantevarde');disp(evNoll);

% evMax = 0;
% loc = 0;
% for k = 1:length(expVal) % letar efter tider d? v?ntev?rdet st?rst
%     if expVal(k) > evMax
%         evMax = expVal(k);
%         loc = k; % sparar indexv?rdet
%         evHeat = expVal(k); % sparar v?ntev?rdet vid uppv?rmningen
%     end
% end

evMax = 0;
loc = 0;
for k = 1:length(expValReal) % letar efter tider d? v?ntev?rdet st?rst
    if expValReal(k) > evMax
        evMax = expValReal(k);
        loc = k; % sparar indexv?rdet
        evHeat = expValReal(k); % sparar v?ntev?rdet vid uppv?rmningen
    end
end

list(p,5) = evMax;
disp('uppvarmning');disp(evHeat);

% evHeatReal = EnergyBin_e.data(loc,floor(evHeat)); % tar fram riktiga energiv?rdet
% srchValRaw = evHeatReal*10; % multiplicerar med 10 f?r att leta efter acceleration
srchValRaw = evHeat*ggrList(o); % multiplicerar med 10 f?r att leta efter acceleration
srchVal = find(EnergyBin_e.data(loc,:) > srchValRaw, 1); % letar efter index f?r motsvarande v?rde i energybin


[accVal, accLoc] = max(fdist_e.data(:,srchVal,:,:)); % hittar de index d?r densiteten ?r som h?gst f?r 10*v?ntev?rdet vid uppv?rmning
truVal = max(accVal(:)); % f?r att f? ett enskilt v?rde

%indVal = find(accVal == truVal); % hittar index f?r maxv?rdet
timeOfAcc = accLoc(accVal == truVal); % sparar tiden d?r maximala accelerationen sker


rankLoc = find(EnergyBin_e.data(loc,:) >= evMax);

% rank = truVal./fdist_e.data(loc, floor(evHeat), 1,1);

rank = truVal./mean(fdist_e.data(loc, rankLoc(1), :,:), 'all');

list(p, 2) = rank;

disp(rank);disp('ranking')

s2 = num2str(p);
R = num2str(1000*rank);
txtHead = ['Num:' s2 ' ' 'Ranking:' R ' ' 'EV:' num2str(evNoll)];
% first plot

specrec = struct('t',[],'f',[],'p',[]);
specrec.t = Eespectomni1.time.epochUnix;
specrec.f = EnergyBin_e.data;
specrec.p = Eespectomni1.data;

FS1 = 16;
FS2 = 0.8*FS1;
nplot = 4; % Number of plots in the ladder plot
h = irf_plot(nplot);


hca = irf_panel('B1');
irf_plot(hca,B1gseall);
title(txtHead,'parent',hca)
irf_pl_mark(hca,specrec.t(timeOfAcc), 'k')
irf_pl_mark(hca,specrec.t(loc), 'r')
irf_legend(hca,{'|B|';'B_{x}';'B_{y}';'B_{z}'},[1.02 0.75],'FontSize',FS2);
%irf_legend(hca,{['win = ' num2str(smoothwinBsec) ' s']},[0.02 0.97],'FontSize',FS2);
hca.FontSize = FS1;
l1 = ylabel(hca,'B_{GSE} (nT)','FontSize',FS1);
set(l1,'interpreter','tex');
hca.FontSize = FS1;
% strng = 'bilder/b1gseall.png';
% 
% set(gcf,'paperpositionmode','auto')
% print(strng, '-dpng')

% second plot

hca = irf_panel('ne1');
irf_plot(hca,ne1);
irf_pl_mark(hca,specrec.t(timeOfAcc), 'k')
irf_pl_mark(hca,specrec.t(loc), 'r')
ylabel(hca,'n_e (cm^{-3})','FontSize',FS1);
hca.FontSize = FS1;
% strng = 'bilder/ne1.png';
% 
% set(gcf,'paperpositionmode','auto')
% print(strng, '-dpng')

%third

hca = irf_panel('E1');
irf_plot(hca,E1gseall);
irf_pl_mark(hca,specrec.t(timeOfAcc), 'k')
irf_pl_mark(hca,specrec.t(loc), 'r')
irf_legend(hca,{'E_{x}';'E_{y}';'E_{z}'},[1.02 0.75],'FontSize',FS2);
hca.FontSize = FS1;
l1 = ylabel(hca,'E_{GSE} (mV/m)','FontSize',FS1);
set(l1,'interpreter','tex');
hca.FontSize = FS1;
% strng = 'bilder/E1.png';
% 
% set(gcf,'paperpositionmode','auto')
% print(strng, '-dpng')

% fourth

hca = irf_panel('newfigure');
irf_spectrogram(hca,specrec)
irf_pl_mark(hca,specrec.t(timeOfAcc), 'k')
irf_pl_mark(hca,specrec.t(loc), 'r')
set(gca,'yscale','log');
ylabel('E (eV)')
colorbar

% expectVal = interp1(expVal, Energybin_e.data, 1:32)



hold on
irf_plot([fdist_e.time.epochUnix expVal])


irf_plot_axis_align
irf_zoom(h,'x',tint_brstData);
%irf_zoom(h,'y');
irf_timeaxis(h);


s1 = 'bilder/';
s3 = num2str(ggrList(o));
angStr = num2str(ang);

ev = num2str(expValReal(1));
evM = num2str(evMax);
evL = num2str(expValReal(timeOfAcc));


% strng = [s1 R '_' 'num' s2 '.png'];
strng = [s1 'num' s2 '_' R '_' angStr '.png'];

set(gcf,'paperpositionmode','auto')
print(strng, '-dpng')
% 
% disp('plot')

set(gcf,'paperpositionmode','auto')
print('bilder/expvalreal', '-dpng')
% s1 = 'bilder/';
% s2 = 'test1';
% s3 = '.png';
% strng = strcat(s1,s2,s3);
% 
% set(gcf,'paperpositionmode','auto')
% print(strng, '-dpng')
end
end

list = sortrows(list, 2);
writematrix(list);


%%




% 
% v = irf_get_data_omni(tint_brstData+[-600 600],'Bx,By,Bz', 'omni_min'); % 1st column time (seconds since 1970), 2nd column velocity
% v(v==9999) = NaN;                % set NaN where no data
% irf_plot(v)
% ylabel('Solar wind speed [km/s]')
% legend({'x', 'y', 'z'});
% xlabel('')