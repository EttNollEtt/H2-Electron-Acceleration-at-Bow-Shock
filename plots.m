%%%
load('shockList.mat')
kB = 1.38*1e-23;
mu0 = 4*pi*1e-7;
mp = 1.67e-27;
me = 9.11e-31;
qe = 1.602e-19;
DegRad = pi/180;

startTime = Sdwb.Data(3,4,1,1);
stopTime = Sdwb.Data(4,4,1,1);
startTime_brst = irf_time(startTime,'epoch>utc');
stopTime_brst = irf_time(stopTime,'epoch>utc');
tint_brstData = irf.tint(startTime_brst,stopTime_brst);
mms.db_init('local_file_db','/data/mms');
fdist_e = mms.db_get_ts('mms1_fpi_brst_l2_des-dist','mms1_des_dist_brst',tint_brstData);
EnergyBin_e = mms.db_get_ts('mms1_fpi_brst_l2_des-dist','mms1_des_energy_brst',tint_brstData); %electron energy


data = importdata('list.txt');
ang = data(end:-1:end-143, 3);

ang = rmmissing(ang);
xvals = 1:length(ang);
% for i=1:length(ang)
%     if ang(i) == ISNAN
%         ang2(i) = [];
%     else
%         pass
%     end
% end

c = polyfit(xvals, ang, 1);
y_est = polyval(c,xvals);
err = 20*ones(1,length(ang));

figure(1);
plot(xvals, y_est, 'k');
hold on
stem(xvals, ang, 'r');
% hold on
% er = errorbar(xvals,ang,err);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none'; 
legend('Fitted line')
ylabel('Angle in degrees')
xlabel('Event number in the ranked list')
set(gcf,'paperpositionmode','auto')
print('angPlot.png', '-dpng')

% arr = zeros(1,length(fdist_e.data(10,:,1,1)));
% for k = 1:length(fdist_e.data(10,:,1,1))
%     arr(k) = mean(fdist_e.data(10,k,:,:), 'all');
% end
% 
% figure(2);
% loglog(EnergyBin_e.data(1500,:), fdist_e.data(1000,:,1,1));
% xlabel('Energy (eV)')
% ylabel('s^{3}cm^{-6}');

set(gcf,'paperpositionmode','auto')
print('distfunc.png', '-dpng')
        