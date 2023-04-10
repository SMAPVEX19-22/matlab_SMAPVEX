close all; clear
%%
global HOMEDIR
%%
domainid   = 'MA';

satellite  = 'SMAP' %'SMOS'; % 'SMAP' %

refSM      = '5cm'; %'10cm'; %'V5'; %  

firstmonth = 5;
lastmonth  = 11;

savefigyes = false; %true;

figdir     = 'figs\';

curdate    = datestr(now, 'yyyymmdd');

%%
dnumA = datenum(2019,  5, 10);
dnumB = datenum(2022, 11,  1);

%%
for w = 1:4
    seasons_1(w) = datenum(2018+w,  5,  1);
    seasons_2(w) = datenum(2018+w, 10, 31);
end

%% Load SMAP and SMOS data over the site
% FROM: [HOMEDIR]_2019 - SMAPVEX19\05 - analysis\compare_SMOS_to_SMAP_and_save.m
fsatpath = [HOMEDIR, '_2019 - SMAPVEX19\05 - analysis\matlab_SMAPVEX\'];
fsatname = sprintf('SMAP_and_SMOS_TB_%s_30_50_v230407.mat', domainid);

load([fsatpath, fsatname]);

% smap = 
%   struct with fields:
% 
%          tbv: [1820×1 double]
%          tbh: [1820×1 double]
%           ad: [1820×1 char]
%         dnum: [1820×1 double]
%      sm_scah: [1820×1 double]
%      sm_scav: [1820×1 double]
%       sm_dca: [1820×1 double]
%     tau_scah: [1820×1 double]
%     tau_scav: [1820×1 double]
%      tau_dca: [1820×1 double]
%        tsurf: [1820×1 double]
% 
% smos = 
%   struct with fields:
% 
%         tbv: [2178×1 double]
%         tbh: [2178×1 double]
%        dnum: [2178×1 double]
%      mininc: 30
%      maxinc: 50
%     anglesN: 60
%         p_v: [0.9726 6.1318]
%         p_h: [0.9474 10.5859]
%          ad: [2178×1 char]

switch satellite
    case 'SMOS'
        sat = smos;
        AM  = 'A';
    case 'SMAP'
        sat = smap;
        AM  = 'D';
end

dvec = datevec(sat.dnum);

%% Load TempNet data

tnetpath  = [HOMEDIR, '_2019 - SMAPVEX19\06 - DATA\TempNet\'];

switch domainid
    case 'MA'
        tnetfilename = 'MA_Data_V0_20221207_20221207';%'MA_Data_V0b_20211119_20211201';%'MA_Data_v0_20200325'; % from read_and_save_TempNet_loggers.m
    case 'MB'
        tnetfilename = 'MB_Data_V0_20221130_20221206';%'MB_Data_V0_20211105_20211201';%'MB_Data_v0_20200325'; % from read_and_save_TempNet_loggers.m        
end

tn = load([tnetpath, tnetfilename]);
% -> 'domainid', 'stationids', 'dnum', 'SM', 'Ts', 'RDC', 'IDC', 'Ta'

%%
tnet_dnum = tn.dnum;

tnet_dataSM = tn.SM;

tnet_dataTs = tn.Ts;

tnet_dataTa = tn.Ta;

tnet_dataSM(tnet_dataTs < 4) = NaN;

tnet_sm5v = squeeze(tnet_dataSM(:,1,:));
tnet_sm05 = squeeze(tnet_dataSM(:,2,:));
tnet_sm10 = squeeze(tnet_dataSM(:,3,:));

tnet_Ts5v = squeeze(tnet_dataTs(:,1,:));
tnet_Ts05 = squeeze(tnet_dataTs(:,2,:));
tnet_Ts10 = squeeze(tnet_dataTs(:,3,:));

tnet_Tair = squeeze(tnet_dataTa(:,1,:));

%% check soil temp
% for m = 1:length(tnet_Ts5v(1,:))
%     figure
%     plot(tnet_dnum, tnet_Ts5v(:,m), '.');
%     title(sprintf('%.0f', tn.stationids(m)));
% end

% for m = 1:length(tnet_Ts5v(1,:))
%     figure
%     plot(tnet_dnum, tnet_sm5v(:,m), '.');
%     title(sprintf('%.0f', tn.stationids(m)));
% end

%% bad SM stations
% for m = 1:length(tnet_sm5v(1,:))
%     figure
%     plot(tnet_dnum, tnet_sm5v(:,m), '.');
%     title(sprintf('%.0f', tn.stationids(m)));
% end

switch domainid
    case 'MA'  
        verylittledata = [402,405,409,419,422,425];
        discontinuous  = [413,411,410,407,406];
        badTs5vstations = [404,414,415,416];
    case 'MB'
        verylittledata = [507];
        discontinuous  = [505,502];
        badTs5vstations = [];
end

find(ismember(tn.stationids, [discontinuous, verylittledata, badTs5vstations]))

tnet_sm5v(:,ismember(tn.stationids, [discontinuous, verylittledata])) = NaN;
tnet_sm05(:,ismember(tn.stationids, [discontinuous, verylittledata])) = NaN;
tnet_sm10(:,ismember(tn.stationids, [discontinuous, verylittledata])) = NaN;

tnet_Ts5v(:,ismember(tn.stationids, [discontinuous, verylittledata, badTs5vstations])) = NaN;
tnet_Ts05(:,ismember(tn.stationids, [discontinuous, verylittledata, badTs5vstations])) = NaN;
tnet_Ts10(:,ismember(tn.stationids, [discontinuous, verylittledata, badTs5vstations])) = NaN;

%%
tnet_sm5vm = nanmean_mathworks(tnet_sm5v, 2);
tnet_sm05m = nanmean_mathworks(tnet_sm05, 2);
tnet_sm10m = nanmean_mathworks(tnet_sm10, 2);

tnet_Ts5vm = nanmean_mathworks(tnet_Ts5v, 2);
tnet_Ts05m = nanmean_mathworks(tnet_Ts05, 2);
tnet_Ts10m = nanmean_mathworks(tnet_Ts10, 2);
tnet_Tairm = nanmean_mathworks(tnet_Tair, 2);

%% matchup

bs = sat.dnum < 0; % initialize 

for w = 1:length(seasons_1)
    
    bs = bs | (sat.dnum >= seasons_1(w) & sat.dnum <= seasons_2(w));
    
end

for i = 1:length(sat.dnum)
    
    if bs(i) & sat.ad(i) == AM %dnum(i) >= dnumA & dnum(i) <= dnumB
        
        [minv, ind] = min(abs(sat.dnum(i) - tnet_dnum));
        
        if minv < 0.5
            tnet_sm5v_q(i,1) = tnet_sm5vm(ind);
            tnet_sm05_q(i,1) = tnet_sm05m(ind);
            tnet_sm10_q(i,1) = tnet_sm10m(ind);
            tnet_dnum_q(i,1) = tnet_dnum(ind);
            tnet_Ts5v_q(i,1) = tnet_Ts5vm(ind);
            tnet_Ts05_q(i,1) = tnet_Ts05m(ind);
            tnet_Ts10_q(i,1) = tnet_Ts10m(ind);
            tnet_Tair_q(i,1) = tnet_Tairm(ind);
        else
            tnet_sm5v_q(i,1) = NaN;
            tnet_sm05_q(i,1) = NaN;
            tnet_sm10_q(i,1) = NaN;
            tnet_dnum_q(i,1) = NaN;
            tnet_Ts5v_q(i,1) = NaN;
            tnet_Ts05_q(i,1) = NaN;
            tnet_Ts10_q(i,1) = NaN;
            tnet_Tair_q(i,1) = NaN;
        end
    else
        tnet_sm5v_q(i,1) = NaN;
        tnet_sm05_q(i,1) = NaN;
        tnet_sm10_q(i,1) = NaN;
        tnet_dnum_q(i,1) = NaN;
        tnet_Ts5v_q(i,1) = NaN;
        tnet_Ts05_q(i,1) = NaN;
        tnet_Ts10_q(i,1) = NaN;
        tnet_Tair_q(i,1) = NaN;
    end
end

%%
switch refSM

    case 'V5'
        tnet_ref   = tnet_sm5vm;
        tnet_ref_q = tnet_sm5v_q;
    case '5cm'
        tnet_ref   = tnet_sm05m;
        tnet_ref_q = tnet_sm05_q;
    case '10cm'
        tnet_ref   = tnet_sm10m;
        tnet_ref_q = tnet_sm10_q;
end
%%

% ev = tbv./tsurf;
% eh = tbh./tsurf;

ev = sat.tbv./(tnet_Ts5v_q+273.15);
eh = sat.tbh./(tnet_Ts5v_q+273.15);

FF = 1;%0.7;

P = 1;%0.6;

Teff_F = (P*tnet_Ts5v_q+(1-P)*tnet_Tair_q + 273.15);

Teff_S = tnet_Ts5v_q + 273.15;

ev = sat.tbv./(FF*Teff_F + (1-FF)*Teff_S);
eh = sat.tbh./(FF*Teff_F + (1-FF)*Teff_S);

%% outliers
%
% ind1 = find(1-ev_MA < 0.08 & tnet_m5v_q > 0.4)
% ind2 = find(1-ev_MA < 0.07 & tnet_m5v_q > 0.28)
% tnet_m5v_q([ind1, ind2]) = NaN;
% ev_MA([ind1, ind2]) = NaN;

%%
rv = 1-ev;
rh = 1-eh;

%%

datelim = dvec(:,2) < firstmonth | dvec(:,2) > lastmonth;

bnan5v  = isnan(tnet_sm5v_q) | isnan(ev) | isnan(eh) | datelim;
bnan05  = isnan(tnet_sm05_q) | isnan(ev) | isnan(eh) | datelim;
bnan10  = isnan(tnet_sm10_q) | isnan(ev) | isnan(eh) | datelim;

bnanref = isnan(tnet_ref_q) | isnan(ev) | isnan(eh) | datelim;

R_m5v_V = corr(tnet_sm5v_q(~bnan5v), 1-ev(~bnan5v))
R_m05_V = corr(tnet_sm05_q(~bnan05), 1-ev(~bnan05))
R_m10_V = corr(tnet_sm10_q(~bnan10), 1-ev(~bnan10))
R_m5v_H = corr(tnet_sm5v_q(~bnan5v), 1-eh(~bnan5v))
R_m05_H = corr(tnet_sm05_q(~bnan05), 1-eh(~bnan05))
R_m10_H = corr(tnet_sm10_q(~bnan10), 1-eh(~bnan10))

R_ref_V = corr(tnet_ref_q(~bnanref), 1-ev(~bnanref))
R_ref_H = corr(tnet_ref_q(~bnanref), 1-eh(~bnanref))

p_m5v_V = polyfit(tnet_sm5v_q(~bnan5v), 1-ev(~bnan5v), 1);
p_m5v_H = polyfit(tnet_sm5v_q(~bnan5v), 1-eh(~bnan5v), 1);

p_ref_V = polyfit(tnet_ref_q(~bnanref), 1-ev(~bnanref), 1);
p_ref_H = polyfit(tnet_ref_q(~bnanref), 1-eh(~bnanref), 1);

%%

figure; hold all
plot(tnet_dnum, tnet_sm5v, '.');
plot(tnet_dnum, tnet_sm5vm, 'k.', 'markersize', 20);
zoom on; grid on;
datetick

figure; hold all
plot(tnet_dnum, tnet_Ts10, '.');
plot(tnet_dnum, tnet_Ts10m, 'k.', 'markersize', 20);
zoom on; grid on;
datetick
%%
figure; hold all
set(gca, 'fontsize', 13);
title([domainid,' TempNet']);
plot(tnet_dnum, tnet_sm5vm, '-', 'linewidth', 3);
plot(tnet_dnum, tnet_sm05m, '-', 'linewidth', 3);
plot(tnet_dnum, tnet_sm10m, '-', 'linewidth', 3);
set(gca, 'xlim', [dnumA, dnumB])
set(gca, 'ylim', [0 0.55]);
ylabel('Uncalibrated SM');
legend('0-5', '5', '10');
zoom on; grid on;
datetick('keeplimits');
%% 

fntsz = 13;
mrksz = 8;

switch domainid
    case 'MA'
       ylimh1 = 0.07;
       ylimh2 = 0.16;    
       ylimv1 = 0.04;
       ylimv2 = 0.13;        
    case 'MB'
       ylimh1 = 0.07;
       ylimh2 = 0.18;    
       ylimv1 = 0.04;
       ylimv2 = 0.14;                
end

%% Reflectivity V-pol and SM
figure; hold all;
set(gcf, 'position', [68         870        3208         420]);
ax1 = gca;
set(ax1, 'ylim', [0 0.55])
set(ax1, 'xlim', [dnumA dnumB])
ax2 = axes('Position', get(ax1, 'Position'), 'XAxisLocation', 'bottom', 'YAxisLocation', 'right', ...
    'Color', 'none', 'XColor', 'k', 'YColor', 'k');
set(ax1, 'Color', 'none');
set(ax2, 'xlim', [dnumA dnumB])
% set(ax2, 'ylim', [0.06 0.10]);
% set(ax2, 'ylim', [0.08 0.14]);
set(ax1, 'fontsize', fntsz);
set(ax2, 'fontsize', fntsz);
set(ax1, 'xtick', []);
set(ax1, 'xticklabel', []);
set(ax2, 'ylim', [ylimv1 ylimv2]);
hold all
plot(ax1, tnet_dnum, tnet_ref, '-');
plot(ax1, tnet_dnum_q, tnet_ref_q, 'k.', 'markersize', 30);
% plot(ax2, dnum_MA(~isnan(rv_MA)), rv_MA(~isnan(rv_MA)), '--.', 'color', colo('r0'), 'markersize', mrksz);
plot(ax2, sat.dnum(~isnan(rh)), rv(~isnan(rv)), '--o', 'color', colo('r0'), 'linewidth', 2, 'markersize', mrksz);
title(['\rm', domainid, ' TempNet (', refSM, ') & ',satellite, ' Reflectivity']);
ylabel(ax1, 'SM [m^3/m^3]')
ylabel('r [-]')
lg = legend(sprintf('r_V (%s)', domainid), 'location', 'northeast');
% set(lg, 'position', [lg.pos]
% set(gca, 'xlim', [datenum(2018, 3, 30), datenum(2018, 10, 30)]);
datetick('keeplimits');
zoom on; grid on;

if savefigyes; saveas(gcf, [figdir, sprintf('%s_rv_and_SM_timeseries_%s_v%s.png', satellite, domainid, curdate)]); end

%% Reflectivity H-pol and SM
figure; hold all;
set(gcf, 'position', [68         355        3208         420]);
ax1 = gca;
set(ax1, 'ylim', [0 0.55])
set(ax1, 'xlim', [dnumA dnumB])
ax2 = axes('Position', get(ax1, 'Position'), 'XAxisLocation', 'bottom', 'YAxisLocation', 'right', ...
    'Color', 'none', 'XColor', 'k', 'YColor', 'k');
set(ax1, 'Color', 'none');
set(ax2, 'xlim', [dnumA dnumB])
% set(ax2, 'ylim', [0.06 0.10]);
% set(ax2, 'ylim', [0.08 0.14]);
set(ax1, 'fontsize', fntsz);
set(ax2, 'fontsize', fntsz);
set(ax1, 'xtick', []);
set(ax1, 'xticklabel', []);
set(ax2, 'ylim', [ylimh1 ylimh2]);
hold all
plot(ax1, tnet_dnum, tnet_ref, '-');
plot(ax1, tnet_dnum_q, tnet_ref_q, 'k.', 'markersize', 30);
% plot(ax2, dnum_MA(~isnan(rv_MA)), rv_MA(~isnan(rv_MA)), '--.', 'color', colo('r0'), 'markersize', mrksz);
plot(ax2, sat.dnum(~isnan(rh)), rh(~isnan(rh)), '--o', 'color', colo('r0'), 'linewidth', 2, 'markersize', mrksz);
title(['\rm', domainid, ' TempNet (', refSM, ') & ',satellite, ' Reflectivity']);
ylabel(ax1, 'SM [m^3/m^3]')
ylabel('r [-]')
lg = legend(sprintf('r_H (%s)', domainid), 'location', 'northeast');
% set(lg, 'position', [lg.pos]
% set(gca, 'xlim', [datenum(2018, 3, 30), datenum(2018, 10, 30)]);
datetick('keeplimits');
zoom on; grid on;

if savefigyes; saveas(gcf, [figdir, sprintf('%s_rh_and_SM_timeseries_%s_v%s.png', satellite, domainid, curdate)]); end

%% reflectivity

% tnet_sm5v_max = nanmax(tnet_sm5v_q(~(isnan(rv) | isnan(rh))));
% tnet_sm5v_min = nanmin(tnet_sm5v_q(~(isnan(rv) | isnan(rh))));

tnet_ref_max = nanmax(tnet_ref_q(~bnan5v));
tnet_ref_min = nanmin(tnet_ref_q(~bnan5v));

rv_reg_max = polyval(p_m5v_V, tnet_ref_max);
rv_reg_min = polyval(p_m5v_V, tnet_ref_min);
rh_reg_max = polyval(p_m5v_H, tnet_ref_max);
rh_reg_min = polyval(p_m5v_H, tnet_ref_min);

d_rv_reg = rv_reg_max - rv_reg_min;
d_rh_reg = rh_reg_max - rh_reg_min;

figure; set(gcf, 'position', [566.6000    9.8000  493.6000  860.0000]);%[1255         435         494         860]);
subplot(2,1,1); hold all;
set(gca, 'fontsize', 13)
title([satellite, ' Reflectivity vs. ', domainid, ' TempNet (', refSM, ')'])
plot(tnet_ref_q(~bnan5v), rv(~bnan5v), 'o');
plot(tnet_ref_min, rv_reg_min, 'k+', 'markersize', 12, 'linewidth', 2);
plot(tnet_ref_max, rv_reg_max, 'k+', 'markersize', 12, 'linewidth', 2);
plot([0 0.5], polyval(p_m5v_V, [0 0.5]), '-');
text(0.35, 0.07, [sprintf('R=%.2f\n', R_ref_V), '\Delta', sprintf('r=%.3f\nSM: [%.2f %.2f]\n', d_rv_reg, tnet_ref_min, tnet_ref_max)], 'fontsize', 13);
text(0.03, 0.15, sprintf('N=%.0f\n(FF=%.1f; P=%.1f)', sum(~bnan5v), FF, P), 'fontsize', 13);
set(gca, 'ylim', [0.05 0.17]);
set(gca, 'xlim', [0 0.5]);
ylabel(sprintf('%s V-pol reflectivity [-]', satellite));
xlabel('Uncalibrated SM [m^3/m^3]');
zoom on; grid on; grid minor;

subplot(2,1,2); hold all;
set(gca, 'fontsize', 13)
plot(tnet_ref_q(~bnan5v), 1-eh(~bnan5v), 'o')
plot(tnet_ref_min, rh_reg_min, 'k+', 'markersize', 12, 'linewidth', 2);
plot(tnet_ref_max, rh_reg_max, 'k+', 'markersize', 12, 'linewidth', 2);
plot([0 0.5], polyval(p_m5v_H, [0 0.5]), '-');
% text(0.35, 0.07, sprintf('R=%.2f', R_m5v_H), 'fontsize', 13);
text(0.35, 0.07, [sprintf('R=%.2f\n', R_ref_H), '\Delta', sprintf('r=%.3f\nSM: [%.2f %.2f]\n', d_rh_reg, tnet_ref_min, tnet_ref_max)], 'fontsize', 13);
text(0.03, 0.15, sprintf('N=%.0f\n(FF=%.1f; P=%.1f)', sum(~bnan5v), FF, P), 'fontsize', 13);
set(gca, 'ylim', [0.05 0.17]);
set(gca, 'xlim', [0 0.5]);
ylabel(sprintf('%s H-pol reflectivity [-]', satellite));
xlabel('Uncalibrated SM [m^3/m^3]');
zoom on; grid on; grid minor;

% saveas(gcf, [figdir, sprintf('SMAP_r_and_insitu_5v_scatter_%s_v%s.png', domainid, curdate)]);

%%
figure; hold all;
plot(1-ev, tnet_sm10_q, 'o')








