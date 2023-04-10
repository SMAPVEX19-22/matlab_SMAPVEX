close all; clear;
%%
global HOMEDIR
%%

site       = 'MA'; %''MB'; % BERMS';%

savefigyes = 0;

savefile   = 0;

curdate    = datestr(now, 'yymmdd');

%%
% mininc  = 20;
% maxinc  = 50;
% anglesN = 90;

mininc  = 30;
maxinc  = 50;
anglesN = 60;

% mininc  = 35;
% maxinc  = 45;
% anglesN = 20;
% %
% mininc  = 37.5;
% maxinc  = 42.5;
% anglesN = 10;
%%

switch site

    case {'MA', 'MB'}

        sitename = 'MA';

        % FROM: /validation/tools/valtools/sandbox/SMAPVEX19/L2SMPE_over_SMAPVEX19_area_timeseries.m
        smappath = [HOMEDIR, '_2019 - SMAPVEX19\06 - DATA\SMAP\'];

        smapfile = 'SMAP_L2SMPE_SMAPVEX19_area_R1829__20190101_20221207.mat';

        x = load([smappath, smapfile]);

        c0MA = 1156;
        r0MA = 263;
        c0MB = 1140;
        r0MB = 270;

        num_files = length(x.smap_files);

        lamin = [];
        lamax = [];
        lomin = [];
        lomax = [];

        k = 0;
        for n = 1:num_files

            %             if strcmp(x.smap_AD{n}, 'D')

            k = k + 1;
            dnumall{k} = datenum(x.smap_yr{n}, x.smap_month{n}, x.smap_day{n}, x.smap_hour{n}, x.smap_minute{n}, 0);

            lamin = min([lamin, min(x.smap_lat{n})]);
            lamax = max([lamax, max(x.smap_lat{n})]);
            lomin = min([lomin, min(x.smap_lon{n})]);
            lomax = max([lomax, max(x.smap_lon{n})]);

            switch site
                case 'MA'
                    ind = x.smap_cols{n} == c0MA & x.smap_rows{n} == r0MA;
                    [i1, i2] = ind2sub(size(x.smap_cols{n}), find(ind));
                case 'MB'
                    ind = x.smap_cols{n} == c0MB & x.smap_rows{n} == r0MB;
                    [i1, i2] = ind2sub(size(x.smap_cols{n}), find(ind));
            end

            if any(ind(:)) & ~isnan(dnumall{k}(i2))
                s.smap_AD(k,1)     = x.smap_AD{n};
                s.smap_scah(k,1)   = x.smap_scah{n}(ind);
                s.smap_scav(k,1)   = x.smap_scav{n}(ind);
                s.smap_dca(k,1)    = x.smap_dca{n}(ind);
                s.smap_tauh(k,1)   = x.tau_scah{n}(ind);
                s.smap_tauv(k,1)   = x.tau_scav{n}(ind);
                s.smap_taudca(k,1) = x.tau_dca{n}(ind);
                s.smap_tbv(k,1)    = x.smap_tbv{n}(ind);
                s.smap_tbh(k,1)    = x.smap_tbh{n}(ind);
                s.smap_st(k,1)     = x.smap_st{n}(ind);
                s.smap_doy(k,1)    = x.smap_doy{n}(i2);
                s.smap_dnum(k,1)   = dnumall{k}(i2);
            else
                s.smap_AD(k,1) = ' ';
                s.smap_scah(k,1)   = NaN;
                s.smap_scav(k,1)   = NaN;
                s.smap_dca(k,1)    = NaN;
                s.smap_tauh(k,1)   = NaN;
                s.smap_tauv(k,1)   = NaN;
                s.smap_taudca(k,1) = NaN;
                s.smap_tbv(k,1)    = NaN;
                s.smap_tbh(k,1)    = NaN;
                s.smap_st(k,1)     = NaN;
                s.smap_doy(k,1)    = NaN;
                s.smap_dnum(k,1)   = NaN;
            end

            %             end

        end

        % SMOS files
        if strcmp(site, 'MA')

            smospath = [HOMEDIR, '_2019 - SMAPVEX19\06 - DATA\SMOS\FromPhilippe\'];

            smosfile = 'GPrdXE_Massachusetts_v700_20180101T094831_20221231T231907.mat';
            dgg      = 'D207855';

        elseif strcmp(site, 'MB')

            % because of interpolation
            anglesN = maxinc - mininc - 2;

            smospath = [HOMEDIR, '_2019 - SMAPVEX19\06 - DATA\SMOS\'];

            smosfile = 'GPrdXE_New_York_v700_20180101T094831_20221231T231907_interpolated_v230407.mat';
            dgg      = 'Dinterp';

        end
        %         scah(scah < 0) = NaN;
        %         scav(scav < 0) = NaN;
        %         dca(dca < 0)   = NaN;

    case 'BERMS'

        sitename = 'BERMS';

        smappath = [HOMEDIR, '_2021 - SMAPVEX boreal\03 - DATA\SMAP\'];

        smapfile = 'SMAP_L2SMPE_BERMS_pixel_R18290_20180101_20221231.mat';

        s = load([smappath, smapfile]);

        % SMOS files
        smospath = [HOMEDIR, '_2019 - SMAPVEX19\06 - DATA\SMOS\FromPhilippe\'];

        smosfile = 'GPrdXE_BERMS_v700_20180101T094831_20221231T231907.mat';
        dgg      = 'D139040';

end
%% Load SMOS data

load([smospath, smosfile]);

smos_dnum = Gprd{1}.(dgg).vd;
smos_ad   = Gprd{1}.(dgg).asc;

%%

dnum_1 = max([min(smos_dnum), min(s.smap_dnum)]);
dnum_2 = min([max(smos_dnum), max(s.smap_dnum)]);


%% filter SMOS data
for n = 1:length(Gprd{1}.(dgg).residuals)

    inc{n}  = Gprd{1}.(dgg).residuals{n}(:,7);

    TBV{n}  = Gprd{1}.(dgg).residuals{n}(:,15);
    TBH{n}  = Gprd{1}.(dgg).residuals{n}(:,14);

    bang{n} = inc{n} > mininc & inc{n} < maxinc;

    if sum(bang{n}) >= anglesN & all(TBH{n}(bang{n}) < 330) & all(TBH{n}(bang{n}) > 0) & all(TBV{n}(bang{n}) < 330) & all(TBV{n}(bang{n}) > 0) ...
            & smos_dnum(n) > dnum_1 & smos_dnum(n) < dnum_2

        if maxinc - mininc > 7
            pv(n,:)    = polyfit(inc{n}(bang{n}), TBV{n}(bang{n}), 1);
            ph(n,:)    = polyfit(inc{n}(bang{n}), TBH{n}(bang{n}), 1);

            TBV40(n,1) = polyval(pv(n, :), 40);
            TBH40(n,1) = polyval(ph(n, :), 40);
        else
            TBV40(n,1)  = nanmean(TBV{n}(bang{n}));
            TBH40(n,1)  = nanmean(TBH{n}(bang{n}));
        end

        Nangles1(n,1) = sum(bang{n});

    else

        pv(n,1:2)       = NaN;
        ph(n,1:2)       = NaN;

        TBV40(n,1)    = NaN;
        TBH40(n,1)    = NaN;

        Nangles1(n,1) = NaN;
    end

end

%%
N_valid_smos  = sum(~isnan(TBV40));
N_valid_smap  = sum(~isnan(s.smap_tbv));

N_mean_angles = nanmean(Nangles1);

%%
TBV40_am           = TBV40(smos_ad == 1);
TBH40_am           = TBH40(smos_ad == 1);
TBV40_pm           = TBV40(smos_ad == 0);
TBH40_pm           = TBH40(smos_ad == 0);

smap_tbv_am        = s.smap_tbv(s.smap_AD == 'D');
smap_tbh_am        = s.smap_tbh(s.smap_AD == 'D');
smap_tbv_pm        = s.smap_tbv(s.smap_AD == 'A');
smap_tbh_pm        = s.smap_tbh(s.smap_AD == 'A');

dnum_smos_local_am = smos_dnum(smos_ad == 1) - 6/24;
dnum_smap_local_am = s.smap_dnum(s.smap_AD == 'D') - 6/24;

dnum_smos_local_pm = smos_dnum(smos_ad == 0) - 6/24;
dnum_smap_local_pm = s.smap_dnum(s.smap_AD == 'A') - 6/24;

dnum_smos_day_am   = floor(dnum_smos_local_am);
dnum_smap_day_am   = floor(dnum_smap_local_am);
dnum_smos_day_pm   = floor(dnum_smos_local_pm);
dnum_smap_day_pm   = floor(dnum_smap_local_pm);

%% matchup SMAP and SMOS

dnums = [floor(min(smos_dnum)):floor(max(smos_dnum))];

smos_tbv_comp = NaN(length(dnums), 2);
smos_tbh_comp = NaN(length(dnums), 2);
smap_tbv_comp = NaN(length(dnums), 2);
smap_tbh_comp = NaN(length(dnums), 2);

% AM

for ij = 1:length(dnums)

    bm_smos_am = dnums(ij) == dnum_smos_day_am;
    bm_smap_am = dnums(ij) == dnum_smap_day_am;

    if sum(bm_smos_am) > 0 && sum(bm_smap_am) > 0 && ~isnan(TBV40_am(bm_smos_am)) && ~isnan(TBH40_am(bm_smos_am)) && ~isnan(smap_tbv_am(bm_smap_am)) && ~isnan(smap_tbh_am(bm_smap_am))

        smos_tbv_comp(ij,1) = TBV40_am(bm_smos_am);
        smos_tbh_comp(ij,1) = TBH40_am(bm_smos_am);

        smap_tbv_comp(ij,1) = smap_tbv_am(bm_smap_am);
        smap_tbh_comp(ij,1) = smap_tbh_am(bm_smap_am);

    end

    bm_smos_pm = dnums(ij) == dnum_smos_day_pm;
    bm_smap_pm = dnums(ij) == dnum_smap_day_pm;

    if sum(bm_smos_pm) > 0 && sum(bm_smap_pm) > 0 && ~isnan(TBV40_pm(bm_smos_pm)) && ~isnan(TBH40_pm(bm_smos_pm)) && ~isnan(smap_tbv_pm(bm_smap_pm)) && ~isnan(smap_tbh_pm(bm_smap_pm))

        smos_tbv_comp(ij,2) = TBV40_pm(bm_smos_pm);
        smos_tbh_comp(ij,2) = TBH40_pm(bm_smos_pm);

        smap_tbv_comp(ij,2) = smap_tbv_pm(bm_smap_pm);
        smap_tbh_comp(ij,2) = smap_tbh_pm(bm_smap_pm);

    end

end

%% comparison metrics: AM and PM separately
md_v    = nanmean(smos_tbv_comp-smap_tbv_comp);
md_h    = nanmean(smos_tbh_comp-smap_tbh_comp);

stdev_v = nanstd(smap_tbv_comp-smos_tbv_comp);
stdev_h = nanstd(smap_tbh_comp-smos_tbh_comp);

rmsd_v  = sqrt(nanmean((smap_tbv_comp-smos_tbv_comp).^2));
rmsd_h  = sqrt(nanmean((smap_tbh_comp-smos_tbh_comp).^2));

bnan_v = isnan(smap_tbv_comp) | isnan(smos_tbv_comp);
bnan_h = isnan(smap_tbh_comp) | isnan(smos_tbh_comp);

for m = 1:2
    R_v(:,m) = corr(smap_tbv_comp(~bnan_v(:,m),m), smos_tbv_comp(~bnan_v(:,m),m));
    R_h(:,m) = corr(smap_tbh_comp(~bnan_h(:,m),m), smos_tbh_comp(~bnan_h(:,m),m));

    p_v_ss(m,:) = polyfit(smap_tbv_comp(~bnan_v(:,m),m), smos_tbv_comp(~bnan_v(:,m),m), 1);
    p_h_ss(m,:) = polyfit(smap_tbh_comp(~bnan_h(:,m),m), smos_tbh_comp(~bnan_h(:,m),m), 1);
end

N_v = sum(~bnan_v);
N_h = sum(~bnan_h);

%% comparison metrics: AM and PM combined

md_v_all    = nanmean(smos_tbv_comp(:) - smap_tbv_comp(:));
md_h_all    = nanmean(smos_tbh_comp(:) - smap_tbh_comp(:));

stdev_v_all = nanstd(smap_tbv_comp(:) - smos_tbv_comp(:));
stdev_h_all = nanstd(smap_tbh_comp(:) - smos_tbh_comp(:));

rmsd_v_all  = sqrt(nanmean((smap_tbv_comp(:) - smos_tbv_comp(:)).^2));
rmsd_h_all  = sqrt(nanmean((smap_tbh_comp(:) - smos_tbh_comp(:)).^2));

R_v_all = corr([smap_tbv_comp(~bnan_v(:,1),1); smap_tbv_comp(~bnan_v(:,2),2)], [smos_tbv_comp(~bnan_v(:,1),1); smos_tbv_comp(~bnan_v(:,2),2)]);
R_h_all = corr([smap_tbh_comp(~bnan_h(:,1),1); smap_tbh_comp(~bnan_h(:,2),2)], [smos_tbh_comp(~bnan_h(:,1),1); smos_tbh_comp(~bnan_h(:,2),2)]);

p_v_all = polyfit([smap_tbv_comp(~bnan_v(:,1),1); smap_tbv_comp(~bnan_v(:,2),2)], [smos_tbv_comp(~bnan_v(:,1),1); smos_tbv_comp(~bnan_v(:,2),2)], 1)
p_h_all = polyfit([smap_tbh_comp(~bnan_h(:,1),1); smap_tbh_comp(~bnan_h(:,2),2)], [smos_tbh_comp(~bnan_h(:,1),1); smos_tbh_comp(~bnan_h(:,2),2)], 1)

N_v_all = sum(~bnan_v(:));
N_h_all = sum(~bnan_h(:));

%% Adjust SMOS 

smos_tbv_comp_adj = (smos_tbv_comp - p_v_all(2))/p_v_all(1);
smos_tbh_comp_adj = (smos_tbh_comp - p_h_all(2))/p_h_all(1);

md_v_adj   = nanmean(smos_tbv_comp_adj(:) - smap_tbv_comp(:))
md_h_adj   = nanmean(smos_tbh_comp_adj(:) - smap_tbh_comp(:))



%% Save SMAP and SMOS parameters
% SMAP parameters
smap.tbv      = s.smap_tbv;
smap.tbh      = s.smap_tbh;
smap.ad       = s.smap_AD; 
smap.dnum     = s.smap_dnum;
smap.doy      = s.smap_doy;
smap.sm_scah  = s.smap_scah;
smap.sm_scav  = s.smap_scav;
smap.sm_dca   = s.smap_dca;
smap.tau_scah = s.smap_tauh;
smap.tau_scav = s.smap_tauv;
smap.tau_dca  = s.smap_taudca;
smap.tsurf    = s.smap_st;      
            
% SMOS parameters
smos.tbv     = TBV40;
smos.tbh     = TBH40;
smos.dnum    = smos_dnum;
smos.mininc  = mininc;
smos.maxinc  = maxinc;
smos.anglesN = anglesN;
smos.p_v     = p_v_all;
smos.p_h     = p_h_all;

% reformat ASC/DESC
smos.ad(smos_ad == 1) = 'A';
smos.ad(smos_ad == 0) = 'D';
smos.ad               =  smos.ad';

% resolve DOY
smosyrs  = datevec(smos.dnum);
smos.doy = smos.dnum - datenum(smosyrs(:,1),1,0); % Jan 1 is day 1


fname = sprintf('SMAP_and_SMOS_TB_%s_%.0f_%.0f_v%s', site, mininc, maxinc, curdate);

if savefile
    save(fname, 'smap', 'smos');
end

%%
sel = 6;

figure; hold all
plot(inc{sel}, TBV{sel}, 's')
plot(inc{sel}, TBH{sel}, 'o')
plot([mininc:maxinc], polyval(pv(sel, :), [mininc:maxinc]));
plot([mininc:maxinc], polyval(ph(sel, :), [mininc:maxinc]));
zoom on; grid on;

%%

figure; hold all
plot(smos_dnum(smos_ad == 1)-smos_dnum(1), TBV40(smos_ad == 1), '+')
% plot(dnum(ad == 1)-dnum(1), TBH40(ad == 1), 'o')
plot(s.smap_dnum(s.smap_AD == 'D')-smos_dnum(1), s.smap_tbv(s.smap_AD == 'D'), 'x')
legend('SMOS', 'SMAP')
% datetick;
zoom on; grid on;

%%
m = 1;
figure; hold all
plot(dnums-smos_dnum(1), smos_tbv_comp(:,m), '+')
plot(dnums-smos_dnum(1), smap_tbv_comp(:,m), 'x')
plot(dnums-smos_dnum(1), smos_tbv_comp_adj(:,m), 'o')
legend('SMOS', 'SMAP')
% datetick;
zoom on; grid on;

%%
minP = 215;
maxP = 285;

fsz  = 14;

AMPM = {'AM', 'PM'};

figure; hold all;
set(gcf, 'Position', 1.0e+03*[0.1810    0.1370    1.1584    1.1176]);
for m = 1:2

    subplot(2,2,m*2-1); hold all;

    set(gca,'fontsize', fsz);
    plot([0,300], [0, 300]);
    scatter(smap_tbv_comp(:,m), smos_tbv_comp(:,m), 'MarkerFaceColor','black', 'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.2);
    plot([minP maxP], polyval(p_v_ss(m,:), [minP maxP]), 'r')
    xlim([minP maxP])
    ylim([minP maxP])
    text(minP+2, 282, sprintf('Pairs N=%.0f', N_v(m)), 'fontsize', fsz);
    text(270, minP+10, sprintf('STDev=%.2f K\nMeanD=%.2f K\nRMSD=%.2f K\nR=%.2f', stdev_v(m), md_v(m), rmsd_v(m), R_v(m)), 'fontsize', fsz);
    text(250, 280, sprintf('y = %.3fx + %.1f', p_v_ss(m,1),  p_v_ss(m,2)), 'color', 'r', 'FontSize', 12);
    xlabel('SMAP [K]')
    ylabel('SMOS [K]')
    title(sprintf('V-pol (%s)', AMPM{m}))
    axis square
    zoom on; grid on;


    %
    subplot(2,2,m*2); hold all;
    set(gca,'fontsize', fsz);
    plot([0,300], [0, 300]);
    scatter(smap_tbh_comp(:,m), smos_tbh_comp(:,m), 'MarkerFaceColor','black', 'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.2);
    plot([minP maxP], polyval(p_h_ss(m,:), [minP maxP]), 'r')
    xlim([minP maxP])
    ylim([minP maxP])
    text(minP+2, 282, sprintf('Pairs N=%.0f', N_v(m)), 'fontsize', fsz);
    text(270, minP+10, sprintf('STDev=%.2f K\nMeanD=%.2f K\nRMSD=%.2f K\nR=%.2f', stdev_h(m), md_h(m), rmsd_h(m), R_h(m)), 'fontsize', fsz);
    text(250, 280, sprintf('y = %.3fx + %.1f', p_h_ss(m,1),  p_h_ss(m,2)), 'color', 'r', 'FontSize', 12);
    xlabel('SMAP [K]')
    ylabel('SMOS [K]')
    title(sprintf('H-pol (%s)', AMPM{m}))
    axis square
    zoom on; grid on;
    if m == 1
        text(210, 294, {sprintf('%s - SMOS Incidence Angle Range: %.1f - %.1f with minimum of %.0f measurements each (mean %.0f)', sitename, mininc, maxinc, anglesN, N_mean_angles), ...
            sprintf('Total Valid N: SMOS=%.0f; SMAP=%.0f (%s - %s)', N_valid_smos, N_valid_smap, datestr(dnum_1, 'mmm yyyy'), datestr(dnum_2, 'mmm yyyy'))}, 'fontsize', fsz, 'HorizontalAlignment','center');
    end

end

savename = sprintf('SMAP_SMOS_scatter_plots_%s_%.0f_%.0f_%.0f.png', site, anglesN, mininc, maxinc);

if savefigyes; saveas(gcf, savename); end

%% AM and PM together
figure; hold all
set(gcf, 'Position', 1.0e+03*[0.1810    0.1370    1.1584    0.6]);

subplot(1,2,1); hold all;
set(gca,'fontsize', fsz);
plot([0,300], [0, 300]);
scatter([smap_tbv_comp(:,1); smap_tbv_comp(:,2)], [smos_tbv_comp(:,1); smos_tbv_comp(:,2)], 'MarkerFaceColor','black', 'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.2);
plot([minP maxP], polyval(p_v_all, [minP maxP]), 'r')
xlim([minP maxP])
ylim([minP maxP])
text(minP+2, 282, sprintf('Pairs N=%.0f', N_h_all), 'fontsize', fsz);
text(270, minP+10, sprintf('STDev=%.2f K\nMeanD=%.2f K\nRMSD=%.2f K\nR=%.2f', stdev_h_all, md_h_all, rmsd_h_all, R_h_all), 'fontsize', fsz);
text(250, 280, sprintf('y = %.3fx + %.1f', p_v_all(1),  p_v_all(2)), 'color', 'r', 'FontSize', 12);
xlabel('SMAP [K]')
ylabel('SMOS [K]')
title(sprintf('V-pol'))
axis square
zoom on; grid on;


subplot(1,2,2); hold all;
set(gca,'fontsize', fsz);
plot([0,300], [0, 300]);
scatter([smap_tbh_comp(:,1); smap_tbh_comp(:,2)], [smos_tbh_comp(:,1); smos_tbh_comp(:,2)], 'MarkerFaceColor','black', 'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.2);
plot([minP maxP], polyval(p_h_all, [minP maxP]), 'r')
xlim([minP maxP])
ylim([minP maxP])
text(minP+2, 282, sprintf('Pairs N=%.0f', N_v_all), 'fontsize', fsz);
text(270, minP+10, sprintf('STDev=%.2f K\nMeanD=%.2f K\nRMSD=%.2f K\nR=%.2f', stdev_v_all, md_v_all, rmsd_v_all, R_v_all), 'fontsize', fsz);
text(250, 280, sprintf('y = %.3fx + %.1f', p_h_all(1),  p_h_all(2)), 'color', 'r', 'FontSize', 12);
xlabel('SMAP [K]')
ylabel('SMOS [K]')
title(sprintf('H-pol'))
axis square
zoom on; grid on;
text(210, 294, {sprintf('%s - SMOS Incidence Angle Range: %.1f - %.1f with minimum of %.0f measurements each (mean %.0f)', sitename, mininc, maxinc, anglesN, N_mean_angles), ...
    sprintf('Total Valid N: SMOS=%.0f; SMAP=%.0f (%s - %s)', N_valid_smos, N_valid_smap, datestr(dnum_1, 'mmm yyyy'), datestr(dnum_2, 'mmm yyyy'))}, 'fontsize', fsz, 'HorizontalAlignment','center');
