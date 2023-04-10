close all; clear;
%%
global HOMEDIR

%%
curdate    = datestr(now, 'yymmdd');

savepath = [HOMEDIR, '_2019 - SMAPVEX19\06 - DATA\SMOS\'];

%% Load SMOS data over the Millbrook site
smospath = [HOMEDIR, '_2019 - SMAPVEX19\06 - DATA\SMOS\FromPhilippe\'];

smosfile = 'GPrdXE_New_York_v700_20180101T094831_20221231T231907';

load([smospath, smosfile, '.mat']);

dgg      = {'D211947', 'D211948', 'D212460'};

%% Read the data

for m = 1:3

    N(m) = length(Gprd{1}.(dgg{m}).residuals);

    for n = 1:N(m)

        ssid{m}{n} = Gprd{1}.(dgg{m}).residuals{n}(:,11);

        inc{m}{n}  = Gprd{1}.(dgg{m}).residuals{n}(:,7);

        % rounded incidence angle for matching for the uniform ramge
        incq{m}{n} = round(inc{m}{n});

        TBV{m}{n}  = Gprd{1}.(dgg{m}).residuals{n}(:,15);
        TBH{m}{n}  = Gprd{1}.(dgg{m}).residuals{n}(:,14);

    end

    dnum{m}    = Gprd{1}.(dgg{m}).vd;

    % rounded timestamp for matching the time of the different DGGs
    dnumq{m}   = round(dnum{m}*10)/10;

    ad{m}      = Gprd{1}.(dgg{m}).asc;
end

%%

inc0 = [0:1:60]';

M = length(inc0);

for n = 1:N(1)

    i12 = find(dnumq{1}(n) == dnumq{2});
    i13 = find(dnumq{1}(n) == dnumq{3});

    if ~isempty(i12) && ~isempty(i13)

        for v = 1:length(inc0)

            i1 = incq{1}{n}   == inc0(v);
            i2 = incq{2}{i12} == inc0(v);
            i3 = incq{3}{i13} == inc0(v);

            TBVint{n}(v,1) = NaN;
            TBHint{n}(v,1) = NaN;

            if sum(i1) > 0 && sum(i2) > 0 && sum(i3) > 0

                tbv1             = TBV{1}{n}(i1);
                tbv2             = TBV{2}{i12}(i2);
                tbv3             = TBV{3}{i13}(i3);

                tbh1             = TBH{1}{n}(i1);
                tbh2             = TBH{2}{i12}(i2);
                tbh3             = TBH{3}{i13}(i3);

                tbv1(tbv1 > 330 | tbv1 < 50) = [];
                tbv2(tbv2 > 330 | tbv2 < 50) = [];
                tbv3(tbv3 > 330 | tbv3 < 50) = [];

                tbh1(tbh1 > 330 | tbh1 < 50) = [];
                tbh2(tbh2 > 330 | tbh2 < 50) = [];
                tbh3(tbh3 > 330 | tbh3 < 50) = [];

                if ~isempty(tbv1) && ~isempty(tbv2) && ~isempty(tbv3) && ~isempty(tbh1) && ~isempty(tbh2) && ~isempty(tbh3)

                    TBVint{n}(v,1) = (nanmean(tbv1) + nanmean(tbv2) + nanmean(tbv3))/3;
                    TBHint{n}(v,1) = (nanmean(tbh1) + nanmean(tbh2) + nanmean(tbh3))/3;

                end

            end

        end

    else

        TBVint{n} = NaN(M,1);
        TBHint{n} = NaN(M,1);

    end

end

%%

clear Gprd;

%%

dggint = 'Dinterp';

Gprd{1}.(dggint).vd  = dnum{1};
Gprd{1}.(dggint).asc = ad{1};

for n = 1:N(1)

    bnan  = isnan(TBVint{n}) | isnan(TBHint{n});

    Gprd{1}.(dggint).residuals{n}(:,7) = inc0(~bnan);

    Gprd{1}.(dggint).residuals{n}(:,15) = TBVint{n}(~bnan);
    Gprd{1}.(dggint).residuals{n}(:,14) = TBHint{n}(~bnan);

end

%%
save([savepath, smosfile, '_interpolated_v', curdate, '.mat'], "Gprd");


