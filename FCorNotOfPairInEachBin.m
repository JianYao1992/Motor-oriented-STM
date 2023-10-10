%% Test significant functional coupling (FC) between simultaneously recorded neurons and determine the direction

clear; clc; close all;

%% Assignment
Group = 'ChR2Group';
prefix = 'TestFC';
TarBinID = {'Baseline','Sample','1stSecofDelay','2ndSecOfDelay','3rdSecOfDelay','4thSecOfDelay','Response'};
CountCriteria = 1000;
to_plot = false;
to_save = true;

%%  FC or not for all neuronal pairs over the time course
for iBin = 1:1:numel(TarBinID)
    BinRange = ['_' num2str(iBin-2) '_' num2str(iBin-1) '_1ms'];
    CurrBin = TarBinID{iBin};
    bin_range = [iBin-2,iBin-1];
    disp(['/// Test FC for each neuronal pair in ' CurrBin ' period ///']);
    Data = [];
    MatFiles = dir(['*' BinRange '*GoNoGo*.mat']);
    Data = load(MatFiles.name);
    % Student t-test
    stats = cell(0);
    thresh = norminv(0.995); % Bonferroni correction for multiple comparison
    for sidx = 1:size(Data.Sums,1) % session
        fprintf('%d of %d\n',sidx, size(Data.Sums,1));
        MiceID = regexp(Data.Sums{sidx,2},'Condition');
        MiceID = str2num(Data.Sums{sidx,2}(MiceID+11:MiceID+12));
        DayID = regexp(Data.Sums{sidx,2},'Day');
        DayID = str2num(Data.Sums{sidx,2}(DayID+3:DayID+4));
        xc_s1 = Data.Sums{sidx,6};
        xshuf_s1 = Data.Sums{sidx,7};
        xc_s2 = Data.Sums{sidx,8};
        xshuf_s2 = Data.Sums{sidx,9};
        for si = 1:(size(xc_s1.xcorr,1)-1)
            su1id = str2double(xc_s1.label{si,1});
            if ismember(su1id,Data.Sums{sidx,3})
                su1 = 'Sustained';
            elseif ismember(su1id,Data.Sums{sidx,4})
                su1 = 'Transient';
            else
                su1 = 'Nonmemory';
            end
            for sj = (si+1):size(xc_s1.xcorr,2)
                su2id = str2double(xc_s1.label{sj,1});
                if ismember(su2id,Data.Sums{sidx,3})
                    su2 = 'Sustained';
                elseif ismember(su2id,Data.Sums{sidx,4})
                    su2 = 'Transient';
                else
                    su2 = 'Nonmemory';
                end
                totalCounts1 = nansum(squeeze(xc_s1.xcorr(si,sj,:)));
                if ~isempty(xc_s2)
                    totalCounts2 = nansum(squeeze(xc_s2.xcorr(si,sj,:)));
                end
                %% Put into 'stats' if totalCount is below criteria, jumping to next iteration directly
                if (isempty(xc_s2) && totalCounts1 < CountCriteria) || (~isempty(xc_s2) && totalCounts1 < CountCriteria && totalCounts2 < CountCriteria)
                    onepair = struct();
                    onepair.fileidx = Data.Sums{sidx,1};
                    onepair.su1_label_idx = si;
                    onepair.su2_label_idx = sj;
                    onepair.totalcounts1 = totalCounts1;
                    onepair.su1_clusterid = su1id;
                    onepair.su2_clusterid = su2id;
                    onepair.wf_stats_su1 = xc_s1.label{si,2}; % cluster_id in both brain areas, FR, vollay-peak, fwhm
                    onepair.wf_stats_su2 = xc_s1.label{sj,2};
                    onepair.wf_su1 = xc_s1.label{si,3};
                    onepair.wf_su2 = xc_s1.label{sj,3};
                    onepair.prefered_sample_su1 = xc_s1.label{si,4};
                    onepair.prefered_sample_su2 = xc_s1.label{sj,4};
                    onepair.FRs1trials_su1 = xc_s1.label{si,5}; % averaged FR in context 1 of su1
                    onepair.FRs2trials_su1 = xc_s1.label{si,6}; % averaged FR in context 2 of su1
                    onepair.FRs1trials_su2 = xc_s1.label{sj,5}; % averaged FR in context 1 of su2
                    onepair.FRs2trials_su2 = xc_s1.label{sj,6}; % averaged FR in context 2 of su2
                    onepair.reg_su1 = xc_s1.label{si,7};
                    onepair.reg_su2 = xc_s1.label{sj,7};
                    onepair.s1_trials = xc_s1.cfg.trials;
                    if ~isempty(xc_s2)
                        onepair.s2_trials = xc_s2.cfg.trials;
                        onepair.totalcounts2 = totalCounts2;
                    else
                        onepair.s2_trials = 0;
                    end
                    stats{end+1} = onepair;
                    clear onepair
                    continue
                end
                
                %% Asymmetry index (AI)
                % asymmetry index for context 1
                if totalCounts1 >= CountCriteria
                    hists1 = squeeze(xc_s1.xcorr(si,sj,:));
                    shufs1 = squeeze(xshuf_s1.shiftpredictor(si,sj,:));
                    diffs1 = hists1 - smooth(shufs1);
                    diffs1(99:102) = 0;
                    stds1 = std(shufs1);
                    scores1 = diffs1(96:105)./stds1;
                    % any score > thresh
                    if any(scores1 > thresh)
                        bincounts1 = diffs1(96:105);
                        bincounts1(scores1<=thresh) = 0;
                        sumdiffs1 = (sum(bincounts1(1:5)) - sum(bincounts1(6:10)));
                        if sumdiffs1 == 0
                            AIs1 = 0;
                        else
                            AIs1 = sumdiffs1/(sum(bincounts1(1:5))+sum(bincounts1(6:end)));
                        end
                    else
                        AIs1 = 0;
                    end
                end
                % asymmetry index for context 2
                if ~isempty(xc_s2) && totalCounts2 >= CountCriteria
                    hists2 = squeeze(xc_s2.xcorr(si,sj,:));
                    shufs2 = squeeze(xshuf_s2.shiftpredictor(si,sj,:));
                    diffs2 = hists2-smooth(shufs2);
                    diffs2(99:102) = 0;
                    stds2 = std(shufs2);
                    scores2 = diffs2(96:105)./stds2;
                    % any score > thresh
                    if any(scores2 > thresh)
                        bincounts2 = diffs2(96:105);
                        bincounts2(scores2<=thresh) = 0;
                        sumdiffs2 = (sum(bincounts2(1:5))-sum(bincounts2(6:10)));
                        if sumdiffs2 == 0
                            AIs2 = 0;
                        else
                            AIs2 = sumdiffs2/(sum(bincounts2(1:5))+sum(bincounts2(6:end)));
                        end
                    else
                        AIs2 = 0;
                    end
                else
                    AIs2 = 0;
                end
                
                %% Put AI into 'stats'
                onepair = struct();
                onepair.fileidx = sidx;
                onepair.su1_label_idx = si;
                onepair.su2_label_idx = sj;
                onepair.su1_sel_type = su1;
                onepair.su2_sel_type = su2;
                onepair.totalcounts1 = totalCounts1;
                onepair.su1_clusterid = su1id;
                onepair.su2_clusterid = su2id;
                onepair.wf_stats_su1 = xc_s1.label{si,2}; % cluster_id in both areas, FR, vollay-peak, fwhm
                onepair.wf_stats_su2 = xc_s1.label{sj,2};
                onepair.wf_su1 = xc_s1.label{si,3};
                onepair.wf_su2 = xc_s1.label{sj,3};
                onepair.prefered_sample_su1 = xc_s1.label{si,4};
                onepair.prefered_sample_su2 = xc_s1.label{sj,4};
                onepair.FRs1trials_su1 = xc_s1.label{si,5}; % averaged FR in context 1 of su1
                onepair.FRs2trials_su1 = xc_s1.label{si,6}; % averaged FR in context 2 of su1
                onepair.FRs1trials_su2 = xc_s1.label{sj,5}; % averaged FR in context 1 of su2
                onepair.FRs2trials_su2 = xc_s1.label{sj,6}; % averaged FR in context 2 of su2
                onepair.reg_su1 = xc_s1.label{si,7};
                onepair.reg_su2 = xc_s1.label{sj,7};
                onepair.s1_trials = xc_s1.cfg.trials;
                if totalCounts1 >= CountCriteria % context 1
                    onepair.hists1 = hists1;
                    onepair.shufs1 = shufs1;
                    onepair.diffs1 = diffs1;
                    onepair.s1_peak_significant = any(scores1>thresh);
                    onepair.AIs1 = AIs1;
                end
                if ~isempty(xc_s2) && totalCounts2 >= CountCriteria % context 2
                    onepair.hists2 = hists2;
                    onepair.shufs2 = shufs2;
                    onepair.diffs2 = diffs2;
                    onepair.s2_peak_significant = any(scores2>thresh);
                    onepair.AIs2 = AIs2;
                end
                if ~isempty(xc_s2)
                    onepair.totalcounts2 = totalCounts2;
                    onepair.s2_trials = xc_s2.cfg.trials;
                else
                    onepair.s2_trials = 0;
                end
                stats{end+1} = onepair;
                
                %% Stem plot
                if to_plot && ((isfield(onepair,'AIs1') && abs(onepair.AIs1)>=0.9) || (isfield(onepair,'AIs2') && abs(onepair.AIs2)>=0.9)) && nnz(scores1>thresh)>=2  && onepair.prefered_sample_su1(currbin+1)>0 && onepair.prefered_sample_su2(currbin+1)>0 && onepair.prefered_sample_su1(currbin+1)==onepair.prefered_sample_su2(currbin+1)
                    fh = figure('Color','w','Position',[100,100,600,800]);
                    if isfield(onepair,'AIs1') && abs(onepair.AIs1)>=0.9 % context 1
                        % spike-pair count for real data in context 1
                        subplot(3,2,1); hold on
                        stem(xc_s1.time(50:51)*1000,onepair.hists1(50:51),'Marker','none','LineWidth',15,'Color',[0.8,0.8,0.8])
                        stem(xc_s1.time([46:49,52:55])*1000,onepair.hists1([46:49,52:55]),'Marker','none','LineWidth',15)
                        title('cross-correlogram')
                        xlabel('time (ms)')
                        ylabel('spike-pair count')
                        % spike-pair count for trial-shifted data in context 1
                        subplot(3,2,3); hold on
                        stem(xc_s1.time(50:51)*1000,onepair.shufs1(50:51),'Marker','none','LineWidth',15,'Color',[0.8,0.8,0.8])
                        stem(xc_s1.time([46:49,52:55])*1000,onepair.shufs1([46:49,52:55]),'Marker','none','LineWidth',15)
                        ssf = smooth(onepair.shufs1);
                        plot(xc_s1.time(46:55)*1000,ssf(46:55),'-r','LineWidth',1.5)
                        title('shift-predictor')
                        xlabel('time (ms)')
                        ylabel('spike-pair count')
                        % shift-predictor subtracted cross correlogram (SSCC) in context 1
                        subplot(3,2,5); hold on
                        sigbin = find((onepair.diffs1(46:55)./stds1)>thresh);
                        insigbin=find((onepair.diffs1(46:55)./stds1)<=thresh);
                        stem(xc_s1.time(insigbin+45)*1000,onepair.diffs1(45+insigbin)./stds1,'Marker','none','LineWidth',15,'Color',[0.8,0.8,0.8]);
                        stem(xc_s1.time(sigbin+45)*1000,onepair.diffs1(45+sigbin)./stds1,'Marker','none','LineWidth',15)
                        yline(thresh,'--r');
                        title('significant peak')
                        xlabel('time (ms)')
                        ylabel('z-score')
                    elseif isfield(onepair,'AIs2') && abs(onepair.AIs2)>=0.9 % context 2
                        % spike-pair count for real data in context 2
                        subplot(3,2,2); hold on
                        stem(xc_s2.time(50:51)*1000,onepair.hists2(50:51),'Marker','none','LineWidth',15,'Color',[0.8,0.8,0.8])
                        stem(xc_s2.time([46:49,52:55])*1000,onepair.hists2([46:49,52:55]),'Marker','none','LineWidth',15)
                        title('cross-correlogram')
                        xlabel('time (ms)')
                        ylabel('spike-pair count')
                        % spike-pair count for trial-shifted data in context 2
                        subplot(3,2,4); hold on
                        stem(xc_s2.time(50:51)*1000,onepair.shufs2(50:51),'Marker','none','LineWidth',15,'Color',[0.8,0.8,0.8])
                        stem(xc_s2.time([46:49,52:55])*1000,onepair.shufs2([46:49,52:55]),'Marker','none','LineWidth',15)
                        ssf = smooth(onepair.shufs2);
                        plot(xc_s2.time(46:55)*1000,ssf(46:55),'-r','LineWidth',1.5)
                        title('shift-predictor')
                        xlabel('time (ms)')
                        ylabel('spike-pair count')
                        % SSCC in context 2
                        subplot(3,2,6); hold on
                        sigbin = find((onepair.diffs2(46:55)./stds2)>thresh);
                        insigbin = find((onepair.diffs2(46:55)./stds2)<=thresh);
                        stem(xc_s2.time(insigbin+45)*1000,onepair.diffs2(45+insigbin)./stds2,'Marker','none','LineWidth',15,'Color',[0.8,0.8,0.8]);
                        stem(xc_s2.time(sigbin+45)*1000,onepair.diffs2(45+sigbin)./stds2,'Marker','none','LineWidth',15)
                        yline(thresh,'--r');
                        title('significant peak')
                        xlabel('time (ms)')
                        ylabel('z-score')
                    end
                    text(0,0.3,sprintf('%d, %d',onepair.su1_clusterid,onepair.su2_clusterid));
                    print(sprintf('xcorr_showcase_%d_%d_%d.pdf',sidx,onepair.su1_clusterid,onepair.su2_clusterid),'-dpdf','-painters');
                    print(sprintf('xcorr_showcase_%d_%d_%d.png',sidx,onepair.su1_clusterid,onepair.su2_clusterid),'-dpng','-painters');
                    close(fh);
                end
                clear onepair
            end
        end
    end
    if to_save
        fprintf('%s_XCORR_stats_%d_%d_1msbin_%s.mat\n',prefix,bin_range(1),bin_range(2),Group);
        save(sprintf('%s_XCORR_stats_%d_%d_1msbin_GoNoGo_%s.mat',prefix,bin_range(1),bin_range(2),Group),'stats','bin_range','-v7.3')
    end
end



