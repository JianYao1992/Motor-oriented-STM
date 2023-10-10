%% Proportion of following memory neurons in FC

clear; clc; close all;

%% Assignment
Group = 'CtrlGroup';
Target = 'Within aAIC';
DelayPeriod = [3 4 5 6];
BeforeDelay = 2;
Context = 'Combine'; % 'hit': FC pair in hit trials; 'CR': in CR trials; 'Combine': combine pairs in hit and CR trials

%% Load FC result
load(['FCofNeurons_PropertyIndividual_1msbin' Target '_' Group '.mat']);
FcPairsEachSecond = FcPairOfEachSecond(:,DelayPeriod);

%% Preference of following neuron /// leading neurons are memory or non-memory neurons ///
PropOfMemoryNeuron = zeros(1,16);
TarFcPairs = cell(1,length(DelayPeriod));
FuPref_MemoryLu = [];
FuPref_NonmemLu = [];
if strcmp(Target,'mPFC-aAIC')
    for DelayBin = 1:length(DelayPeriod)
        switch Context
            case 'hit'
                TarFcPairs{DelayBin} = FcPairsEachSecond{DelayBin}.s1;
            case 'CR'
                TarFcPairs{DelayBin} = FcPairsEachSecond{DelayBin}.s2;
            case 'Combine'
                TarFcPairs{DelayBin} = [FcPairsEachSecond{DelayBin}.s1 FcPairsEachSecond{DelayBin}.s2];
        end
    end
    if strcmp(Context,'Combine')
        for DelayBin = 1:length(TarFcPairs)
            PairID = [];
            for iPair = 1:length(TarFcPairs{DelayBin})
                PairID = [PairID TarFcPairs{DelayBin}{iPair}.PairID];
            end
            [C,ia,~] = unique(PairID,'stable');
            TarFcPairs{DelayBin} = TarFcPairs{DelayBin}(:,ia);
        end
    end
    % preference of following neuron
    for DelayBin = 1:length(DelayPeriod)
        for iPair = 1:length(TarFcPairs{DelayBin})
            if TarFcPairs{DelayBin}{iPair}.prefered_sample_su2(BeforeDelay+DelayBin) > 0
                FuPref_MemoryLu(end+1) = TarFcPairs{DelayBin}{iPair}.prefered_sample_su1(BeforeDelay+DelayBin);
            else
                FuPref_NonmemLu(end+1) = TarFcPairs{DelayBin}{iPair}.prefered_sample_su1(BeforeDelay+DelayBin);
            end
        end
    end
elseif strcmp(Target,'aAIC-mPFC')
    for DelayBin = 1:length(DelayPeriod)
        switch Context
            case 'hit'
                TarFcPairs{DelayBin} = FcPairsEachSecond{DelayBin}.s1;
            case 'CR'
                TarFcPairs{DelayBin} = FcPairsEachSecond{DelayBin}.s2;
            case 'Combine'
                TarFcPairs{DelayBin} = [FcPairsEachSecond{DelayBin}.s1 FcPairsEachSecond{DelayBin}.s2];
        end
    end
    if strcmp(Context,'Combine')
        for DelayBin = 1:length(TarFcPairs)
            PairID = [];
            for iPair = 1:length(TarFcPairs{DelayBin})
                PairID = [PairID TarFcPairs{DelayBin}{iPair}.PairID];
            end
            [C,ia,~] = unique(PairID,'stable');
            TarFcPairs{DelayBin} = TarFcPairs{DelayBin}(:,ia);
        end
    end
    % preference of following neuron
    for DelayBin = 1:length(DelayPeriod)
        for iPair = 1:length(TarFcPairs{DelayBin})
            if TarFcPairs{DelayBin}{iPair}.prefered_sample_su1(BeforeDelay+DelayBin) > 0
                FuPref_MemoryLu(end+1) = TarFcPairs{DelayBin}{iPair}.prefered_sample_su2(BeforeDelay+DelayBin);
            else
                FuPref_NonmemLu(end+1) = TarFcPairs{DelayBin}{iPair}.prefered_sample_su2(BeforeDelay+DelayBin);
            end
        end
    end
else
    PrefEachBin = cell(1,length(DelayPeriod));
    PairID = cell(1,length(DelayPeriod));
    % construct target data
    for DelayBin = 1:length(DelayPeriod)
        switch Context
            case 'hit'
                for iPair = 1:length(FcPairsEachSecond{DelayBin}.s1)
                    if FcPairsEachSecond{DelayBin}.s1{iPair}.AIs1 > 0
                        PrefEachBin{DelayBin} = [PrefEachBin{DelayBin}; horzcat(FcPairsEachSecond{DelayBin}.s1{iPair}.prefered_sample_su1(DelayPeriod),FcPairsEachSecond{DelayBin}.s1{iPair}.prefered_sample_su2(DelayPeriod))];
                    else
                        PrefEachBin{DelayBin} = [PrefEachBin{DelayBin}; horzcat(FcPairsEachSecond{DelayBin}.s1{iPair}.prefered_sample_su2(DelayPeriod),FcPairsEachSecond{DelayBin}.s1{iPair}.prefered_sample_su1(DelayPeriod))];
                    end
                end
            case 'CR'
                for iPair = 1:length(FcPairsEachSecond{DelayBin}.s2)
                    if FcPairsEachSecond{DelayBin}.s2{iPair}.AIs2 > 0
                        PrefEachBin{DelayBin} = [PrefEachBin{DelayBin}; horzcat(FcPairsEachSecond{DelayBin}.s2{iPair}.prefered_sample_su1(DelayPeriod),FcPairsEachSecond{DelayBin}.s2{iPair}.prefered_sample_su2(DelayPeriod))];
                    else
                        PrefEachBin{DelayBin} = [PrefEachBin{DelayBin}; horzcat(FcPairsEachSecond{DelayBin}.s2{iPair}.prefered_sample_su2(DelayPeriod),FcPairsEachSecond{DelayBin}.s2{iPair}.prefered_sample_su1(DelayPeriod))];
                    end
                end
            case 'Combine'
                % Go trials
                for iPair = 1:length(FcPairsEachSecond{DelayBin}.s1)
                    if FcPairsEachSecond{DelayBin}.s1{iPair}.AIs1 > 0
                        tempPairID = FcPairsEachSecond{DelayBin}.s1{iPair}.PairID;
                        if ~ismember(tempPairID,PairID{DelayBin})
                            PairID{DelayBin} = [PairID{DelayBin} tempPairID];
                            PrefEachBin{DelayBin} = [PrefEachBin{DelayBin}; horzcat(FcPairsEachSecond{DelayBin}.s1{iPair}.prefered_sample_su1(DelayPeriod),FcPairsEachSecond{DelayBin}.s1{iPair}.prefered_sample_su2(DelayPeriod))];
                        end
                    else
                        tempPairID = FcPairsEachSecond{DelayBin}.s1{iPair}.PairID + length(TargetPairID);
                        if ~ismember(tempPairID,PairID{DelayBin})
                            PairID{DelayBin} = [PairID{DelayBin} tempPairID];
                            PrefEachBin{DelayBin} = [PrefEachBin{DelayBin}; horzcat(FcPairsEachSecond{DelayBin}.s1{iPair}.prefered_sample_su2(DelayPeriod),FcPairsEachSecond{DelayBin}.s1{iPair}.prefered_sample_su1(DelayPeriod))];
                        end
                    end
                end
                % NoGo trials
                for iPair = 1:length(FcPairsEachSecond{DelayBin}.s2)
                    if FcPairsEachSecond{DelayBin}.s2{iPair}.AIs2 > 0
                        tempPairID = FcPairsEachSecond{DelayBin}.s2{iPair}.PairID;
                        if ~ismember(tempPairID,PairID{DelayBin})
                            PairID{DelayBin} = [PairID{DelayBin} tempPairID];
                            PrefEachBin{DelayBin} = [PrefEachBin{DelayBin}; horzcat(FcPairsEachSecond{DelayBin}.s2{iPair}.prefered_sample_su1(DelayPeriod),FcPairsEachSecond{DelayBin}.s2{iPair}.prefered_sample_su2(DelayPeriod))];
                        end
                    else
                        tempPairID = FcPairsEachSecond{DelayBin}.s2{iPair}.PairID + length(TargetPairID);
                        if ~ismember(tempPairID,PairID{DelayBin})
                            PairID{DelayBin} = [PairID{DelayBin} tempPairID];
                            PrefEachBin{DelayBin} = [PrefEachBin{DelayBin}; horzcat(FcPairsEachSecond{DelayBin}.s2{iPair}.prefered_sample_su2(DelayPeriod),FcPairsEachSecond{DelayBin}.s2{iPair}.prefered_sample_su1(DelayPeriod))];
                        end
                    end
                end
        end
    end
    % preference of following neuron
    for DelayBin = 1:length(DelayPeriod)
        for iPair = 1:size(PrefEachBin{DelayBin},1)
            if PrefEachBin{DelayBin}(iPair,DelayBin) > 0
                FuPref_MemoryLu(end+1) = PrefEachBin{DelayBin}(iPair,DelayBin+length(DelayPeriod));
            else
                FuPref_NonmemLu(end+1) = PrefEachBin{DelayBin}(iPair,DelayBin+length(DelayPeriod));
            end
        end
    end
end

%% Chi-square test
PropOfMemoryNeuron = [nnz(FuPref_MemoryLu),numel(FuPref_MemoryLu),nnz(FuPref_NonmemLu),numel(FuPref_NonmemLu)];
[~,~,p] = crosstab(1:PropOfMemoryNeuron(2)+PropOfMemoryNeuron(4)>PropOfMemoryNeuron(2),[1:PropOfMemoryNeuron(2)>PropOfMemoryNeuron(1),1:PropOfMemoryNeuron(4)>PropOfMemoryNeuron(3)]);

%% Save result
save(sprintf('Proportion of following memory neurons in FC in %s trials_%s_%s_1msbin',Context,Target,Group),'PropOfMemoryNeuron','-v7.3');
close all;

%% Bar plot
figure('Color','w','Position',[100,100,235,260])
hold on
bar(1.1,PropOfMemoryNeuron(1)/PropOfMemoryNeuron(2),0.6,'EdgeColor','k','FaceColor','k');
bar(1.9,PropOfMemoryNeuron(3)/PropOfMemoryNeuron(4),0.6,'EdgeColor','k','FaceColor','w');
title(['p = ' num2str(p)]);
xlim([0.6 2.4])
set(gca,'YTick',0:0.2:1,'ylim',[0 1],'XTick',zeros(1,0));
ylabel('Proportion of memory F.U. (%)');
set(gcf,'Renderer','Painter'); saveas(gcf,['Proportion of following memory neuron in FC in ' Context ' trials_' Target '_' Group '_1msbin'],'fig'); close;