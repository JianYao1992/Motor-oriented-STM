%% Heatmap of FC between memory neurons

clear; clc; close all;

%% Assignment
Target = 'Local+Cross-region'; % 'Within mPFC','Within aAIC', 'mPFC-aAIC'
Group = 'CtrlGroup';
Neurontype = 'Within memory neurons';
WhichTrials = 'hit trials';
BinRange = 1:1:4;
if strcmp(Target,'Within mPFC')
    tarreg1 = 1; tarreg2 = 1;
elseif strcmp(Target,'Within aAIC')
    tarreg1 = 2; tarreg2 = 2;
elseif strcmp(Target,'mPFC-aAIC') || strcmp(Target,'Local+Cross-region')
    tarreg1 = 1; tarreg2 = 2;
end
TotalBinNum = 7;

%% FC pattern
FCinEachBin_mPFC = [];
FCinEachBin_aAIC = [];
FCinEachBin_Proj = [];
FCinEachBin_Tar = [];
for bin = 1:length(BinRange)
    load(['FCofNeurons_PropertyPopulation_1msbin_' num2str(BinRange(bin)) '_' num2str(BinRange(bin)+1) '_' Group '.mat']);
    for iPair = 1:size(conn_chain_go,1) % hit
        if conn_chain_go(iPair,3) < 0
            temp = conn_chain_go(iPair,1); conn_chain_go(iPair,1) = conn_chain_go(iPair,2); conn_chain_go(iPair,2) = temp;
            conn_chain_go(iPair,3) = -1*conn_chain_go(iPair,3);
            temp = pref_chain_go(iPair,1:7); pref_chain_go(iPair,1:7) = pref_chain_go(iPair,8:14); pref_chain_go(iPair,8:14) = temp;
            temp = reg_chain_go(iPair,1); reg_chain_go(iPair,1) = reg_chain_go(iPair,2); reg_chain_go(iPair,2) = temp;
            reg_chain_go(iPair,3) = -1*reg_chain_go(iPair,3);
        end
    end
    for iPair = 1:size(conn_chain_nogo,1) % CR
        if conn_chain_nogo(iPair,3) < 0
            temp = conn_chain_nogo(iPair,1); conn_chain_nogo(iPair,1) = conn_chain_nogo(iPair,2); conn_chain_nogo(iPair,2) = temp;
            conn_chain_nogo(iPair,3) = -1*conn_chain_nogo(iPair,3);
            temp = pref_chain_nogo(iPair,1:7); pref_chain_nogo(iPair,1:7) = pref_chain_nogo(iPair,8:14); pref_chain_nogo(iPair,8:14) = temp;
            temp = reg_chain_nogo(iPair,1); reg_chain_nogo(iPair,1) = reg_chain_nogo(iPair,2); reg_chain_nogo(iPair,2) = temp;
            reg_chain_nogo(iPair,3) = -1*reg_chain_nogo(iPair,3);
        end
    end
    if strcmp(Neurontype,'Within memory neurons')
        reg_chain_mPFC_Hit = reg_chain_go(reg_chain_go(:,1)==1 & reg_chain_go(:,2)==1 & ~all(pref_chain_go(:,3:6)==0,2) & ~all(pref_chain_go(:,10:13)==0,2),:); % mPFC
        conn_chain_mPFC_Hit = conn_chain_go(reg_chain_go(:,1)==1 & reg_chain_go(:,2)==1 & ~all(pref_chain_go(:,3:6)==0,2) & ~all(pref_chain_go(:,10:13)==0,2),:);
        pref_chain_mPFC_Hit = pref_chain_go(reg_chain_go(:,1)==1 & reg_chain_go(:,2)==1 & ~all(pref_chain_go(:,3:6)==0,2) & ~all(pref_chain_go(:,10:13)==0,2),:);
        mouse_chain_mPFC_Hit = mouse_chain_go(reg_chain_go(:,1)==1 & reg_chain_go(:,2)==1 & ~all(pref_chain_go(:,3:6)==0,2) & ~all(pref_chain_go(:,10:13)==0,2),:);
        reg_chain_aAIC_Hit = reg_chain_go(reg_chain_go(:,1)==2 & reg_chain_go(:,2)==2 & ~all(pref_chain_go(:,3:6)==0,2) & ~all(pref_chain_go(:,10:13)==0,2),:); % aAIC
        conn_chain_aAIC_Hit = conn_chain_go(reg_chain_go(:,1)==2 & reg_chain_go(:,2)==2 & ~all(pref_chain_go(:,3:6)==0,2) & ~all(pref_chain_go(:,10:13)==0,2),:);
        pref_chain_aAIC_Hit = pref_chain_go(reg_chain_go(:,1)==2 & reg_chain_go(:,2)==2 & ~all(pref_chain_go(:,3:6)==0,2) & ~all(pref_chain_go(:,10:13)==0,2),:);
        mouse_chain_aAIC_Hit = mouse_chain_go(reg_chain_go(:,1)==2 & reg_chain_go(:,2)==2 & ~all(pref_chain_go(:,3:6)==0,2) & ~all(pref_chain_go(:,10:13)==0,2),:);
        reg_chain_Proj_Hit = reg_chain_go(reg_chain_go(:,1)==1 & reg_chain_go(:,2)==2 & ~all(pref_chain_go(:,3:6)==0,2) & ~all(pref_chain_go(:,10:13)==0,2),:); % mPFC-aAIC
        conn_chain_Proj_Hit = conn_chain_go(reg_chain_go(:,1)==1 & reg_chain_go(:,2)==2 & ~all(pref_chain_go(:,3:6)==0,2) & ~all(pref_chain_go(:,10:13)==0,2),:);
        pref_chain_Proj_Hit = pref_chain_go(reg_chain_go(:,1)==1 & reg_chain_go(:,2)==2 & ~all(pref_chain_go(:,3:6)==0,2) & ~all(pref_chain_go(:,10:13)==0,2),:);
        mouse_chain_Proj_Hit = mouse_chain_go(reg_chain_go(:,1)==1 & reg_chain_go(:,2)==2 & ~all(pref_chain_go(:,3:6)==0,2) & ~all(pref_chain_go(:,10:13)==0,2),:);
    elseif strcmp(Neurontype,'all')
        if strcmp(Target,'Within aAIC') || strcmp(Target,'Within mPFC') || strcmp(Target,'mPFC-aAIC')
            reg_chain_Tar_Hit = reg_chain_go(reg_chain_go(:,1) == tarreg1 & reg_chain_go(:,2) == tarreg2,:);
            conn_chain_Tar_Hit = conn_chain_go(reg_chain_go(:,1) == tarreg1 & reg_chain_go(:,2) == tarreg2,:);
            pref_chain_Tar_Hit = pref_chain_go(reg_chain_go(:,1) == tarreg1 & reg_chain_go(:,2) == tarreg2,:);
            mouse_chain_Tar_Hit = mouse_chain_go(reg_chain_go(:,1) == tarreg1 & reg_chain_go(:,2) == tarreg2,:);
        elseif strcmp(Target,'Local+Cross-region')
            % mPFC
            reg_chain_mPFC_Hit = reg_chain_go(reg_chain_go(:,1) == 1 & reg_chain_go(:,2) == 1,:);
            conn_chain_mPFC_Hit = conn_chain_go(reg_chain_go(:,1) == 1 & reg_chain_go(:,2) == 1,:);
            pref_chain_mPFC_Hit = pref_chain_go(reg_chain_go(:,1) == 1 & reg_chain_go(:,2) == 1,:);
            mouse_chain_mPFC_Hit = mouse_chain_go(reg_chain_go(:,1) == 1 & reg_chain_go(:,2) == 1,:);
            % aAIC
            reg_chain_aAIC_Hit = reg_chain_go(reg_chain_go(:,1) == 2 & reg_chain_go(:,2) == 2,:);
            conn_chain_aAIC_Hit = conn_chain_go(reg_chain_go(:,1) == 2 & reg_chain_go(:,2) == 2,:);
            pref_chain_aAIC_Hit = pref_chain_go(reg_chain_go(:,1) == 2 & reg_chain_go(:,2) == 2,:);
            mouse_chain_aAIC_Hit = mouse_chain_go(reg_chain_go(:,1) == 2 & reg_chain_go(:,2) == 2,:);
            % mPFC-aAIC
            reg_chain_Proj_Hit = reg_chain_go(reg_chain_go(:,1) == tarreg1 & reg_chain_go(:,2) == tarreg2,:);
            conn_chain_Proj_Hit = conn_chain_go(reg_chain_go(:,1) == tarreg1 & reg_chain_go(:,2) == tarreg2,:);
            pref_chain_Proj_Hit = pref_chain_go(reg_chain_go(:,1) == tarreg1 & reg_chain_go(:,2) == tarreg2,:);
            mouse_chain_Proj_Hit = mouse_chain_go(reg_chain_go(:,1) == tarreg1 & reg_chain_go(:,2) == tarreg2,:);
        end
    end
    selecbin = bin;
    if strcmp(WhichTrials,'hit trials')
        % mPFC
        if exist('conn_chain_mPFC_Hit')
            conn_chain_mPFC = conn_chain_mPFC_Hit;
            pref_chain_mPFC = pref_chain_mPFC_Hit;
            mouse_chain_mPFC = mouse_chain_mPFC_Hit;
        end
        % aAIC
        if exist('conn_chain_aAIC_Hit')
            conn_chain_aAIC = conn_chain_aAIC_Hit;
            pref_chain_aAIC = pref_chain_aAIC_Hit;
            mouse_chain_aAIC = mouse_chain_aAIC_Hit;
        end
        % mPFC-aAIC
        if exist('conn_chain_Proj_Hit')
            conn_chain_Proj = conn_chain_Proj_Hit;
            pref_chain_Proj = pref_chain_Proj_Hit;
            mouse_chain_Proj = mouse_chain_Proj_Hit;
        end
        % target
        if exist('reg_chain_Tar_Hit')
            conn_chain_Tar = conn_chain_Tar_Hit;
            pref_chain_Tar = pref_chain_Tar_Hit;
            mouse_chain_Tar = mouse_chain_Tar_Hit;
        end
    else
        % mPFC
        if exist('conn_chain_mPFC_CR')
            conn_chain_mPFC = conn_chain_mPFC_CR;
            pref_chain_mPFC = pref_chain_mPFC_CR;
            mouse_chain_mPFC = mouse_chain_mPFC_CR;
        end
        % aAIC
        if exist('conn_chain_aAIC_CR')
            conn_chain_aAIC = conn_chain_aAIC_CR;
            pref_chain_aAIC = pref_chain_aAIC_CR;
            mouse_chain_aAIC = mouse_chain_aAIC_CR;
        end
        % mPFC-aAIC
        if exist('conn_chain_Proj_CR')
            conn_chain_Proj = conn_chain_proj_CR;
            pref_chain_Proj = pref_chain_proj_CR;
            mouse_chain_Proj = mouse_chain_proj_CR;
        end
        % target
        if exist('conn_chain_Tar_CR')
            conn_chain_Tar = conn_chain_tar_CR;
            pref_chain_Tar = pref_chain_tar_CR;
            mouse_chain_Tar = mouse_chain_tar_CR;
        end
    end
    % mPFC
    if exist('conn_chain_mPFC')
        tempCouplinginEachBin_mPFC = ConstructFcMatrix(conn_chain_mPFC,pref_chain_mPFC,mouse_chain_mPFC,selecbin,TotalBinNum,FCinEachBin_mPFC);
        FCinEachBin_mPFC = tempCouplinginEachBin_mPFC;
    end
    % aAIC
    if exist('conn_chain_aAIC')
        tempCouplinginEachBin_aAIC = ConstructFcMatrix(conn_chain_aAIC,pref_chain_aAIC,mouse_chain_aAIC,selecbin,TotalBinNum,FCinEachBin_aAIC);
        FCinEachBin_aAIC = tempCouplinginEachBin_aAIC;
    end
    % mPFC-aAIC
    if exist('conn_chain_Proj')
        tempCouplinginEachBin_proj = ConstructFcMatrix(conn_chain_Proj,pref_chain_Proj,mouse_chain_Proj,selecbin,TotalBinNum,FCinEachBin_Proj);
        FCinEachBin_Proj = tempCouplinginEachBin_proj;
    end
    % target
    if exist('conn_chain_Tar')
        tempCouplinginEachBin_Tar = ConstructFcMatrix(conn_chain_Tar,pref_chain_Tar,mouse_chain_Tar,selecbin,TotalBinNum,FCinEachBin_Tar);
        FCinEachBin_Tar = tempCouplinginEachBin_Tar;
    end
end
if strcmp(Target,'Local+Cross-region')
    FCinEachBin = vertcat(FCinEachBin_mPFC, FCinEachBin_Proj, FCinEachBin_aAIC);
    save(sprintf('FcTemporalPattern_1msbin_%s_%s_%s_%s.mat',Target,Neurontype,Group,WhichTrials),'FCinEachBin','FCinEachBin_mPFC','FCinEachBin_aAIC','FCinEachBin_Proj','-v7.3');
elseif strcmp(Target,'Within aAIC') || strcmp(Target,'Within mPFC') || strcmp(Target,'mPFC-aAIC')
    save(sprintf('FcTemporalPattern_1msbin_%s_%s_%s_%s.mat',Target,Neurontype,Group,WhichTrials),'FCinEachBin_Tar','-v7.3');
end

%% Heatmap plot
figure('Color','w','Position',[200,200,400,700])
% mPFC
if exist('conn_chain_mPFC')
    sumI_mPFC = SortFcPairs(FCinEachBin_mPFC);
end
% aAIC
if exist('conn_chain_aAIC')
    sumI_aAIC = SortFcPairs(FCinEachBin_aAIC);
end
% mPFC-to-aAIC
if exist('conn_chain_Proj')
    sumI_Proj = SortFcPairs(FCinEachBin_Proj);
end
% target
if exist('conn_chain_Tar')
    sumI_Tar = SortFcPairs(FCinEachBin_Tar);
end
if exist('conn_chain_mPFC')
    FCinEachBin_mPFC = FCinEachBin_mPFC(sumI_mPFC,3:end-1); % mPFC
    FCinEachBin_Proj = FCinEachBin_Proj(sumI_Proj,3:end-1); % mPFC-aAIC
    FCinEachBin_aAIC = FCinEachBin_aAIC(sumI_aAIC,3:end-1); % aAIC
    SortedFCinEachBin = vertcat(FCinEachBin_mPFC,FCinEachBin_Proj,FCinEachBin_aAIC);
    imagesc(SortedFCinEachBin);
    title(['mPFC' num2str(size(FCinEachBin_mPFC,1)) 'mPFC-to-aAIC' num2str(size(FCinEachBin_Proj,1)) 'aAIC' num2str(size(FCinEachBin_aAIC,1))]);
elseif exist('conn_chain_Tar')
    imagesc(FCinEachBin_Tar(sumI_Tar,3:end-1));
    title([Target num2str(size(FCinEachBin_mPFC,1))]);
end
cmap = [1 1 1;22/255 73/255 157/255;231/255 31/255 25/255];
colormap(cmap)
ylabel('Pairs ID')
xlabel('Time in delay period (s)')
set(gca,'XTick',0.5:1:4.5,'XTickLabel',{'0','1','2','3','4'}); hold on;
ylim([0.5,max(ylim())]);
set(gcf,'Renderer','Painter'); saveas(gcf,[Neurontype '-' Target '-' Group '-' WhichTrials '_1msbin'],'fig'); close;





