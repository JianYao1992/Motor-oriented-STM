function [UnitsIDwithPutaSelec,BinsIDwithPutaSelec]= GetPutativeSigUnitID(FRinS1,FRinS2,EleBinNum,BaseDura,SampOdorLen,DelayDuration,TestOdorLen,ShownBaseDura,ShownRespWindowDura,IsForHitCR)

p = []; FDR = []; 
UnitsIDwithPutaSelec = [];
BinsIDwithPutaSelec = cell(1,size(FRinS1,2));

for i = 1:size(FRinS1,2) % neuron
    fprintf('Analysing %dth neuron of %d neurons\n',i,size(FRinS1,2))
    for j = 1:floor(size(FRinS1{1,1},2)/EleBinNum)
        p(i,j) =  ranksum(mean(FRinS1{i}(:,2+(j-1)*EleBinNum:1+j*EleBinNum),2),mean(FRinS2{i}(:,2+(j-1)*EleBinNum:1+j*EleBinNum),2));
    end
    FDR(i,:) = p(i,:)*(ShownBaseDura + SampOdorLen + DelayDuration + TestOdorLen + ShownRespWindowDura)*10/EleBinNum;
    tempBinsIDwithPutaSelec = [];
    for j = 1:floor(size(FRinS1{1,1},2)/EleBinNum)
        if FDR(i,j) <= 0.05
            tempBinsIDwithPutaSelec = [tempBinsIDwithPutaSelec j];
        end
    end
    BinsIDwithPutaSelec{i} = tempBinsIDwithPutaSelec;
end

if IsForHitCR == 1
    for i = 1:size(FDR,1) % neuron
        for j = (BaseDura-ShownBaseDura)*10/EleBinNum+1:(BaseDura + SampOdorLen + DelayDuration + TestOdorLen + ShownRespWindowDura)*10/EleBinNum
            if FDR(i,j) <= 0.05
                if isempty(find(UnitsIDwithPutaSelec == i))
                    UnitsIDwithPutaSelec = [UnitsIDwithPutaSelec; i];
                end
            end
        end
    end
else
    for i = 1:size(FDR,1) % neuron
        for j = BaseDura + SampOdorLen + 1:BaseDura + SampOdorLen + DelayDuration
            if FDR(i,j) <= 0.05
                if isempty(find(UnitsIDwithPutaSelec == i))
                    UnitsIDwithPutaSelec = [UnitsIDwithPutaSelec; i];
                end
            end
        end
    end
end
