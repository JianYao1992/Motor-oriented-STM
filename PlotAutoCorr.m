%% Spike autocorrelogram

function y = PlotAutoCorr(NewData)
y = [];
Units = NewData(:,1:3);
Boundary = 0.0195;
NeuronID = unique(NewData(:,2));
NeuronID(NeuronID==0) = [];
for iNeuron = 1:length(NeuronID)
    tempSpikeTime = Units(Units(:,2)==NeuronID(iNeuron),3); % targeted neuron
    IntervalNum = length(-Boundary:0.001:Boundary)+1;
    tempAutoCorrlogram = zeros(1,IntervalNum);
    for iSpike=1:length(tempSpikeTime) % each spike
        if length(tempSpikeTime) <= 200
            temp3 = tempSpikeTime;
        else
            if iSpike <= 100
                temp3 = tempSpikeTime(1:iSpike+100);
            elseif iSpike >= length(tempSpikeTime)-100
                temp3 = tempSpikeTime(iSpike-100:length(tempSpikeTime));
            else
                temp3 = tempSpikeTime(iSpike-100:iSpike+100);
            end
        end
        temp3 = setdiff(temp3,tempSpikeTime(iSpike));
        for itr4 = -Boundary:0.001:Boundary
            if isempty(find(temp3+itr4-0.0005<tempSpikeTime(iSpike) & temp3+itr4+0.0005>tempSpikeTime(iSpike),1))
            else
                tempAutoCorrlogram(floor(itr4*1000+floor(IntervalNum/2)+1))=...
                    tempAutoCorrlogram(floor(itr4*1000+floor(IntervalNum/2)+1))+length(find(temp3+itr4-0.0005<tempSpikeTime(iSpike) & temp3+itr4+0.0005>tempSpikeTime(iSpike), 1)); 
            end
        end
    end
    y = [y; tempAutoCorrlogram];
end