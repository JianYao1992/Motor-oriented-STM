function Ps = GetDataMatrixForMixedRepeatedAnova(Data1, Data2, name)

between_group = [];
Mice_ID = [];
within_group = [];
performance = [];
DayNum = length(Data1);
for i = 1:length(Data1{1}) % control group
    between_group = [between_group; 1*ones(DayNum,1)];
    Mice_ID = [Mice_ID; i*ones(DayNum,1)];
    for j = 1:DayNum % day
        within_group = [within_group; j];
        performance = [performance; Data1{j}(i)];
    end
end

for i = 1:length(Data2{1}) % experimental group
    between_group = [between_group; 2*ones(DayNum,1)];
    Mice_ID = [Mice_ID; (i+length(Data1{1}))*ones(DayNum,1)];
    for j=1:DayNum % day
        within_group = [within_group; j];
        performance = [performance; Data2{j}(i)];
    end
end
Frame = [performance between_group within_group Mice_ID];
[SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(Frame,0,[name '_stat']);