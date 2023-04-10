%% Calculate the decision variable(DV) of each trial
% DV_X(i) = X(i)*(mean(X)(j~=i)-mean(Y));
% DV_Y(i) = Y(i)*(mean(X)-mean(Y)(j~=i));

function [DVX,DVY] = DecisionVariableCalculation(XAveragedFR,YAveragedFR)

%% DV of X trials
if size(XAveragedFR,1) > 1
   XAveragedFR = XAveragedFR'; 
end
DVX = zeros(1,length(XAveragedFR));
for iTrial = 1:length(XAveragedFR)
    % FR of individual trial
    Xi = XAveragedFR(1,iTrial);
    % averaged FR of remaining trials
    NewXAveragedFR = XAveragedFR;
    NewXAveragedFR(iTrial) = [];
    MeanX = mean(NewXAveragedFR);
    MeanY = mean(YAveragedFR);
    % decision variable
    DVX(1,iTrial) = Xi*(MeanX-MeanY);
end

%% DV of Y trials
if size(YAveragedFR,1)>1
    YAveragedFR = YAveragedFR';
end
DVY = zeros(1,length(YAveragedFR));
for iTrial=1:length(YAveragedFR)
    % FR of individual trial
    Yi = YAveragedFR(1,iTrial);
    % averaged FR of remaining trials
    NewYAveragedFR = YAveragedFR;
    NewYAveragedFR(iTrial) = [];
    MeanX = mean(XAveragedFR);
    MeanY = mean(NewYAveragedFR);
    % decision variable
    DVY(1,iTrial) = Yi*(MeanX-MeanY);
end