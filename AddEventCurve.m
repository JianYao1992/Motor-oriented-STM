function AddEventCurve(OdorMN,DelayMN,ResponseMN,WaterMN,YMax,YMin,LineNum)

if nargin < 7
    LineNum = 6;
end
plot([0 0],[YMin YMax],'k--','linewidth',1)
hold on
plot([OdorMN OdorMN],[YMin YMax],'k--','linewidth',1)
if LineNum > 2
    plot([OdorMN+DelayMN OdorMN+DelayMN],[YMin YMax],'k--','linewidth',1)
    plot([OdorMN+DelayMN+ResponseMN OdorMN+DelayMN+ResponseMN],[YMin YMax],'k--','linewidth',1)
end
if LineNum > 4
    plot([OdorMN+DelayMN+ResponseMN+WaterMN OdorMN+DelayMN+ResponseMN+WaterMN],[YMin YMax],'k--','linewidth',1)
end
