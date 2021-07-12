function [ CNR_MA_Hbo, CNR_MA_Hbr, CNR_MI_Hbo, CNR_MI_Hbr ] = CNR( data, BasData, MIData, MAData )
%Idx apenas para canais longos 
idx_long = find(data.probe.link.detector < 8);

%Calcular o CNR quer para MA como MI (f�rmula)
t_CNR_MA = (mean(MAData,1)- mean(BasData,1)) ./ (sqrt(std(MAData,1) + std(BasData,1)));
t_CNR_MI = (mean(MIData,1)- mean(BasData,1)) ./ (sqrt(std(MIData,1) + std(BasData,1)));

%Considerar apenas os canais longos
CNR_MA = t_CNR_MA(:,idx_long);
CNR_MI = t_CNR_MI(:, idx_long);

%Dividir em canais Hbo e Hbr
CNR_MA_Hbo = CNR_MA(1:2:end);
CNR_MA_Hbr = CNR_MA(2:2:end);

CNR_MI_Hbo = CNR_MI(1:2:end);
CNR_MI_Hbr = CNR_MI(2:2:end);

%% Plot CNR por canal
%Execu��o Motora
figure;
subplot(1,2,1)
plot(CNR_MA_Hbo,'o');
ylim([-10, size(CNR_MA_Hbo,2)])
title('CNR da condi��o execu��o motora para canais HbO');
ylabel('CNR = (\mu MA - \mu Bas) / sqrt(\sigma MA + \sigma Bas)')
xlabel('Canal')

subplot(1,2,2)
plot(CNR_MA_Hbr, 'o');
ylim([-10, size(CNR_MA_Hbr,2)])
title('CNR da condi��o execu��o motora para canais HbR');
ylabel('CNR = (\mu MA - \mu Bas) / sqrt(\sigma MA + \sigma Bas)')
xlabel('Canal')

%Imagina��o Motora
figure;
subplot(1,2,1)
plot(CNR_MI_Hbo, 'o');
ylim([-10, size(CNR_MI_Hbo,2)])
title('CNR da condicao imagina��o motora para canais HbO');
ylabel('CNR = (\mu MI - \mu Bas) / sqrt(\sigma MI + \sigma Bas)')
xlabel('Canal')

subplot(1,2,2)
plot(CNR_MI_Hbr, 'o');
ylim([-10, size(CNR_MI_Hbr,2)])
title('CNR da condicao imagina��o motora para canais HbR');
ylabel('CNR = (\mu MI - \mu Bas) / sqrt(\sigma MI + \sigma Bas)')
xlabel('Canal')

end
