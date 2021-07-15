%% Baseado no artigo: Gagnon, L., Cooper, R. J., Yücel, M. A., Perdue, K. L., Greve, D. N., & Boas, D. A. (2012). 
%Short separation channel location impacts the performance of short channel regression in NIRS. NeuroImage, 59(3), 2518–2528.
function [RegressedData]=ShortSeparationRegression(hb)

UnregressedChannelPlot=hb.data(:,1);
CanalCurto_indice=find(hb.probe.link.detector >= 8);
CanalLongo_indice=find(hb.probe.link.detector < 8);
RegressedData=hb;
for i=1:size(CanalLongo_indice,1)
     [value, IdxMaisProximo] = min(abs(CanalLongo_indice(i)-CanalCurto_indice));
     
    %Calcular alpha por canal longo
    %vetor coluna - não a distâncias, mas sim os dados
    %alpha=S'.L % (S'.S)
     
    t_calc1 = hb.data(:,CanalCurto_indice(IdxMaisProximo))' * hb.data(:, CanalLongo_indice(i));
    t_calc2 = hb.data(:,CanalCurto_indice(IdxMaisProximo))' *  hb.data(:,CanalCurto_indice(IdxMaisProximo));
    
    alpha(i) = t_calc1/ t_calc2;
    
    %L'=L-alpha.L
    
    RegressedData.data(:,CanalLongo_indice(i)) = hb.data(:,CanalLongo_indice(i)) - alpha(i)* hb.data(:,CanalCurto_indice(IdxMaisProximo));
    
end
RegressedChannelPlot=RegressedData.data(:,1);

%% Representação gráfica
% Referente ao par fonte-detetor S1-D1, sinal HbO)
figure;
subplot(2,1,1);
plot(UnregressedChannelPlot); 
ylabel('Amplitude');
xlabel('Timepoints');
title('Sinal de um canal sem remoção do ruido fisiológico') 
hold on;

subplot(2,1,2);
plot(RegressedChannelPlot);
ylabel('Amplitude');
xlabel('Timepoints');
title('Sinal do mesmo canal com remoçao do ruido fisiológico') 
hold on;
