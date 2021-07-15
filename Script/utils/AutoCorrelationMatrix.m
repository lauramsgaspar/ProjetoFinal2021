function [ diagonalCorrelacao, diagonalCorrelacaoLongo, mediaAutoCorrLongo, stdAutoCorrLongo] = AutoCorrelationMatrix(fnrisdata)
%Para cada canal dividir em Hbo e Hbr 
CanalHbo = fnrisdata.data(:, 1:2:end);
CanalHbr = fnrisdata.data (:, 2:2:end);

figure;
MatrizCorrelacao=corr(CanalHbo, CanalHbr);
%imagesc-Display image with scaled colors
imagesc(MatrizCorrelacao);
title('Matriz de Autocorrelacao HbO e HbR')
xlabel('Dados dos Canais HbO')
ylabel('Dados dos Canais HbR')

diagonalCorrelacao = diag(MatrizCorrelacao);
%Calcular média e std da diagonal
mediaAutoCorr = mean(diagonalCorrelacao);
stdAutoCorr = std(diagonalCorrelacao);

%Apenas com os canais longos
CanalLongo_indice=find(fnrisdata.probe.link.detector < 8);
DadosCanalLongo=fnrisdata.data(:, CanalLongo_indice);
CanalLongoHbo = DadosCanalLongo (:, 1:2:end);
CanalLongoHbr = DadosCanalLongo (:, 2:2:end);

MatrizCorrelacaoCanalLongo = corr(CanalLongoHbo,CanalLongoHbr);
diagonalCorrelacaoLongo = diag(MatrizCorrelacaoCanalLongo);
mediaAutoCorrLongo = mean(diagonalCorrelacaoLongo);
stdAutoCorrLongo = std(diagonalCorrelacaoLongo);

figure;
%imagesc-Display image with scaled colors
imagesc(MatrizCorrelacaoCanalLongo);
title('Matriz de Autocorrelacao HbO e HbR')
xlabel('Dados dos Canais Longo HbO','FontSize',12)
ylabel('Dados dos Canais Longo HbR','FontSize',12)
end