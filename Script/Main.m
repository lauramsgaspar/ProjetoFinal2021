%% Est�gio Curricular 
%% Laura Gaspar - a21280722@isec.pt
%A pipeline de an�lise de dados, pr�-processamento de dados e an�lise estat�stica do impacto do pr�-processamento
%de dados fNIRS em tarefas motoras (BAS, MI, MA), desenvolvido no �mbito deste est�gio, � apresentada seguidamente:

%% Load Data
clear all
close all
addpath('utils')
%nirs_file =nirs.io.loadDirectory('C:\Users\HP\Documents\GitHub\fNIRS_ISEC2021'); 
nirs_file =nirs.io.loadDirectory('C:\Users\User\Documents\GitHub\fNIRS_ISEC2021\Nirs data');  

%% Tratamento de dados
j=nirs.modules.Resample(); 
j.Fs=10;
nirs_file = j.run(nirs_file); 

%Renomear as condi��es
j=nirs.modules.RenameStims(); 
j.listOfChanges={
      'stim_channel21' 'Baseline';
      'stim_channel22' 'ImaginacaoMotora';
      'stim_channel23' 'ExecucaoMotora'  
      };
nirs_file = j.run(nirs_file); %executar

%Alterar a dura��o dos blocos de condi��o para 30 segundos
%nirs_file=nirs.viz.StimUtil(nirs_file);
duracao_bloco=30;
nirs_file=nirs.design.change_stimulus_duration(nirs_file,[],duracao_bloco);
nirs_file=j.run(nirs_file);

%% Avalia��o da qualidade do sinal fNIRS (CV, SNR, Matriz Autocorrel��o entre HbR/HbO, CNR)
%% Excluir canais com base no valor de Vari�ncia (CV(%)) e r�cio Sinal-Ruido (SNR)
CVThreshold=0.15;
SNRThreshold=20;
[idx] = DetetarMausCanaisWL (nirs_file, CVThreshold, SNRThreshold); 

NirsData_CR = nirs_file;  %NirsData after Channel Rejection
NirsData_CR.data(:,idx)=[];

%% Transformar dados de intensidade de luz para dados de hemoglobina aplicando a densidade �tica e a Lei Modificada de Beer Lambert
j = nirs.modules.OpticalDensity(  );
j = nirs.modules.BeerLambertLaw( j );
hb=j.run(NirsData_CR);
Fs=hb.Fs;

%% Representa��o gr�fica do antes e depois da Lei Modificada de Beer Lambert
% referente ao par fonte-detetor: S1-D1, sinal HbO
figure;
subplot(2,1,1)
plot(nirs_file.data(:,1));
title('Dados em formato de comprimento de onda');
ylabel('Amplitude');
xlabel('Timepoints');
subplot(2,1,2)
plot(hb.data(:,1));
title('Dados de hemoglobina');
ylabel('Amplitude'); %concentra��es de hemoglobina
xlabel('Timepoints');

%% GLM - Dados iniciais
j=nirs.modules.GLM(); 
SubjStats=j.run(hb); 
SubjStats.draw('beta', [-2 2], 'p < 0.05')
ValoresBetaInicial = SubjStats.beta; % Valores de beta guardados numa coluna [156x1]

%Nota: esta coluna contem hbo, hbr, canais longos, canais curtos e as 3 condi��es

%% Agrupar dados fNIRS por condi��o 

% BAS - antes do processamento
NRofChannels = size(hb.probe.link, 1);
NrOfOnsetB = size(hb.stimulus.values{1}.onset,1);
Baseline_Vector = zeros(duracao_bloco*hb.Fs,NRofChannels,NrOfOnsetB-1);

for i=1 : NrOfOnsetB-1
    Idx_iB= hb.stimulus.values{1}.onset(i);
    Idx_fB= hb.stimulus.values{1}.onset(i) + duracao_bloco*hb.Fs;
    Baseline_Vector(:,:,i) =hb.data(Idx_iB:((Idx_fB)-1),:);
end
%Matlab reshape
Baseline_Data = reshape(Baseline_Vector, [(NrOfOnsetB-1)*size(Baseline_Vector,1), NRofChannels]);


% MI - antes do processamento
NrOfOnsetMI = size(hb.stimulus.values{2}.onset,1);
MotorImagery_Vector = zeros(duracao_bloco*hb.Fs,NRofChannels, NrOfOnsetMI);
for i =1 : NrOfOnsetMI
    Idx_iMI=hb.stimulus.values{2}.onset(i);
    Idx_fMI=hb.stimulus.values{2}.onset(i)  + duracao_bloco*hb.Fs;
    MotorImagery_Vector(:,:,i) =hb.data(Idx_iMI:((Idx_fMI)-1),:);
end
%Matlab reshape
MotorImagery_Data = reshape(MotorImagery_Vector, [NrOfOnsetMI*size(MotorImagery_Vector,1), NRofChannels]);

% MA - antes do processamento
NrOfOnsetMA = size(hb.stimulus.values{3}.onset,1);
MotorAction_Vector = zeros(duracao_bloco*hb.Fs,NRofChannels,NrOfOnsetMA);
for i =1 : NrOfOnsetMA
    Idx_iMA=hb.stimulus.values{3}.onset(i);
    Idx_fMA=hb.stimulus.values{3}.onset(i)  + duracao_bloco*hb.Fs;
    MotorAction_Vector(:,:,i) =hb.data(Idx_iMA:((Idx_fMA)-1),:);
end
%Matlab reshape
MotorAction_Data = reshape(MotorAction_Vector, [NrOfOnsetMA*size(MotorAction_Vector,1), NRofChannels]);

%% Diferenciar dados por condi��o para canais HbO e HbR - Antes do processamento
%BAS
BaselineData_HBO = Baseline_Data(:,1:2:end);
BaselineData_HBR = Baseline_Data(:,2:2:end);
%MI
MIData_HBO = MotorImagery_Data(:,1:2:end);
MIData_HBR = MotorImagery_Data(:,2:2:end);
%MA
MAData_HBO=MotorAction_Data(:,1:2:end);
MAData_HBR=MotorAction_Data(:,2:2:end);

%% Boxplot dados ap�s aplica��o da Lei Modificada de Beer-Lambert (figura 1 HbO e figura 2 HbR)

%Figura 1
CombinedBAS_HBO = reshape(BaselineData_HBO, [(size(BaselineData_HBO,1)* size(BaselineData_HBO,2)), 1]);
CombinedMI_HBO = reshape(MIData_HBO, [(size(MIData_HBO,1)* size(MIData_HBO,2)), 1]);
CombinedMA_HBO = reshape(MAData_HBO, [(size(MAData_HBO,1)* size(MAData_HBO,2)), 1]);

% xhbo = [CombinedBAS_HBO; CombinedMI_HBO; CombinedMA_HBO];
% ghbo = [zeros(length(CombinedBAS_HBO), 1); ones(length(CombinedMI_HBO), 1); 2*ones(length(CombinedMA_HBO), 1)];
% boxplot(xhbo, ghbo)
% title('Combina��o das condi��es para canais Hbo');
% xlabel('Condi��es');
% 
%Figura 2
CombinedBAS_HBR = reshape(BaselineData_HBR, [(size(BaselineData_HBR,1)* size(BaselineData_HBR,2)), 1]);
CombinedMI_HBR = reshape(MIData_HBR, [(size(MIData_HBR,1)* size(MIData_HBR,2)), 1]);
CombinedMA_HBR = reshape(MAData_HBR, [(size(MAData_HBR,1)* size(MAData_HBR,2)), 1]);

% xhbr = [CombinedBAS_HBR; CombinedMI_HBR; CombinedMA_HBR];
% ghbr = [zeros(length(CombinedBAS_HBR), 1); ones(length(CombinedMI_HBR), 1); 2*ones(length(CombinedMA_HBR), 1)];
% boxplot(xhbr, ghbr)
% title('Combina��o das condi��es para canais Hbo');
% xlabel('Condi��es');

%% Contrast-to-noise ration (CNR) - Antes do processamento
[CNR_MA_Hbo, CNR_MA_Hbr, CNR_MI_hbo, CNR_MI_Hbr ] = CNR( hb, Baseline_Data, MotorImagery_Data, MotorAction_Data);

%% Matriz de Autocorrela��o HbO e HbR - Antes do processamento (m�trica mais especifica para o fNIRS)
[diagonalCorrelacaoAntes, diagonalCorrelacaoLongoAntes, mediaAutoCorrLongoAntes, stdAutoCorrLongoAntes] = AutoCorrelationMatrix(hb);

%% T�nicas de Pre-processamento (Filtro passa-banda butterworth de 3�ordem, TDDR, SSR)
%% Filtragem de dados - Filtro passa-banda Butterworth de 3�ordem
UnfilteredSignal=hb.data; 
sample_rate=hb.Fs;
%ordem
FilterOrder=3;
%Frequencia de corte
FreqInterval1=0.01;
FreqInterval2=0.09;

[FilteredSignal, Filtered1, BandPassFilteredSignal] = filterButterworth(sample_rate, UnfilteredSignal, FreqInterval1, FreqInterval2,FilterOrder);

%% An�lise espectograma
% referente ao par fonte-detetor S2-D3, sinal HbO
figure;
subplot(2,1,1)
spectrogram(UnfilteredSignal(:,11));
title('Espectograma do canal 11 N�o Filtrado')
% spectrogram(Filtered1(:,11))
% title('Espectograma do canal 11 depois da aplica��o do Filtro Butterworth Passa-Alto de 3�ordem com Fc=0.01 Hz ')
% figure;
% spectrogram(FilteredSignal(:,11))
% title('Espectograma do canal 11 depois da aplica��o do Filtro Butterworth Passa-Baixo de 3� ordem com Fc=0.09 Hz ')
% figure;
subplot(2,1,2)
spectrogram(BandPassFilteredSignal(:,11))
title('Espectograma do canal 11 aplicando o Filtro Butterworth Passa-Banda de 3� ordem com wn=[0.01, 0.09] Hz ')

%% Corre��o do ru�do de movimento 
% Foi implementado o algoritmo de corre��o Temporal Derivative Distribution
% Repair (TDDR)
MotionCorrected = TDDR(BandPassFilteredSignal, hb.Fs);

%Representa��o gr�fica
%referente ao par fonte-detetor S1-D1, sinal HbO
figure;
subplot(3,1,1)
plot(hb.data(:,1));
title('Sinal sem qualquer tipo de processamento(canal 1)')
xlabel('Timepoints');
ylabel('Amplitude');
subplot(3,1,2)
plot(BandPassFilteredSignal(:,1));
title('Sinal antes da aplica��o do TDDR(canal 1)');
xlabel('Timepoints');
ylabel('Amplitude');
subplot(3,1,3)
plot(MotionCorrected(:,1));
title('Sinal ap�s aplica��o do TDDR(canal 1)');
xlabel('Timepoints');
ylabel('Amplitude');
ylim([-100 100])

%% Remo��o de ru�do extracerebral 
% Recorreu-se ao m�todo Short Separation Regression - SSR
hb.data = MotionCorrected;
[RegressedData] = ShortSeparationRegression(hb);

%% Autocorrela��o entre HbR e HbO - Ap�s processamento
% NOTA: O DATASET DE INTERESSE AGORA � (RegressedData)

FinalHb = RegressedData;
[diagonalCorrelacaoDepois , diagonalCorrelacaoLongoDepois, mediaAutoCorrLongoDepois, stdAutoCorrLongoDepois] = AutoCorrelationMatrix(FinalHb);

boxplot([diagonalCorrelacaoLongoAntes, diagonalCorrelacaoLongoDepois])
title('Autocorrela��o com canais longos antes e ap�s o Processamento')

%% Agrupamento de dados - Ap�s processamento

%BAS - Ap�s processamento
Baseline_Vector2 = zeros(duracao_bloco*FinalHb.Fs,NRofChannels,NrOfOnsetB-1);

for i=1 : NrOfOnsetB-1
    Idx_iB= FinalHb.stimulus.values{1}.onset(i);
    Idx_fB= FinalHb.stimulus.values{1}.onset(i) + duracao_bloco*FinalHb.Fs;
    Baseline_Vector2(:,:,i) =FinalHb.data(Idx_iB:((Idx_fB)-1),:);
end

%Matlab reshape
Baseline_DataFinal = reshape(Baseline_Vector2, [(NrOfOnsetB-1)*size(Baseline_Vector2,1), NRofChannels]);


%MI - Ap�s processamento
MotorImagery_Vector2 = zeros(duracao_bloco*FinalHb.Fs,NRofChannels, NrOfOnsetMI);
for i =1 : NrOfOnsetMI
    Idx_iMI=FinalHb.stimulus.values{2}.onset(i);
    Idx_fMI=FinalHb.stimulus.values{2}.onset(i)  + duracao_bloco*FinalHb.Fs;
    MotorImagery_Vector2(:,:,i) =FinalHb.data(Idx_iMI:((Idx_fMI)-1),:);
end
%Matlab reshape
MotorImagery_DataFinal = reshape(MotorImagery_Vector2, [NrOfOnsetMI*size(MotorImagery_Vector2,1), NRofChannels]);


%MA - Ap�s processamento
MotorAction_Vector2 = zeros(duracao_bloco*hb.Fs,NRofChannels,NrOfOnsetMA);
for i =1 : NrOfOnsetMA
    Idx_iMA=FinalHb.stimulus.values{3}.onset(i);
    Idx_fMA=FinalHb.stimulus.values{3}.onset(i)  + duracao_bloco*FinalHb.Fs;
    MotorAction_Vector2(:,:,i) =FinalHb.data(Idx_iMA:((Idx_fMA)-1),:);
end
%Matlab reshape
MotorAction_DataFinal = reshape(MotorAction_Vector2, [NrOfOnsetMA*size(MotorAction_Vector2,1), NRofChannels]);

%% Diferenciar dados por condi��o para canais HbO e HbR - Ap�s processamento

Baseline_Final_HBO = Baseline_DataFinal(:,1:2:end);
Baseline_Final_HBR = Baseline_DataFinal(:,2:2:end);

MIData_Final_HBO = MotorImagery_DataFinal(:,1:2:end);
MIData_Final_HBR = MotorImagery_DataFinal(:,2:2:end);

MAData_Final_HBO=MotorAction_DataFinal(:,1:2:end);
MAData_Final_HBR=MotorAction_DataFinal(:,2:2:end);

%% Verifica��o do impacto do processamento atrav�s do boxplot

% Antes do processamento
CombinedBAS_HBO = reshape(BaselineData_HBO, [(size(BaselineData_HBO,1)* size(BaselineData_HBO,2)), 1]);
CombinedMI_HBO = reshape(MIData_HBO, [(size(MIData_HBO,1)* size(MIData_HBO,2)), 1]);
CombinedMA_HBO = reshape(MAData_HBO, [(size(MAData_HBO,1)* size(MAData_HBO,2)), 1]);
% Depois do processamento
CombinedBAS_HBO_Pos = reshape(Baseline_Final_HBO, [(size(Baseline_Final_HBO,1)* size(Baseline_Final_HBO,2)), 1]);
CombinedMI_HBO_Pos = reshape(MIData_Final_HBO, [(size(MIData_Final_HBO,1)* size(MIData_Final_HBO,2)), 1]);
CombinedMA_HBO_Pos = reshape(MAData_Final_HBO, [(size(MAData_Final_HBO,1)* size(MAData_Final_HBO,2)), 1]);

xhbo = [CombinedBAS_HBO; CombinedMI_HBO; CombinedMA_HBO; CombinedBAS_HBO_Pos; CombinedMI_HBO_Pos; CombinedMA_HBO_Pos];
ghbo = [zeros(length(CombinedBAS_HBO), 1); ones(length(CombinedMI_HBO), 1); 2*ones(length(CombinedMA_HBO), 1);3*ones(length(CombinedBAS_HBO_Pos), 1);4*ones(length(CombinedMI_HBO_Pos), 1);5*ones(length(CombinedMA_HBO_Pos), 1)];
boxplot(xhbo, ghbo)
title('Impacto do processamento para dados de HbO');
xlabel('Condi��es');

%% Contrast-to-noise ratio (CNR) - Ap�s processamento
[CNR_MA_Hbo, CNR_MA_Hbr, CNR_MI_hbo, CNR_MI_Hbr ] = CNR( hb, Baseline_DataFinal, MotorImagery_DataFinal , MotorAction_DataFinal );

%% M�dia de bloco - Block Average
% Canal de interesse possivel(S2D3) = 6;
% MI
%Inicial
Channel11MI_inicial = MIData_HBO (:,6);
NrRepetitions = 4;
t_BlockAverage11MI_inicial = reshape(Channel11MI_inicial, [hb.Fs*duracao_bloco,NrRepetitions]);
BlockAverage11MI_inicial = mean(t_BlockAverage11MI_inicial,2);

%Final
Channel11MI_final = MIData_Final_HBO (:,6);
t_BlockAverage11MI_final = reshape(Channel11MI_final, [hb.Fs*duracao_bloco,NrRepetitions]);
BlockAverage11MI_final = mean(t_BlockAverage11MI_final,2);

% Representa��o gr�fica
figure;
subplot(1,2,1)
plot(BlockAverage11MI_inicial)
xlabel('Timepoints');
ylabel('Amplitude');
title('M�dia de bloco da condi��o MI antes do processamento (canal 11)'); 
subplot(1,2,2)
plot(BlockAverage11MI_final)
xlabel('Timepoints');
ylabel('Amplitude');
title('M�dia de bloco da condi��o MI depois do processamento (canal 11)'); 

%MA
%Inicial
Channel11MA_inicial = MAData_HBO (:,6);
NrRepetitions = 4;
t_BlockAverage11MA_inicial = reshape(Channel11MA_inicial, [hb.Fs*duracao_bloco,NrRepetitions]);
BlockAverage11MA_inicial = mean(t_BlockAverage11MA_inicial,2);

%Final
Channel11MA_final = MAData_Final_HBO (:,6);
t_BlockAverage11MA_final = reshape(Channel11MA_final, [hb.Fs*duracao_bloco,NrRepetitions]);
BlockAverage11MA_final = mean(t_BlockAverage11MA_final,2);

% Representa��o gr�fica
figure;
subplot(1,2,1)
plot(BlockAverage11MA_inicial)
xlabel('Timepoints');
ylabel('Amplitude');
title('M�dia de bloco da condi��o MA antes do processamento (canal 11)'); 
subplot(1,2,2)
plot(BlockAverage11MA_final)
xlabel('Timepoints');
ylabel('Amplitude');
title('M�dia de bloco da condi��o MA depois do processamento (canal 11)'); 

%% GLM - Dados finais
%Agora sim, vamos correr o GLM
j=nirs.modules.GLM(); 
ActivacaoFinal=j.run(FinalHb); 
ActivacaoFinal.draw('beta', [-2 2], 'p < 0.05') %mapas de ativa��o

%% An�lise de betas
%Remover os canais curtos, considerar 52 links por condi��o
%A variavel 'SubjectStats.variables'est� organizada por condi��o 52 links
%cada condi��o pela seguinte ordem 'Baseline', 'MI', 'MA'

IndiceBetaBaseline = find(SubjStats.variables.detector(1:size(hb.probe.link,1))<8);
%Mesma estrutura, portanto basta adicionar o numero de canais a cada indice
IndiceBetaMI = size(hb.probe.link,1) + IndiceBetaBaseline;
IndiceBetaMA = size(hb.probe.link,1) + IndiceBetaMI;

BetasBaselineInicio = SubjStats.beta(IndiceBetaBaseline);
BetasMIInicio = SubjStats.beta(IndiceBetaMI);
BetasMAInicio = SubjStats.beta(IndiceBetaMA);

BetasBaselineFinal = ActivacaoFinal.beta(IndiceBetaBaseline);
BetasMIFinal = ActivacaoFinal.beta(IndiceBetaMI);
BetasMAFinal = ActivacaoFinal.beta(IndiceBetaMA);


% Teste de normalidade das amostras dos Beta por condi�ao (Antes e Depois) 
[h1k, p1k]= lillietest(BetasBaselineInicio)
[h2k, p2k]= lillietest(BetasMIInicio)
[h3k, p3k]= lillietest(BetasMAInicio)

[h4k, p4k]= lillietest(BetasBaselineFinal)
[h5k, p5k]= lillietest(BetasMIFinal)
[h6k, p6k]= lillietest(BetasMAFinal)

% Testes Estatisticos sobre a Distribui��o dos Betas
[p1, h1] = ranksum(BetasBaselineFinal,BetasBaselineInicio);
[p2, h2] = ranksum(BetasMIFinal,BetasMIInicio);
[p3, h3] = ranksum(BetasMAFinal,BetasMAInicio);

%Se quisermos ver por tipo de sinal
%HbO
[p4, h4] = ranksum(BetasBaselineFinal(1:2:end),BetasBaselineInicio(1:2:end))
[p5, h5] = ranksum(BetasMIFinal(1:2:end),BetasMIInicio(1:2:end))
[p6, h6] = ranksum(BetasMAFinal(1:2:end),BetasMAInicio(1:2:end))

%HbR
[p7, p7] = ranksum(BetasBaselineFinal(2:2:end),BetasBaselineInicio(2:2:end))
[p8, p8] = ranksum(BetasMIFinal(2:2:end),BetasMIInicio(2:2:end))
[p9, p9] = ranksum(BetasMAFinal(2:2:end),BetasMAInicio(2:2:end))

[p1, h1, stats1] = ranksum(BetasBaselineFinal,BetasBaselineInicio);
[p2, h2, stats2] = ranksum(BetasMIFinal,BetasMIInicio);
[p3, h3, stats3] = ranksum(BetasMAFinal,BetasMAInicio);

%Se quisermos ver por tipo de sinal
%Hbo
[p4, h4, stats4] = ranksum(BetasBaselineFinal(1:2:end),BetasBaselineInicio(1:2:end))
[p5, h5, stats5] = ranksum(BetasMIFinal(1:2:end),BetasMIInicio(1:2:end))
[p6, h6, stats6] = ranksum(BetasMAFinal(1:2:end),BetasMAInicio(1:2:end))

%Hbr
[p7, p7, stats7] = ranksum(BetasBaselineFinal(2:2:end),BetasBaselineInicio(2:2:end))
[p8, p8, stats8] = ranksum(BetasMIFinal(2:2:end),BetasMIInicio(2:2:end))
[p9, p9, stats9] = ranksum(BetasMAFinal(2:2:end),BetasMAInicio(2:2:end))

%% An�lise de Contraste (MI>MA)
SubjStats.conditions

cImaginacaoSuperiorExecucao = [0 -1 +1]; % Motor Imagery > Motor Action

MISuperiorMA_inicio = SubjStats.ttest(cImaginacaoSuperiorExecucao);
MISuperiorMA_final = ActivacaoFinal.ttest(cImaginacaoSuperiorExecucao);

close all;  

MISuperiorMA_inicio.draw('beta', [-2 2], 'p < 0.05')
MISuperiorMA_final.draw('beta', [-2 2], 'p < 0.05')

%Analise de betas : contraste MI>MA

IndiceContraste = find(MISuperiorMA_inicio.variables.detector<8);

BetaContrasteInicio = MISuperiorMA_inicio.beta(IndiceContraste);
BetaContrasteFinal = MISuperiorMA_final.beta(IndiceContraste);

%Testar a normalidade dos betas de contraste
[h1b, p1b]= lillietest(BetaContrasteInicio);
[h2b, p2b]= lillietest(BetaContrasteFinal);

%Se quisermos ver por tipo de sinal
%Hbo
[p10, h10, stats10] = ranksum(BetaContrasteInicio(1:2:end),BetaContrasteFinal(1:2:end));
%Hbr
[p11, h11, stats11] = ranksum(BetaContrasteInicio(2:2:end),BetaContrasteFinal(2:2:end));

%% Sorting betas
% Top 3 canais mais ativos/valores de beta maiores Por condi��o
%% Sorting betas - Dados n�o processados
[BetaBasAntes , t_indexBasAntes]= sort(SubjStats.beta(1:2:size(hb.probe.link,1)),'descend');
indexBasAntes=2*t_indexBasAntes -1;
top3_indexBasAntes=(indexBasAntes(1:10));

[BetaMIAntes, t_indexMIPos]= sort(SubjStats.beta(size(hb.probe.link,1)+1:2:2*size(hb.probe.link,1)),'descend');
indexMIAntes= size(hb.probe.link,1) + 2*t_indexMIPos -1;
top3_indexMIAntes=(indexMIAntes(1:10));

[BetaMAAntes , t_indexMAAntes]= sort(SubjStats.beta(2*size(hb.probe.link,1)+1:2:3*size(hb.probe.link,1)),'descend');
indexMAAntes= 2*size(hb.probe.link,1) + 2*t_indexMAAntes -1;
top3_indexMAAntes=(indexMAAntes(1:10));

%% Sorting betas - Dados processados
%Top 3 canais mais ativos/valores de beta maiores Por condi��o
%Pos Processamento
[BetaBasPos , t_indexBasPos]= sort(ActivacaoFinal.beta(1:2:size(hb.probe.link,1)),'descend');
indexBasPos=2*t_indexBasPos -1;
top3_indexBasPos=(indexBasPos(1:10));

[BetaMIPos , t_indexMIPos]= sort(ActivacaoFinal.beta(size(hb.probe.link,1)+1:2:2*size(hb.probe.link,1)),'descend');
indexMIPos= size(hb.probe.link,1) + 2*t_indexMIPos -1;
top3_indexMIPos=(indexMIPos(1:10));

[BetaMAPos , t_indexMAPos]= sort(ActivacaoFinal.beta(2*size(hb.probe.link,1)+1:2:3*size(hb.probe.link,1)),'descend');
indexMAPos= 2*size(hb.probe.link,1) + 2*t_indexMAPos -1;
top3_indexMAPos=(indexMAPos(1:10));

%% Sorting betas - Contraste MI>MA
% Dados N�o Processados
[BetaContrasteAntes , t_indexContrasteAntes]= sort(MISuperiorMA_inicio.beta(1:2:end),'descend');
indexContrasteAntes=2*t_indexContrasteAntes -1;
top3_indexContrasteAntes=(indexContrasteAntes(1:10));

%Dados Processados
[BetaContrastePos , t_indexContrastePos]= sort(MISuperiorMA_final.beta(1:2:end),'descend');
indexContrastePos=2*t_indexContrastePos -1;
top3_indexContrastePos=(indexContrastePos(1:10));
