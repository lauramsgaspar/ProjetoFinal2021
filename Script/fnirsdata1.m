%% Tarefas
%% Load Data
clear all
close all
addpath('utils')
nirs_file =nirs.io.loadDirectory('C:\Users\HP\Desktop\fnirsdata');  
%% a)defenir a frequencia de amostragem
j=nirs.modules.Resample(); 
j.Fs=10;
nirs_file = j.run(nirs_file); 

%% b)Substituir o nome das nossas condições
j=nirs.modules.RenameStims(); 
j.listOfChanges={
      'stim_channel21' 'Baseline';
      'stim_channel22' 'ImaginacaoMotora';
      'stim_channel23' 'ExecucaoMotora'  
      };
nirs_file = j.run(nirs_file); %executar

%% c) mudar as condições para 30 segundos
%nirs_file=nirs.viz.StimUtil(nirs_file);
duracao_bloco=30;
nirs_file=nirs.design.change_stimulus_duration(nirs_file,[],duracao_bloco);
nirs_file=j.run(nirs_file);

%% Excluir canais com base no valor de Variância e racio Sinal-Ruido
CVThreshold=0.15;
SNRTheshold=20;
[idx] = DetetarMausCanaisWL (nirs_file, CVThreshold, SNRTheshold);

CleanedData = nirs_file;
CleanedData.data(:,idx)=[];

%% Agora vamos transformar dados de intensidade de luz para dados de hemoglobina através LBLM
j = nirs.modules.OpticalDensity(  );
j = nirs.modules.BeerLambertLaw( j );
hb=j.run(nirs_file);
Fs=hb.Fs;

%% Dados fNIRS por Condição 
%% Condição Baseline  

%Numero de vezes que a condição aparece
NrOfOnsetB = size(hb.stimulus.values{1}.onset,1);
Baseline_Vector = zeros(duracao_bloco*hb.Fs,NrOfOnsetB);

ChannelOfInterest = 25;

for i=1 : NrOfOnsetB
    Idx_iB= hb.stimulus.values{1}.onset(i);
    Idx_fB= hb.stimulus.values{1}.onset(i) + duracao_bloco*hb.Fs;
    Baseline_Vector(:,i) =hb.data(Idx_iB:((Idx_fB)-1), ChannelOfInterest);
end
%Matlab reshape
Baseline_Data = reshape(Baseline_Vector, [NrOfOnsetB*size(Baseline_Vector,1), 1]);

 %% João P. fez assim?
%Numero de vezes que a condição aparece
NrOfOnsetB = size(hb.stimulus.values{1}.onset,1);
Baseline_Vector = zeros(duracao_bloco*hb.Fs,NrOfOnsetB);

for i=1 : NrOfOnsetB-1
    [value,Idx_iB]= min(abs(hb.time - hb.stimulus.values{1}.onset(i)));
    [value2,Idx_fB]= min(abs(hb.time - hb.stimulus.values{1}.onset(i) + duracao_bloco*hb.Fs));
    Baseline_Vector(:,:,i) =hb.data(Idx_iB:Idx_fB,:);
end
%Matlab reshape
Baseline_Data = reshape(Baseline_Vector, [NrOfOnsetB*size(Baseline_Vector,1), 1]);




% figure;
% plot(Baseline_Vector(:,1))
% title('Os primeiros 30 segundos da condição de Baseline ')
% xlabel('Timepoints');
% 
% figure;
% plot(Baseline_Data)
% title('Condição de Baseline ')
% xlabel('Timepoints');

%% Imaginação Motora 

%Numero de vezes que a condição aparece
NrOfOnsetMI = size(hb.stimulus.values{2}.onset,1);
MotorImagery_Vector = zeros(duracao_bloco*hb.Fs,NrOfOnsetMI);
for i =1 : NrOfOnsetMI
    Idx_iMI=hb.stimulus.values{2}.onset(i);
    Idx_fMI=hb.stimulus.values{2}.onset(i)  + duracao_bloco*hb.Fs;
    MotorImagery_Vector(:,i) =hb.data(Idx_iMI:((Idx_fMI)-1), ChannelOfInterest);
end
%Matlab reshape
MotorImagery_Data = reshape(MotorImagery_Vector, [NrOfOnsetMI*size(MotorImagery_Vector,1), 1]);

%% Ação Motora 

%Numero de vezes que a condição aparece
NrOfOnsetMA = size(hb.stimulus.values{3}.onset,1);
MotorAction_Vector = zeros(duracao_bloco*hb.Fs,NrOfOnsetMA);
for i =1 : NrOfOnsetMA
    Idx_iMA=hb.stimulus.values{3}.onset(i);
    Idx_fMA=hb.stimulus.values{3}.onset(i)  + duracao_bloco*hb.Fs;
    MotorAction_Vector(:,i) =hb.data(Idx_iMA:((Idx_fMA)-1), ChannelOfInterest);
end
%Matlab reshape
MotorAction_Data = reshape(MotorAction_Vector, [NrOfOnsetMA*size(MotorAction_Vector,1), 1]);

figure;
x = [Baseline_Data; MotorImagery_Data; MotorAction_Data];
g = [zeros(length(Baseline_Data), 1); ones(length(MotorImagery_Data), 1); 2*ones(length(MotorAction_Data), 1)];
boxplot(x, g)
title('Boxplot das 3 condicoes')

%% To do:  
% %% Diferenciar dados por condição para canais Hbo e Hbr 
% CanalHbo_Baseline_Data=Baseline_Data(1:2:end,:);
% CanalHbr_Baseline_Data=Baseline_Data(2:2:end,:);
% 
% CanalHbo_MotorImagery_Data=MotorImagery_Data(1:2:end,:);
% CanalHbr_MotorImagery_Data=MotorImagery_Data(2:2:end,:);
% 
% CanalHbo_MotorAction_Data=MotorAction_Data(1:2:end,:);
% CanalHbr_MotorAction_Data=MotorAction_Data(2:2:end,:);
% 
% %% Boxplot dados apos beer-lambert (figura 1 Hbo e figura 2 Hbr)
% xhbo = [CanalHbo_Baseline_Data; CanalHbo_MotorImagery_Data; CanalHbo_MotorAction_Data];
% ghbo = [zeros(length(CanalHbo_Baseline_Data), 1); ones(length(CanalHbo_MotorImagery_Data), 1); 2*ones(length(CanalHbo_MotorAction_Data), 1)];
% boxplot(xhbo, ghbo)
% 
% xhbr = [CanalHbr_Baseline_Data; CanalHbr_MotorImagery_Data; CanalHbr_MotorAction_Data];
% ghbr = [zeros(length(CanalHbr_Baseline_Data), 1); ones(length(CanalHbr_MotorImagery_Data), 1); 2*ones(length(CanalHbr_MotorAction_Data), 1)];
% boxplot(xhbr, ghbr)

%% Matriz de Anticorrelacao Hbo e Hbr - mais especifica para o fnirs

[diagonalCorrelacao, diagonalCorrelacaoLongo] = AutoCorrelationMatrix(hb);


%% Short Separation Regression 

[CleanedData] = ShortSeparationRegression(hb);

%% Verificar autocorrelação apos SSR

[diagonalCorrelacaoSSR , diagonalCorrelacaoLongoSSR] = AutoCorrelationMatrix(CleanedData);

boxplot([diagonalCorrelacaoLongo, diagonalCorrelacaoLongoSSR])

title('Autocorrelação com canais longos antes e após o SSR')
%% Filtro passa baixo - experimentar

% hbFiltered=CleanedData;
% tam=size(CleanedData.data,2);
% for f=1:tam
%     hbFiltered.data(:,f)=lowpass(CleanedData.data(:,f),0.2, 10); %ver para R2016a
% end

%% Filtro Butterworth de 3ªordem

UnfilteredSignal=CleanedData.data; 
sample_rate=CleanedData.Fs;
FilterOrder=3;
%Frequencia de corte
FreqInterval1=0.01;
FreqInterval2=0.09;

[FilteredSignal, Filtered1, BandPassFilteredSignal] = filterButterworth(sample_rate, UnfilteredSignal, FreqInterval1, FreqInterval2,FilterOrder);

%% Spectogram Analysis

figure;
spectrogram(UnfilteredSignal)
title('Sinal Não Filtrado')
figure;
spectrogram(Filtered1)
title('Depois do Filtro Butterworth Passa-Alto com Fc=0.01 Hz ')
figure;
spectrogram(FilteredSignal)
title('Depois do Filtro Butterworth Passa-Baixo com Fc=0.09 Hz ')
figure;
spectrogram(BandPassFilteredSignal)
title('Depois do Filtro Butterworth Passa-Banda com wn=[0.01 - 0.09] Hz ')

%% TDDR - aplicar em dados já filtrados
[signal_corrected] = TDDR(BandPassFilteredSignal , sample_rate );

%usando Toolbox
nirs.modules.TDDR

%% Ver métricas apos SSR : Autocorr e Boxplot dados por condição por sinal --> (Bas1 Bas2 ; Mi1 Mi2; MA1 Ma2);

%Condição Baseline após SSR
NrOfOnsetB2 = size(CleanedData.stimulus.values{1}.onset,1);
Baseline_Vector2 = zeros(duracao_bloco*CleanedData.Fs,NrOfOnsetB2);

for i=1 : NrOfOnsetB
    Idx_iB2=CleanedData.stimulus.values{1}.onset(i);
    Idx_fB2=CleanedData.stimulus.values{1}.onset(i)  + duracao_bloco*CleanedData.Fs;
    Baseline_Vector2(:,i) =CleanedData.data(Idx_iB2:((Idx_fB2)-1), ChannelOfInterest);
end
Baseline_Data2 = reshape(Baseline_Vector2, [NrOfOnsetB2*size(Baseline_Vector2,1), 1]);

%Condição Imaginação Motora após SSR
NrOfOnsetMI2 = size(CleanedData.stimulus.values{2}.onset,1);
MotorImagery_Vector2 = zeros(duracao_bloco*CleanedData.Fs,NrOfOnsetMI2);

for i =1 : NrOfOnsetMI2
    Idx_iMI2=CleanedData.stimulus.values{2}.onset(i);
    Idx_fMI2=CleanedData.stimulus.values{2}.onset(i)  + duracao_bloco*CleanedData.Fs;
    MotorImagery_Vector2(:,i) = CleanedData.data(Idx_iMI2:((Idx_fMI2)-1), ChannelOfInterest);
end
MotorImagery_Data2 = reshape(MotorImagery_Vector2, [NrOfOnsetMI2*size(MotorImagery_Vector2,1), 1]);

%Condição Execução Motora após SSR
NrOfOnsetMA2 = size(CleanedData.stimulus.values{3}.onset,1);
MotorAction_Vector2 = zeros(duracao_bloco*CleanedData.Fs,NrOfOnsetMA2);

for i =1 : NrOfOnsetMA2
    Idx_iMA2=CleanedData.stimulus.values{3}.onset(i);
    Idx_fMA2=CleanedData.stimulus.values{3}.onset(i)  + duracao_bloco*CleanedData.Fs;
    MotorAction_Vector2(:,i) =CleanedData.data(Idx_iMA2:((Idx_fMA2)-1), ChannelOfInterest);
end
MotorAction_Data2 = reshape(MotorAction_Vector2, [NrOfOnsetMA2*size(MotorAction_Vector2,1), 1]);

%Boxplot dados por condiçao com antes e depois do SSR
xantesdepois = [Baseline_Data; Baseline_Data2; MotorImagery_Data; MotorImagery_Data2; MotorAction_Data; MotorAction_Data2];
gantesdepois = [zeros(length(Baseline_Data), 1); ones(length(Baseline_Data2), 1); 2*ones(length(MotorImagery_Data), 1); 3*ones(length(MotorImagery_Data2), 1); 4*ones(length(MotorAction_Data), 1); 5*ones(length(MotorAction_Data2), 1)];
boxplot(xantesdepois, gantesdepois)
title('Boxplot dos dados Hbo por condiçao com o antes e o depois do SSR')


%% d)Correr o GLM e ver o plot dos dados do canal com maior valor de beta para a condição 'Motor Action'.
%Agora sim, vamos correr o GLM
j=nirs.modules.GLM(); 
SubjStats=j.run(hb); 

%Tentativa com erros de variaveis mas correta 

% SubjStats(1).draw('beta');
% %criei uma variável para guardar os valores beta de todas as condiçoes
% beta=SubjStats(1).beta; 
% tabela=SubjStats(1).table;
% canal=hb(1).probe.link;
% maiorbeta=0;
% canalmaiorbeta=0;
% indice=0;
% %através de um ciclo for selecionei as condições correspondentes à
% %condição2 -imaginação motora
% for i=53:104
%     if(maiorbeta<beta(i))
%         maiorbeta=beta(i);
%         indice=i;
%     source=tabela(indice,1);
%     detetor=tabela(indice,2);
%     tipo=tabela(indice,3);
%     end
% end
% numerocanais=size(hb.data);
% for m=1:(numerocanais(:,2))
%     for n=1:3
%         if((table2array(source)==table2array(canal(m,1))) && (table2array(detetor)==(table2array(canal(m,2))))&& (strcmp(char(tipo), char(canal(m,3)))))
%            canalmaiorbeta=m;
%            %hb.draw(m);
%            figure;
%           plot(hb.data(:,canalmaiorbeta));
%         end
%     end
% end

%Resolução do João P.
% [max_beta1 , t_index1]= max(SubjStats.beta(1:2:end)); %apenas consideramos o hbo
% index1 =2*t_index1 -1; %estávamos indexando apenas metade do vetor


%% To do : média, desvio padrão, boxplot das tres condiçoes
% 
%Condição 1 - Baseline
 
% mediaB=mean(Baseline_Data); %média de um vetor
% desvioB=std(Baseline_Data); %Desvio Padrão de um vetor
% BoxPlotB=boxplot(Baseline_Data); %criar boxplot 
% boxplot(Baseline_Vector(:,[1,2,3,4,5,6,7,8,9,10])); %box plot de quando aparece a 1º,2º,...da 1 cond
% boxplot(Baseline_Data);
% 
% %Condição 2 - Motor Imaginary
% mean(MotorImagery_Vector); %média de um vetor
% std(MotorImagery_Vector); %Desvio Padrão de um vetor
% boxplot(MotorImagery_Vector); %criar boxplot 
% boxplot(MotorImagery_Vector(:,1)); %box plot de quando aparece a 1ºa 1 cond
%  
% %Condição 3- Motor Action
% mean(MotorAction_Vector); %média de um vetor
% std(MotorAction_Vector); %Desvio Padrão de um vetor
% boxplot(MotorAction_Vector); %criar boxplot 
% boxplot(MotorAction_Vector(:,1)); %box plot de quando aparece a 1ºa 1 cond

