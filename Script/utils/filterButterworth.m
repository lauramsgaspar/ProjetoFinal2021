function [FilteredSignal, Filtered1 ,BandPassFilteredSignal] = filterButterworth(sample_rate, UnfilteredSignal, FreqInterval1, FreqInterval2,FilterOrder)
%[b,a]=butter(n,Wn,ftype)

%Normalizar as frequências
wn1=FreqInterval1/(sample_rate/2);
wn2=FreqInterval2/(sample_rate/2);

%Filtro Passa-Alto
[b1 a1]=butter(FilterOrder, wn1, 'high'); 
Filtered1=filter(b1,a1,UnfilteredSignal);

%Filtro Passa-Baixo
[b2 a2]=butter(FilterOrder, wn2, 'low'); 
FilteredSignal=filter(b2,a2,Filtered1);

%h1=fvtool(b1,a1);
%h2=fvtool(b2,a2);


%Comparar com Filtro Passa-Banda (permite a passagem das frequências de uma certa
%faixa e rejeita as frequências fora dessa faixa)
%Filtro Passa-Banda
wn=[wn1 wn2];
[b a]=butter(FilterOrder, wn, 'bandpass'); 
h=fvtool(b,a);
BandPassFilteredSignal=filter(b,a,UnfilteredSignal);

%% Plot 

figure;
subplot(4,1,1)
plot(UnfilteredSignal(:,11))
title('Dados Não Filtrados de um canal à escolha')
xlabel('Timepoints');
ylabel('Amplitude');
subplot(4,1,2)
plot(Filtered1(:,11))
title('Aplicação do Filtro Butterworth Passa-Alto de 3ªordem com Fc=0.01 Hz do mesmo canal')
xlabel('Timepoints');
ylabel('Amplitude');
subplot(4,1,3)
plot(FilteredSignal(:,11))
title('Aplicação do Filtro Butterworth Passa-Baixo de 3ªordem com Fc=0.09 Hz do mesmo canal')
xlabel('Timepoints');
ylabel('Amplitude');
subplot(4,1,4);
plot(BandPassFilteredSignal(:,11))
title('Aplicação direta do Filtro Butterworth Passa-Banda de 3ªordem com wn=[0.01, 0.09] Hz do mesmo canal')
xlabel('Timepoints');
ylabel('Amplitude');
end

