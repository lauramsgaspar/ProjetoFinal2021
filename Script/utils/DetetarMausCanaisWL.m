function [ind_excluir ] = DetetarMausCanaisWL( nirs_file, CVThreshold, SNRThreshold )

ind_excluir=[];
for i=1:size(nirs_file.probe.link,1)
    media(i)=mean(nirs_file.data(:,i));
    desviopadrao(i)=std(nirs_file.data(:,i));
    CV(i)=(desviopadrao(i)/media(i));
    SNR(i)= 20*(log10(media(i)/desviopadrao(i)));
    if ((CV(i)> CVThreshold) || (SNR(i)<SNRThreshold))
        t_idx=i;
        ind_excluir(i)=t_idx;
        %newnirs_file.data(:,ind_excluir)=[];
    else
        continue;
    end
end

%% Representação gráfico CV(%) and SNR por canal
figure;
subplot(1,2,1)
plot(CV*100)
title('Coeficiente de Variação(%) por Canal')
ylim([0 20])
xlim([0 size(nirs_file.probe.link,1)])
xlabel('Canal');
ylabel('CV = \sigma / \mu *100')
y=(CVThreshold*100);
line([0,52],[y,y],'Color','red','LineStyle','--');
hold on;

subplot(1,2,2)
plot(SNR)
title('Signal to Noise por Canal')
ylim([0 50])
xlim([0 size(nirs_file.probe.link,1)])
xlabel('Canal');
ylabel('SNR = 20*log10( \mu / \sigma )')
y=SNRThreshold;
line([0,52],[y,y],'Color','red','LineStyle','--');
hold on;

end

