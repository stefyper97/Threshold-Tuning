# Threshold-Tuning
Code to fine tuning the threshold on the ALPIDE 

Nel file logbook_Thrscan.txt sono contenute le informazioni sui test fatti: nella prima colonna va inserito il numero dell'HIC testato, nella seconda i valore di Vcasn, nella terza il valore di Ithr e nella quarta il Path alla cartella thrscan_analysis_output.
I file sono ordinati in modo da poter mettere le acquisizioni più recenti nelle prime righe del file, inoltre per semplicità in questo momento sono inserite in modo che il parametro che varia sia crescente: 
es:
HIC        Vcasn          Ithr          Path
27          50                40           HIC_027/ITHR_40-70/20211007_131809_thrscan_HIC_027_dark_VCASN_50_ITHR_40/thrscan_analysis_output/
27          50                45           HIC_027/ITHR_40-70/20211007_134258_thrscan_HIC_027_dark_VCASN_50_ITHR_45/thrscan_analysis_output/

eccetera, eccetera.
Nel file Threshold_ithr_vcasn.C sono implementate due macro: Threshold(const char* directory) che per ogni acquisizione crea un vector di TH1D* dove, per ciascun chip, viene creato un istogramma per le soglie e uno per il noise entrambi interpolati con una gaussiana. Tutti gli istogrammi e i canvas vengono salvati in un ROOTfile chiamato Threshold_noise.root creato all'interno della cartella thrscan_analysis_output. In quella cartella inoltre viene creato un file txt che contiene i risultati del fit della soglia e del rumore per ciascun chip (fit_data.txt). i dati memorizzati sono: ChipID  Soglia  Sigma(della gaussiana)  Rumore   Sigma(della gaussiana del noise):
es:
ChipID      Threshold       Sigma(Thr)       Noise           Sigma(Noise)
0x70         225.677           22.439             4.55856       1.23757
0x71         242.221          19.2319            4.65675       1.23704

eccetera. eccetera.

Il "main" della macro Threshold_ithr_vcasn.C è la funzione Threshold_ithr_vcasn(const char *dirFile) che prende in ingresso il file logbook_Thrscan.txt e legge riga per riga il file eseguendo la funzione Threshold nel giusto Path, e salvando le variabili in un TTree che sarà memorizzato in un TFile chiamato Threshold_Parameters.root. Il TTree contiene: i valori di Ithr e Vcasn, il numero dell'HIC testato, il Path, 4 array contenenti la soglia, la dispersione della soglia, il rumore e la dispersione del rumore per ciascun chip (array di 10 elementi), valor medio della soglia sull'HIC e rispettiva deviazione standard, questo per ogni acquisizione contenuta in logbook_Thrscan.txt.

Lanciando i comandi:
root -l
.L Threshold_ithr_vcasn.C+
Threshold_ithr_vcasn("logbook_Thrscan.txt");

si ottiene di far girare la macro Threshold su tutte le cartelle inserite nel file logbook_Thrscan.txt e di avere le variabili memorizzate in un TTree all'interno del file Threshold_Parameters.root (creato nella cartella dove si esegue tutto quanto).
