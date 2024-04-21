# Space-time-based source separation

Corso: Sound Analysis, Synthesis and Processing

Autori: Riccardo Pomarico, Alice Marazzi.

Il progetto è stato realizzato in Matlab e si concentra sulla separazione delle sorgenti acustiche utilizzando due segnali misti acquisiti da microfoni ideali distanti 9 cm l'uno dall'altro. Ogni segnale contiene una miscela di 3 segnali vocali provenienti da sorgenti posizionate ad angoli specifici rispetto alla direzione dei microfoni. L'obiettivo è separare queste sorgenti utilizzando un mascheramento binario applicato alla trasformata di Fourier a breve termine (STFT) di uno dei segnali.
Per ottenere le maschere binarie, abbiamo estratto un vettore di caratteristiche per ogni posizione tempo-frequenza della STFT e applicato un algoritmo di clustering con 3 cluster. Le maschere risultanti sono state quindi utilizzate per stimare i segnali sorgente.
Sono stati generati e confrontati diversi grafici, inclusi spettrogrammi logaritmici dei segnali misti e dei segnali sorgente veri e stimati, oltre a maschere binarie e plot di densità delle caratteristiche estratte.
