# installazione e caricamento dei pacchetti

install.packages("FactoMineR", type="binary")
install.packages("factoextra")
#install.packages("clustertend")
install.packages("clValid")
install.packages("hopkins")
install.packages("mclust")
install.packages("mixtools")
install.packages(c("dplyr", "ggplot2", "plotly", "caret", "car", "corrplot", "tidyverse"))

#install.packages("tibble", dependencies = TRUE)
#library(tibble)

library(FactoMineR)# data visualization
library(cluster) # applicazione algoritmi di clustering
library(NbClust) # determinare il numero di clusters ottimale
#library(clustertend) # per l'assessment sulla tendenza alla clusterizzazione - deprecated - sostituito da 'hopkins'
library(hopkins)
library(clValid) # per confrontare gli algoritmi di clustering tra loro
library(factoextra) # matrice di dissimilarità
library(fpc) #per le statistiche del clustering come l'indice di Dunn o altri indici
library(mvtnorm)
library(mclust) # Gaussian Mixture Model
library(mixtools)
library(MASS)

library(ggplot2)
library(corrplot)
library(plotly)
library(dplyr)
library(caret)
library(car)
library(GGally)
library(tidyverse)

# importazione e pulizia preliminare (null, NA, ...)

data <- read.csv("penguins_size.csv")
glimpse(data)

# Individuazione eventuali valori NA

valori_null <- is.na(data)
data_null <- data[rowSums(valori_null)>0, ] 
data <- na.omit(data) #eliminazione valori Na

# Individuazione eventuali duplicati

duplicati <- data[duplicated(data) | duplicated(data, fromLast = TRUE), ] # non vengono rilevati duplicati

# pulizia dell'environment
rm(duplicati)
rm(data_null)
rm(valori_null)

# osservazione preliminare del dataset: visualizzazione delle distribuzioni sugli attributi

# grafico a torta per sesso degli esemplari di pinguino
frequenze <- table(data$sex)
pie(frequenze, labels = c("Maschi", "Femmine"), col = c("#77DD77", "#F7D451"))
table(data$sex) # esiste una osservazione con sesso non specificato che viene rimossa
data <- subset(data, sex %in% c("MALE", "FEMALE"))
table(data$sex) # quasi 50/50

# grafico a torta per specie degli esemplari di pinguino
frequenze <- table(data$species)
pie(frequenze, labels = c("Adelie", "Chinstrap", "Gentoo"), col = c("#99D9EA", "violet", "#FF7F50"))
table(data$species) # esiste una osservazione con sesso non specificato che viene rimossa

# grafico a torta per isola degli esemplari di pinguino
frequenze <- table(data$island)
pie(frequenze, labels = c("Adelie", "Chinstrap", "Gentoo"), col = c("#99D9EA", "violet", "#FF7F50"))
table(data$species) # esiste una osservazione con sesso non specificato che viene rimossa

# istogramma della distribuzione per peso del pinguino 

hist(data$body_mass_g, 
     breaks = 80, # Numero di barre
     main = "",
     xlab = "Peso [g]",
     ylab = "Esemplari", 
     col = "#99D9EA")

# istogramma della lunghezza dell'ala
hist(data$flipper_length_mm, col= "#F7D451", main = "", xlab="Lunghezza ala", ylab="Esemplari")

# ulteriore grafico per visualizzare una serie di altre cose tutte insieme:
options(repr.plot.width=12, repr.plot.height=8)
ggpairs(data[,-c(1,2)],aes(color=sex))+theme_bw(base_size = 16) #color-coded for sex  # è un correlation plot

# calcolo della correlazione tra gli attributi:

encode_labels <- function(column) {
  unique_labels <- unique(column)
  labels_to_numbers <- as.numeric(factor(column, levels = unique_labels))
  return(labels_to_numbers)
}

# l'attributo relativo al sesso viene trasformato in numerico (MALE,FEMALE)=(0,1)
data$sex<-as.numeric(as.character(factor(data$sex, c("MALE", "FEMALE"), labels=c(0,1)))) 
data$species <- as.numeric(as.factor(data$species))
data$island <- as.numeric(as.factor(data$island))

cormatrix <- cor(data)
View(cormatrix)

heatmap(cormatrix, 
        symm = TRUE,  # simmetrica
        col = colorRampPalette(c("#99D9EA", "white", "#FF7F50"))(100), 
        main = "Heatmap della correlazione tra attributi")

# il peso dell'esemplare risulta correlato, intuitivamente, alla lunghezza delle ali e alle dimensioni del becco (sia spessore e, in misura maggiore, lunghezza)
# altra correlazione osservabile esiste tra le due dimensioni del becco e lunghezza delle ali. qui si possono effettuare delle considerazioni...
# È noto che la correlazione indichi solo una relazione statistica tra due variabili, e non implichi necessariamente un rapporto di causa-effetto.
# In questo caso, il peso potrebbe essere un fattore confondente che influenza sia la lunghezza delle ali che la lunghezza/spessore del becco, creando una correlazione virtuale tra queste ultime due variabili.

# l'isola di provenienza è correlata con la specie di appartenenza e con le dimensioni
# la specie a sua volta è correlata alle dimensioni dell'esemplare (in termini di peso, lunghezza del becco e lunghezza delle ali, quest'ultima in particolare)

# infine il sesso e peso dell'animale mostrano un'ulteriore correlazione, suggerendo che gli esemplari di peso maggiore siano generalmente maschi.

# con l'eliminazione della variabile relativa alla specie (variabile sulla quale si vuole fare clustering) e quella relativa all'isolad i provenienza (a causa delle correlazioni alte in media che ha con gli altri attributi)
# restano per fare clustering un totale di 5 variabili. il numero è considerato accettabile senza necessità di PCA

remove(encode_labels)
remove(frequenze)
remove(cormatrix)

###################################################################################################################################################

# 1) assessment della tendenza del dataset alla clusterizzazione (ed eventuale preprocessamento/standardizzazione)
# Effettuato con visualizzazione attraverso grafici, statistica di Hopkins e algoritmo VAT (matrice di dissimilarità)

# Un problema significativo nell'analisi dei cluster è che i metodi di clustering restituiranno un risultato anche se i dati non presentano una struttura a cluster ben definita.
#In altre parole, applicare un metodo di clustering a un dataset qualsiasi porterà a una divisione dei dati, perché è per questo che è stato progettato, ma non è detto che tale divisione abbia alcun senso.
# per questo si fa un assessment preliminare dell'idoneità del dataset alla clusterizzazione

backup<-data # creazione di un dataset di backup 
df<-scale(data[,-c(1,2)]) # standardizzazione perché le variabili sono espresse in unità di misura differenti

# ---------------------- #
#     analisi visuale    #
# ---------------------- #

pairs(df, gap =0, pch=16)

# già da questi scatterplot preliminari si individua la presenza di due cluster per lo meno
# per verificare se l'intuizione è corretta, si utilizza rapidamente il kmeans senza prestare troppa attenzione ai dettagli

set.seed(123)
km.temp <- kmeans(df, 2) # prova con 2 cluster che sono quelli che si intuiscono visivamente
cl1 <- km.temp$cluster
pairs(df, gap=0, pch=cl1, col=c( "#77DD77", "#F7D451")[cl1])

# conferma ricevuta: dando in input all'algoritmo k=2 aleatoriamente, nello spazio i due cluster individuati sono visibilmente distinguibili.
# Il kmeans verrà ripetuto più avanti, assieme ad altri metodi, con l'adeguato numero di cluster k in input

remove(km.temp)
remove(cl1)

# prima di procedere alla computazione dellas tatistica di hopkins e algoritmo VAT, una breve parentesi doverosa,
# che verrà ripresa più avanti nel dettaglio, per visualizzare una rappresentazione non più deterministica,
# come quella appena offerta dal kmeans, ma probabilistica, attraverso l'uso dei mixture models.

# ---------------------- #
#      parentesi GMM     #
# ---------------------- #

#GMM_df <- df[,-5]  # senza la variabile sul sesso (anche se è stata trasformata in numerica, continua ad essere categoriale in linea di principio)

# Stima della densità di distribuzione col GMM ad elementi finiti
plot(densityMclust(df), what = "density", type = "hdr", data = df, points.cex = 0.5) #linee con sfumature diverse e punti nello spazio

# a livello grafico, anche dal punto di vista probabilistico continuano ad essere ben distinguibili insiemi aggregati di osservazioni.
# si procede con misure più puntuali della tendenza a clusterizzazione con statistica di Hopkins e algoritmo VAT

# ---------------------- #
#  statistica di Hopkins #
# ---------------------- #

# Valore vicino a 1: Forte aggregazione (cluster ben definiti)
# Valore vicino a 0: Scarsa aggregazione (cluster non ben definiti)

hopkins(df, m=nrow(df)-1) # m=nrows(df)-1 excludes self-comparisons and gives a fairer representation of the average intra-cluster distance. This allows for a more accurate assessment of cluster separation through the Hopkins statistic.
# clustertrend package is deprecated - downloading and using "hopkins" package
# SE DA QUESTO ERRORE: "unused argument (m = nrow(df) - 1)" ATTIVARE LA LIBRARY HOPKINS DI NUOVO E RI-ESEGUIRE. MUST BE MASKED BY SOMETHING

# valore H = 0.997, molto alto: indica una forte aggregazione dei dati all'interno di un cluster. I punti all'interno del cluster sono molto simili tra loro e ben distinti dai punti degli altri cluster.

# ---------------------- #
#     algoritmo VAT      #
# ---------------------- #

# algoritmo VAT (visual assessment of cluster tendency) - usa la matrice di dissimilarità (con distanza euclidea) e "ordina" i quadrati perché siano mostrati graficamente in output all'algoritmo

# calcolo della matrice di dissimilarità
dis_matrix <- dist(df, method = "euclid") # numero più alto = maggiore dissimilarità
round(as.matrix(dis_matrix)[1:6, 1:6],digits = 2)

fviz_dist(dis_matrix, show_labels = FALSE) +
  labs(title = "Ordered Dissimilarity Image (of ODM)") #factoextra package + ggplot2

# ricordando che il rosso indica alta similarità e il blu alta dissimilarità, anche qui si inviduano visualmente dei cluster come nel grafico precedente

# l'algoritmo VAT detecta la tendenza al clustering in forma visuale contando il numero di quadrati rossi lungo la diagonale
# e l'immagine generata dall'algoritmo VAT mette in evidenza l'esistenza di alcuni cluster

# si conclude, a valle di tutte le valutazioni espresse in questa sezione, che il dataset si presta a clusterizzazione e si procede con l'implementazione di algoritmi di clustering.

###############################################################################################################################

##########################
#      CLUSTERING        #
##########################

###############################################################################################################################

# CLAUDIA SARA GIORGIO!
# NELL'R MARKDOWN, PRIMA DI INIZIARE CON LA PRIMA TECNICA, QUI BISOGNA METTERE TUTTA LA TEORIA.
# DESCRIVERE L'ELBOW, IL SILHOUETTE, INIDICE DI DUNN, TUTTI I PASSAGGI RIPETUTI PER OGNI ALGORITMO, PERO' TUTTI QUI ALL'INIZIO

# Clustering model is a notion used to signify what kind of clusters we are trying to identify. The four most common models of clustering methods are hierarchical clustering, k-means clustering, model-based clustering, and density-based clustering:
# - Hierarchical clustering. It creates a hierarchy of clusters, and presents the hierarchy in a dendrogram. This method does not require the number of clusters to be specified at the beginning. Distance connectivity between observations is the measure.
# - k-means clustering. It is also referred to as flat clustering. Unlike hierarchical clustering, it does not create a hierarchy of clusters, and it requires the number of clusters as an input. However, its performance is faster than hierarchical clustering. Distance from mean value of each observation/cluster is the measure.
# - k-medoids
# - Model-based clustering (or Distribution models). Both hierarchical clustering and k-means clustering use a heuristic approach to construct clusters, and do not rely on a formal model. Model-based clustering assumes a data model and applies an EM algorithm to find the most likely model components and the number of clusters. Significance of statistical distribution of variables in the dataset is the measure.
# - Density-based clustering. It constructs clusters in regard to the density measurement. Clusters in this method have a higher density than the remainder of the dataset. Density in data space is the measure.

###############################################################################################################################

# 2) HIERARCHICAL CLUSTERING

# in questo caso specifico !!!NON!!! bisogna stabilire il numero di cluster a priori

# Qual è il più adatto per il caso preso in esame? Divisivo o Agglomerativo?

# Clustering gerarchico agglomerativo:
# - Più adatto a dati con cluster compatti e ben separati.
# - Più efficiente per dataset di grandi dimensioni.
# - I dendrogrammi facilitano l'interpretazione.

# Clustering gerarchico divisivo:
# - Più adatto a dati con cluster irregolari o sovrapposti.
# - Meno sensibile al rumore.
# - Può essere più complesso da interpretare.

# Il risultato della statistica di Hopkins suggerirebbe cluster ben separati, eppure come verrà mostrato più avanti esiste sovrapposizione.
# Questa ragione, e la maggiore tolleranza al rumore dei dati dell'algoritmo diana(), incoraggiano l'utilizzo del clustering divisivo, anche perchè si tratta di un dataset relativamente piccolo ove non è necessaria grande potenza computazionale.

# Divisive Analysis Clustering

res.diana <- diana(dis_matrix) # non vengono considerati altri parametri perché viene data in input una matrice di dissimilarità calcolata precedentemente (metodo di determinazione delle distanze e standardizzazione già effettuati)

#si effettua un taglio nel dendogramma (scelta reiterata di k=6) e si analizzano i cluster che ne fuoriescono
fviz_dend(res.diana, k = 6,
          cex = 0.5, # label size
          k_colors = c("#99D9EA", "#77DD77", "violet", "#F7D451","#FF7F50", "grey"),
          color_labels_by_k = TRUE, # colori per cluster
          rect = TRUE # aggiunta di rettangoli tratteggiati intorno a ciascun cluster
) 

data$diana <- cutree(res.diana, k = 6)
data%>%
  group_by(diana)%>%
  summarise(C_length_mean= mean(culmen_length_mm),
            C_depth_mean= mean(culmen_depth_mm),
            F_length_mean= mean(flipper_length_mm),
            B_mass_mean= mean(body_mass_g),
            Sex=(mean(sex)))

# commento sui cluster: come si evince graficamente, non contengono tutti lo stesso numero di osservazioni, con i primi due che contengono molte piu osservazioni degli ultimi due
# il cluster 4, popolato dai pinguini con le ali più lunghe e i pesi maggiori (le due variabili sono correlate) è a prevalenza maschile (media della variabile sex=0). Si tratta di uno dei 2 gruppi meno popolosi.
# tra i cluster più numerosi invece ce n'è uno che si colloca dal lato opposto (il cluster 2):valore minimo di media dei pesi e ali più corte, a prevalenza femminile stavolta. Rispetto ai pinguini del cluster precedente, questi hanno becchi meno lunghi, ma più spessi.

# altro grafico per visualizzare i cluster nello spazio
fviz_cluster(eclust(df, FUNcluster="diana", k=6, hc_metric="euclidean", hc_method="complete"), dfs, geom = "point")

# result validation

# ---------------------- #
#  distanza cofenetica   #
# ---------------------- #

# verifica del dendogramma con la distanza cofenetica, una misura di similarità utilizzata nel clustering
# Quantifica la distanza tra due punti dati nello spazio dei dati originali in base alla loro distanza nello spazio dei cluster.
# In altre parole dice quanto la distanza tra due punti in un dendrogramma rifletta la loro distanza effettiva nei dati originali.

res.coph <- cophenetic(res.diana) # calcolo della distanza cofenetica
cor(dis_matrix, res.coph) # calcolo della correlazione tra la distanza della matrice di dissimilarità e la distanza cofenetica
# 0.83 la correlazione è alta: bene

# ---------------------- #
#   coeff. Silhouette    #
# ---------------------- #

# Misurano quanto bene ogni punto sia assegnato al suo cluster rispetto ad altri cluster.
# Valori positivi: Indicano una buona assegnazione dei punti ai cluster. Più il valore è vicino a 1, migliore è la separazione del punto dal suo cluster rispetto agli altri cluster.
# Valori negativi: Indicano che il punto potrebbe essere meglio assegnato a un altro cluster. Più il valore è vicino a -1, peggiore è la sua assegnazione.
# Valore 0: Indica che il punto è equidistante da due o più cluster.
# Tuttavia, il valore del coefficiente è sensibile al metodo di calcolo della distanza e alla forma dei cluster.

# silhouette plot (con numero di cluster k=5)
fviz_silhouette(eclust(df, FUNcluster="diana", k=6, hc_metric="euclidean", hc_method="complete"))

# In console compaiono i coefficienti per ogni cluster. Sono tutti prossimi a 0.5, in particolare la media tra tutti i coefficienti è: 0,51, positivo. questo indica una discreta affidabilità del clustering effettuato dall'hierarchical divisive.
# graficamente si osserva come si distribuiscano le silhouette rispetto alla linea tratteggiata in corrispondenza della media, indicando come ad avere una buona assegnazione al cluster 'corretto' sono sopratutto le osservazioni appartenenti ai cluster 4 e 5.
# Il cluster numero 1, che è anche il più popoloso, sembrerebbe essere quello col maggior numero di assegnazioni "sbagliate", o per meglio dire potenzilametne migliorabili. D'altronde si tratta del cluster con minor valore di coefficiente silhouette.

# ---------------------- #
#    indice di Dunn      #
# ---------------------- #

# Il suo scopo è identificare cluster che siano sia compatti (i punti all'interno del cluster sono vicini tra loro) sia ben separati (i punti tra cluster sono distanti).
# Si calcola come il rapporto tra la distanza minima inter-cluster e la distanza massima intra-cluster:
# Compattezza: Considera la distanza massima tra due punti qualsiasi all'interno di un cluster. Rappresenta la dispersione dei punti all'interno di un cluster.
# Separazione: Considera la distanza minima tra due punti qualsiasi in cluster diversi. Rappresenta l'isolamento di ciascun cluster dagli altri.

#Interpretazione del valore:

# 0: I cluster sono completamente sovrapposti.
# 0,5: I cluster sono ben separati, ma c'è un po' di sovrapposizione.
# 1: I cluster sono completamente separati e non c'è sovrapposizione.

diana.stats <- cluster.stats(dis_matrix, data$diana)
diana.stats$dunn

# l'indice è pari a 0.18 approssimativamente. Contraddice in parte quanto appena derivato dal metodo silhouette,
# tuttavia c'è da tenere presente che in corrispondenza di cluster di dimensioni diverse (anche se non drasticamente, è il caso in esame), l'indice di dunn potrebbe fornire un'indicazione errata.

# ---------------------- #
# confronto con la specie#
# ---------------------- #

table(data$diana, backup$species)

x<-backup
x$cluster <- as.factor(data$diana)
options(repr.plot.width = 15, repr.plot.height = 7)

ggplot(x, aes(flipper_length_mm, body_mass_g)) + 
  geom_point(aes(color=cluster), size=3) +
  facet_wrap(~ species, ncol=3) 

remove(x)

# oltre ai risultati appena ottenuti, viene mostrato anche il grafico relativo all'esecuzione del clustering con k=4:

data$diana4 <- cutree(res.diana, k = 4)
data%>%
  group_by(diana4)%>%
  summarise(C_length_mean= mean(culmen_length_mm),
            C_depth_mean= mean(culmen_depth_mm),
            F_length_mean= mean(flipper_length_mm),
            B_mass_mean= mean(body_mass_g),
            Sex=(mean(sex)))

x<-backup
x$cluster <- as.factor(data$diana4)
options(repr.plot.width = 15, repr.plot.height = 7)
ggplot(x, aes(flipper_length_mm, body_mass_g)) + 
  geom_point(aes(color=cluster), size=3) +
  facet_wrap(~ species, ncol=3) 

# viene mantenuto il k=6 soprattutto a seguito del confronto tra l'assegnazione delle osservazioni effettuata dal clustering gerarchico divisivo e l'etichetta nota sulla specie (nel dataset originale)
# si noti come il grafico generato con numero di cluster pari a 4 mostri una più scarsa separazione tra i cluster, confrontato con l'anteriore che invece mostrava maggiore precisione e consentiva una distinzione più netta dei pinguini nelle tre specie

remove(x)
data<-data[,-9]

# si procede con k-means e k-medoids per valutare se possano dare risutati migliori rispetto al clustering gerarchico.

###############################################################################################################################

# 2) K-MEANS

# L'algoritmo K-Means è un metodo di clustering iterativo che raggruppa i dati in base alla similarità. Inizia con centroidi casuali e iterativamente assegna i punti dati ai cluster più vicini, ricalcolando i centroidi in base alle osservazioni assegnate. Questo processo si ripete fino a quando la configurazione dei cluster converge.Centroidi iniziali (centri dei cluster): Vengono scelti casualmente dei centroidi iniziali per rappresentare il centro di ogni cluster.
# 1) Assegnazione dei dati ai cluster: Ogni osservazione (punto dati) viene assegnata al cluster con il centroide più vicino, utilizzando una misura di distanza (come la distanza euclidea).
# 2) Aggiornamento dei centroidi: I centroidi vengono ricalcolati come media (centroide) di tutte le osservazioni che appartengono a ciascun cluster.
# 3) Iterazione: I passaggi 2 e 3 vengono ripetuti più volte. L'obiettivo è minimizzare la varianza totale all'interno di ogni cluster. Questo processo iterativo continua finché non ci sono più cambiamenti significativi nella posizione dei centroidi, oppure fino a raggiungere un livello di tolleranza predefinito.

# A differenza del caso precedente, qui invece BISOGNA stabilire il numero di clusters a priori. Per far ciò, si puo' ricorrere a diverse metodologie.

# A seguire sono presentati i risultati dell'analisi condotta con i metodi: Elbow, Silhouette e GAP Statistic

# ------------------------ #
#           ELBOW          #
# ------------------------ #

fviz_nbclust(df, kmeans, nstart=25, method="wss") + 
  geom_vline(xintercept = 6, linetype = 2) +
  labs(title= "Elbow method per il numero ottimale di clusters k")

# in questo caso il metodo elbow da risultati un po' ambigui perché graficamente non si individua un "gomito" brusco nella curva
# si sceglie K=6 perché è il punto a partire del quale ulteriori clusterizzazioni non portino grandi miglioramenti dal putno di vista del wss (whitin-cluster sum of squares)

# ------------------------ #
#         SILHOUETTE       #
# ------------------------ #

fviz_nbclust(df, kmeans, method="silhouette") +
  geom_vline(xintercept = 4, linetype = 2) +
  labs(title= "Silhouette method per il numero ottimale di clusters k")

# per k = 6 si ha il valore medio del coefficiente di silhouette più alto, in questo punto i cluster sono più compatti e meglio separati
# da k=6 in poi, all'aumentare di k, l'andamento della curva silhouette è quasi strettametne decrescente
# il numero di cluster = 6 è considerato il numero ottimale di cluster per il dataset.

# ------------------------ #
#       GAP Statistic      #
# ------------------------ #

fviz_nbclust(df, kmeans, nstart =25, method = "gap_stat", nboot = 50) +
  labs(title= "Gap Statistic method per il numero ottimale di clusters k")

# all'aumentare di k, anche la gap statistic continua a crescere, con un unico piccolo punto di inflessione in corrispondenza di k=6, che verrà preso come numero ideale di cluster

# ulteriore valutazione con il pacchetto NbClust che determina il numero ottimale di cluster in base al metodo utilizzato e il tipo di distanza tra le osservazioni (se euclidea, manhattan, canberra etc.)

# ------------------------ #
#           NbClust        #
# ------------------------ #

nb <- NbClust(data = df, diss = NULL, distance = "euclidean",
              min.nc = 2, max.nc = 10, method = "kmeans")

remove(nb)
# con il kmeans per esempio, per la regola della maggioranza, il numero di cluster ottimale identificato è pari a 4 con distanza euclidea. 

# risultati a confronto:

# elbow: 6 cluster
# silhouette: 4 cluster
# gap statistic: 6 cluster
# NbClust: 4 cluster

# viene effettuato il kmeans con k=6 

set.seed(123)
res.km <- eclust(df, "kmeans", k=6, nstart=50, graph=FALSE) # 6 è il numero di cluster scelto tra elbow, silhouette e gap statistic

data$kmeans <- res.km$cluster # aggiunta dell'ultima colonna "kmeans" al dataframe iniziale (NON standardizzato perché serve osservare i dati originali)
data%>%
  group_by(kmeans)%>%
  summarise(C_length_mean= mean(culmen_length_mm),
            C_depth_mean= mean(culmen_depth_mm),
            F_length_mean= mean(flipper_length_mm),
            B_mass_mean= mean(body_mass_g),
            Sex=(mean(sex)),
            Count=(length(kmeans))) #aggiunta di un'ultima colonna con il numero di osservazioni appartenenti a ciascun cluster

# commento sui cluster appena ottenuti:
# il cluster più popoloso(78 osservazioni), il numero 5, si compone sopratttutto esemplari femmine (media di "sex" 1), di peso relativamente basso (il più basso in media tra tutti i cluster, insieme alla lunghezza dell'ala e la lunghezza del becco - queste variabili sono correlate per cui non ci sorprende)
# oppure il cluster con meno esemplari (36 osservazioni), il numero 1, che contiene soprattutto esemplari maschi, con dimensioni del becco maggiori rispetto a tutti gli altri cluster (sia spessore che lunghezza).
# gli esemplari con peso in media maggiore si trovano invece nel cluster numero 6, dove ci sono esemplari maschi in maggioranza, con la lunghezza dell'ala massima del dataset (correlata al peso).

fviz_cluster(res.km, geom="point", ellipse.type="norm",
                          palette=c("#99D9EA", "#77DD77", "violet", "#F7D451","#FF7F50", "orange"), 
                          ggtheme=theme_minimal())

# result validation

# ---------------------- #
#   coeff. Silhouette    #
# ---------------------- #

KMEANS1_S <- fviz_silhouette(eclust(df, FUNcluster="kmeans", k=6, hc_metric="euclidean", hc_method="complete", graph=FALSE))
KMEANS1_S

# ---------------------- #
#    indice di Dunn      #
# ---------------------- #

km.stats <- cluster.stats(dis_matrix, data$kmeans)
km.stats$dunn

# l'indice di dunn è equivalente a quello dell'algoritmo diana().

# ---------------------- #
# confronto con la specie#
# ---------------------- #

table(res.km$cluster, backup$species)

x<-backup
x$cluster <- as.factor(res.km$cluster)
options(repr.plot.width = 15, repr.plot.height = 7)

ggplot(x, aes(flipper_length_mm, body_mass_g)) + 
  geom_point(aes(color=cluster), size=3) +
  facet_wrap(~ species, ncol=3) 

remove(x)

###############################################################################################################################

# 2) K-MEDOIDS

# A differenza del k-means, il k-medoids usa centroidi "reali": non calcolati come medie, ma scelti tra i punti effettivamente presenti nel dataset come centroide. Questo lo rende più resistente a valori anomali e rumore.
# Altra differenza col k-means è che non si basa sulla distanza euclidea (linea retta), ma utilizza una misura di dissimilarità più generale, come la distanza di Manhattan o la distanza di Chebyshev.
# Offre infine maggiore flessibilità nella scelta del numero di cluster (k) e nella selezione dei medoids iniziali. Viene utilizzato l'algoritmo PAM (partitioning around medoids) agevolmente nelle sezioni a seguire, visto che si tratta di un dataset con meno di 2000 osservazioni 
# (tuttavia, si tene presente che una bassa dimensionalità potrebbe portare ad una scelta dei medoidi equivalente a quella dei centroidi appena effettuata con il k-means, portando di conseguenza agli stessi risultati. Di fatto, si evincerà alla fine dell'implementazione che ciò si verifica)

# Per ottenere il numero ottimale di clusters k si implementano gli stessi metodi adoperati fino ad ora: Elbow, Silhouette e GAP Statistic

# ------------------------ #
#           ELBOW          #
# ------------------------ #

fviz_nbclust(df, pam, nstart=25, method="wss") + 
  geom_vline(xintercept = 6, linetype = 2) +
  labs(title= "Elbow method per il numero ottimale di clusters k")

# in questo caso il metodo elbow da risultati un po' ambigui perché graficamente non si individua un "gomito" brusco nella curva
# si sceglie K=6 perché è il punto a partire del quale ulteriori clusterizzazioni non portino grandi miglioramenti dal putno di vista del wss (whitin-cluster sum of squares)

# ------------------------ #
#         SILHOUETTE       #
# ------------------------ #

#library(cluster)
fviz_nbclust(df, pam, method="silhouette") +
  geom_vline(xintercept = 6, linetype = 2) +
  labs(title= "Silhouette method per il numero ottimale di clusters k")

# per k = 6 si ha il valore medio del coefficiente di silhouette più alto, in questo punto i cluster sono più compatti e meglio separati
# il numero di cluster corrispondente a questo punto è considerato il numero ottimale di cluster per il dataset.

# ------------------------ #
#       GAP Statistic      #
# ------------------------ #

fviz_nbclust(df, pam, nstart =5, method = "gap_stat", nboot = 5) +
  labs(title= "Gap Statistic method per il numero ottimale di clusters k") # questo impiegherà del tempo

# unico punto di inflessione in corrispondenza di k=6, che verrà preso come numero ideale di cluster

# ------------------------ #
#           NbClust        #
# ------------------------ #

nb <- NbClust(data = df, diss = NULL, distance = "euclidean",
              min.nc = 2, max.nc = 10, method = "kmeans")
remove(nb)

# risultati a confronto:

# elbow: 6 cluster
# silhouette: 6 cluster
# gap statistic: 6 cluster
# nbClust: 4

# si decide di scegliere k= 6 come numero di clusters

res.pam <- pam(df, 6)

# si aggiunge questa clusterizzazione data dal PAM al dataset df standardizzato a partire dall'originale (si aggiunge una colonna "kmedoids")

data$kmedoids <- res.pam$cluster # aggiunta dell'ultima colonna "kmedoids" al dataframe iniziale (NON standardizzato perché serve osservare i dati originali)
data%>%
  group_by(kmedoids)%>%
  summarise(C_length_mean= mean(culmen_length_mm),
            C_depth_mean= mean(culmen_depth_mm),
            F_length_mean= mean(flipper_length_mm),
            B_mass_mean= mean(body_mass_g),
            Sex=(mean(sex)),
            Count=(length(kmedoids))) # aggiunta di una colonna col conto del numero di osservazioni 

# i cluster ottenuti visivamente

fviz_cluster(res.pam, geom="point", ellipse.type="norm",
             palette=c("#99D9EA", "#77DD77", "violet", "#F7D451","#FF7F50", "orange"),
             ggtheme=theme_minimal())

# Sia da un'analisi dei cluster ottenuti che, banalmente, visivamente osservando il grafico si evince come il kmeans e il kmedoids abbiano dato, per k=6, risultati quasi uguali.
# Il fatto che i risultati siano pressoché uguali non è anomalo in determinate condizioni, corrispondenti a quelle del caso in esame: 

# in dataset "piccoli", la scelta di medoids vs centroidi ha generalmente un impatto minore, cosa che rende i risultati più simili.
# Altro punto importante è la scelta del numero di cluster k=6. Il fatto che entrambi gli algoritmi usino lo stesso k contribuisce a farli convergere verso la stessa soluzione.
# Distribuzione uniforme: Se i dati sono uniformemente distribuiti nello spazio delle caratteristiche, entrambi gli algoritmi potrebbero assegnare i punti ai cluster in modo simile.
# Infine, ad influire è la specificità del dataset. La struttura di quest'ultimo potrebbe coincidere con le proprietà che k-means e k-medoids cercano di ottimizzare, portandoli a convergere verso la stessa soluzione.

# result validation

# ---------------------- #
#   coeff. Silhouette    #
# ---------------------- #

fviz_silhouette(eclust(df, FUNcluster="pam", k=4, hc_metric="euclidean", hc_method="complete"))

# qui la media tra tutti i coefficienti è 0,51. (kmeans era 0.45, diana 0.51) questo indica una discreta affidabilità del clustering effettuato dal kmedodis.

# ---------------------- #
#    indice di Dunn      #
# ---------------------- #

kmed.stats <- cluster.stats(dis_matrix, data$kmedoids)
kmed.stats$dunn

# peggiore dei casi precedenti, anche qui l'indice di dunn non è particolarmente positivo, valore piccolo pari a 0,17

# ---------------------- #
# confronto con la specie#
# ---------------------- #

table(res.pam$cluster, backup$species)

x<-backup
x$cluster <- as.factor(res.pam$cluster)
options(repr.plot.width = 15, repr.plot.height = 7)
ggplot(x, aes(flipper_length_mm, body_mass_g)) + 
  geom_point(aes(color=cluster), size=3) +
  facet_wrap(~ species, ncol=3) 

remove(x)

# il risultato ottenuto confrontando l'esito del k-medoids con le etichette di cui si dispone nel dataaset originale non differisce dal caso del k-means con k=6.

###############################################################################################################################

# +++++++++++++++++++++++++++++++++++++++ #
# Confronto tra i tre algoritmi adoperati #
# +++++++++++++++++++++++++++++++++++++++ #

# innanzitutto si mettono a confronto connettività, indice di dunn e metodo silhouette, con un confronto interno.

# il connecivity value, che non era stato descritto fino ad ora
# Un valore di connettività vicino a 0 suggerisce un cluster ben separato, in cui la maggior parte dei punti dati ha vicini dello stesso cluster tra i loro k vicini più prossimi.
#Al contrario, un valore di connettività più alto (lontano da 0) indica un cluster scarsamente separato, in cui i punti dati potrebbero avere vicini provenienti da altri cluster tra i loro k vicini più prossimi.

clmethods <- c("diana","kmeans","pam")
internal_val <- clValid(df, nClust = 4:6, clMethods = clmethods, validation = "internal")
summary(internal_val)

# ciò che deriva da questo confronto è inequivocabile, a riportare i migliori valori per i tre indici è il metodo gerarchico divisivo. In particolare, per numero di cluster pari a 6 (silhouette) e 4 (dunn)
# Un indice di Dunn e Silhouette elevati suggeriscono cluster ben separati e compatti da una prospettiva gerarchica.
# per complementare questo confronto si ricorre ulteriormente alle misure di stabilità

## Misure di stabilità

#Le misure di stabilità sono una variante speciale delle misure interne che valutano la stabilità di un risultato di clustering confrontandolo con i cluster ottenuti rimuovendo una colonna alla volta.
# Queste misure includono la proporzione media di non sovrapposizione (APN), la distanza media (AD), la distanza media tra medie (ADM) e la misura di merito (FOM).

#APN (Average Proportion of Non-overlap): Misura la percentuale media di osservazioni non assegnate allo stesso cluster sia nella condizione originale che in quella con una colonna rimossa.
#AD (Average Distance): Calcola la distanza media tra osservazioni assegnate allo stesso cluster in entrambi i casi (originale e con colonna rimossa).
#ADM (Average Distance between Means): Misura la distanza media tra i centri dei cluster per le osservazioni assegnate allo stesso cluster in entrambi i casi.
# FOM (Figure of Merit): Calcola la varianza intra-cluster media della colonna eliminata, considerando il clustering basato sulle colonne rimanenti (non eliminate).

#Considerazioni generali:
#  La media viene calcolata su tutte le colonne eliminate.
#Tutte le misure dovrebbero essere minimizzate per indicare una maggiore stabilità.
#In sintesi, queste misure di stabilità aiutano a valutare quanto i risultati di clustering rimangono affidabili anche se si eliminano singole colonne dal dataset. Misure più basse indicano una maggiore stabilità del clustering.

stab <- clValid(df, nClust = 4:6, clMethods = clmethods, validation = "stability")
summary(stab)

# KMeans: FOM elevato indica una buona concordanza generale tra diverse soluzioni di assegnazione di clustering 
# KMedoids: APN, AD e ADM elevati suggeriscono cluster ben separati, basandosi sulle distanze tra i medoidi. Tuttavia,l'algoritmo PAM è sensibile all'inizializzazione, e questi due valori sono relativi alla specifica inizializzazione del caso preso in esame.

# In generale, L'indice di Dunn, il metodo Silhouette, FOM, APN, AD e ADM misurano tutti aspetti differenti della qualità del clustering. Non esiste una risposta definitiva univoca su quale algoritmo tra i tre messi a confronto sia il più adatto al dataset esamianto, almeno non basansosi esclusivamente su queste metriche di valutazione senza informazioni addizionali.
# Infatti, Un algoritmo potrebbe avere un valore alto per una metrica e un valore più basso per un'altra. Spesso c'è un compromesso tra compattezza (similarità all'interno del cluster) e separazione (differenza tra cluster).

# Ad essere determinante nella scelta, è quale sia l'obiettivo dell'analisi: con le operazioni di clustering su questo dataset l'intenzione era quella di identificare l'eventuale esistenza di clusters quanto più definiti, che enfatizzassero le differenze tra i diversi esemplari di pinguini nel dataset per raggiungere una classificazione in specie diverse.
# Si valuta dunque che il kmedoids effettuato con k=6 sia quello più attendibile per il raggiungimento degli obiettivi relativi a questo studio.

remove(clmethods)
remove(internal_val)
remove(stab)

##########################################################################################################################################################

##################
#####  PCA   #####
##################

# PCA sull'algoritmo di clustering più adeguato al dataset per capire se c'è un miglioramento 

res.pca <- PCA(scale(backup[,-c(1,2)]), graph = FALSE)
print(res.pca)

# Calcolo delle componenti principali

cp <- get_eigenvalue(res.pca)
cp

# vengono mostrate le componenti principali individuate con la PCA. 
# visti i risultati ottenuti, si decide di scegliere solo le prime due componenti principali per due motivi:
# il primo è che le ultime tre non hanno autovalore superiore a 1, il secondo è che, combinandole, la varianza percentuale cumulata arriva a un valore soddisfacente di 85% della rappresentazione dello spazio di osservazioni originale

# Anche dal punto di vista grafico si può raggiungere la stessa conclusione:

# Scree Plot: un metodo visuale per determinare il numero di componenti principali
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

# a seguire, vengono mostrate le correlazioni tra le variabili originali e le componenti principali, con alcuni grafici

# correlation plot
var <- get_pca_var(res.pca)  
var$coord

fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#99D9EA", "#F7D451", "#FF7F50")
)

round(var$contrib,digits = 3)
corrplot(var$contrib, is.corr=FALSE)

# si evince ancora più chiaramente come le prime due componenti siano già sufficienti a "coprire" lo spazio originale.
# la seconda in particolare è rappresentativa in particolare dello spessore del becco e del sesso. la prima, anche se in maniera meno marcata, è influenzata quasi da tutti gli attributi originali

# nello specifico:
# Autovalore 1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# Autovalore 2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)

# creazione nuovo dataset
df_pca <- as.data.frame(res.pca$ind$coord)

df_pca<-df_pca[,-c(3,4,5)]
names(df_pca)[1] <- "Dim1"
names(df_pca)[2] <- "Dim2"

remove(var)
remove(cp)

# osservazione preliminare del nuovo dataset:
pairs(df_pca, gap =0, pch=16)

# già da questi scatterplot preliminari si individua la presenza di due cluster per lo meno
# per verificare se l'intuizione è corretta, si utilizza rapidamente il kmeans senza prestare troppa attenzione ai dettagli

set.seed(123)
km.temp <- kmeans(df_pca, 4) # prova con 4 cluster che sono quelli che si intuiscono visivamente
cl1 <- km.temp$cluster
pairs(df_pca, gap=0, pch=cl1, col=c( "#77DD77", "#F7D451", "violet", "lightblue")[cl1])

# conferma ricevuta: dando in input all'algoritmo k=4 aleatoriamente, nello spazio i due cluster individuati sono visibilmente distinguibili.
# Il kmeans verrà ripetuto più avanti, assieme ad altri metodi, con l'adeguato numero di cluster k in input

remove(km.temp)
remove(cl1)

###############################################################################################################################
###############################################################################################################################
###############################################################################################################################

# applicazione del diana() al nuovo dataset a seguito della PCA

# calcolo della nuova matrice di dissimilarità
dis_matrix_pca <- dist(df_pca, method = "euclid") # numero più alto = maggiore dissimilarità
round(as.matrix(dis_matrix_pca)[1:6, 1:6],digits = 2)

res.diana_pca <- diana(dis_matrix_pca) 

#si effettua un taglio nel dendogramma (scelta reiterata di k=6) e si analizzano i cluster che ne fuoriescono
fviz_dend(res.diana_pca, k = 6,
          cex = 0.5, # label size
          k_colors = c("#99D9EA", "#77DD77", "violet","#FF7F50", "#F7D451", "orange"),
          color_labels_by_k = TRUE, # colori per cluster
          rect = TRUE # aggiunta di rettangoli tratteggiati intorno a ciascun cluster
) 

df_pca$diana_pca <- cutree(res.diana_pca, k = 6)
df_pca%>%
  group_by(diana_pca)%>%
  summarise(Dim1=mean(Dim1),
            Dim2=mean(Dim2))

fviz_cluster(eclust(df_pca, FUNcluster="diana", k=6, hc_metric="euclidean", hc_method="complete"), dfs, geom = "point")

# result validation

# ---------------------- #
#  distanza cofenetica   #
# ---------------------- #

res.coph_pca <- cophenetic(res.diana_pca) 
cor(dis_matrix_pca, res.coph_pca) 

# ---------------------- #
#   coeff. Silhouette    #
# ---------------------- #

# silhouette plot (con numero di cluster k=6)
fviz_silhouette(eclust(df_pca, FUNcluster="diana", k=6, hc_metric="euclidean", hc_method="complete"))

# ---------------------- #
#    indice di Dunn      #
# ---------------------- #

diana.stats_pca <- cluster.stats(dis_matrix_pca, df_pca$diana_pca)
diana.stats_pca$dunn

# ---------------------- #
# confronto con la specie#
# ---------------------- #

table(df_pca$diana_pca, backup$species)

x <- cbind(df_pca[,-3], data$species)
x$cluster <- as.factor(df_pca$diana_pca)
options(repr.plot.width = 15, repr.plot.height = 7)

ggplot(x, aes(Dim1, Dim2)) + 
  geom_point(aes(color=cluster), size=3) +
  facet_wrap(~ data$species, ncol=3) 

remove(x)

##############################################################################################################################
# applicazione del kmeans() al nuovo dataset a seguito della PCA
##############################################################################################################################

df_pca<-df_pca[,-3]

# ------------------------ #
#           ELBOW          #
# ------------------------ #

fviz_nbclust(df_pca, kmeans, nstart=25, method="wss") + 
  geom_vline(xintercept = 4, linetype = 2) +
  labs(title= "Elbow method per il numero ottimale di clusters k")

# ------------------------ #
#         SILHOUETTE       #
# ------------------------ #

fviz_nbclust(df_pca, kmeans, method="silhouette") +
  geom_vline(xintercept = 4, linetype = 2) +
  labs(title= "Silhouette method per il numero ottimale di clusters k")

# ------------------------ #
#       GAP Statistic      #
# ------------------------ #

fviz_nbclust(df_pca, kmeans, nstart =25, method = "gap_stat", nboot = 50) +
  labs(title= "Gap Statistic method per il numero ottimale di clusters k")

# ------------------------ #
#          NbClust         #
# ------------------------ #

nb <- NbClust(data = df_pca, diss = NULL, distance = "euclidean",
              min.nc = 2, max.nc = 10, method = "kmeans")

# risultati a confronto:

# elbow: 4 cluster
# silhouette: 4 cluster
# gap statistic: 6 cluster
# NbClust: 4 cluster

set.seed(123)
res.km_pca <- eclust(df_pca, "kmeans", k=4, nstart=50, graph=FALSE)

df_pca$kmeans_pca <- res.km_pca$cluster # aggiunta dell'ultima colonna "kmeans_pca" al dataframe iniziale (NON standardizzato perché serve osservare i dati originali)
df_pca%>%
  group_by(kmeans_pca)%>%
  summarise(Dim1=mean(Dim1),
            Dim2=mean(Dim2))

fviz_cluster(res.km_pca, geom="point", ellipse.type="norm",
                    palette=c("#99D9EA", "#77DD77", "violet", "#F7D451"), 
                    ggtheme=theme_minimal())

# result validation

# ---------------------- #
#   coeff. Silhouette    #
# ---------------------- #

KMEANS2_S <- fviz_silhouette(eclust(df_pca, FUNcluster="kmeans", k=4, hc_metric="euclidean", hc_method="complete", graph=FALSE))
KMEANS2_S

# ---------------------- #
#    indice di Dunn      #
# ---------------------- #

km.stats_pca <- cluster.stats(dis_matrix_pca, df_pca$kmeans_pca)
km.stats_pca$dunn

# l'indice di dunn è equivalente a quello dell'algoritmo diana().

# ---------------------- #
# confronto con la specie#
# ---------------------- #

table(df_pca$kmeans_pca, backup$species)

x <- cbind(df_pca[,-3], data$species)
x$cluster <- as.factor(df_pca$kmeans_pca)
options(repr.plot.width = 15, repr.plot.height = 7)

ggplot(x, aes(Dim1, Dim2)) + 
  geom_point(aes(color=cluster), size=3) +
  facet_wrap(~ data$species, ncol=3) 

remove(x)

###############################################################################################################################

# +++++++++++++++++++++++++++++++++++++++ #
# Confronto tra i 2   algoritmi adoperati #
# +++++++++++++++++++++++++++++++++++++++ #

clmethods_pca <- c("diana","kmeans")
internal_val <- clValid(df_pca, nClust = 4:6, clMethods = clmethods_pca, validation = "internal")
summary(internal_val)

# si esegue diana() un'altra volta con quattro cluster:

# taglio nel dendogramma a k=4
fviz_dend(res.diana_pca, k = 4,
          cex = 0.5, # label size
          k_colors = c("#99D9EA", "#77DD77", "violet","#FF7F50"),
          color_labels_by_k = TRUE, # colori per cluster
          rect = TRUE # aggiunta di rettangoli tratteggiati intorno a ciascun cluster
) 

df_pca$diana_pca <- cutree(res.diana_pca, k = 4)
df_pca%>%
  group_by(diana_pca)%>%
  summarise(Dim1=mean(Dim1),
            Dim2=mean(Dim2))

fviz_cluster(eclust(df_pca, FUNcluster="diana", k=4, hc_metric="euclidean", hc_method="complete"), dfs, geom = "point")

# result validation
# ---------------------- #
#  distanza cofenetica   #
# ---------------------- #
res.coph_pca <- cophenetic(res.diana_pca) 
cor(dis_matrix_pca, res.coph_pca) 
# ---------------------- #
#   coeff. Silhouette    #
# ---------------------- #
# silhouette plot (con numero di cluster k=4)
fviz_silhouette(eclust(df_pca, FUNcluster="diana", k=4, hc_metric="euclidean", hc_method="complete"))
# ---------------------- #
#    indice di Dunn      #
# ---------------------- #
diana.stats_pca <- cluster.stats(dis_matrix_pca, df_pca$diana_pca)
diana.stats_pca$dunn
