
##########################
#     MIXTURE MODELS     #
##########################

colnames(df)
colnames(backup)

# rimozione della colonna 'sex' prima di dare il dataset in pasto ai mixture models
df_mm <- as.data.frame(df[,-5])
# df_mm <- as.data.frame(scale(backup[,-c(1,7)]))

glimpse(df_mm)

# nel dataset che verrà utilizzato mancano anche le informazioni su specie e isola di appartenenza

# codice per includerle invece:

#df_mm <- df[,-7] # rimozione colonna 'sex'
#df_mm$species <- as.numeric(as.factor(df_mm$species))
#df_mm$island <- as.numeric(as.factor(df_mm$island))
#df_mm <- scale(df_mm)
#df_mm <- df_mm[,-1]
#df_mm <- cbind(df_mm, backup[1])
#df_mm$species <- as.numeric(as.factor(df_mm$species))
#colnames(df_mm)[6] <- "Classes"

# prima di procedere si osservino le distribuzioni approssimate delle variabili del dataset

source("EM2_plot.r")
source("mixt2_univ.r")

em2_weight <- EM2_plot(df_mm$body_mass_g, "Distribuzione delle osservazioni per peso")     
em2_culmen_length <- EM2_plot(df_mm$culmen_length_mm, "Distribuzione della lunghezza del becco")  
em2_culmen_depth <- EM2_plot(df_mm$culmen_depth_mm, "Distribuzione dello spessore del becco")  
em2_flipper_length <- EM2_plot(df_mm$flipper_length, "Distribuzione della lunghezza delle ali") 

####################################################################################################################

# Applicazione dei mixture models

multiv_model <- Mclust(df_mm)
summary(multiv_model, parameters=TRUE)

# Risultato:
# Mclust VEE (ellipsoidal, equal shape and orientation) model with 4 components: 

# log-likelihood   n df       BIC       ICL
# -1126.211      333 32 -2438.283 -2477.667

# Clustering table:
#  1   2   3   4 
#  149  65  61  58

plot(multiv_model, what = "BIC")

# Il mixture model Mclust individua autonomamente il metodo più efficace per clusterizzare le osservazioni nel VEE
# con quattro componenti. Il grafico del BIC conferma graficamente questa scelta in corrispondenza del 4 sull'asse delle ascisse.

plot(multiv_model, what = "density", type = "hdr", data = df_mm, points.cex = 0.5) 

# type = "hdr" permette di ottenere un grafico di densità: più è definita la curva al perimetro della forma, più è densa la distribuzione dei dati.

plot(multiv_model, what="classification", points.cex = 0.5)

# in questo grafico i punti nel piano sono classificati con colori e forme diverse sulla base del cluster al quale il modello li ha assegnati

plot (multiv_model, what = "uncertainty", points.cex =26)

# quest'ultimo grafico permette di visualizzare l'incertezza associata al modello adottato

###############################################################################################################################

# valutazione della bontà del clustering col metodo silhouette

sil <- silhouette(multiv_model$classification, dist(df_mm))
summary(sil)
plot(sil)

###############################################################################################################################

library(stats) #capire se serve ancora

# Get the log-likelihood values for each data point
log_likelihoods <- matrix(predict(multiv_model)$z, ncol = multiv_model$G)

# Calculate the anomaly scores based on log-likelihoods
anomaly_scores <- apply(log_likelihoods, 1, max)

# Set a threshold to classify anomalies
threshold <- quantile(anomaly_scores, 0.95)

# Identify anomalies based on the threshold
anomalies <- df_mm[anomaly_scores > threshold, ]

# Visualize the results
plot(df_mm, pch = 19, col = ifelse(anomaly_scores > threshold, "#FF7F50", "#99D9EA"),
     main = "Anomaly Detection using Gaussian Mixture Model")
#Nel contesto dell'Anomaly Detection con GMM, i punteggi di anomalia vengono calcolati per ciascuna osservazione. 
#Questi punteggi rappresentano quanto un'osservazione si discosta dal comportamento generale del dataset.

###############################################################################################################################

# RIDUZIONE DELLA DIMENSIONALITÀ

# rimuovere / spiegare la diversità nello scopo con l'LDA

# Riduzione della Dimensionalità con Modelli di Miscela Finita Gaussiana
#I modelli di miscela finita (GMM) possono gestire dati ad alta dimensionalità, ovvero con molte caratteristiche. Le tecniche di riduzione della dimensionalità, come la funzione MclustDR, possono essere utili per:

#Visualizzare i dati in uno spazio a dimensionalità inferiore (ad esempio, grafici 2D o 3D) per una più facile interpretazione. Questo permette di osservare come i dati si raggruppano e di identificare eventuali pattern che potrebbero essere nascosti in una dimensionalità elevata.
#Identificare le caratteristiche più importanti che contribuiscono alla struttura dei cluster. La riduzione di dimensionalità può aiutare a capire quali caratteristiche sono più significative per distinguere i diversi gruppi presenti nei dati. Questo può aiutare a concentrarsi sulle variabili chiave per l'analisi e l'interpretazione dei risultati.
# 
red_m_model <- MclustDR(multiv_model)
summary(red_m_model)

plot(red_m_model,what="pairs")
plot(red_m_model,what="boundaries", ngrid=200)
plot(red_m_model, what = "density", type = "persp")

set.seed(123)

red_m_model <- red_m_model$dir[,-3]

df_em <- mvnormalmixEM(red_m_model, k = 2)

plot(df_em, density = TRUE, cex.axis = 1.4, cex.lab = 1.5, cex.main = 1.5, main2 = "Clustering LDA", ylab2 = "LDA 1", xlab2 = "LDA 2")

class_mixt_em <- as.numeric(apply(df_em$posterior, MARGIN=1, FUN=which.max))
class_mixt_em
penguin_km <- kmeans(df_mm, centers=2)
class_km <- penguin_km$cluster
df_post <- cbind(df_em$posterior,class_mixt_em,class_km)
colnames(df_post) <- c("probG1", "probG2","class_mixt", "class_km")
head(df_post)
table(class_mixt_em,class_km)

###############################################################################################################################

# LDA

# RIDUZIONE DELLA DIMENSIONALITÀ - Riduzione della Dimensionalità con Modelli di Miscela Finita Gaussiana
#I modelli di miscela finita (GMM) possono gestire dati ad alta dimensionalità, ovvero con molte caratteristiche. Le tecniche di riduzione della dimensionalità, come la funzione MclustDR, la LDA, possono essere utili per:

#Visualizzare i dati in uno spazio a dimensionalità inferiore (ad esempio, grafici 2D o 3D) per una più facile interpretazione. Questo permette di osservare come i dati si raggruppano e di identificare eventuali pattern che potrebbero essere nascosti in una dimensionalità elevata.
#Identificare le caratteristiche più importanti che contribuiscono alla struttura dei cluster. La riduzione di dimensionalità può aiutare a capire quali caratteristiche sono più significative per distinguere i diversi gruppi presenti nei dati. Questo può aiutare a concentrarsi sulle variabili chiave per l'analisi e l'interpretazione dei risultati.

# creazione di un dataset con tutte le informazioni standardizzate, meno la colonna della specie che resta numerica

df_mm <- backup
df_mm$species <- as.numeric(as.factor(df_mm$species))
df_mm$island <- as.numeric(as.factor(df_mm$island))
# df_mm$sex <- as.numeric(as.factor(df_mm$sex))
df_mm <- scale(df_mm)
df_mm <- df_mm[,-1]
df_mm <- cbind(backup[1],df_mm)
df_mm$species <- as.numeric(as.factor(df_mm$species))

# as data.frame 

df_mm <- as.data.frame(df_mm)
glimpse(df_mm)

# applicazione LDA
df_lda <- lda(species ~ ., data=df_mm)
df_lda

# probabilità a posteriori per ogni osservazione del dataset - viene stampata dopo con la classe finale

df_lda.values <- predict(df_lda)

# Scatterplots of the Discriminant Functions

plot(df_lda.values$x[,1],df_lda.values$x[,2],  cex=0.8,pch=20, xlab="LDA axis 1", ylab="LDA axis 2") # make a scatterplot
text(df_lda.values$x[,1],df_lda.values$x[,2],df_lda.values$class,cex=0.7,pos=4,col="red") # add labels
plot(df_lda.values$x[,1],df_lda.values$x[,2],cex=0.8, col=df_lda.values$class, pch=20, xlab="LDA axis 1", ylab="LDA axis 2") # make a scatterplot

# Dallo scatterplot delle prime due discriminanti, i tre cluster che rappresentano le tre specie sono graficamente distinguibili
#far vedere le prime righe, per far vedere che i valori maggiori corrispondono alla classe 1 ecc

# mostra per ogni osservazione del dataset contemporaneamente la probabilità a posteriori di appartenenza a ciascuna delle tre classi e c'è l'ultima colonna col "verdetto" su quale sia la classe alla quale viene assegnato
lda_result <- cbind(as.data.frame(df_lda.values$posterior),as.data.frame(df_lda.values$class))
lda_result

# dettaglio esclusivamente sulla classe per ogni osservazione/unità
df_lda.values$class

# matrice di confusione
confMat.LDA<-confusionMatrix((as.factor(df_mm$species)), df_lda.values$class)$table
confMat.LDA

# errore percentuale di erroneità della classificazione
delta.LDA=(confMat.LDA[1,2]+confMat.LDA[2,1])/nrow(df_mm)*100 
delta.LDA

# accuratezza in percentuale della previsione effettuata dall'LDA
accuracy.LDA<-sum(diag(confMat.LDA))/(nrow(df_mm))*100
accuracy.LDA

#
# analysis on training and test set
#

set.seed(123)

train_penguins_rows <- createDataPartition(backup$species, p = 0.7, 
                                  list = FALSE, 
                                  times = 1)
train_penguins <-backup[train_penguins_rows,]
ntrain_rows<-nrow(train_penguins)
ntrain_rows
table(train_penguins$species)
table(train_penguins$species)/ntrain_rows*100
round(table(train_penguins$species)/ntrain_rows*100,2)

test_penguins<-backup[-train_penguins_rows,]
ntest_rows<-nrow(test_penguins)
ntest_rows
table(test_penguins$species)
table(test_penguins$species)/ntest_rows*100 
round(table(test_penguins$species)/ntest_rows*100,2)

#
# analysis on the training set
#

lda.fit_train=lda(species ~ ., data=train_penguins)  
lda.predict_train=predict(lda.fit_train,train_penguins)
pred.class_train=lda.predict_train$class
confMat.LDA_train=table(train_penguins$species,pred.class_train)
confMat.LDA_train
accuracy.LDA_train<-sum(diag(confMat.LDA_train))/ntrain_rows*100
accuracy.LDA_train

#
# then fit the model on the test set
#

lda.predict_test=predict(lda.fit_train,test_penguins)
pred.class_test=lda.predict_test$class
confMat.LDA_test=table(test_penguins$species,pred.class_test)
confMat.LDA_test
accuracy.LDA_test<-sum(diag(confMat.LDA_test))/ntest_rows*100
accuracy.LDA_test
accuracy.LDA_train
