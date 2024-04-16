##########################
#     MIXTURE MODELS     #
##########################

# df_mm <- backup[,-7]
df_mm<- df[,-5]
#df_mm$species <- as.numeric(as.factor(df_mm$species))
#df_mm$island <- as.numeric(as.factor(df_mm$island))
#df_mm <- scale(df_mm)
#df_mm <- df_mm[,-1]

#df_mm <- cbind(df_mm, backup[1])
#df_mm$species <- as.numeric(as.factor(df_mm$species))
#colnames(df_mm)[6] <- "Classes"

multiv_model <- Mclust(df_mm)
summary(multiv_model, parameters=TRUE)
plot(multiv_model, what = "BIC")

plot(multiv_model, what = "density", type = "hdr", data = df_mm, points.cex = 0.5) 
plot(multiv_model, what="classification", points.cex = 0.5)
plot (multiv_model, what = "uncertainty", points.cex =26)

###############################################################################################################################

# CLAUDIA SARA GIORGIO VALUTARE SE RIMUOVERE

# RIDUZIONE DELLA DIMENSIONALITÀ - Riduzione della Dimensionalità con Modelli di Miscela Finita Gaussiana
#I modelli di miscela finita (GMM) possono gestire dati ad alta dimensionalità, ovvero con molte caratteristiche. Le tecniche di riduzione della dimensionalità, come la funzione MclustDR, possono essere utili per:

#Visualizzare i dati in uno spazio a dimensionalità inferiore (ad esempio, grafici 2D o 3D) per una più facile interpretazione. Questo permette di osservare come i dati si raggruppano e di identificare eventuali pattern che potrebbero essere nascosti in una dimensionalità elevata.
#Identificare le caratteristiche più importanti che contribuiscono alla struttura dei cluster. La riduzione di dimensionalità può aiutare a capire quali caratteristiche sono più significative per distinguere i diversi gruppi presenti nei dati. Questo può aiutare a concentrarsi sulle variabili chiave per l'analisi e l'interpretazione dei risultati.

red_m_model <- MclustDR(multiv_model)
summary(red_m_model)

plot(red_m_model,what="pairs")
plot(red_m_model,what="boundaries", ngrid=200)
plot(red_m_model, what = "density", type = "persp")

###############################################################################################################################

# CLAUDIA SARA GIORGIO VALUTARE SE RIMUOVERE

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

###############################################################################################################################

# LDA

df_mm <- backup[,-7]
df_mm$species <- as.numeric(as.factor(df_mm$species))
df_mm$island <- as.numeric(as.factor(df_mm$island))
df_mm <- scale(df_mm)
df_mm <- df_mm[,-1]

df_mm <- cbind(df_mm, backup[1])
df_mm$species <- as.numeric(as.factor(df_mm$species))
#colnames(df_mm)[6] <- "Classes"

df.lda <- lda(species ~ ., data=df_mm)
df.lda

# This means that the first discriminant function is a linear combination of the 
# variables: −0.403∗Alcohol+0.165∗Malic+...+−0.003∗Proline. 
# for convenience, the value for each discriminant function 
# (eg. the first discriminant function) are scaled so that their mean value 
# is zero and its variance is one.
#
# The “proportion of trace” that is printed when you type “wine.lda” 
# (the variable returned by the lda() function) is the percentage 
# separation achieved by each discriminant function. 
# For example, for the wine data we get the same 
# values as just calculated (68.75% and 31.25%).
#
df.lda.values <- predict(df.lda)
df.lda.values

# Scatterplots of the Discriminant Functions
# We can obtain a scatterplot of the best two discriminant functions, with the data 
# points labelled by cultivar
plot(df.lda.values$x)
     
##########################################################
data <- as.data.frame(df_mm)
source("C:/Users/Giodra/Downloads/EM2_plot.r")
source("C:/Users/Giodra/Downloads/mixt2_univ.r")
em2_weight <- EM2_plot(data$body_mass_g, "Penguin Weight distribution")     
em2_culmen_length <- EM2_plot(data$culmen_length_mm, "Culmen Length distribution")  
em2_culmen_depth <- EM2_plot(data$culmen_depth_mm, "Culmen Depth distribution")  
em2_flipper_length <- EM2_plot(data$flipper_length, "Flipper Length distribution") 













