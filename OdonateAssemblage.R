#########################################################################
###################        Working directory    #########################
#########################################################################

setwd("E:/McIntyre_Lab/Thesis/Data/OdonateCommunity")

##########################################################################
########################## Packages used #################################
##########################################################################
library(BiodiversityR)
library(labdsv)
library(MASS)
library(MVA)
library(vegan)
library(plyr)
library(betapart)
library(ggplot2)
library(ggrepel)
library(ggdendro)
library(grid)
library(corrplot)
library(lattice)
library(picante)
library(stats)
library(sars)
library(dplyr)
library(optpart)
library(cluster)
library(geosphere)
library(tidyverse)
library(caret)
library(mvtnorm)

#########################################################################
########################     Data used       ############################
#########################################################################

#Data that includes adults seen at each site
odes <- read.csv("odes-withadult.csv", row.names = 1, header = T)
View(odes) 
odes_matrix <- as.matrix(odes)

#Data for nymphs sampled at each site
odes.comm <- read.csv("nmds_odes.csv", row.names = 1, header = T)
View(odes.comm)

ode.matrix <- as.matrix(odes.comm)

ode.dist <- vegdist(odes, method = "bray", binary = T)

sites <- read.csv("oa_sites.csv", row.names = 1, header = T)

#Sites types without EAFB03 (no environmental data)
sites.nmds <- read.csv("oa_sites-nmds.csv", row.names = 1, header = T)
View(sites.nmds)

#Environmental data for each site
env <- read.csv("nmds_env.t.csv", row.names = 1, header = T)
env_scale <- scale(env, center = F, scale = T)
env_matrix <- as.matrix(env_scale)
View(env_matrix)

#Environmental data (without the stream comp) for correlation analyses
corr.env <- read.csv("corr_env.trans.csv", row.names = 1, header = T)
corr.env.s <- scale(corr.env, center = F, scale = T)
View(corr.env.s)

#########################################################################
#####################    Diversity indices     ##########################
#########################################################################

#Species richness per site
ode.rich <- specnumber(odes)
ode.rich

#Richness estimate
sp_est <- specpool(odes) 
sp_est

#Examine the richness across steephead and non-steephead sites
ode.sp <- specpool(odes, pool = sites$type)
ode.sp

#Richness and estimates for steephead sites
ode.sp["Steephead",]

#Richness and estimates for non-steephead sites
ode.sp["Non-Steephead",]

#If you want a particular richness metric for all stream types
ode.sp[,c("boot", "boot.se")]

#species accumulation curves for estimates of total species richness based on 
#iterative sampling with confidence intervals:
(ode.sp = poolaccum(odes))
plot(ode.sp)

#To plot just one of these metrics on their own:
ode.est <- poolaccum(odes)
ode.est
plot(ode.est, col="black", 
     strip=function(..., bg) strip.default(..., bg="gray90"), 
     display=c("boot", "S"))

#Determine whether there is a difference in richness by stream type:
oderich_aov <- t.test(ode.rich ~ type, data = sites)
summary(oderich_aov)
oderich_aov

#Which stream type was the most diverse
#Plot the average number of species at each steephead vs non-steephead site
boxplot(specnumber(odes) ~ type, data = sites, xlab= "Stream type", 
        ylab = "# of species")

#Calculate the shannon's diversity index for each site
diversity(odes[-1], index = "shannon")

sites.2 <- read.csv("oa_sites.csv", row.names = 1, stringsAsFactors = T)

#Shannon's diversity index for each stream type
diversitycomp(odes, y = sites.2, factor1 = "type", index ="Shannon",
              method = "pooled")

#Which stream type had individuals the most evenly distributed over space
diversitycomp(odes, y=sites.2, factor1="type", index="Jevenness", 
              method = "pooled")
gamma <- specpool(odes)
gamma

#Which stream type had the highest alpha diversity
alpha <- tapply(specnumber(odes), sites, FUN= mean)
alpha

beta <- (gamma$Species/alpha) - 1

#How similar were steephead sites to non-steephead sites (beta diversity)
ode <- rowsum(odes, group= sites$type)
betadiver(ode, method="w")



#Is beta diversity significantly different between steepheads and non-steepheads?
ode.b <- betadiver(odes, method = "w")
ode.ano <- anosim(ode.b, sites$type)
ode.ano


#Check if any sites are outliers with respect to ode richness
occupancy <- read.csv("nmds_odes.csv", header=TRUE)
boxplot(specnumber(occupancy), xlab = "Sites 1-23", ylab = "# of species")

ano <- anosim(ode.matrix, sites.nmds$type, distance = "bray", permutations = 9999)
ano

#species occurrence patterns
abuocc(odes)


#########################################################################
######################    Mantel Test  ##################################
#########################################################################

### Canopy ###
canopy <- env_matrix[,1]
View(canopy)
dist.ode <- vegdist(ode.matrix, method = "bray", binary = T)
dist.can <- dist(canopy, method = "euclidean")
ode.can <- mantel(dist.ode, dist.can, method = "spearman", permutations = 9999, na.rm = T)
ode.can
### p= 0.0561

### Temperature ###
temp <- env_matrix[,2]
View(temp)
dist.ode <- vegdist(ode.matrix, method = "bray", binary = T)
dist.temp <- dist(temp, method = "euclidean")
ode.temp <- mantel(dist.ode, dist.temp, method = "spearman", permutations = 9999, na.rm = T)
ode.temp

### Velocity ###
v <- env_matrix[,3]
View(v)
dist.ode <- vegdist(ode.matrix, method = "bray", binary = T)
dist.v <- dist(v, method = "euclidean")
ode.v <- mantel(dist.ode, dist.v, method = "spearman", permutations = 9999, na.rm = T)
ode.v

### pH ###
ph <- env_matrix[,4]
View(ph)
dist.ode <- vegdist(ode.matrix, method = "bray", binary = T)
dist.ph <- dist(ph, method = "euclidean")
ode.ph <- mantel(dist.ode, dist.ph, method = "spearman", permutations = 9999, na.rm = T)
ode.ph

### DO ###
do <- env_matrix[,5]
View(do)
dist.ode <- vegdist(ode.matrix, method = "bray", binary = T)
dist.do <- dist(do, method = "euclidean")
ode.do <- mantel(dist.ode, dist.do, method = "spearman", permutations = 9999, na.rm = T)
ode.do

### TDS ###
tds <- env_matrix[,6]
View(tds)
dist.ode <- vegdist(ode.matrix, method = "bray", binary = T)
dist.tds <- dist(tds, method = "euclidean")
ode.tds <- mantel(dist.ode, dist.tds, method = "spearman", permutations = 9999, na.rm = T)
ode.tds

### Width ###
wid <- env_matrix[,7]
View(wid)
dist.ode <- vegdist(ode.matrix, method = "bray", binary = T)
dist.wid <- dist(wid, method = "euclidean")
ode.wid <- mantel(dist.ode, dist.wid, method = "spearman", permutations = 9999, na.rm = T)
ode.wid
# p=0.0566

### Turbidity ###
turb <- env_matrix[,8]
View(turb)
dist.ode <- vegdist(ode.matrix, method = "bray", binary = T)
dist.turb <- dist(turb, method = "euclidean")
ode.turb <- mantel(dist.ode, dist.turb, method = "spearman", permutations = 9999, na.rm = T)
ode.turb

### Depth ###
dep <- env_matrix[,9]
View(dep)
dist.ode <- vegdist(ode.matrix, method = "bray", binary = T)
dist.dep <- dist(dep, method = "euclidean")
ode.dep <- mantel(dist.ode, dist.dep, method = "spearman", permutations = 9999, na.rm = T)
ode.dep

### Sand ###
sand <- env_matrix[,10]
View(sand)
dist.ode <- vegdist(ode.matrix, method = "bray", binary = T)
dist.sand <- dist(sand, method = "euclidean")
ode.sand <- mantel(dist.ode, dist.sand, method = "spearman", permutations = 9999, na.rm = T)
ode.sand

### Gravel ###
grv <- env_matrix[,11]
View(grv)
dist.ode <- vegdist(ode.matrix, method = "bray", binary = T)
dist.grv <- dist(grv, method = "euclidean")
ode.grv <- mantel(dist.ode, dist.grv, method = "spearman", permutations = 9999, na.rm = T)
ode.grv

### Cobble ###
cob <- env_matrix[,12]
View(cob)
dist.ode <- vegdist(ode.matrix, method = "bray", binary = T)
dist.cob <- dist(cob, method = "euclidean")
ode.cob <- mantel(dist.ode, dist.cob, method = "spearman", permutations = 9999, na.rm = T)
ode.cob

### Boulders ###
bo <- env_matrix[,13]
View(bo)
dist.ode <- vegdist(ode.matrix, method = "bray", binary = T)
dist.bo <- dist(bo, method = "euclidean")
ode.bo <- mantel(dist.ode, dist.bo, method = "spearman", permutations = 9999, na.rm = T)
ode.bo

### fpom ###
fpom <- env_matrix[,14]
View(fpom)
dist.ode <- vegdist(ode.matrix, method = "bray", binary = T)
dist.fpom <- dist(fpom, method = "euclidean")
ode.fpom <- mantel(dist.ode, dist.fpom, method = "spearman", permutations = 9999, na.rm = T)
ode.fpom

### cpom ###
cpom <- env_matrix[,15]
View(cpom)
dist.ode <- vegdist(ode.matrix, method = "bray", binary = T)
dist.cpom <- dist(cpom, method = "euclidean")
ode.cpom <- mantel(dist.ode, dist.cpom, method = "spearman", permutations = 9999, na.rm = T)
ode.cpom

### Bryophytes ###
bry <- env_matrix[,16]
View(bry)
dist.ode <- vegdist(ode.matrix, method = "bray", binary = T)
dist.bry <- dist(bry, method = "euclidean")
ode.bry <- mantel(dist.ode, dist.bry, method = "spearman", permutations = 9999, na.rm = T)
ode.bry

### Small wood ###
sw <- env_matrix[,17]
View(sw)
dist.ode <- vegdist(ode.matrix, method = "bray", binary = T)
dist.sw <- dist(sw, method = "euclidean")
ode.sw <- mantel(dist.ode, dist.sw, method = "spearman", permutations = 9999, na.rm = T)
ode.sw

### Large wood ###
lw <- env_matrix[,18]
View(lw)
dist.ode <- vegdist(ode.matrix, method = "bray", binary = T)
dist.lw <- dist(lw, method = "euclidean")
ode.lw <- mantel(dist.ode, dist.lw, method = "spearman", permutations = 9999, na.rm = T)
ode.lw

### Cumulative effect ###
dist.env <- dist(env_matrix, method = "euclidean")
ode.env <- mantel(dist.ode, dist.env, method = "spearman", permutations = 9999, na.rm = T)
ode.env

#########################################################################
########################         NMDS         ###########################
#########################################################################

set.seed(123)
meta.ode <- metaMDS(ode.matrix, distance = "bray", k= 3, 
                    binary = T, weakties = F, autotransform = F)
meta.ode
plot(meta.ode)
data.scores <- as.data.frame(scores(meta.ode)$sites)
data.scores$type <- sites.nmds$type

for(i in 2:5) print(metaMDS(ode.matrix, distance = "bray", k=i,
                            binary = T, weakties = F, autotransform = F, trace = F)$stress*100)
#3 dimensions is best

stressplot(meta.ode)

## Plot the site x species NMDS in ggplot ##
nmds.plot <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 4, aes(shape = type, color = type)) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", color = "Type", y = "NMDS2", shape = "Type")  + 
  scale_colour_manual(values = c("#009E73", "#E69F00"))

nmds.plot


## envfit to add environmetal data to the NMDS plot ##

#Correlation analyses to decide if any metrics should be removed
env_cor <- corr.env
env_cor <- as.matrix(corr.env)
env_cor <- na.omit(env_cor)
cor_mat <- cor(env_cor, method = 'pearson')
corrplot(cor_mat, method = "number", type = "lower", diag = F, tl.col = "black")
#remove spc

en <- envfit(meta.ode, env_matrix, permutations = 999, na.rm = T)
en

en_coord <- as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)

nmds.env <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = data.scores, size = 4, aes(shape = type, color = type)) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord, linewidth =1, alpha = 0.5, color = "grey30") +
  geom_text_repel(data = en_coord, aes(x = NMDS1, y = NMDS2), color = "grey30", 
                  fontface = "bold", label = row.names(en_coord)) +
  annotate(geom = "label", x = -1, y = 1.25, size = 5,
           label = paste("Stress: ", round(meta.ode$stress, digits = 3))) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", color = "Type", y = "NMDS2", shape = "Type")  + 
  scale_colour_manual(values = c("#009E73", "#E69F00"))
nmds.env

scores(en, "vectors")
#########################################################################
########################     Cluster Analysis     #######################
#########################################################################

#### Abitoic clusters ####
env_clust <- hclust(dist.env, method = "average")
ec <- as.dendrogram(env_clust)
plot(as.dendrogram(env_clust), ylab = "Euclidean dissimilarity")
rect.hclust(env_clust, k = 2, border = 3:4)

#r coefficient for how well the analysis fits my data
env.conph <- cophenetic(env_clust)
cor(dist.env, env.conph, method = "spearman")

env.clust.cut <- cutree(env_clust, k=2)
plot(env_clust, labels = as.character(env.clust.cut))

#Fusion levels
par(mfrow=c(2,2))

summary(env_clust)
plot(env_clust$height, nrow(env_matrix):2, type = "S", main = "Fusion levels - Chord - Average",
     ylab = "k (number of clusters)", xlab = "h (node height)", col = "grey")
text(env_clust$height, nrow(env_matrix):2, nrow(env_matrix):2, col = "red", cex = 0.8)

#Select optimum number of clusters using Silhouette widths
asw <- numeric(nrow(env_matrix))

for (k in 2:(nrow(env_matrix)-1)) {
  sil <- silhouette(cutree(env_clust, k=k), dist.env)
  asw[k] <- summary(sil)$avg.width
}

k.best <- which.max(asw)

plot(1:nrow(env_matrix), asw, type = "h", main = "Silhouette- optimal number of clusters, Average linkage",
     xlab = "k (number of groups)", ylab = "Average silhouette width")
axis(1, k.best, paste("optimum", k.best, sep = "\n"), col = "red", font = 2, 
     col.axis = "red")
points(k.best, max(asw), pch=16, col="red", cex=1.5)
cat("", "Silhouette-optimal number of clusters k=", k.best, "\n", "with an average
    silhouette width of", max(asw), "\n")


envpam <- pam(dist.env, k=2)
attributes(envpam)

envpam$medoids
envpam$clustering
plot(envpam)

clusplot(envpam, labels = 3, cex = 0.65)


#### Assemblage clusters ####
clust.dist <- vegdist(odes, method = "bray", binary = T)

#Hierarchical clustering, average linkage
ode_clust <- hclust(clust.dist, method = "average")
plot(as.dendrogram(ode_clust), ylab = "Bray-Curtis' dissimilarity")
rect.hclust(env_clust, k = 3)


comm.coph <- cophenetic(ode_clust)

cor(clust.dist, comm.coph) #pretty good fit

odes.clust.cut <- cutree(ode_clust, k=3)
plot(ode_clust, labels = as.character(odes.clust.cut)) 


asw.s <- numeric(nrow(odes_matrix))

for (k in 2:(nrow(odes_matrix)-1)) {
  sils <- silhouette(cutree(ode_clust, k=k), clust.dist)
  asw.s[k] <- summary(sils)$avg.width
}

k.best.s <- which.max(asw.s)

plot(1:nrow(odes_matrix), asw.s, type = "h", main = "Silhouette- optimal number of clusters, Average linkage",
     xlab = "k (number of groups)", ylab = "Average silhouette width")
axis(1, k.best.s, paste("optimum", k.best.s, sep = "\n"), col = "red", font = 2, 
     col.axis = "red")
points(k.best.s, max(asw.s), pch=16, col="red", cex=1.5)
cat("", "Silhouette-optimal number of clusters k=", k.best.s, "\n", "with an average
    silhouette width of", max(asw.s), "\n")


#Complete linkage
odes.clust.comp <- hclust(clust.dist, method = "complete")
plot(odes.clust.comp, ylab = "Morisitas-Horn dissimilarity")

comp.comm.coph <- cophenetic(odes.clust.comp)

cor(clust.dist, comp.comm.coph) #not a good fit

#Non-hierarchical
demopam <- pam(clust.dist, k=3)
attributes(demopam)

demopam$medoids
demopam$clustering
plot(demopam)

clusplot(demopam, labels = 3, cex = 0.65)

#Pres/abs
ode.beta <- betadiver(odes, method = "w")
oc <- hclust(ode.beta)
plot(oc)

ode.clus <- pam(ode.beta, k=2)
plot(ode.clus)
clusplot(ode.clus, labels= 3, cex=0.65)

#########################################################################
#######################     Dis. Func. Analysis     #####################
#########################################################################

#Bring in data
envdfa <- read.csv("env_dfa.csv", row.names = 1, header = T)
envdfa

#Split the data into training (80%) and test set (20%)
set.seed(1234)

train.ind <- envdfa$type %>%
  createDataPartition(p = 0.8, list = F)

train.data <- envdfa[train.ind,]
test.data <- envdfa[-train.ind,]

#Estimate preprocessing parameters
preproc.pm <- train.data %>%
  preProcess(method = c("center", "scale"))

preproc.pm

#Transform data using the estimated parameters
train.trans <- preproc.pm %>% predict(train.data)
test.trans <- preproc.pm %>% predict(test.data)

#Fit the model
model <- lda(type~., data = train.trans)
#Warning: variables are collinear

#Make predictions
predictions <- model %>% predict(test.trans)
predictions$x

#Model accuracy
mean(predictions$class==test.trans$type)

test.data$type <- as.factor(test.data$type)
model2 <- lda(type~., data = train.data)
predicted <- predict(model, newdata = test.data)
confusionMatrix(predicted$class, test.data$type)

model2
predicted

### Nancy's way ###
model.mass <- lda(type~., envdfa)
plot(model.mass)
par(mar=c(1,1,1,1))
model.mass
predict(model.mass)
train <- sort(sample(1:120,60))
table(type[train])
train

### Graph the output ###



var_covar <- matrix(data = c(1.5, 0.4, 0.4, 1.5), nrow = 2)

Xplus1 <- rmvnorm(400, mean = c(3,3), sigma = var_covar)

Xminus1 <- rmvnorm(600, mean = c(3, 3), sigma = var_covar) 

Y_samples <- c(rep(1, 400), rep(-1, 600))

dataset <- as.data.frame(cbind(rbind(Xplus1, Xminus1), Y_samples))
colnames(dataset) <- c("X1", "X2", "Y")
dataset$Y <- as.character(dataset$Y)

ggplot(data = dataset) + geom_point(aes(X1, X2, color = Y))


set.seed(1)
sample <- sample(c(T, F), nrow(envdfa), replace = T, prob = c(0.7, 0.3))
train <- envdfa[sample, ]
test <- envdfa[!sample, ]
model3 <- lda(type~., data = train)
model3
predictor <- predict(model3, test)
names(predictor)
mean(predictor$class==test$type)
predictor
lda_plot <- cbind(train, predict(model3)$x)
colnames
ggplot(lda_plot) + geom_point(col = type)
#####################################################################
######        Stacked Bar Plot with Colors and Legend       #########
#####################################################################

sites.graph <- c("RES01", "RES01", "RES03", "RES03", "RES04", "RES04","RES05",
                 "RES05", "WB01", "WB01", "TRSF01", "TRSF01", "BB01", "BB01", 
                 "CCR01", "CCR01","LCWMA01", "LCWMA01", "LCWMA02","LCWMA02", 
                 "LCWMA03","LCWMA03", "JB01", "JB01","PB01","PB01", "EAFB01", 
                 "EAFB01", "EAFB02", "EAFB02", "EAFB03","EAFB03","BRSF01", 
                 "BRSF01", "BC01", "BC01", "LTSF01","LTSF01", "LTSF02", "LTSF02",
                 "ARBP01", "ARBP01","TSP01", "TSP01", "TSP02", "TSP02","CB01", "CB01")
# repeat each site twice, once for nymphal observations and another for adult
sites.graph

#Make a data frame for species richness at each site
ode.comm <- data.frame(Site = sites.graph, sp.rich = c(3,0,3,0,5,3,2,0,8,0,5,1,2,
                                                       0,2,0,6,3,1,4,10,4,5,4,1,
                                                       4,8,5,4,0,0,3,7,0,6,1,12,
                                                       0,3,2,11,11,5,3,1,1,9,3))
#Add a column to the data frame for life stage for each value of richness at a site
ode.comm$stage <- c("Nymph", "Adult", "Nymph", "Adult", "Nymph", "Adult", 
                    "Nymph", "Adult", "Nymph", "Adult", "Nymph", "Adult", 
                    "Nymph", "Adult", "Nymph", "Adult", "Nymph", "Adult",  
                    "Nymph", "Adult", "Nymph", "Adult", "Nymph", "Adult", 
                    "Nymph", "Adult", "Nymph", "Adult", "Nymph", "Adult",
                    "Nymph", "Adult", "Nymph", "Adult","Nymph", "Adult",
                    "Nymph", "Adult","Nymph", "Adult","Nymph", "Adult",
                    "Nymph", "Adult","Nymph", "Adult","Nymph", "Adult")
# alternate adult nymph so that each site has an adult and a nymph column
ode.comm


ode.comm$Site <- factor(ode.comm$Site, levels = c("RES01", "RES03", "RES04",
                                                  "RES05","WB01", "TRSF01", 
                                                  "BB01", "CCR01","LCWMA01", 
                                                  "LCWMA02","LCWMA03","JB01",
                                                  "PB01","EAFB01", "EAFB02", 
                                                  "EAFB03","BRSF01","BC01", 
                                                  "LTSF01", "LTSF02","ARBP01",
                                                  "TSP01","TSP02","CB01"), ordered = T)

# Sort by site and life stage (adult or nymph)
odecomm_sort <- arrange(ode.comm, Site, stage)
head(odecomm_sort)

# Calculate the cumulative sum of sp.rich for each site
ode_cumsum <- ddply(odecomm_sort, "Site", transform, tot.rich=cumsum(sp.rich))

# change the data source in ggplot to this new data frame

#Stacked bar graph of species richness per site broken down to number of species
#documented as nymphs and number of species documented as adults
ode.plot <- ggplot(data = ode_cumsum, aes(x=Site, y=sp.rich, label=sp.rich, fill=stage)) + 
  geom_bar(stat = "identity") + 
  geom_text(data=subset(ode.comm, sp.rich != 0), aes(y=sp.rich), 
            position = position_stack(vjust = 0.5), color="white", size=3.5) +
  scale_fill_brewer(palette = "Paired") +
  labs(x = "Site", y = "Species richness", fill = "Life stage") +
  theme_classic() +
  theme(panel.grid= element_blank(), 
        axis.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16),
        axis.text = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10,
                                 color = "black")) +
  scale_y_continuous(name = "Species richness", limits = c(0,25))
ode.plot



