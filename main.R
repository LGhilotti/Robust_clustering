######################################################################################
### code implementing the algorithm in "A reweighting approach to robust clustering" #
######################################################################################
library(tclust)
library(plyr)
library(dplyr)
library(readr)
library("readxl")

set.seed(1234)

# data
#data = iris[,1:3]
#data=as.matrix(data)
#nonstand_data = read.csv("wine-clustering.csv") - 3 classes

#nonstand_data = read_excel("Date_Fruit_Datasets.xlsx") - 7 classes
#class = nonstand_data["Class"]
#data = nonstand_data[1:34]

############### Pumpkins ##################################
#nonstand_data <- read_excel("Pumpkin_Seeds_Dataset.xlsx")
#class = nonstand_data["Class"]
#to_remove <- c("Class", "Convex_Area")
#nonstand_data_noclass = nonstand_data[, !(names(nonstand_data) %in% to_remove)]

#data = scale(nonstand_data_noclass)
#########################################################

#nonstand_data <- read_excel("Raisin_Dataset.xlsx")
#class = nonstand_data["Class"]
#to_remove <- c("ConvexArea", "Perimeter", "Class")
#nonstand_data_noclass = nonstand_data[, !(names(nonstand_data) %in% to_remove)]

#data = (nonstand_data_noclass)

nonstand_data_na <- read_csv("water_potability.csv")
nonstand_data <- nonstand_data_na %>% drop_na()
class = nonstand_data["Potability"]
to_remove <- c("Potability")
nonstand_data_noclass = nonstand_data[, !(names(nonstand_data) %in% to_remove)]

data = (nonstand_data_noclass)

# PCA
#pca <- prcomp(data, center = TRUE,scale. = TRUE)

#num_pc <- which(cumsum(pca$sdev^2)/sum(pca$sdev^2)>0.95)[1]
#data <- pca$x[,1:num_pc]

n = nrow(data)
p = ncol(data)

# number of clusters, known in advance
k = 2

# initial trimming level
alpha0 = 0.3

# number of iterations L
L = 50

# final trimming level
alphaL = 0.05
eps = (alpha0 - alphaL) /L

# 1) initialization
init_est <- tclust(data, k = k, alpha = alpha0, nstart= 50, iter.max = 50, restr = "eigen", restr.fact = 12, 
        equal.weights = FALSE)

mus <- init_est$centers
sigmas <- alply(init_est$cov,3)
proportions <- c(init_est$size/(n*(1-alpha0)) , 0) # k+1 elements, the last one indicates the contamination

# 2) Reweighting process
for (l in 1:L){
  
  alphal = alpha0 - l*eps
  
  # 2.1) Update proportions
  # for each observation, save the minimum Mahalanobis distance from the k centers and the index of that center
  
  # ma_d contains the Mahalanobis distance for each datum from each center (n x k matrix)
  ma_dist <- matrix(nrow=n,ncol=k)
  for (j in 1:k){
    ma_dist[,j] <- mahalanobis(data,mus[,j],sigmas[[j]])
  }
  # D is a n-vector containing the minimum Mahalanobis for each datum
  D <- apply(ma_dist, 1, FUN = min)
  # z is a n-vector containing the index of the closest cluster
  z <- apply(ma_dist, 1, FUN = which.min)
  
  D_quant <- quantile(D, probs = 1-alphal ) # even if not a datum, it is close
  
  # set A
  A <- which(D <= D_quant)
  # set B
  B <- which(D <= qchisq(1-alphaL, df=p))
  # intersection
  AB <- intersect(A,B)
  
  # H is a list of k elements, each one corresponding to H_i
  H <- lapply((1:k), function(j) intersect(AB, which(z==j)))
  n0 <- sum(lengths(H))
  
  # update the contamination level
  proportions[k+1] <- 1 - length(B)/n
  # update the cluster weights
  proportions[1:k] <- lengths(H)/n0*(1-proportions[k+1])
  
  # 2.2) Update locations and scatters
  mus <- sapply(H, function(h) colMeans(data[h,]))
  sigmas <- lapply(H, function(h) cov(data[h,]))
  
  beta <- n0/(n*(1-proportions[k+1]))
  if (beta < 1)
    sigmas <- lapply(sigmas,"*", beta/pchisq(qchisq(1-beta,df=p), df=p+2))
   
}

# 3) Output of the algorithm
ma_dist <- matrix(nrow=n,ncol=k)
for (j in 1:k){
  ma_dist[,j] <- mahalanobis(data,mus[,j],sigmas[[j]])
}
# D is a n-vector containing the minimum Mahalanobis for each datum
D <- apply(ma_dist, 1, FUN = min)
# z is a n-vector containing the index of the closest cluster
z <- apply(ma_dist, 1, FUN = which.min)

# set B
B <- which(D <= qchisq(1-alphaL, df=p))

# H is a list of k elements, each one corresponding to H_i
H <- lapply((1:k), function(j) intersect(B, which(z==j)))

# Trimmed observations (proportion alphaL)
trimmed <- setdiff(1:n, B)

########### Water
# report the not standardized data summaries
means <- sapply(H, function(h) colMeans(data[h,]))
vec_class <- pull(class,Potability)
table(z[B], vec_class[B])
########################## Pumpkins
# report the not standardized data summaries
means <- sapply(H, function(h) colMeans(data[h,]))
vec_class <- pull(class,Class)
table(z[B], vec_class[B])

# for each outlier, compute the distance from the means
dist <- data[trimmed,]-

# for pumpkins
cols <- character(nrow(data))
cols[] <- "black"

cols[class == "Çerçevelik"] <- "blue"
cols[ class == "Ürgüp Sivrisi"] <- "red"
cols[trimmed] <- "black"
dev.new()
pairs(data,col=cols)
dev.off()

color = c(rep("black",length(trimmed)), rep("red",k))
points = rbind(data[trimmed,],t(means))
dev.new()
pairs(points, col=color)

km <- kmeans(data, centers = t(init_est$centers))
km_cen <- km$centers

table(km$cluster, vec_class)
table(z[B], vec_class[B])
table(init_est$cluster[which(init_est$cluster > 0)], vec_class[which(init_est$cluster > 0)])

color = c(rep("black",k), rep("red",k))
points = rbind(km_cen,t(means))
dev.new()
pairs(points, col=color)