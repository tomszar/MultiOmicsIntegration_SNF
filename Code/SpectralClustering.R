#### FOLLOWING FUNCTIONS WERE MODIFIED FROM THE SNFtool library
spectralClust <- function(affinity, K, type=3) {
    # Implements spectral clustering on given affinity matrix into K clusters.
    #
    # Args:
    #   affinity: Affinity matrix (size NxN) to perform clustering on
    #   K: Number of clusters 
    #   type (default 3): Used to speciy the type of spectral clustering
    #
    # Returns:
    #   labels: A vector of length N assigning a label 1:K to each sample
 
    d <- rowSums(affinity)
    d[d == 0] <- .Machine$double.eps
    D <- diag(d)
    L <- D - affinity
 
    if (type == 1) {
        NL <- L
 
    } else if (type == 2) {
        Di <- diag(1 / d)
        NL <- Di %*% L
 
    } else if(type == 3) {
        Di <- diag(1 / sqrt(d))
        NL <- Di %*% L %*% Di
    }
 
    eig <- eigen(NL)
    res <- sort(abs(eig$values),index.return = TRUE)
    U <- eig$vectors[,res$ix[1:K]]
    normalize <- function(x) x / sqrt(sum(x^2))
 
    if (type == 3) {
        U <- t(apply(U,1,normalize))
    }
 
    eigDiscrete <- .discretisation(U)
    eigDiscrete <- eigDiscrete$discrete
    labels <- apply(eigDiscrete,1,which.max)
    return(list(labels = labels,
                eigen_values = eig$values,
                eigen_vectors = eig$vectors)) }

.discretisation <- function(eigenVectors) {
  
  normalize <- function(x) x / sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors,1,normalize))
  
  n = nrow(eigenVectors)
  k = ncol(eigenVectors)
  
  R = matrix(0,k,k)
  R[,1] = t(eigenVectors[round(n/2),])
  
  mini <- function(x) {
    i = which(x == min(x))
    return(i[1])
  }
  
  c = matrix(0,n,1)
  for (j in 2:k) {
    c = c + abs(eigenVectors %*% matrix(R[,j-1],k,1))
    i = mini(c)
    R[,j] = t(eigenVectors[i,])
  }
  
  lastObjectiveValue = 0
  for (i in 1:20) {
    eigenDiscrete = .discretisationEigenVectorData(eigenVectors %*% R)
    
    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    U = svde[['u']]
    V = svde[['v']]
    S = svde[['d']]
    
    NcutValue = 2 * (n-sum(S))
    if(abs(NcutValue - lastObjectiveValue) < .Machine$double.eps) 
      break
    
    lastObjectiveValue = NcutValue
    R = V %*% t(U)
    
  }
  
  return(list(discrete=eigenDiscrete,continuous =eigenVectors))
}
 
.discretisationEigenVectorData <- function(eigenVector) {
  
  Y = matrix(0,nrow(eigenVector),ncol(eigenVector))
  maxi <- function(x) {
    i = which(x == max(x))
    return(i[1])
  }
  j = apply(eigenVector,1,maxi)
  Y[cbind(1:nrow(eigenVector),j)] = 1
  
  return(Y)
  
}

plot_tsne <- function(W, labels){
    library(Rtsne)
    dimred <- Rtsne(as.dist(0.5 - W), is_distance = TRUE)$Y
    colnames(dimred) <- paste("tSNE", seq_len(ncol(dimred)), sep = "")
    
    g <- ggplot(as.data.frame(dimred), aes(tSNE1, tSNE2)) + 
            geom_point(aes(colour=labels)) +
            theme_bw()
    if(is.numeric(labels)){
        g_color <- scale_color_viridis_c()
        g <- g + g_color
    }
    print(g)
}