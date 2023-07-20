# Original code for the gamma-metric https://github.com/Pierre-Michel/gammaMetric

# Ensemble des fonctions permettant le calcul de la gamma-metric
if(!require(corpcor)){
  install.packages("corpcor")
  library(corpcor)
}

## Fonction calculant la distance entre k1 et k2 ------------------------------------------------------
d_k1k2 <- function(mu1, mu2, eigen.values1, eigen.values2, eigen.vectors1, eigen.vectors2){
  
  # Vecteurs mu1mu2
  mu1mu2 <- mu2 - mu1
  
  # Coordonnées du vecteur mu1mu2 dans la base de l'ellipse
  mu1_coord <- solve(eigen.vectors1) %*% as.matrix(normalize(mu1mu2))
  mu2_coord <- solve(eigen.vectors2) %*% as.matrix(normalize(mu1mu2))
  
  # Calcul de alpha (facteur de normalisation)
  alpha_k1k2 <- sqrt(sum(eigen.values1)) + sqrt(sum(eigen.values2))
  
  # Calcul distances du centre au bord de l'ellipse (dans la direction mu1mu2)
  dk1 <- 1/sqrt(sum((mu1_coord^2)/(abs(eigen.values1)), na.rm = TRUE))
  dk2 <- 1/sqrt(sum((mu2_coord^2)/(abs(eigen.values2)), na.rm = TRUE))
  
  # Calcul des distances du centre au bord de l'ellipse (dans la direction mu1mu2)
  d_k1k2 <- (1/alpha_k1k2)*(norm(mu1mu2) - (dk1 + dk2))
  
  return(d_k1k2)
}

## Fonction utilitaire : renvoie le vecteur divisé par la somme des coefficients ----------------------
divsum <- function(x){
  return(x/sum(x))
}

## Fonction qui renvoie une ellipse (crédits : J.Ravaglia) --------------------------------------------
ellipse <- function(n, l1, l2, axe, center){
  angles <- seq(0, 2*pi, length.out = n)
  w <- atan2(axe[2], axe[1])
  X <- l1*cos(angles)
  Y <- l2*sin(angles)
  x <- center[1] + X*cos(w) - Y*sin(w)
  y <- center[2] + X*sin(w) - Y*cos(w)
  return(cbind(x,y))
}

## Fonction qui calcul la gamma-metric ---------------------------------------------------------------
gammaMetric <- function(X, class, covEstim = c("empiric","shrink",'mixte'), plot = F){
  X <- as.matrix(X)
  
  # Nombre de lignes
  n <- nrow(X)
  
  # Liste des k différentes classes
  muL <- list()
  eigen.valuesL <- list()
  eigen.vectorsL <- list()
  k <- length(unique(class))
  length(muL) = length(eigen.valuesL) = length(eigen.vectorsL) = k
  for(i in 1:k){
    XX <- X[class == sort(unique(class))[i],]
    XX <- apply(as.matrix(XX), 2, as.numeric)
    
    # Barycentre
    mu <- apply(as.matrix(XX), 2, mean, na.rm = TRUE)
    
    # Matrice de variance-covariance
    if(covEstim == "empiric"){
      WkM <- cov(XX, use = "complete.obs")
    }else if(covEstim == "shrink"){
      WkM <- corpcor::cov.shrink(XX, verbose = FALSE)
    }else{
      if(nrow(XX)<ncol(XX)){
        WkM <- corpcor::cov.shrink(XX, verbose = FALSE)
      }else{
        WkM <- cov(XX, use = 'complete.obs')
      }
    }
    
    # Les valeurs propres
    eigen.values <- eigen(WkM)$values
    
    # Les vecteurs propres
    eigen.vectors <- eigen(WkM)$vectors
    
    muL[[i]] <- mu
    eigen.valuesL[[i]] <- eigen.values
    eigen.vectorsL[[i]] <- eigen.vectors
    
  }
  
  # Calcul de la gamma-metric
  gamma <- 0
  for(k1 in 1:(k-1)){
    for(k2 in (k1+1):k){
      gamma <- gamma + d_k1k2(muL[[k1]], muL[[k2]], eigen.valuesL[[k1]], eigen.valuesL[[k2]], eigen.vectorsL[[k1]], eigen.vectorsL[[k2]])
    }
  }
  
  # On affiche les deux nuages de points
  if(plot == TRUE){
    plot(X[,1], X[,2], xlab = "x1", ylab = "x2", col = class, main = paste("gamma = ", round(gamma,2)))
    for(i in 1:k){
      lines(x = c(muL[[i]][1], muL[[i]][1] + (eigen.vectorsL[[i]][1,1]*eigen.valuesL[[i]][1])),
            y = c(muL[[i]][2], muL[[i]][2] + (eigen.vectorsL[[i]][2,1]*eigen.valuesL[[i]][1])), col = "blue")
      lines(x = c(muL[[i]][1], muL[[i]][1] + (eigen.vectorsL[[i]][1,2]*eigen.valuesL[[i]][2])),
            y = c(muL[[i]][2], muL[[i]][2] + (eigen.vectorsL[[i]][2,2]*eigen.valuesL[[i]][2])), col = "blue")
      ellipse <- ellipse(n = 200, l1 = eigen.values[[i]][1], l2 = eigen.valuesL[[i]][2], axe = eigen.vectorsL[[i]][,1], center = muL[[i]])
      points(ellipse, t = "l", col = "green")
      
      for(k2 in (i+1):k){
        if(k2 > k)break
        lines(x = c(muL[[i]][1], muL[[k2]][1]), y = c(muL[[i]][2], muL[[k2]][2]))
      }
    }
  }
  return(gamma)
}

## Fonction norme d'un vecteur ------------------------------------------------------------------------
norm <- function(x){
  return(sqrt(sum(x^2)))
}

## Fonction qui normalise un vecteur ------------------------------------------------------------------
normalize <- function(x){
  return(x/norm(x))
}
