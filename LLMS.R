LLMS <- function(X                                    ,
                 algo        = "R-kNN"                ,
                 k           = 0.2                    ,
                 theta       = seq(0,2,length.out=100),
                 m           = 5                      ,
                 q           = 0.2                    ,
                 batch.size  = 0.2                    ,
                 graph.model = FALSE){

    

  
# ---------------------------------------------------------------------------------------------------------------#
#                                                                                                                #
#                                 Local Linear Manifold Selection (LLMS)                                         #
#                                                                                                                #
#                                                                                                                #
#   This function implements the three versions of the Local Linear Manifold Selection estimator, as well as     #
#   the graphical representation. In the first case, the output is a vector of selection probabilities, of       #
#   the same length as the number of variables. In the second case, the output is a matrix of selection          #
#   probabilities, that represents the edges of the graph.                                                       #
#                                                                                                                # 
# ---------------------------------------------------------------------------------------------------------------#

  
  
  
  
## -------------------------------------
## Check the input of the functions and
## throw error if needed
X <- data.matrix(X)
index_not_numeric = which(is.na(as.numeric(X)))
if(length(index_not_numeric)>0){
  stop(paste("Input data X is not numeric at index position(s): ", index_not_numeric))
}
n = nrow(X)
p = ncol(X)
if(any(c(k,m,q,batch.size)<0)){
  stop("ERROR: No parameter can be negative.")
}
if(k>=1){
  stop("ERROR: The parameter 'k' cannot be larger than or equal to 1.")
}
if(m>p){
  stop("ERROR: The parameter 'm' cannot be larger than the number of columns of X.")
}
if(q>=1){
  stop("ERROR: The parameter 'q' cannot be larger than or equal to 1.")
}
if(any(theta<0)){
  stop("ERROR: The penalty parameter 'theta' cannot be negative.")
}
if(!identical(algo,"vanilla") & !identical(algo,"R-kNN") & !identical(algo,"q-kNN")){
  message("ERROR: The parameter 'algo' is mispecified, using default R-kNN instead.")
  algo = "R-kNN"
}
if(!is.logical(graph.model)){
  message("ERROR: The parameter 'graph.model' is mispecified, using default FALSE instead.")
  graph.model = FALSE
}

## -------------------------------------
## Change k into its integer counterpart 
k = floor(k*n)

## ----------------------------
## Length of the penalty vector
nb_theta = length(theta)
  
## ----------------------
## EigenThresholding (ET)
ET <- function(R,theta){
  eig <- eigen(R,symmetric=TRUE)
  V <- eig$vectors
  L <- eig$values
  tol = 1e-10
  L[abs(L)<=tol]= abs(L[abs(L)<=tol])
  L[L<.Machine$double.eps]=0
  return(list("V"=V,"L"=L))
}

## ---------------
## q-kNN algorithm
##  First it computes all pairwise distances
##  Then, each k-neighbourhood is selected
##  as the q-neighbourhoods of the unvisited 
##  q neighbours of a given point i
qkNN <- function(X,i,k,q){
  n.X <- dim(X)[1]
  p   <- dim(X)[2]
  '%ni%' <- Negate('%in%')
  J1 = matrix(1,p,1)
  J2 = matrix(1,1,n.X)
  distances <- sqrt(((X^2)%*%J1)%*%J2 + t(((X^2)%*%J1)%*%J2) - 2*(X%*%t(X)))
  Neighbourhood <- apply(distances,2,order)
  X_nn <- array(0,c(k,p,length(i)))
  for(j in 1:length(i)){
    nn = i[j]
    idx <- NULL
    while(length(idx) < k ){
      NN1 <- NULL
      for(j2 in 1:length(nn)){
        NN0 <- Neighbourhood[,nn[j2]]
        unseen <- which(NN0%ni%idx)[1:q]
        if(!(unseen[1]>q)){
          NN0 <- NN0[unseen]
          NN1 <- c(NN1,NN0)
        }
      }
      if(length(NN1)==0){
        message("WARNING : ALL POINTS VISITED")
        for(j2 in 1:length(nn)){
          NN0 <- Neighbourhood[,nn[j2]]
          unseen <- which(NN0%ni%idx)[1]
            NN0 <- NN0[unseen]
            NN1 <- c(NN1,NN0)
        }
      }
      nn <- unique(as.vector(NN1))
      if(length(idx)+length(nn)>k){
        nn <- nn[1:(k-length(idx))]
      }
      idx = c(idx,nn)
    }
    X_nn[1:length(idx),,j] <- X[idx,]
  }
return(X_nn)
}
    
## --------------------------------------------------------
## Indices of the random batch on which the LLMS is applied
N <- sample(n,floor(batch.size*n),replace=TRUE)

## --------------
## Pre-allocation
P <- matrix(0,p,nb_theta)
if(graph.model){
  G <- array(0,c(p,p,nb_theta))
}else{
  G <- NULL
}
count <- matrix(0,nb_theta,1)
    
## -----------------------------------
## Scaled data for computing distances
xscaled <- scale(X)
  
## ----------------------  
## Initial Subset Weights
if(algo == "R-kNN"){
  w <- rep(1/p,p)
}

## --------------------------------------
## If the algorithm is chosen to be q-kNN
## All the neighbourhoods are constructed 
## once, to speedup compuations.
if(algo == "q-kNN"){
  X_qknn <- qkNN(xscaled,N,k,floor(q*k))
}

## -------------------
## Progess bar display
pb_i <- txtProgressBar(min = 0, max = floor(batch.size*n), style = 3,width=100)

## ----------------------
## Loop over Observations
for(i in 1:floor(batch.size*n)){
  
  ## -----------------------
  ## Indent the progress bar
  setTxtProgressBar(pb_i, i)
  
  ## ----------------------------------------------------------
  ## Choice of the algorithm for constructing the neighbourhood
  
  if(algo == "vanilla"){
    ## --------------------------------------
    ##                kNN
    ## It uses all variables to compute the 
    ## distances.
    idx <- order(sqrt(colSums((t(xscaled) - xscaled[N[i],])^2)))[2:k]
    Xlocal <- X[idx,]
    
  }else if(algo == "R-kNN"){
    ## --------------------------------------
    ##                R-kNN
    ## It choses a random subset of variables
    ## to compute the distances.
    subset_var <- sample(p,m,prob=w,replace=TRUE)
    idx <- order(sqrt(colSums((t(xscaled[,subset_var]) - xscaled[N[i],subset_var])^2)))[2:k]
    Xlocal <- X[idx,]
    
  }else if(algo == "q-kNN"){
    ## --------------------------------------
    ##                q-kNN
    ## Already computed neighbourhoods
    ## See lines 108-114 of this file
    Xlocal <- X_qknn[,,i]
  }
  
  ## ----------------------------------------
  ## Correlation matrix of the neighbourhood.
  ## To cope with small non-linearities
  ## it uses the spearman correlation.
  ## It also handles the missing values
  ## by using only complete observations.
  R <- cor(Xlocal,method="spearman",use="complete.obs")
  
  ## ---------------------------------------------
  ## Apply the ET as the linear selection operator
  R.tilde <- ET(R)
  V <- R.tilde$V
  L <- R.tilde$L
  VL <- abs(V%*%diag(L))
  # Loop over the penalty parameter
  for(index_theta in 1:nb_theta){
    VL[VL < theta[index_theta]] = 0
    S = VL%*%t(VL)
    S <- S - diag(diag(S))
    if(graph.model){
      G[,,index_theta] = G[,,index_theta] + 1*(S!=0)
    }
    S = colSums(S)
    P[,index_theta] = P[,index_theta] + 1*(S!=0)
    # Skip the neighbourhoods with empty selection
    if(any(S!=0)){
      count[index_theta] = count[index_theta] + 1
    }
  }
  
  ## -------------------------
  ## Update the subset weights
  if(algo == "R-kNN"){
    w = w + (1/floor(batch.size*n))*apply(P,1,mean)
  }
    
}

## ----------------------
## Close the progress bar
close(pb_i)

## --------------------------------------------
## Divide per number of non-zero neighbourhoods
count[count==0]=1
for(index_theta in 1:nb_theta){
  P[,index_theta] = P[,index_theta]/count[index_theta] 
  if(graph.model){
    G[,,index_theta] = G[,,index_theta]/count[index_theta] 
  }
}
    
## ---------------------------
## Compute the selection paths
Paths <- t(P)

## ---------------    
## Optimal Penalty
## The value of theta that 
## maximises the variance 
## of the paths
theta_star <- which.max(apply(Paths,1,var))
P <- P[,theta_star]


## ----------------------
## Output of the function
return(
    list(
      "P"         = P,
      "Paths"     = Paths,
      "Graph"     = G,
      "theta_opt" = theta_star
        )
      )

## ------------------------  
## End of the LLMS function  
}