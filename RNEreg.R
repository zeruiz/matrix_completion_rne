library(tidyverse)

BiDistance <- function(dat){
  dat[,1] <- as.numeric(factor(dat[,1]))
  dat[,2] <- as.numeric(factor(dat[,2]))
  U <- unique(dat[,1])
  I <- unique(dat[,2])
  colnames(dat) <- c("U", "I", "R")
  Ud <- data.frame("U"=U, "0<d<=1"=0, "1<d<=5"=0, "d>5"=0, check.names = F)
  Id <- data.frame("I"=I, "0<d<=1"=0, "1<d<=5"=0, "d>5"=0, check.names = F) 
  for (u in U) {
    iu <- dat[which(dat[,1]==u),2]
    nghbr <- dat[which(dat[,2] %in% iu),]
    r_mat <- as.matrix(spread(nghbr, I, R))
    du <- sweep(r_mat, 2, r_mat[which(r_mat[,1]==u),])
    if (dim(r_mat)[2] >= 3 & dim(r_mat)[1] >= 2){
      du <- apply(du[,-1], 1, function(x){sqrt(sum(x^2, na.rm=T))/sum(!is.na(x))})
    }else{
      du <- sqrt(du[,-1]**2)
    }
    Ud[which(Ud$U==u),2:4] <- c(sum(du>0 & du<=1), sum(du>1 & du<=5), sum(du>5))
    cat("Ud Finsihed",u, "\n")
  }
  for (i in I) {
    ui <- dat[which(dat[,2]==i),1]
    nghbr <- dat[which(dat[,1] %in% ui),]
    r_mat <- as.matrix(spread(nghbr, U, R))
    di <- sweep(r_mat, 2, r_mat[which(r_mat[,1]==i),])
    if (dim(r_mat)[2] >= 3 & dim(r_mat)[1] >= 2){
      di <- apply(di[,-1], 1, function(x){sqrt(sum(x^2, na.rm=T))/sum(!is.na(x))})
    }else{
      di <- sqrt(di[,-1]**2)
    }
    Id[which(Id$I==i),2:4] <- c(sum(di>0 & di<=1), sum(di>1 & di<=5), sum(di>5))
    cat("Id Finsihed",i, "\n")
  }
  list(Ud=Ud, Id=Id)
}

SparseDistance <- function(dat, idx, flag){
  colnames(dat) <- c("U","I","R")
  Du <- list()
  Di <- list()
  if (flag == "U"){
    iu <- dat[which(dat$U==idx),]$I
    nghbr <- dat[which(dat$I %in% iu),]
    Uj <- sort(unique(nghbr$U))
    r_mat <- as.matrix(spread(nghbr, I, R))
    du <- sweep(r_mat, 2, r_mat[which(r_mat[,1]==idx),])
    if (dim(r_mat)[2] >= 3 & dim(r_mat)[1] >= 2){
      du <- apply(du[,-1], 1, function(x){sqrt(sum(x^2, na.rm=T))/sum(!is.na(x))})
    }else{
      du <- sqrt(du[,-1]**2)
    }
    return(data.frame("Ui"=idx, "Uj"=Uj, "du"=du))
  }
  if (flag == "I"){
    ui <- dat[which(dat$I==idx),]$U
    nghbr <- dat[which(dat$U %in% ui),]
    Ij <- sort(unique(nghbr$I))
    r_mat <- as.matrix(spread(nghbr, U, R))
    di <- sweep(r_mat, 2, r_mat[which(r_mat[,1]==idx),])
    if (dim(r_mat)[2] >= 3 & dim(r_mat)[1] >= 2){
      di <- apply(di[,-1], 1, function(x){sqrt(sum(x^2, na.rm=T))/sum(!is.na(x))})
    }else{
      di <- sqrt(di[,-1]**2)
    }
    return(data.frame("Ii"=idx, "Ij"=Ij, "di"=di))
  }
}

BiNeighborDistance <- function(train, test){
  colnames(train) <- c("U", "I", "R")
  colnames(test) <- c("U", "I", "R")
  # Bi-normalization
  Rubar <- train %>% group_by(U) %>% summarise(rubar=mean(R), .groups='drop') # .groups is experimental
  Ribar <- train %>% group_by(I) %>% summarise(ribar=mean(R), .groups='drop')
  train <- left_join(train, Rubar, by="U") %>% 
    left_join(., Ribar, by="I") %>% 
    group_by(U, I) %>% 
    mutate(R=R-mean(c(rubar,ribar), na.rm=T)) %>% ungroup %>% 
    select(U,I,R)
  test <- left_join(test, Rubar, by="U") %>% 
    left_join(., Ribar, by="I") %>% 
    group_by(U, I) %>% 
    mutate(R=R-mean(c(rubar,ribar), na.rm=T)) %>% ungroup %>% 
    select(U,I,R)
  
  Utst <- unique(test$U)
  Itst <- unique(test$I)
  Utrn <- unique(train$U)
  Itrn <- unique(train$I)
  Du <- list()
  Di <- list()
  for (u in Utst) {
    if (u %in% Utrn){
      Du[[u]] <- SparseDistance(train, u, "U")
    }else{
      Du[[u]] <- data.frame("Ui"=NA, "Uj"=NA, "du"=NA)
    }
  }
  for (i in Itst) {
    if (i %in% Itrn){
      Di[[i]] <- SparseDistance(train, i, "I")
    }else{
      Di[[i]] <- data.frame("Ii"=NA, "Ij"=NA, "di"=NA)
    }
  }
  D <- list()
  Y <- list()
  for (k in 1:dim(test)[1]){
    u <- test$U[k]
    i <- test$I[k]
    R_u_trgt <- train[which(train$U==u),]
    u_nghbr <- train %>% filter(I %in% R_u_trgt$I) %>% select(U) %>% unique
    R_u_nghbr <- train %>% filter(U %in% u_nghbr$U) %>% 
      select(U, I, R)
    R_i_trgt <- train[which(train$I==i),]
    i_nghbr <- train %>% filter(U %in% R_i_trgt$U) %>% select(I) %>% unique
    R_i_nghbr <- train %>% filter(I %in% i_nghbr$I) %>% 
      select(U,I,R)
    Nui <- union(R_u_nghbr, R_i_nghbr)
    M <- Nui %>% 
      left_join(Du[[u]], by=c("U"="Uj")) %>% 
      left_join(Di[[i]], by=c("I"="Ij"))
    if (length(which(M$U==u & M$I==i))>0){
      M <- M[-which(M$U==u & M$I==i),]
    }
    D[[k]] <- M %>% select(du, di) %>% as.matrix()
    Y[[k]] <- M %>% select(R) %>% as.matrix()
    if(k %% 1000 == 0){cat("Finsihed", k, "\n")}
  }
  return(list(D=D, Y=Y))
}

BiPairwiseDistance <- function(train, test){
  D <- list()
  Y <- list()
  tag <- list()
  Rubar <- train %>% group_by(userId) %>% summarise(rubar=mean(rating))
  Ribar <- train %>% group_by(movieId) %>% summarise(ribar=mean(rating))
  train <- dplyr::left_join(train, Rubar, by="userId") %>% left_join(., Ribar, by="movieId") %>% 
    group_by(userId, movieId) %>% mutate(rating=rating-mean(c(rubar,ribar), na.rm=T)) %>% ungroup %>% 
    select(userId, movieId, rating)
  test <- dplyr::left_join(test, Rubar, by="userId") %>% left_join(., Ribar, by="movieId") %>% 
    group_by(userId, movieId) %>% mutate(rating=rating-mean(c(rubar,ribar), na.rm=T)) %>% ungroup %>% 
    select(userId, movieId, rating)
  for (k in 1:dim(test)[1]) {
    u <- test$userId[k]
    i <- test$movieId[k]
    R_u_nghbr <- data.frame(userId = numeric(0), movieId = numeric(0), rating = numeric(0))
    R_i_nghbr <- data.frame(userId = numeric(0), movieId = numeric(0), rating = numeric(0))
    D_u <- NULL
    D_i <- NULL
    if (u %in% train$userId){
      R_u_trgt <- train %>% filter(userId==u) %>% select(userId, movieId, rating)
      u_nghbr <- train %>% filter(movieId %in% R_u_trgt$movieId) %>% select(userId) %>% unique
      R_u_nghbr <- train %>% filter(userId %in% u_nghbr$userId) %>% 
        select(userId, movieId, rating)
      R_u <- R_u_nghbr %>% spread(movieId, rating)
      rowidx_u <- R_u$userId
      trgt_u <- which(R_u$userId == u)
      R_u <- as.matrix(R_u[,-1])
      d_u <- sweep(R_u, 2, R_u[trgt_u,]) %>% apply(., 1, function(x){sqrt(sum(x^2, na.rm=T))/sum(!is.na(x))})
      D_u <- data.frame(up = rowidx_u, distance.u=d_u)
    }
    if (i %in% train$movieId){
      R_i_trgt <- train %>% filter(movieId==i) %>% select(movieId, userId, rating)
      i_nghbr <- train %>% filter(userId %in% R_i_trgt$userId) %>% select(movieId) %>% unique
      R_i_nghbr <- train %>% filter(movieId %in% i_nghbr$movieId) %>% 
        select(movieId, userId, rating)
      R_i <- R_i_nghbr %>% spread(userId, rating)
      rowidx_i <- R_i$movieId
      trgt_i <- which(R_i$movieId == i)
      R_i <- as.matrix(R_i[,-1])
      d_i <- sweep(R_i, 2, R_i[trgt_i,]) %>% apply(., 1, function(x){sqrt(sum(x^2, na.rm=T))/sum(!is.na(x))})
      D_i <- data.frame(ip=rowidx_i, distance.i=d_i)
    }
    # Filtering step: no NA and no larger than 2 in distance for both du, di
    R_ui <- dplyr::union(R_u_nghbr, R_i_nghbr)
    
    if (dim(R_ui)[1] == 0){
      cat(k)
      df <- data.frame(userId=numeric(0), movieId=numeric(0), rating=numeric(0),
                       distance.u=numeric(0), distance.i=numeric(0))
    }else{
      if (!purrr::is_empty(D_u)){
        df <- left_join(R_ui, D_u, by=c("userId"="up"))
      }else{
        df <- cbind(R_ui, distance.u=NA)
      }
      if (!purrr::is_empty(D_i)){
        df <- left_join(df, D_i, by=c("movieId"="ip"))
      }else{
        df <- cbind(df, distance.i=NA)
      }
    }
    # compensation for only u or only i
    #df <- df %>% mutate(distance.i = replace(distance.i, distance.i>2, NA),
    #                    distance.u = replace(distance.u, distance.u>2, NA)) %>% 
    #  filter(! (is.na(distance.u) & is.na(distance.i)) )
    
    if (dim(df)[1] <= 1){
      D[[k]] <- matrix(c(NA, NA), nrow=1)
      Y[[k]] <- NA
    }else{
      D[[k]] <- df %>% select(-c(userId, movieId, rating)) %>% as.matrix
      Y[[k]] <- df$rating
      N <- dim(R_ui)[1]
      du_na <- which(is.na(df$distance.u))
      di_na <- which(is.na(df$distance.i))
      dui_na <- sum(du_na %in% di_na)
      tag[[k]] <- cbind(length(du_na)/N, length(di_na)/N, dui_na/N)
    }
    if(k %% 1000 == 0){cat("Finsihed", k, "\n")}
  }
  list(D=D, Y=Y, tag=tag)
}

NWregUniNorm <- function(D, Y, hu, hi, train, test){
  colnames(train) <- c("U","I","R")
  colnames(test) <- c("U","I","R")
  # D: two columns of distances for user and movie
  ytrue <- test$R
  Rubar <- train %>% group_by(U) %>% summarise(rubar=mean(R))
  Ribar <- train %>% group_by(I) %>% summarise(ribar=mean(R))
  
  Ku <- lapply(D, function(x){dnorm(x[,1]/hu)})
  Ki <- lapply(D, function(x){dnorm(x[,2]/hi)})
  rhatu <- mapply(function(X,Y) {sum(X*Y, na.rm=T)/(sum(X,na.rm=T)+1e-8)}, X=Ku, Y=Y)
  rhati <- mapply(function(X,Y) {sum(X*Y, na.rm=T)/(sum(X,na.rm=T)+1e-8)}, X=Ki, Y=Y)
  df <- data.frame( "pred"= apply(cbind(rhatu, rhati), 1, function(x){mean(x, na.rm=T)}),
                    "true"= ytrue)
  
  dfE <- left_join(test, Rubar, by="U") %>% 
    left_join(., Ribar, by="I") %>%
    group_by(U, I) %>% 
    mutate(demean=mean(c(rubar,ribar), na.rm=T)) %>% ungroup %>% 
    select(U, I, demean)
  df$pred <- df$pred + dfE$demean
  df
}

UniHGridSearch <- function(D, Y, train, test){
  res <- expand.grid("hu"=c(0.001, 0.01, 0.1, 0.5, 1, 5, 7),
                     "hi"=c(0.001, 0.01, 0.1, 0.5, 1, 5, 7),
                     stringsAsFactors = FALSE) %>% as.data.frame()
  res <- cbind(res, "rmse"=rep(0, dim(res)[1]))
  for (i in 1:dim(res)[1]){
    model <- NWregUniNorm(D, Y, res[i,1], res[i,2], train, test)
    res[i,3] <- sqrt(mean((model$pred-model$true)^2))
  }
  #idx <- which.min(res[,3])
  #list(hu = res[idx,1], hi = res[idx,2])
  res
}

NWregBiNorm <- function(D, Y, hu, hi, train, test){
  colnames(train) <- c("U","I","R")
  colnames(test) <- c("U","I","R")
  ytrue <- test$R
  K <- lapply(D, function(x){
    # apply(x, 1, function(y){
    #if (is.na(y[1])){ dnorm(y[2]/hi) }
    #else if (is.na(y[2])) {dnorm(y[1]/hu)}
    #else { mvtnorm::dmvnorm(x=c(y[1]/hu, y[2]/hi), mean = c(0, 0))}
    # prod(dnorm(y / c(hu, hi)), na.rm = T)
    # })
    x <- as.matrix(x)
    matrixStats::colProds(dnorm(t(x)/c(hu, hi)), na.rm = T)
  })
  df <- data.frame( "pred"= mapply(function(X,Y) {sum(X*Y, na.rm = T)/(sum(X, na.rm = T)+1e-8)}, X=K, Y=Y),
                    "true"= ytrue)
  df$pred[which(df$pred==0)] <- NA
  Rubar <- train %>% group_by(U) %>% summarise(rubar=mean(R))
  Ribar <- train %>% group_by(I) %>% summarise(ribar=mean(R))
  dfE <- left_join(test, Rubar, by="U") %>% 
    left_join(., Ribar, by="I") %>%
    group_by(U, I) %>% 
    mutate(demean=mean(c(rubar,ribar), na.rm=T)) %>% ungroup %>% 
    select(U, I, demean)
  df$pred <- df$pred + dfE$demean
  df
}

BiHGridSearch <- function(D, Y, train, test){
  res <- expand.grid("hu"=c(0.001, 0.01, 0.1, 0.5, 1, 5, 7),
                     "hi"=c(0.001, 0.01, 0.1, 0.5, 1, 5, 7),
                     stringsAsFactors = FALSE) %>% as.data.frame()
  res <- cbind(res, "rmse"=rep(0, dim(res)[1]))
  for (i in 1:dim(res)[1]){
    model <- NWregBiNorm(D, Y, res[i,1], res[i,2], train, test)
    res[i,3] <- sqrt(mean((model$pred-model$true)^2, na.rm=T))
  }
  # idx <- which.min(res[,3])
  # list(hu = res[idx,1], hi = res[idx,2])
  res
}

NWregDistanceFilter <- function(D, hu, hi){
  Ud <- lapply(D, function(x){unique(na.omit(x[,1]))} ) %>% 
    lapply(., function(x){c(sum(x>0 & x<=1), sum(x>1 & x<=5), sum(x>5))}) %>% 
    do.call(rbind, .) %>% unique %>% as.data.frame()
  Id <- lapply(D, function(x){unique(na.omit(x[,2]))} ) %>% 
    lapply(., function(x){c(sum(x>0 & x<=1), sum(x>1 & x<=5), sum(x>5))}) %>% 
    do.call(rbind, .) %>% unique %>% as.data.frame()
  UHd <- lapply(D, function(x){unique(na.omit(x[,1]/hu))} ) %>% 
    lapply(., function(x){c(sum(x>0 & x<=1), sum(x>1 & x<=5), sum(x>5))}) %>% 
    do.call(rbind, .) %>% unique %>% as.data.frame()
  IHd <- lapply(D, function(x){unique(na.omit(x[,2]/hi))} ) %>% 
    lapply(., function(x){c(sum(x>0 & x<=1), sum(x>1 & x<=5), sum(x>5))}) %>% 
    do.call(rbind, .) %>% unique %>% as.data.frame()
  colnames(Ud) <- c("0<d<=1", "1<d<=5", "d>5")
  colnames(Id) <- c("0<d<=1", "1<d<=5", "d>5")
  colnames(UHd) <- c("0<d<=1", "1<d<=5", "d>5")
  colnames(IHd) <- c("0<d<=1", "1<d<=5", "d>5")
  list(Ud=Ud, UHd=UHd, Id=Id, IHd=IHd)
}

####################################################
## LOOCV on train
####################################################
LOOBiNeighbor <- function(dat){
  colnames(dat) <- c("U", "I", "R")
  # Bi-normalization
  Rubar <- dat %>% group_by(U) %>% summarise(rubar=mean(R), .groups='drop') # .groups is experimental
  Ribar <- dat %>% group_by(I) %>% summarise(ribar=mean(R), .groups='drop')
  dat <- left_join(dat, Rubar, by="U") %>% 
    left_join(., Ribar, by="I") %>% 
    group_by(U, I) %>% 
    mutate(R=R-mean(c(rubar,ribar), na.rm=T)) %>% ungroup %>% 
    select(U,I,R)

  Hgrid <- expand.grid("hu"=c(0.01, 0.1, 0.5, 1, 3, 5),
                      "hi"=c(0.01, 0.1, 0.5, 1, 3, 5),
                      stringsAsFactors = FALSE)
  N <- dim(dat)[1]
  ngrid <- dim(Hgrid)[1]
  rmse1 <- data.frame(matrix(0, nrow=N, ncol=ngrid)) -> rmse2
  for (k in 1:N){
    u <- dat$U[k]
    i <- dat$I[k]
    ytrue <- dat$R[k]
    R_u_nghbr <- data.frame(U = numeric(0), I = numeric(0), R = numeric(0))
    R_i_nghbr <- data.frame(U = numeric(0), I = numeric(0), R = numeric(0))
    D_u <- NULL
    D_i <- NULL
    LOOdat <- dat[-k,]
    if (u %in% LOOdat$U){
      R_u_trgt <- LOOdat %>% filter(U==u) %>% select(U,I,R)
      u_nghbr <- LOOdat %>% filter(I %in% R_u_trgt$I) %>% select(U) %>% unique
      R_u_nghbr <- LOOdat %>% filter(U %in% u_nghbr$U) %>% select(U,I,R)
      R_u <- R_u_nghbr %>% spread(I, R)
      rowidx_u <- R_u$U
      trgt_u <- which(R_u$U == u)
      R_u <- as.matrix(R_u[,-1])
      d_u <- sweep(R_u, 2, R_u[trgt_u,]) %>% apply(., 1, function(x){sqrt(sum(x^2, na.rm=T))/sum(!is.na(x))})
      D_u <- data.frame(up = rowidx_u, du=d_u)
    }
    if (i %in% LOOdat$I){
      R_i_trgt <- LOOdat %>% filter(I==i) %>% select(I,U,R)
      i_nghbr <- LOOdat %>% filter(U %in% R_i_trgt$U) %>% select(I) %>% unique
      R_i_nghbr <- LOOdat %>% filter(I%in% i_nghbr$I) %>% select(I,U,R)
      R_i <- R_i_nghbr %>% spread(U, R)
      rowidx_i <- R_i$I
      trgt_i <- which(R_i$I == i)
      R_i <- as.matrix(R_i[,-1])
      d_i <- sweep(R_i, 2, R_i[trgt_i,]) %>% apply(., 1, function(x){sqrt(sum(x^2, na.rm=T))/sum(!is.na(x))})
      D_i <- data.frame(ip=rowidx_i, di=d_i)
    }
    R_ui <- dplyr::union(R_u_nghbr, R_i_nghbr)
    
    if (dim(R_ui)[1] == 0){
      cat(k)
      df <- data.frame(U=numeric(0), I=numeric(0), R=numeric(0), du=numeric(0), di=numeric(0))
    }else{
      if (!purrr::is_empty(D_u)){
        df <- left_join(R_ui, D_u, by=c("U"="up"))
      }else{
        df <- cbind(R_ui, du=NA)
      }
      if (!purrr::is_empty(D_i)){
        df <- left_join(df, D_i, by=c("I"="ip"))
      }else{
        df <- cbind(df, di=NA)
      }
    }
    # compensation for only u or only i
    df <- df %>% mutate(di = replace(di, di>2, NA),
                        du = replace(du, du>2, NA)) %>% 
      filter(! (is.na(du) & is.na(di)) )
    
    if (dim(df)[1] <= 1){
      rmse1[k,] <- rep(NA, ngrid)
      rmse2[k,] <- rep(NA, ngrid)
    }else{
      D <- df %>% select(-c(U,I,R)) %>% as.matrix
      Y <- df$R
      for (h in 1:ngrid){
        Ku <- dnorm(D[,1]/Hgrid[h,1])
        Ki <- dnorm(D[,2]/Hgrid[h,2])
        rhatu <- sum(Ku*Y, na.rm=T)/(sum(Ku,na.rm=T)+1e-8)
        rhati <- sum(Ki*Y, na.rm=T)/(sum(Ki,na.rm=T)+1e-8)
        rmse1[k,h] <- (mean(c(rhatu,rhati))-ytrue)^2
        
        Ku[is.na(Ku)] <- 1
        Ki[is.na(Ki)] <- 1
        K.degenerated <- Ku*Ki -> K
        K[which(K.degenerated==1)] <- NA
        rmse2[k,h] <- (sum(K*Y, na.rm = T)/(sum(K, na.rm = T)+1e-8)-ytrue)^2
      }
    }
    if(k %% 1000 == 0){cat("Finsihed", k, "\n")}
  }
  return(list(unirmse = rmse1, birmse = rmse2))
}

LOOSmoothCheck <- function(dat){
  colnames(dat) <- c("U", "I", "R")
  # Bi-normalization
  Rubar <- dat %>% group_by(U) %>% summarise(rubar=mean(R), .groups='drop') # .groups is experimental
  Ribar <- dat %>% group_by(I) %>% summarise(ribar=mean(R), .groups='drop')
  dat <- left_join(dat, Rubar, by="U") %>% 
    left_join(., Ribar, by="I") %>% 
    group_by(U, I) %>% 
    mutate(R=R-mean(c(rubar,ribar), na.rm=T)) %>% ungroup %>% 
    select(U,I,R)
  
  N <- dim(dat)[1]
  distcathu <- matrix(0, nrow=N, ncol=5)
  distcathi <- matrix(0, nrow=N, ncol=5)
  for (k in 1:N){
    u <- dat$U[k]
    i <- dat$I[k]
    ytrue <- dat$R[k]
    R_u_nghbr <- data.frame(U = numeric(0), I = numeric(0), R = numeric(0))
    R_i_nghbr <- data.frame(U = numeric(0), I = numeric(0), R = numeric(0))
    D_u <- NULL
    D_i <- NULL
    LOOdat <- dat[-k,]
    if (u %in% LOOdat$U){
      R_u_trgt <- LOOdat %>% filter(U==u) %>% select(U,I,R)
      u_nghbr <- LOOdat %>% filter(I %in% R_u_trgt$I) %>% select(U) %>% unique
      R_u_nghbr <- LOOdat %>% filter(U %in% u_nghbr$U) %>% select(U,I,R)
      R_u <- R_u_nghbr %>% spread(I, R)
      rowidx_u <- R_u$U
      trgt_u <- which(R_u$U == u)
      R_u <- as.matrix(R_u[,-1])
      d_u <- sweep(R_u, 2, R_u[trgt_u,]) %>% apply(., 1, function(x){sqrt(sum(x^2, na.rm=T))/sum(!is.na(x))})
      D_u <- data.frame(up = rowidx_u, du=d_u)
    }
    if (i %in% LOOdat$I){
      R_i_trgt <- LOOdat %>% filter(I==i) %>% select(I,U,R)
      i_nghbr <- LOOdat %>% filter(U %in% R_i_trgt$U) %>% select(I) %>% unique
      R_i_nghbr <- LOOdat %>% filter(I%in% i_nghbr$I) %>% select(I,U,R)
      R_i <- R_i_nghbr %>% spread(U, R)
      rowidx_i <- R_i$I
      trgt_i <- which(R_i$I == i)
      R_i <- as.matrix(R_i[,-1])
      d_i <- sweep(R_i, 2, R_i[trgt_i,]) %>% apply(., 1, function(x){sqrt(sum(x^2, na.rm=T))/sum(!is.na(x))})
      D_i <- data.frame(ip=rowidx_i, di=d_i)
    }
    R_ui <- dplyr::union(R_u_nghbr, R_i_nghbr)
    
    if (dim(R_ui)[1] == 0){
      cat(k)
      df <- data.frame(U=numeric(0), I=numeric(0), R=numeric(0), du=numeric(0), di=numeric(0))
    }else{
      if (!purrr::is_empty(D_u)){
        df <- left_join(R_ui, D_u, by=c("U"="up"))
      }else{
        df <- cbind(R_ui, du=NA)
      }
      if (!purrr::is_empty(D_i)){
        df <- left_join(df, D_i, by=c("I"="ip"))
      }else{
        df <- cbind(df, di=NA)
      }
    }
    # compensation for only u or only i
    #df <- df %>% mutate(di = replace(di, di>2, NA),
    #                    du = replace(du, du>2, NA)) %>% 
    #  filter(! (is.na(du) & is.na(di)) )
    
    if (dim(df)[1] <= 1){
      distcathu[k,] <- rep(NA,5)
      distcathi[k,] <- rep(NA,5)
    }else{
      D <- df %>% select(-c(U,I,R)) %>% as.matrix
      Y <- df$R
      Ku <- D[,1]
      Ki <- D[,2]
      Kuna <- sum(is.na(Ku))
      Kina <- sum(is.na(Ki))
      Ku <- na.omit(Ku)
      Ki <- na.omit(Ki)
      distcathu[k,] <- c(sum(Ku==0), sum(Ku>0 & Ku<=1), sum(Ku>1 & Ku<=3), sum(Ku>3), Kuna)
      distcathi[k,] <- c(sum(Ki==0), sum(Ki>0 & Ki<=1), sum(Ki>1 & Ki<=3), sum(Ki>3), Kina)  
    }
    if(k %% 1000 == 0){cat("Finsihed", distcathu[k,], "\n")}
  }
  return(list(hu = distcathu, hi = distcathi))
}


####
#### New version of multi-level neighboring smooth
####
BiPairwiseDistance <- function(train, test){
  D <- list()
  Y <- list()
  tag <- list()
  Rubar <- train %>% group_by(userId) %>% summarise(rubar=mean(rating))
  Ribar <- train %>% group_by(movieId) %>% summarise(ribar=mean(rating))
  train <- dplyr::left_join(train, Rubar, by="userId") %>% left_join(., Ribar, by="movieId") %>% 
    group_by(userId, movieId) %>% mutate(rating=rating-mean(c(rubar,ribar), na.rm=T)) %>% ungroup %>% 
    select(userId, movieId, rating)
  test <- dplyr::left_join(test, Rubar, by="userId") %>% left_join(., Ribar, by="movieId") %>% 
    group_by(userId, movieId) %>% mutate(rating=rating-mean(c(rubar,ribar), na.rm=T)) %>% ungroup %>% 
    select(userId, movieId, rating)
  for (k in 1:dim(test)[1]) {
    u <- test$userId[k]
    i <- test$movieId[k]
    NIu <- data.frame(userId = numeric(0), movieId = numeric(0), rating = numeric(0))
    NUi <- data.frame(userId = numeric(0), movieId = numeric(0), rating = numeric(0))
    if (u %in% train$userId){
      if (u != test$userId[k-1] || k == 1){
        Iu <- train[which(train$userId==u),]$movieId
        UIu <- train[which(train$movieId %in% Iu),]$userId %>% unique
        NIu <- train[which(train$userId %in% UIu),]
        Mu <- NIu %>% spread(movieId, rating)
        rowidx_u <- Mu$userId
        ur <- as.matrix(Mu[which(Mu$userId==u),-1])
        Mu <- as.matrix(Mu[,-1])
        #d_u <- as.matrix(dist(rbind(ur, Mu)))[1]
        d_u <- apply(Mu, 1, function(x){sqrt(sum((x-ur)^2, na.rm=T)/sum(!is.na(x)))})
        D_u <- data.frame(up = rowidx_u, distance.u=d_u)
        D_up <- D_u
      }else{
        D_u <- D_up
      }
    }else{
      D_u <- NULL
    }
    if (i %in% train$movieId){
      Ui <- train[which(train$movieId==i),]$userId
      IUi <- train[which(train$userId %in% Ui),]$movieId %>% unique
      NUi <- train[which(train$movieId %in% IUi),]
      Mi <- NUi %>% spread(userId, rating)
      rowidx_i <- Mi$movieId
      ir <- as.matrix(Mi[which(Mi$movieId==i),-1])
      Mi <- as.matrix(Mi[,-1])
      d_i <- apply(Mi, 1, function(x){sqrt(sum((x-ir)^2, na.rm=T)/sum(!is.na(x)))})
      D_i <- data.frame(ip=rowidx_i, distance.i=d_i)
    }else{
      D_i <- NULL
    }
    # Filtering step: no NA and no larger than 2 in distance for both du, di
    R_ui <- dplyr::union(NIu, NUi)
    
    if (dim(R_ui)[1] == 0){
      cat(k, ",")
      df <- data.frame(userId=numeric(0), movieId=numeric(0), rating=numeric(0),
                       distance.u=numeric(0), distance.i=numeric(0))
    }else{
      if (!is.null(D_u)){
        df <- left_join(R_ui, D_u, by=c("userId"="up"))
      }else{
        df <- cbind(R_ui, distance.u=NA)
      }
      if (!is.null(D_i)){
        df <- left_join(df, D_i, by=c("movieId"="ip"))
      }else{
        df <- cbind(df, distance.i=NA)
      }
    }
    # compensation for only u or only i
    #df <- df %>% mutate(distance.i = replace(distance.i, distance.i>2, NA),
    #                    distance.u = replace(distance.u, distance.u>2, NA)) %>% 
    #  filter(! (is.na(distance.u) & is.na(distance.i)) )
    
    if (dim(df)[1] <= 1){
      D[[k]] <- matrix(c(NA, NA), nrow=1)
      Y[[k]] <- NA
    }else{
      D[[k]] <- df[,c("distance.u", "distance.i")] %>% as.matrix
      Y[[k]] <- df$rating
    }
    if(k %% 1000 == 0){cat("Finsihed", k, "\n")}
  }
  list(D=D, Y=Y)
}


#########
### Modified on Jan 29, 2022
#########
library(tidyverse)
bi_center <- function(train, test){
  colnames(train) <- c("U", "I", "R")
  colnames(test) <- c("U", "I", "R")
  Rubar <- train %>% group_by(U) %>% summarise(rubar=mean(R))
  Ribar <- train %>% group_by(I) %>% summarise(ribar=mean(R))
  trainc <- dplyr::left_join(train, Rubar, by="U") %>% 
    left_join(., Ribar, by="I") %>% 
    group_by(U, I) %>% 
    mutate(R = R - mean(c(rubar,ribar), na.rm=T)) %>% 
    ungroup %>% 
    select(U, I, R)
  testc <- dplyr::left_join(test, Rubar, by="U") %>% 
    left_join(., Ribar, by="I") %>% 
    group_by(U, I) %>% 
    mutate(R = R - mean(c(rubar,ribar), na.rm=T)) %>% 
    ungroup %>% 
    select(U, I, R)
  return(list(train=trainc, test=testc))
}

pair_distance <- function(dat, byrow){
  colnames(dat) <- c("U", "I", "R")
  if (isTRUE(byrow)){
    M <- tidyr::spread(dat, U, R)
  }else{
    M <- tidyr::spread(dat, I, R)
  }
  D <- as.matrix(dist(M[,-1], method="euclidean"))
  adj_D <- sqrt(D^2/dim(M[,-1])[2])
  D_df <- na.omit(reshape2::melt(adj_D))
  return(D_df)
}

radial_neighbor_profile <- function(dist_u, dist_i, trainc, testc=NULL){
  # NOTE: trainc, dist_u, dist_i have been updated colnames as U, I, R
  # trainc only for bandwidth optimization
  # testc for prediction
  prof <- list()
  if (is.null(testc)){
    L <- dim(trainc)[1]
    target_set <- trainc
  }else{
    L <- dim(testc)[1]
    target_set <- testc
  }
  NIu <- data.frame(U = numeric(0), I = numeric(0), R = numeric(0))
  NUi <- data.frame(U = numeric(0), I = numeric(0), R = numeric(0))
  for (k in 1:L) {
    u <- target_set$U[k]
    i <- target_set$I[k]
    if (u %in% trainc$U){
      if (u != trainc$U[k-1] || k == 1){
        Iu <- trainc[which(trainc$U==u),]$I
        UIu <- trainc[which(trainc$I %in% Iu),]$U %>% unique
        NIu <- trainc[which(trainc$U %in% UIu),]
      }else{
        NIu <- NIu
      }
    }else{
      NIu <- data.frame(U = numeric(0), I = numeric(0), R = numeric(0))
    }
    if (i %in% trainc$I){
      Ui <- trainc[which(trainc$I==i),]$U
      IUi <- trainc[which(trainc$U %in% Ui),]$I %>% unique
      NUi <- trainc[which(trainc$I %in% IUi),]
    }else{
      NUi <- data.frame(U = numeric(0), I = numeric(0), R = numeric(0))
    }
    radial_neighbor <- dplyr::filter(dplyr::union(NIu, NUi), U != u | I != i)
    radial_neighbor_d <- radial_neighbor %>% 
      dplyr::left_join(dist_u[which(dist_u$Var2==u), c(1,3)], 
                       by=c("U"="Var1")) %>% 
      dplyr::left_join(dist_i[which(dist_i$Var2==i), c(1,3)],
                       by=c("I"="Var1"))
    colnames(radial_neighbor_d)[4:5] <- c("du", "di")
    prof[[k]] <- dplyr::filter(radial_neighbor_d, !is.na(du) | !is.na(di))
    # cat("Finished", k, "\n")
  }
  return(prof)
}

NWreg_univ <- function(prof, hu, hi, train, test=NULL){
  colnames(train) <- c("U","I","R")
  Rubar <- train %>% group_by(U) %>% summarise(rubar=mean(R))
  Ribar <- train %>% group_by(I) %>% summarise(ribar=mean(R))
  Ku <- lapply(prof, function(x){dnorm(x$du/hu)})
  Ki <- lapply(prof, function(x){dnorm(x$di/hi)})
  Y <- lapply(prof, function(x){x$R})
  rhatu <- mapply(function(X,Y) {sum(X*Y, na.rm=T)/(sum(X,na.rm=T)+1e-8)}, X=Ku, Y=Y)
  rhati <- mapply(function(X,Y) {sum(X*Y, na.rm=T)/(sum(X,na.rm=T)+1e-8)}, X=Ki, Y=Y)
  if (!is.null(test)){
    colnames(test) <- c("U","I","R")
    df <- data.frame( "pred"= apply(cbind(rhatu, rhati), 1, function(x){mean(x, na.rm=T)}),
                      "true"= test$R)
    df$pred[which(df$pred==0)] <- NA
    dfE <- left_join(test, Rubar, by="U") %>% 
      left_join(., Ribar, by="I") %>%
      group_by(U, I) %>% 
      mutate(demean=mean(c(rubar,ribar), na.rm=T)) %>% ungroup %>% 
      select(U, I, demean)
  }else{
    df <- data.frame( "pred"= apply(cbind(rhatu, rhati), 1, function(x){mean(x, na.rm=T)}),
                      "true"= train$R)
    df$pred[which(df$pred==0)] <- NA
    dfE <- left_join(train, Rubar, by="U") %>% 
      left_join(., Ribar, by="I") %>%
      group_by(U, I) %>% 
      mutate(demean=mean(c(rubar,ribar), na.rm=T)) %>% ungroup %>% 
      select(U, I, demean)
  }
  df$pred <- df$pred + dfE$demean
  return(df)
}

bandwith_grid_search_univ <- function(prof, train){
  res <- expand.grid("hu"=c(0.001, 0.01, 0.1, 0.5, 1, 5, 7),
                     "hi"=c(0.001, 0.01, 0.1, 0.5, 1, 5, 7),
                     stringsAsFactors = FALSE) %>% as.data.frame()
  res <- cbind(res, "rmse"=rep(0, dim(res)[1]))
  for (i in 1:dim(res)[1]){
    model <- NWreg_univ(prof, res[i,1], res[i,2], train)
    res[i,3] <- sqrt(mean((model$pred-model$true)^2, na.rm=T))
  }
  idx <- which.min(res[,3])
  list(hu = res[idx,1], hi = res[idx,2])
  # res
}

NWreg_biv <- function(prof, hu, hi, train, test=NULL){
  colnames(train) <- c("U","I","R")
  Rubar <- train %>% group_by(U) %>% summarise(rubar=mean(R))
  Ribar <- train %>% group_by(I) %>% summarise(ribar=mean(R))
  Y <- lapply(prof, function(x){x$R})
  W1 <- lapply(prof, function(x){dnorm(as.matrix(x[, 4])/hu)})
  W2 <- lapply(prof, function(x){dnorm(as.matrix(x[, 5])/hu)})
  K <- mapply(function(x, y){x * y}, x=W1, y=W2)
  if (!is.null(test)){
    colnames(test) <- c("U","I","R")
    df <- data.frame( "pred"= mapply(function(X,Y) {sum(X*Y, na.rm = T)/(sum(X, na.rm = T)+1e-8)}, X=K, Y=Y),
                      "true"= test$R)
    df$pred[which(df$pred==0)] <- NA
    dfE <- left_join(test, Rubar, by="U") %>% 
      left_join(., Ribar, by="I") %>%
      group_by(U, I) %>% 
      mutate(demean=mean(c(rubar,ribar), na.rm=T)) %>% ungroup %>% 
      select(U, I, demean)
  }else{
    df <- data.frame( "pred"= mapply(function(X,Y) {sum(X*Y, na.rm = T)/(sum(X, na.rm = T)+1e-8)}, X=K, Y=Y),
                      "true"= train$R)
    df$pred[which(df$pred==0)] <- NA
    dfE <- left_join(train, Rubar, by="U") %>% 
      left_join(., Ribar, by="I") %>%
      group_by(U, I) %>% 
      mutate(demean=mean(c(rubar,ribar), na.rm=T)) %>% ungroup %>% 
      select(U, I, demean)
  }
  df$pred <- df$pred + dfE$demean
  return(df)
}

bandwith_grid_search_biv <- function(prof, train){
  res <- expand.grid("hu"=c(0.001, 0.01, 0.1, 0.5, 1, 5, 7),
                     "hi"=c(0.001, 0.01, 0.1, 0.5, 1, 5, 7),
                     stringsAsFactors = FALSE) %>% as.data.frame()
  res <- cbind(res, "rmse"=rep(0, dim(res)[1]))
  for (i in 1:dim(res)[1]){
    model <- NWreg_biv(prof, res[i,1], res[i,2], train)
    res[i,3] <- sqrt(mean((model$pred-model$true)^2, na.rm=T))
  }
  idx <- which.min(res[,3])
  list(hu = res[idx,1], hi = res[idx,2])
  # res
}

# center the data based on the train set
dat_center <- bi_center(mltrn, mltst)
# construct the distances between users and items
dist_u <- pair_distance(dat_center$train, byrow = TRUE)
dist_i <- pair_distance(dat_center$train, byrow = FALSE)
# radial neighbor profiles for train set
profile_train <- radial_neighbor_profile(dist_u, dist_i, dat_center$train)
# radial neighbor profiles for test set based on train
profile_test <- radial_neighbor_profile(dist_u, dist_i, dat_center$train, dat_center$test)
# optimize bandwith over grid search
H_univ <- bandwith_grid_search_univ(profile_train, dat_center$train)
H_biv <- bandwith_grid_search_biv(profile_train, dat_center$train)
# predict on the test set
fit_nw_univ <- NWreg_univ(profile_test, H_univ$hu, H_univ$hi, dat_center$train, dat_center$test)
fit_nw_biv <- NWreg_biv(profile_test, H_univ$hu, H_univ$hi, dat_center$train, dat_center$test)
sqrt(mean((fit_nw_univ$true-fit_nw_univ$pred)^2, na.rm=T))
sqrt(mean((fit_nw_biv$true-fit_nw_biv$pred)^2, na.rm=T))
sum(is.na(fit_nw_univ$pred))/dim(mltst)[1]
sum(is.na(fit_nw_biv$pred))/dim(mltst)[1]

