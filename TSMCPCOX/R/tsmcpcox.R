#' Two stage multiple change point detection in Cox models.
#'
#' This function provides a two-stage procedure for simultaneously detecting multiple
#' change-points in Cox models. In the cutting stage, the change-point problem is
#' converted into a model selection problem so that a modern model selection method
#' can be applied. In the refining stage, the change-points obtained in the cutting
#' stage are finalized via a refining method. The tuning parameter lambda is chosen
#' by BIC.
#'
#'
#' @param Y the right-censored failure time.
#' @param Delta the censoring indicator
#' @param X the n-by-p design matrix.
#' @param Z the covariate where the change points are.
#' @param method the method to be used by  mcp or scad.
#' @param c The number of candidates for the segment length
#' @return tsmcpcox returns an object of class "tsmcpcox".
#' An object of class "tsmcpcox" is a list containing the following components:
#' @return \item{the number of change.points}{estimators of change points.}
#' @return \item{change.points}{estimators of change points.}
#' @return \item{coefficients}{estimators of coefficients}
#' @return \item{m.opt}{the selected segment length}
#' @export tsmcpcox
#' @import ncvreg SIS survival NMOF
#' @seealso plus lars

#' @examples ##  Cox model
#'
#' ##
#' library(TSMCPCOX)
# library(survival)
# library(NMOF)
# library(ncvreg)
# library(SIS)
# set.seed(10123)
# n=1000
# u<-runif(n,0,1)
# X<-cbind(rbinom(n,1,0.5),rnorm(n,0,1))
# Z=sort(runif(n,0,1))
# C<-rexp(n,0.5)  #censoring time
# tau<-c(0.3,0.7) ## true change point location
# beta<-c(1,1) #cofficients in first segment
# b1<-c(-2,0)  #cofficients in second segment
# b2<-c(3,-3)  #cofficients in third segment
# r<-drop(X%*%beta)+drop(X%*%b1)*(Z>tau[1])+drop(X%*%b2)*(Z>tau[2])
# T<--log(u)/exp(r)# set lambda0(t)=1
# Y<-ifelse(T<=C,T,C)
# Delta<-ifelse(T<=C,1,0)
# 
# ##estimate change points
#   tsmcpcox(Y=Y, Delta=Delta,X1=X,Z=Z, method = "MCP")




  tsmcpcox <- function(Y, Delta,X1,X2=NULL,Z, method = c("MCP", "SCAD"), c=25) {
    n <- length(Y)
    p <- dim(X1)[2]
    if (length(p) == 0) {
      p <- 0
    }
    p2<-dim(X2)[2]
    if (length(p2) == 0) {
      p2 <- 0
    }
    len<-seq(5,30,length=c)
    num.td<-NULL
    bicm<-NULL
    loc.td<-NULL
    x1_temp <-  X1[order(Z),]
    x2_temp<- X2[order(Z),]
    y_temp <- Y[order(Z)]
    delta_temp <- Delta[order(Z)]
    z_temp<-sort(Z)
  for(l in 1:c){
    m<-ceiling(len[l]*n^(1/5))
    q <- floor(n/m)
    f<- n-m*(q-1)
    K_temp <- matrix(0, nrow = q, ncol = q, byrow = TRUE)

    for (i in 1:q) K_temp[i, 1:i] <- rep(1, i)
    x1 <- NULL
    x2 <- NULL
    y <- NULL
    delta <- NULL
    z <- NULL
    x1[[1]] <- as.matrix(x1_temp[1:((n - (q - 1) * m)), ])
    y[[1]] <- y_temp[1:((n - (q - 1) * m))]
    delta[[1]] <- delta_temp[1:((n - (q - 1) * m))]
    z[[1]] <- z_temp[1:((n - (q - 1) * m))]
    for (i in 2:q) {
      x1[[i]] <- as.matrix(x1_temp[(n - (q - i + 1) * m + 1):((n - (q -
                                                                     i) * m)), ])
       y[[i]] <- y_temp[(n - (q - i + 1) * m + 1):((n - (q - i) * m))]
       delta[[i]] <- delta_temp[(n - (q - i + 1) * m + 1):((n - (q - i) * m))]
       z[[i]] <- z_temp[(n - (q - i + 1) * m + 1):((n - (q - i) * m))]
    }
    if(length(x2_temp)!=0)
    {
      x2[[1]] <- as.matrix(x2_temp[1:((n - (q - 1) * m)), ])
      for (i in 2:q) {
        x2[[i]] <- as.matrix(x2_temp[(n - (q - i + 1) * m + 1):((n - (q -  i) * m)), ])
      }
    }

   
    X_temp1 <- lapply(1:length(x1), function(j, mat, list) kronecker(mat[j,
                                                                        , drop = FALSE], list[[j]]), mat = K_temp, list = x1)
    Xn <- cbind(do.call("rbind", X_temp1),x2_temp)

    #########################################adlasso######################

    y_surv <- survival::Surv(y_temp,delta_temp)
    object <- tune.fit1(Xn, y_surv, family="cox",q,m,f,
                        penalty = method, tune = "bic")

    coeff<-rep(0,p*q+p2)

    coeff[object$ix] <- object$beta

    adcp.coef.s <- sum(abs(coeff))

    adcp.coef.v.m <- abs(matrix(c(coeff[1:(p*q)]), q, p , byrow = T))

    adcp.coef.m <- c(apply(adcp.coef.v.m, 1, max))

    adcp.cp <- which(adcp.coef.m != 0)

    if (length(adcp.cp) > 1) {
      for (i in 2:length(adcp.cp)) {
        if (adcp.cp[i] - adcp.cp[i - 1] == 1)
          adcp.cp[i] <- 0
      }
    }


    adcp.cp1 <- adcp.cp[adcp.cp > 1 & adcp.cp < q]

    d1 <- length(adcp.cp1)

    if (d1 == 0) {
      adcpcss.cp <- integer(0)

    }

    if (d1 >= 1) {
      # step 4 find change points
      adcpcss.cp <- NULL


      adcp.cp1 <- c(0, adcp.cp1, q + 1)
      for (i in 1:d1) {


        yt <- NULL
        x1t <- NULL
        x2t<-NULL
        deltat<-NULL
        zt<-NULL
        for (k in (adcp.cp1[i + 1] - 1):(adcp.cp1[i + 1])) 
        {
          yt <- c(yt, y[[k]])
          zt<-c(zt,z[[k]])
          deltat<-c(deltat,delta[[k]])
          x1t <- rbind(x1t, as.matrix(x1[[k]]))


        }
        if(length(x2_temp)!=0){
          for (k in (adcp.cp1[i + 1] - 1):(adcp.cp1[i + 1])){
              x2t <- rbind(x2t, as.matrix(x2[[k]]))
          }
        }

        # using profile-likelihood find change point.
        cp <- order(abs(z_temp-css(x1t,x2t, zt,yt,deltat,1.2*m)))[1]
        if (length(cp) == 0)
          next
        if (cp != 0) {
          if (adcp.cp1[i + 1] == 0)
            adcpcss.cp <- c(adcpcss.cp, cp)
          if (adcp.cp1[i + 1] > 0)
            adcpcss.cp <- c(adcpcss.cp, cp)

        }


      }

      if (length(adcpcss.cp) == 0)
        adcpcss.cp <- integer(0) else adcpcss.cp <- adcpcss.cp

    }

    tt <- which(abs(diff(adcpcss.cp)) < 10)
    if (length(tt) > 0)
      adcpcss.cp <- adcpcss.cp[-tt]


    #return(c(length(adcpcss.cp)))
    num.td[l]<-length(adcpcss.cp)
    adcpcss.cp1<-c(0,adcpcss.cp,n)
    loc.td[[l]]<-adcpcss.cp1
    q <- length(adcpcss.cp1)-1
    K_temp <- matrix(0, nrow = q, ncol = q, byrow = TRUE)

    for (i in 1:q) K_temp[i, 1:i] <- rep(1, i)
    x1 <- NULL
    for (i in 1:q) {
      x1[[i]] <- as.matrix(x1_temp[(adcpcss.cp1[i]+1):adcpcss.cp1[i+1], ])
   }

   
    X1_temp1 <- lapply(1:length(x1), function(j, mat, list) kronecker(mat[j,
                                                                        , drop = FALSE], list[[j]]), mat = K_temp, list = x1)
    Xn <- cbind(do.call("rbind", X1_temp1),x2_temp)



    coxfit<-coxph(Surv(Y,Delta)~ Xn-1)
    coeff<-coxfit$coefficients
    likvalue<-coxfit$loglik[2]
    bicm[l]<--2*likvalue+(p*(num.td[l]+1)+p2)*log(n)
  }

  l.opt<-min(which(bicm==min(bicm)))
  m.opt<-ceiling(l.opt*n^(1/5))
  loc.opt<-loc.td[[l.opt]]
  num.opt<-num.td[l.opt]
  adcpcss.cp<-loc.opt[loc.opt>1&loc.opt<n]
  if(length(loc.opt)==(length(tau)+2)){
    adcpcss.cp1<-loc.opt
    q <- length(adcpcss.cp1)-1
    K_temp <- matrix(0, nrow = q, ncol = q, byrow = TRUE)

    for (i in 1:q) K_temp[i, 1:i] <- rep(1, i)

  
    x1 <- NULL
    for (i in 1:q) {
      x1[[i]] <- as.matrix(x1_temp[(adcpcss.cp1[i]+1):adcpcss.cp1[i+1], ])
        }

X_temp1 <- lapply(1:length(x1), function(j, mat, list) kronecker(mat[j,
                                                                        , drop = FALSE], list[[j]]), mat = K_temp, list = x1)
    Xn <- cbind(do.call("rbind", X_temp1),x2_temp)
    y_surv <- survival::Surv(y_temp,delta_temp)
    object <- tune.fit(Xn, y_surv, family="cox",
                       penalty = method, tune = "bic")
    coeff<-rep(0,p*q+p2)

    coeff[object$ix] <-object$beta
    result<-list(
     change.points.num=length(Z[adcpcss.cp]),
     change.points=Z[adcpcss.cp],
     coeff=coeff,
     m=m.opt
     )
    return(result)}else{
      return( result<-list(
        change.points.num=0,
        change.points=NULL,
        coeff=coeff,
        m=m.opt
      ))
    }
}
