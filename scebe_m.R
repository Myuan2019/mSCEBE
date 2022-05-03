scebe_m <- function(mydata,x){
  N <- length(unique(mydata$ID))
  q <- ncol(x)
  p <- 2

  fit <- try(lmer(Conc~Time+(Time|ID),data = mydata),silent = TRUE)

  
  if(inherits(fit,"try-error")) return("base model error") else{
    ind.parm <- ranef(fit)$ID
    ind.parm = ind.parm[order(as.numeric(rownames(ind.parm))), ] 
    
    
    b0hat <- ind.parm$`(Intercept)`
    b1hat <- ind.parm$Time
    
    Bhat <- cbind(b0hat,b1hat)
    
    xc <- scale(x,center = TRUE,scale = FALSE)
    cov.x <- var(x)
    
    ###NEBE

    {
      gammahat <- t(Bhat)%*%xc%*%solve(t(xc)%*%xc)
      
      
      
      df=N-q-1
      
      x_all <- cbind(rep(1,N),x)
      
      H=x_all%*%solve(t(x_all)%*%x_all)%*%t(x_all)
      
      ubsigmahat <- t(Bhat)%*%(diag(N)-H)%*%Bhat/df
      
      vgammahat <- kronecker(solve(t(xc)%*%xc),ubsigmahat)
      
      
      sdmat <- matrix(sqrt(diag(vgammahat)),nrow = p,ncol = q)
      tmat <- gammahat/sdmat
      
      
      pmat <- 2*pt(abs(tmat),df=df,lower.tail = F)
      
      rst.nebe <- cbind(array(gammahat),array(sdmat),array(tmat),array(pmat))
      colnames(rst.nebe)=c("est","sd","t-value","p-value")
      
    }
    

    
    
    ###rmscebe   
    

    {
      gammahat <- t(Bhat)%*%xc%*%solve(t(xc)%*%xc)
      Ip <- diag(p)
      
      Rhat <- matrix(as.numeric(VarCorr(fit)$ID),ncol=2)
      
      nt=numeric(N)
      for (i in 1:N) {
        nt[i]=sum(mydata$ID==i)
      }
      
      Z=cbind(rep(1,nrow(mydata)),mydata$Time)
      
      
      list_OMEGA <- vector(mode = 'list',length = N)
      list_Wxi <- vector(mode = 'list',length = N)
      list_part1 <- vector(mode = 'list',length = N)
      list_part2 <- vector(mode = 'list',length = N)
      list_Ai <- vector(mode = 'list',length = N)
      BB <- matrix(0,nrow = N*p,ncol = p)
      CC <- matrix(0,nrow = p,ncol = nrow(mydata))
      
      for (i in 1:N){
        Gi <- diag(sigma(fit)^2,nt[i])
        Zi <- Z[((1+sum(nt[-(i:N)])):sum(nt[1:i])),]
        sigmai <- Zi%*%Rhat%*%t(Zi)+Gi
        ZGZRi=solve(crossprod(Zi,solve(Gi,Zi))+solve(Rhat))
        Si <- ZGZRi%*%solve(Rhat)
        list_OMEGA[[i]] <- crossprod(Zi,solve(sigmai))%*%Zi
        list_Wxi[[i]] <- kronecker(t(x[i,]),list_OMEGA[[i]])
        list_part1[[i]] <- kronecker(xc[i,]%*%t(x[i,]),Ip-Si)
        list_part2[[i]] <-kronecker(xc[i,],Si)
        list_Ai[[i]] <- ZGZRi%*%crossprod(Zi,solve(Gi))
        BB[(1+p*(i-1)):(p*i),] <- Si
        CC[,((1+sum(nt[-(i:N)])):sum(nt[1:i]))] <- crossprod(Zi,solve(sigmai))
        
      }
      
      sOMEGA <- solve(Reduce('+',list_OMEGA))
      
      Wx <- sOMEGA%*%Reduce('+',list_Wxi)
      
      
      part1 <- Reduce('+',list_part1)
      
      part2 <- Reduce('+',list_part2)
      
      
      SC <- kronecker(solve(t(xc)%*%xc),Ip)%*%(part1+part2%*%Wx)
      
      vgammasc <- solve(SC)%*%array(gammahat)
      
      gammasc <- matrix(vgammasc,nrow = p,ncol = q)
      
      
      AA <- as.matrix(bdiag(list_Ai))
      CC=sOMEGA%*%CC
      DD <- AA+BB%*%CC
      
      Rhat1 <- Rhat-gammasc%*%cov.x%*%t(gammasc)
      
      list_sigmai1 <- vector(mode = 'list',length = N)
      
      for (i in 1:N){
        Gi <- diag(sigma(fit)^2,nt[i])
        Zi <- Z[((1+sum(nt[-(i:N)])):sum(nt[1:i])),]
        list_sigmai1[[i]] <- Zi%*%Rhat1%*%t(Zi)+Gi
      }
      
      SIG1 <- as.matrix(bdiag(list_sigmai1))
      
      k_xc_ip <- kronecker(solve(t(xc)%*%xc)%*%t(xc),Ip)
      vvecg <- k_xc_ip%*%DD
      vvecg <- vvecg%*%SIG1  
      vvecg <- vvecg%*%t(DD)
      vvecg <- vvecg%*%t(k_xc_ip)
      
      Vhatrscebe <- solve(SC)%*%vvecg%*%t(solve(SC))
      sdmatrscebe <- matrix(sqrt(diag(Vhatrscebe)),nrow = p,ncol = q)
      
      tmatrscebe <- gammasc/sdmatrscebe
      
      pmatrscebe <- 2*pnorm(abs(tmatrscebe),lower.tail = F)
      
      
      rst.rmscebe <- cbind(array(gammasc),array(sdmatrscebe),array(tmatrscebe),array(pmatrscebe))
      colnames(rst.rmscebe)=c("est","sd","t-value","p-value")
      
      
      
      
    }
    

    
    ####scebe
    list_sigmai <- vector(mode = 'list',length = N)
    
    for (i in 1:N){
      Gi <- diag(sigma(fit)^2,nt[i])
      Zi <- Z[((1+sum(nt[-(i:N)])):sum(nt[1:i])),]
      list_sigmai[[i]] <- Zi%*%Rhat%*%t(Zi)+Gi
    }
    
    SIG <- as.matrix(bdiag(list_sigmai))
    
    k_xc_ip <- kronecker(solve(t(xc)%*%xc)%*%t(xc),Ip)
    vvecg <- k_xc_ip%*%DD
    vvecg <- vvecg%*%SIG
    vvecg <- vvecg%*%t(DD)
    vvecg <- vvecg%*%t(k_xc_ip)
    
    Vhatscebe <- solve(SC)%*%vvecg%*%t(solve(SC))
    
    sdmatscebe <- matrix(sqrt(diag(Vhatscebe)),nrow = p,ncol = q)
    
    tmatscebe <- gammasc/sdmatscebe
    
    pmatscebe <- 2*pnorm(abs(tmatscebe),lower.tail = F)
    
    
    rst.mscebe <- cbind(array(gammasc),array(sdmatscebe),array(tmatscebe),array(pmatscebe))
    colnames(rst.mscebe)=c("est","sd","t-value","p-value")
    
    rst.all <- list(rst.nebe,rst.mscebe,rst.rmscebe)
    names(rst.all) <- c("NEBE","mSCEBE","rmSCEBE")
    return(rst.all)
  } 
  
}










