# MuMAC (Multi-model Analysis of Cultural Data)
#
# Carlos A. Botero 
# North Carolina State University
# 12 Aug 2014 (C)

####################################################
MuMAC_distance <- function(mydata, response, Predictors, VerticalDependency = 'VerticalDependency',
                  test=F, mySeed = runif(1)*1E5, Nneighbors = 10) {
  require(nlme)
  require(lme4)
  require(ROCR)
  require(MuMIn) 
  require(spdep)
  require(ape)
  require(geiger)
  require(geosphere)
  
  ##################################################################
  # first some house keeping...
  
  # are we dealing with a binary or a continuous dependent variable?
  proceed <- T
  
  if ( eval(parse(text=paste('!is.null(levels(mydata$', response,'))', sep = ''))) ) { 
    is.binary <- T 
    if (class(VerticalDependency) == "phylo") { 
      warning ( 'Sorry: phylogenetic logistic regression is not yet supported by MuMAC')
      proceed <- F 
    }
  }else { is.binary = F }
  
  if ( proceed ) {    
    # are vertical dependencies specified as a phylogeny, a taxonomy, or not specified?
    if ( !is.na(VerticalDependency)[[1]] ) {
      if (class(VerticalDependency) != "phylo") {
        names( mydata ) [which(names(mydata) == VerticalDependency)] <- 'VerticalDependency'    
      }
      else { mytree <- VerticalDependency }
    }
    
    ###############################################################################
    # evaluate potential for cultural diffusion by contact (here we consider the 
    # potential effect of traditions in neighboring groups on each focal society)
    print('Computing spatial dependencies...')
    
    # score for the 10 nearest neighbors - use ENTIRE dataset
    myDependent <- eval(parse(text=paste('mydata$', response, sep='')))   
    
    # xys <- as.matrix(cbind(mydata$longitude, mydata$latitude))
    # nebs <- knn2nb(knearneigh(xys, k=Nneighbors, longlat=TRUE))
    # nebs = nndist(X = xys[,1], Y = xys[,2], k = 10)
    # NumNeigh <- sapply(nebs, length)
    # Adj <- unlist(nebs)
    
    xys = sp::spDists(mydata[,c("longitude", "latitude")] %>% as.matrix,
                      longlat = TRUE)
    row = xys[1,]
    
    distance = Nneighbors
    
    nebs = apply(xys, 1, function(r) which(r < distance & r > 0))
    neighbours_score = lapply(nebs, function(n) mean(myDependent[n])) %>% 
      unlist(.)
    # assume if no neighbours are within distance, then 0
    neighbours_score = ifelse(is.nan(neighbours_score), 0 , neighbours_score)
    
    N <- dim(mydata)[1]  #number of observations
    # adjIDs=array(NA, dim=c(N,max(NumNeigh)))
    # cellNum=0
    # for (j in 1:N){
    #   for (adjNum in 1:NumNeigh[j]){
    #     cellNum=cellNum+1
    #     adjIDs[j,adjNum] <- Adj[cellNum]
    #   }
    # }
    # 
    # ymat=array(NA, dim=c(N,max(NumNeigh)))
    # 
    # for (j in 1:N){
    #   for (i in 1:NumNeigh[j]){
    #     ymat[j,i] <- as.numeric(myDependent[adjIDs[j,i]])
    #   }
    # }
    
    # compute average score for neighbors
    #mydata$NeighborScore <- rowMeans(ymat,na.rm=T) 
    mydata$NeighborScore <- neighbours_score
    ################################################
    if (test==T) {
      # separate into training and test datasets
      set.seed(floor(mySeed)) # set the seed so the sample is repeatable
      Obs2Use <- sample(1:dim(mydata)[1], round(dim(mydata)[1]*0.66))
      
      mydata.train <- mydata[ Obs2Use, ]
      mydata.test <- mydata[ -Obs2Use, ]
    }
    else { 
      mydata.train <<- mydata 
      mydata.test <<- NULL
    }
    
    # assign to global environment
    mydata.train <<- mydata.train
    mydata.test <<- mydata.test
    
    # add to predictors
    Predictors <- c(Predictors, 'NeighborScore')
    
    ##############################################
    # Run all possible model parameterizations  
    
    # define global model
    GlobalMod <- paste(response, '~ 1')
    
    for (i in 1:length(Predictors)) { GlobalMod <- paste (GlobalMod, '+', Predictors[i] ) }
    
    print('Estimating all possible model combinations...')
    
    if ( !is.na(VerticalDependency)[[1]] ) {
      if ( class(VerticalDependency) != "phylo" ) {
        if (is.binary) {
          mymod.VD <- glmer( formula(paste(GlobalMod, '+ (1|VerticalDependency)')), data = mydata.train, family=binomial, na.action = 'na.pass')
        }
        else {
          mymod.VD <- lmer( formula(paste(GlobalMod, '+ (1|VerticalDependency)')), data = mydata.train, REML=F, na.action = 'na.pass')
        }
      }
      else {
        if (is.binary) {
          # not supported yet
        }
        else {
          ThisTree <<- treedata(mytree, mydata.train, warning=F)$phy
          eval(parse(text=paste("mymod.VD <- gls(", GlobalMod, ", correlation=corBrownian(phy=ThisTree), data=mydata.train, na.action = 'na.pass', method = 'ML')")))  
        }
      }
      myRes.VD <- dredge(mymod.VD,subset=TRUE)
    }
    
    # now models without language
    if (is.binary) {
      mymod.NoVD <- glm(formula(GlobalMod), data = mydata.train, family=binomial, na.action = 'na.pass' )
    }
    else {
      if ( class(VerticalDependency) != "phylo" ) {
        mymod.NoVD <- glm(formula(GlobalMod), data = mydata.train, family=gaussian, na.action = 'na.pass' )        
      }
      else {
        eval(parse(text=paste("mymod.NoVD <- gls(", GlobalMod,", data=mydata.train, na.action = 'na.pass', method = 'ML')")))  
      }
    }
    myRes.noVD <- dredge(mymod.NoVD,subset=TRUE)
    
    #######################################################
    # extract AICc values and compute delta AIC and Akaike weights
    if ( !is.na(VerticalDependency)[[1]] ) { aics<-c(myRes.VD$AICc, myRes.noVD$AICc) }
    else { aics <- myRes.noVD$AICc }
    
    daic<-aics-min(aics)
    waic<-exp(-.5*daic)/sum(exp(-.5*daic))
    
    #############################################################
    # Now compute the average model based on Akaike weights
    
    print('[Please be patient]')
    mod.NoVD <- get.models(myRes.noVD, subset=TRUE)
    if ( !is.na(VerticalDependency)[[1]] ) { 
      mod.VD <- get.models(myRes.VD, subset=TRUE)
      avgmod <- model.avg(c(mod.VD,mod.NoVD))
    }
    else { avgmod <- model.avg(mod.NoVD) }
    
    avgmod_summary = summary(avgmod)
    
    MM_average <- as.data.frame(summary(avgmod)$coefmat.full[,c('Estimate','Std. Error')])
    
    MM_average$RVI <- NA
    MM_average$PredictiveAccuracy <- NA
    
    # compute relative variable importance (by adding weights of all models in which 
    # a given predictor appears)
    if ( !is.na(VerticalDependency)[[1]] ) { 
      for (i in Predictors) {
        MM_average[i,'RVI'] <- sum(waic[eval(parse(text=paste("which(!is.na(c(myRes.VD$'", i,"', myRes.noVD$'", i, "') ))", sep=''))) ])
      }  
      MM_average['(Intercept)','RVI'] <- sum(waic[eval(parse(text=paste("which(!is.na(c(myRes.VD$'(Intercept)', myRes.noVD$'(Intercept)') ))", sep=''))) ])
      MM_average['VerticalDependency','RVI'] <- sum(waic[1:dim(myRes.VD)[1]])
    }
    else {
      for (i in Predictors) {
        MM_average[i,'RVI'] <- sum(waic[eval(parse(text=paste("which(!is.na( myRes.noVD$'", i, "' ))", sep=''))) ])
      }  
      MM_average['(Intercept)','RVI'] <- sum(waic[eval(parse(text=paste("which(!is.na( myRes.noVD$'(Intercept)') )", sep=''))) ])
    }
    
    #########################################################
    # output table with model weights
    if ( !is.na(VerticalDependency)[[1]] ) { ModWeights <- as.data.frame(matrix(NA, dim(myRes.VD)[1]*2, 3)) }
    else { ModWeights <- as.data.frame(matrix(NA, dim(myRes.noVD)[1], 3)) }
    names(ModWeights) <- c('Model', 'delta AICc', 'AICc weight')
    
    if ( !is.na(VerticalDependency)[[1]] ) {
      # first output the model with language as a random effect
      for (i in 1:dim(myRes.VD)[1]) {
        ThisModelPredictors <- NULL
        
        for (j in 1:dim(myRes.VD)[2]) {
          if ( !is.na(myRes.VD[i,j]) & names(myRes.VD)[j] %in% c(Predictors, '(Intercept)')) { 
            ThisModelPredictors <- paste(ThisModelPredictors, names(myRes.VD)[j]) 
          }
        }
        
        ModWeights$Model[i] <- paste(ThisModelPredictors, 'VerticalDependency')
        ModWeights$'delta AICc'[i] <- daic[i]
        ModWeights$'AICc weight'[i] <- waic[i]
      }
      
      # now output the same model without language
      for (i in 1:dim(myRes.noVD)[1]) {
        ThisModelPredictors <- NULL
        
        for (j in 1:dim(myRes.noVD)[2]) {
          if ( !is.na(myRes.noVD[i,j]) & names(myRes.noVD)[j] %in% c(Predictors, '(Intercept)') ) { 
            ThisModelPredictors <- paste(ThisModelPredictors, names(myRes.noVD)[j]) 
          }
        }
        
        ModWeights$Model[i + dim(myRes.VD)[1]] <- ThisModelPredictors
        ModWeights$'delta AICc'[i + dim(myRes.VD)[1]] <- daic[i + dim(myRes.VD)[1]]
        ModWeights$'AICc weight'[i + dim(myRes.VD)[1]] <- waic[i + dim(myRes.VD)[1]]
      }
    }
    else {
      # ONLY output models without language
      for (i in 1:dim(myRes.noVD)[1]) {
        ThisModelPredictors <- NULL
        
        for (j in 1:dim(myRes.noVD)[2]) {
          if ( !is.na(myRes.noVD[i,j]) & names(myRes.noVD)[j] %in% c(Predictors, '(Intercept)') ) { 
            ThisModelPredictors <- paste(ThisModelPredictors, names(myRes.noVD)[j]) 
          }
        }
        
        ModWeights$Model[i] <- ThisModelPredictors
        ModWeights$'delta AICc'[i] <- daic[i]
        ModWeights$'AICc weight'[i] <- waic[i]
      }
    }
    
    #########################################################
    if ( test ) {
      # compute Predictive accuracy for average model
      preds=predict(avgmod, mydata.test, type="response")
      
      if (is.binary) {
        eval(parse(text=paste("MM_avg_AUC <- performance(prediction(preds, mydata.test$'", response,"'), 'auc')@y.values[[1]]", sep='')))
      }
      else {
        MM_avg_Rsq <- eval(parse(text=paste('cor.test(preds, mydata.test$', response, ')', sep='')))$estimate[[1]]
      }
      #########################################################
      # compute predictive accuracies for each predictor
      # first the intercept 
      if (is.binary) {
        # (check: this one should be similar to 0.5)
        mymod <- glm(formula(paste(response, '~ 1', sep='')), data = mydata.test, family=binomial)
        preds=predict(mymod, mydata.test, type="response")
        eval(parse(text=paste("MM_average['(Intercept)','PredictiveAccuracy'] <- performance(prediction(preds, mydata.test$'", response,"'), 'auc')@y.values", sep='')))
      }
      else {
        # I'm not computing this because the predictions with just an intercept provide
        # a constant output and you can't compute a correlation under those conditions (as the SD == 0)
#         mymod <- glm(formula(paste(response, '~ 1', sep='')), data = mydata.test, family=gaussian)
#         preds=predict(mymod, mydata.test, type="response")
#         MM_average['(Intercept)','PredictiveAccuracy'] <- eval(parse(text=paste('cor.test(preds, mydata.test$', response, ')', sep='')))$estimate[[1]]
      }
      
      # now the predictors
      for (i in Predictors) {
        i <- gsub('[*]', ':', i)
        if (is.binary) {
          mymod <- glm(formula(paste(response, '~ 1 +', i, sep='')), data = mydata.test, family=binomial)
          preds=predict(mymod, mydata.test, type="response")
          eval(parse(text=paste("MM_average['", i,"','PredictiveAccuracy'] <- performance(prediction(preds, mydata.test$'", response,"'), 'auc')@y.values", sep='')))
        }
        else {
          mymod <- glm(formula(paste(response, '~ 1 +', i, sep='')), data = mydata.test, family=gaussian)
          preds=predict(mymod, mydata.test, type="response")
          MM_average[i,'PredictiveAccuracy'] <- eval(parse(text=paste('cor.test(preds, mydata.test$', response, ')', sep='')))$estimate[[1]]
        }    
      }
      
      # and finally... the vertical dependency
      if ( !is.na(VerticalDependency)[[1]] ) {
        if ( class(VerticalDependency) != "phylo" ) {
          if (is.binary) {
            mymod <- glmer(formula(paste(response, '~ 1 + (1|VerticalDependency)', sep='')), data = mydata.test, family=binomial)
            preds=predict(mymod, mydata.test, type="response")
            eval(parse(text=paste("MM_average['VerticalDependency','PredictiveAccuracy'] <- performance(prediction(preds, mydata.test$'", response,"'), 'auc')@y.values", sep='')))
          }
          else {
            mymod <- lmer(formula(paste(response, '~ 1 + (1|VerticalDependency)', sep='')), data = mydata.test, REML = F)
            preds=predict(mymod, mydata.test, type="response")
            MM_average['VerticalDependency','PredictiveAccuracy'] <- eval(parse(text=paste('cor.test(preds, mydata.test$', response, ')', sep='')))$estimate[[1]]
          }
        }
        else {
          # Same here: Can't compute a correlation when one variable is a constant
#           ThisTree <<- treedata(mytree,mydata.test, warning=F)$phy
#           mymod <- gls(formula(paste(response, '~ 1', sep='')), correlation=corBrownian(phy=ThisTree), data=mydata.test, na.action = 'na.pass', method = 'ML')  
#           preds=predict(mymod, mydata.test, type="response")
#           MM_average['VerticalDependency','PredictiveAccuracy'] <- eval(parse(text=paste('cor.test(preds, mydata.test$', response, ')', sep='')))$estimate[[1]]
        }
      }
    }
    
    if (class(VerticalDependency) == "phylo") { rm(list = "ThisTree", pos = ".GlobalEnv") }
    #########################################################
    # output results
    myOutput <- list()
    myOutput$ModWeights <- ModWeights[sort(ModWeights$'delta AICc',index.return=T)$ix,]
    myOutput$MM_average <- MM_average[sort(row.names(MM_average),index.return=T)$ix,]
    
    if ( test ) {
      if (is.binary) { myOutput$MM_avg_AUC <- MM_avg_AUC }
      else { myOutput$MM_avg_Rsq <- MM_avg_Rsq }
    }
    
    rm(list = c("mydata.train", "mydata.test"), pos = ".GlobalEnv")
    
    return(myOutput)
  }
}