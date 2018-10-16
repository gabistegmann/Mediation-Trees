MedTrees <- function(totalXY ,
                     bc_effect,
                     a_effect,
                        data,
                        rPartFormula,
                        group= ~ id,
                        control = rpart.control()){
  
  groupingName = attr(terms(nlme::splitFormula(group,'~')[[1]]),"term.labels")
  groupingFactor = data[,names(data)==groupingName]
  
  #    terms = attr(terms(lmeFormula),"term.labels")
  # continuous = !is.factor(data[,names(data)==terms[1]])
  
  ### The 3 subfunctions necessary for rpart to work.
  # The evaluation function.
  # Called once per node:
  # returns a list of two variables: a label for the node
  # and a deviance value for the node.  The deviance is
  # of length one, equal to 0 if the node is perfect/ unsplittable
  # larger equals worse
  evaluation <- function(y, wt, parms){
    
    
    dataNode=parms[groupingFactor%in%y,]
    fit = lm(totalXY, data = dataNode)
    fitt=1
    list(label=fitt,deviance=sum((fit$residuals)^2))
    
   
  }
  
  
  # The split function, where the work occurs.  This is used to decide
  # on the optimal split for a given covariate.
  # Called once per split candidate per node
  ### If the covariate, x, is continuous:
  # x variable is ordered
  # y is provided in the sort order of x
  # returns two vectors of length (n-1)
  #      goodness: goodness of the split, with larger numbers better. 0=no split
  #      direction: -1 = send y < cutpoint to the left
  #                  1 = send y < cutpoint to the right
  #
  ### If x is non-continuous
  # x is a set of integers (NOT the original values) defining the
  # groups for an unordered predictor.
  # Again, return goodness and direction
  #      direction: a vector of length m (number of groups), which is the applied
  #                 ordering for the unordered categories.  This is done so that
  #                 m-1 splits are performed, instead of all possible splits
  #      goodness: m-1 values, same idea as before
  
  ### pass in the dataset through the parms variable, with subj as y
  split <- function(y, wt, x, parms, continuous) {
    print(paste("splitting:", length(unique(x)), "values"))
    dev = vector()
    xUnique = unique(x)
    
    #  rootDev = lme(lmeFormula, data = parms[groupingFactor %in%
    #                                           y, ], random = randomFormula, correlation = R,
    #                na.action = na.omit)$logLik
   dataNode=parms[groupingFactor%in%y,]
    fit = lm(totalXY, data = dataNode)
    rootDev = -sum((fit$residuals)^2)
    

    if (continuous) {
      for (i in xUnique) {
        yLeft = y[x <= i]
        yRight = y[x > i]
        if (length(yLeft) < control$minbucket || length(yRight) <
            control$minbucket) {
          dev = c(dev, 0)
        }
        else {
          dataL = parms[groupingFactor %in%
                                 yLeft, ]
          modelLeft = lm(totalXY, data = dataL)
          dataR = parms[groupingFactor %in%
                                 yRight, ]
          modelRight = lm(totalXY, data = dataR)
          
          
          dev = c(dev, -sum((modelLeft$residuals)^2, (modelRight$residuals)^2))
          
        }
      }
      good = rep(0, length(x))
      for (i in 1:length(xUnique)) {
        good[x == xUnique[i]] = dev[i]
      }
      good = good[1:(length(good) - 1)]
      list(goodness = good + abs(rootDev) * (good != 0) *
             2, direction = rep(-1, length(good)))
    }
    else {
      order = rep(0, length(xUnique))
      # response = parms[, names(parms) == responseName]  #  -          ----- problem
      # for (i in 1:length(xUnique)) {
      #   order[i] = mean(response[x == xUnique[i]], na.rm = TRUE)
      # }
      dir = sort(order, index.return = TRUE)$ix
      for (i in 1:(length(dir) - 1)) {
        yLeft = y[x %in% dir[1:i]]
        yRight = y[x %in% dir[(i + 1):length(dir)]]
        if (length(yLeft) < control$minbucket || length(yRight) <
            control$minbucket) {
          dev = c(dev, 0)
        }
        else {
          
          dataL = parms[groupingFactor %in%
                          yLeft, ]
          modelLeft = lm(totalXY, data = dataL)
          dataR = parms[groupingFactor %in%
                          yRight, ]
          modelRight = lm(totalXY, data = dataR)
          
          
          dev = c(dev, -sum((modelLeft$residuals)^2, (modelRight$residuals)^2)) 
        }
      }
      list(goodness = dev + abs(rootDev) * (dev != 0) *
             2, direction = dir)
    }
  }
  # The init function.  This is used, to the best of my knowledge, to initialize the process.
  # summary is used to fill print the report summary(model), and text is used to add text to
  # the plot of the tree.
  initialize <- function(y,offset,parms=0,wt){
    
  
    list(
      y=y,
      parms=parms,
      numresp=1,
      numy=1,
      summary=function(yval,dev,wt,ylevel,digits){paste("deviance (-2logLik)",format(signif(dev),3),"fitt",signif(yval,2))},
      text= function(yval,dev,wt,ylevel,digits,n,use.n){
        if(!use.n){paste("m:",format(signif(yval,1)))}
        else{paste("n:",n)}
      }
    )
    
    
  }
  model <- list()
  model.rpart = rpart(paste(groupingName,c(rPartFormula)),
                      method=list(eval=evaluation,
                      split=split,
                      init=initialize),
                      control=control,
                      data=data,
                      parms=data)
  model$rpart_out <- model.rpart
  
  
  model$leaf_node <- model.rpart$where
  summary_XM = summary_XY = summary_XMY =  list()
  for (j in 1:length(table(model.rpart$where))) {
    id <- names(table(model.rpart$where))[j] == model.rpart$where
    
    
      model_XM <- lm(a_effect, data = data[id, ])
      s_XM = summary(model_XM)
      
      summary_XM[[as.numeric(names(table(model.rpart$where)))[j]]] <- s_XM
      
      model_XY <- lm(totalXY, data = data[id, ])
      s_XY = summary(model_XY)
      
      summary_XY[[as.numeric(names(table(model.rpart$where)))[j]]] <- s_XY
      
      model_XMY <- lm(bc_effect, data = data[id, ])
      s_XMY = summary(model_XMY)
      
      summary_XMY[[as.numeric(names(table(model.rpart$where)))[j]]] <- s_XMY
      

  }
  
  model$summary_XM <- summary_XM
  model$summary_XY <- summary_XY
  model$summary_XMY <- summary_XMY
  
  model$data = data
  model$groupingName = groupingName
  model$rPartFormula = rPartFormula
  model$call <- match.call()
  class(model) <- "MedTree"
  
  model
  }
 