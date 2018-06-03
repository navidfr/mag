###### GENERIC Clustering code


###############################
###############################
############################### SIMULATION and secondary   SHAPE    ADD ID    MAG.VAR                           

cluster.sim<-function(edep, efrq, shape, p=0.1, n=10000){
  
  edep.amp <- edep/p
  lambda.alt <- efrq * edep.amp
  lambda.ref <- edep.amp - lambda.alt
  
  scale.ref <- lambda.ref/shape
  scale.alt <- lambda.alt/shape  
  
  #number of copies amplified from each unique molecule       K=0.005 => shape=200
  #keep shape=200,  lambda= shape*scale       Shape=K=1/k    scale=lambda*k    scale=lambda/shape  dispersion= k=0.005     gamma(1/k,lambda*k)
  ref.amp <- rgamma(n, shape=shape, scale=scale.ref)             #lambda=90       Total # of reads    
  alt.amp <- rgamma(n, shape=shape, scale=scale.alt)        #lambda=10        Total # of altreads.1 (for 1/10 freq)
  
  ref.cnt <- sapply(ref.amp*p, rpois, n=1)
  alt.cnt <- sapply(alt.amp*p, rpois, n=1)
  dep <- alt.cnt + ref.cnt   
  frq <- alt.cnt/dep  
  
  clust.out<-cbind(ref.cnt, alt.cnt, dep, frq)
  return(clust.out)
}

get.var <-function(x, x.nb, s, weight=1) {
  a <- mean(x.nb[as.logical(x)]);
  b <- var(x.nb[as.logical(x)])*weight;
  c <- sum(as.logical(x));
  d <- s;
  return(c(a, b, c, d));
}

global.s1 <<- 1.04;
global.s2 <<- 1.05;

######find shape with this:
#normal is the germline data
fitdistr(x = Normal.data.MC231$GB.ref.dp, densfun = "Gamma")
library(MASS)


#####feb 28 update: 
#add id





#MAY 20th
#does all the removing of 0 and 1 s


#new functions


mbhc.prepdata<- function(arg.data){
  
  ### arg data should have the frequencies in rows and columns for samples. for 3 samples, 100 mutation: 3 columns,100rows. 
  arg.data <- as.matrix(arg.data)
  arg.data <- round(arg.data,3)
  arg.data <- ifelse(arg.data<1e-3, 1e-3, arg.data)
  arg.data <- ifelse(arg.data> 1-(1e-3), 1-(1e-3), arg.data)
  
  prep.data<- data.frame(arg.data, row.names=NULL)
  prep.data$ID<- c(1:dim(prep.data)[1])
  return(prep.data)
}




###############################
###############################
############################### OPTIM                                           


### may 16th   range and variance
optim.obj <- function(par1,exp.depth=exp.depth, y, weight=-0.0001){     
  
  R <-dbeta(y, shape1=par1, shape2=exp.depth-par1, log=T) #log should be F- because if Log=T we could have positive and negative values. We should add log
  #print(R)
  #R[which(R < 1e-10)] <- 0;
  #if(length(which(R == 0)) > 0) {
  #return(1000);
  #}
  # print(R)
  #	print(paste(y, shapes))
  # print(mean(R))
  v <- sd(y)
  range<- max(y)-min(y)
  #print(v)
  #print('v')
  #print(v)
  
  #print('range')
  #print(range)
  
  # print('v,meanR,logR, expWei')
  v <- ifelse(v<0.005, 0.005, v)
  range<- ifelse(range<0.005,0.005,range)
  
  #R <- R/v
  R <- R/(range*v)
  #	print(mean(R))
  #R <- log(R)
  # print('mean')
  #  print(mean(R))
  #print(exp(5*(weight+0.0001)))
  results <- -mean(R)/exp(5*(weight+0.0001))
  #print(paste(results, 'end'))
  return(results)  
}  #MAY 6th version. Uses range- only one par.




###############################
###############################
############################### FITS   AND WEIGHTS                                    





#####MAY 6th version
fit.variants.weight <- function(variants){   #to calculate the initial matrix
  weight <- variants[length(variants)]
  variants <- variants[-length(variants)]
  #fit<-optim(par1, y=variants, fn=optim.obj, weight=weight,method="Brent",lower=1.5, upper=1000)#,  control=list(factr=1e-3, maxit=0, trace=0)) #, lower=1.5, upper=1000, method='L-BFGS-B'
  fit1<-optimize(y=variants,weight=weight,exp.depth=exp.depth ,f=optim.obj, interval = c(1,199),tol=0.001)#,lower=1, upper=200, control=list(factr=1e-3, maxit=200, trace=0)) #, lower=1.5, upper=1000, method='L-BFGS-B'
  #fit2<- optim(par=par1,y=variants, weight=weight, fn=optim.obj)
  #fit3<- DEoptim(optim.obj,exp.depth=exp.depth ,y=variants, weight=weight, 
  #              lower = 1,upper=200, DEoptim.control(trace=F, reltol=10e-1, itermax=15) )
  
  #global.s1 <- fit$par[1]
  #global.s2 <- fit$par[2]
  #print('vars')
  #print(variants)
  #return(rbind(c(fit1$objective, fit1$minimum),c(fit2$value, fit2$par),c(fit3$optim$bestval, fit3$optim$bestmem)))
  return(c(fit1$object,fit1$minimum))#,c(fit3$optim$bestval, fit3$optim$bestmem)))
}



fit.variants <- function(variants){   #to calculate the initial matrix
  return(fit.variants.weight(c(variants, -0.0001)));
}

fit.elements <- function(element.1, element.2, arg.data){   
  data.1 <- arg.data[as.logical(element.1)]
  data.2 <- arg.data[as.logical(element.2)]
  data.merge <- c(data.1, data.2)
  fit <- fit.variants(data.merge)
  return(fit)
}




### cleaned up MAY 20 ---- changed ss to arg.data
fit.elements.weight <- function(element.1, element.2, arg.data, arg.data.other){   
  data.1 <- arg.data[as.logical(element.1)]
  data.2 <- arg.data[as.logical(element.2)]
  weight <- get.weights.within(element.1 + element.2, arg.data.other)
  data.merge <- c(data.1, data.2, 10*weight)
  fit <- fit.variants.weight(data.merge)
  print(data.merge)
  return(fit)
}

fit.elements.weight.reduce <- function(element.1, element.2, arg.data, arg.data.other){   
  data.1 <- arg.data[as.logical(element.1)]
  data.2 <- arg.data[as.logical(element.2)]
  weight <- get.weights.within(element.1 + element.2, arg.data.other)
  data.merge <- c(data.1, data.2, 10*weight)
  fit <- fit.variants.weight(data.merge)
  #print(data.merge)
  return(fit)
}









##### CLEANED UP: MAY 20th
get.weights.within <- function(x1, arg.data) {   #arg.data should be the remaining samples
  weights <- c()
  #print(arg.data)
  if(sum(x1)>1){  
    weights <- apply(arg.data, 2, function(x){ return(max(dist(x[as.logical(x1)]))) })
  } else{
    weights <- 0
  }  
  #print(weis)
  return(max(weights))
}

# x1 is matrix of all remaining clusters
# x2 is the last cluster.
# arg.data should be the remaining samples
get.weights.between <- function(x1, x2, arg.data) {  
  xx <- x1;
  if(is.vector(x1)) {
    xx <- t(as.data.frame(x1));
  }
  weights <- apply(xx, 1, function(x) { x <- x + x2; return(get.weights.within(x, arg.data=arg.data)) } );
  return(weights)
}







###############################
###############################
############################### MAG SINGLE                                     



#### NEW OPTIM may 16th
mbhc.single <- function(arg.data, sampleNum=1){     
  if(is.null(arg.data$ID)){ return(print('Run mag.prep on the data'))}
  time0 <- as.numeric(Sys.time())
  arg.data<- round(arg.data,3)
  arg.data[,-dim(arg.data)[2]]<- ifelse(arg.data[,-dim(arg.data)[2]] < 1e-3, 1e-3, arg.data[,-dim(arg.data)[2]])
  arg.data[,-dim(arg.data)[2]]<- ifelse(arg.data[,-dim(arg.data)[2]] > 1-1e-3, 1-1e-3, arg.data[,-dim(arg.data)[2]])
  
  x.nb <- arg.data;
  if(!is.null(ncol(arg.data)) && ncol(arg.data) > 0) {
    x.nb <- arg.data[,sampleNum]
  }
  params<-c();  x.values<-c();  x.steps<-c();  row.loop<-c(); x.step<-c()
  x.elements.list<-list()
  x.elements<-data.frame(diag(length(x.nb)))
  x.el.first<- data.frame(diag(length(x.nb)))
  
  x.el.first.okay<- apply(x.el.first, 1, function(x) paste(x, collapse = '')) # april 15
  
  sort.x<-sort(arg.data[,sampleNum])
  #sort.x<- ifelse(sort.x < 1e-9, 1e-9, sort.x)
  #sort.x<- ifelse(sort.x > 1-1e-4, 1-1e-3, sort.x)
  
  
  #x.loop<-sort.x   
  sort.arg.data<- arg.data[order(arg.data[,sampleNum]),]
  
  names(x.elements)<-sort.arg.data$ID
  names(x.el.first)<- sort.arg.data$ID
  
  l<-length(sort.x)
  mat.loop<-matrix(10^5, ncol=l, nrow=l)
  par.1.loop<-matrix(0, ncol=l, nrow=l)
  #par.2.loop<-matrix(0, ncol=l, nrow=l)
  
  s.pairs <- cbind(sort.x[1:(l-1)], sort.x[2:l])  
  #print(s.pairs)
  s.likls <- t(apply(s.pairs, 1, fit.variants))  #calculated likelihoods on the pairs
  
  #mat.loop<-diag(0, nrow=l, ncol=l)
  for(i in 2:l){
    mat.loop[i-1,i] <- s.likls[i-1, 1]      ### the s.likls is the likelihood between i-1)th element and the i th. 
    par.1.loop[i-1,i] <- s.likls[i-1, 2]
    # par.2.loop[i-1,i] <- s.likls[i-1, 3]
  }  
  
  time1 <- as.numeric(Sys.time())
  var.mean<-c() 
  s<-1
  while(nrow(x.elements)>= 2){
    # print("INJA")
    # print(s)
    #print(mat.loop)
    
    #print(min(mat.loop))
    x.values[s]<-min(mat.loop)
    el1<-which(mat.loop==min(mat.loop), arr.ind = T)[1,1]
    el2<-which(mat.loop==min(mat.loop), arr.ind = T)[1,2]
    el.rm<-c(el1,el2)
    # print("Remove ina")
    # print(el.rm)
    params<-rbind(params,par.1.loop[el1,el2])
    
    #what happenes at each step
    x.step <- x.elements[el.rm[1],] + x.elements[el.rm[2],]   ###this is very clever! I dont need to update the xloop.
    ###only use x.elements rows for each cluster!
    x.steps<-rbind(x.steps, x.step)
    #print(x.steps)   
    x.elements<-rbind(x.elements[-el.rm,], x.step)
    x.elements.list[[s]]<-x.elements
    
    
    #### april 14 add "okay"
    temp<- t(apply(x.elements, 1, get.var, x.nb=sort.x, s=s ))          ####for oc.v9 I changed x.nb=x.nb to x.nb=
    #####
    #okay<- apply(x.elements, 1, function(x) nrow(merge(t(x), x.el.first)))               ################# APRIL 14 2018   takeeeesss a long time
    ####################### faster version? :
    okay<- apply(x.elements, 1, function(x) sum(paste(x, collapse = '')%in%x.el.first.okay))
    weights<- rep(0, length(okay))  ##### just to have weight  APRIL 15
    
    temp<- cbind(temp,weights ,okay)
    
    var.mean<- rbind(var.mean, temp)
    
    
    
    
    #####UPDATE MATRIX HERE    
    mat.loop<-mat.loop[-el.rm, -el.rm]      
    par.1.loop<-par.1.loop[-el.rm, -el.rm]
    #par.2.loop<-par.2.loop[-el.rm, -el.rm]    
    
    mat.column.update<-c()
    par1.column.update<-c()
    #par2.column.update<-c()
    
    max<-max(which(x.step%in%1))
    min<-min(which(x.step%in%1))
    #		print(paste("min", min))
    nei<-c()       
    if(min == 1 & max != l){
      neiIND<-which(x.elements[,max+1]==1)
      nei<-x.elements[neiIND,]
      #			print("if1")
    }else if(max == l & min != 1){
      neiIND<-which(x.elements[,min-1]==1)
      nei<-x.elements[neiIND,]
      #			print("if2")
    }else if( max!=l & min!=1){
      nei1<-which(x.elements[,min-1]==1)
      nei2<-which(x.elements[,max+1]==1)
      neiIND<-c(nei1,nei2)
      nei<-x.elements[neiIND,]
      #			print("if3")
    }
    
    if(length(nei) > 0){
      x.last <- x.elements[nrow(x.elements), ]
      #print('sortx')
      #print(sort.x)
      #print(x.last)
      #print(nei)
      #print('done')
      temp <- apply(nei, 1, fit.elements, arg.data=sort.x, element.2=x.last)
      
      ###try to only look at cluster neighbours 
      ###instead of dim i will use sqrt of length. Because mat.loop is always square matrix   WORKS!!! SEP12
      l2<-sqrt(length(mat.loop))
      mat.column.update<-rep(10^5,l2)
      mat.column.update[neiIND]<-temp[1,]
      par1.column.update<-rep(0, l2)
      par1.column.update[neiIND]<-temp[2,]
      #par2.column.update<-rep(0, l2)
      #par2.column.update[neiIND]<-temp[3,]
      
      mat.loop<-cbind(mat.loop, mat.column.update )
      mat.loop<-rbind(mat.loop, rep(10^5,dim(mat.loop)[2]))    #just to keep mat.loop square. 
      par.1.loop<-cbind(par.1.loop, par1.column.update)
      par.1.loop<-rbind(par.1.loop, rep(0, dim(par.1.loop)[2]))
      #par.2.loop<-cbind(par.2.loop, par2.column.update)
      #par.2.loop<-rbind(par.2.loop, rep(0, dim(par.2.loop)[2]))
    }
    
    s<-s+1 
  }
  
  #var.mean <-cbind(var.mean, 1/(var.mean[,3]))
  colnames(var.mean)<-c("mean", "var","NumberOfPoints", "step", "weights", 'okay')
  
  time2 <- as.numeric(Sys.time())
  cat("while took ", (time2 - time1), " seconds.\n")
  cat("total took ", (time2 - time0), " seconds.\n")
  
  result <- list(x.el=x.elements.list, params=params, var.mean=var.mean, data=arg.data, data.sorted=sort.arg.data,freq.s= sort.x, time=list(time0,time1, time2))
  return(result)
}



###############################
###############################
############################### find.cut.off                                  


### april 9. works with added ids to mag single.                     ### modified the code so the mag.single does.not output data.frame
mbhc.find.cut.off <- function(mag.output){
  library(dplyr)
  time0 <- as.numeric(Sys.time())
  var.mean <- mag.output$var.mean
  var.mean <- data.frame(var.mean, row.names = NULL)    
  var.mean[is.na(var.mean)] <- 0
  var.mean$var <- round(var.mean$var, 4)
  vars<- var.mean$var[var.mean$var>0]
  var.bound<- quantile(vars)[4]*1.5
  
  
  flag <- FALSE
  m <- c()
  str <- c()  #clone structure line from temp var.mean matrix
  str.f <- c()
  cutstr <- c()
  cutstr2 <- c()  
  
  
  i <- max(var.mean$step)
  allPoints <- mag.output$x.el[[i]]
  ok.points <- allPoints - allPoints  ##### To generate clean slate
  n <- length(allPoints)
  #	print(allPoints)
  id.ok <-c() #clusters with okay variance! 
  temp.mv <- c();
  temp.points <- c();
  temp.step <- c();
  temp.point.step <- c();
  while(flag == FALSE){
    #print(flag)
    print(i); flush.console();
    v <- i
    
    ############    
    ## I need to use x.el to get the breaks, because using StepI would cause problem 
    ## in identical clusters and mutations.
    ## TO FIX THIS I JUST KEPT MY ORIGINAL METHOD - 
    ## I just duplicated the rows in the case of identical mutations
    ## because there should be onlytwo rows at each step 
    ## and having only 1 row means identical clusters. 
    ############
    stepI <- var.mean[var.mean$step==i,]
    stepI.x.el <- mag.output$x.el[[i]] 
    stepI.x.el.step <- stepI.x.el[nrow(stepI.x.el), ]  #Only look at the most recent cluster
    
    stepII <- var.mean[var.mean$step==i-1, ]  ## previous step
    #stepII$var<-stepII$var,)
    stepII.x.el <- mag.output$x.el[[i-1]]
    print("II")
    #print(stepII.x.el)
    
    #if( !(id.sum%in%id.ok)){    ###this means if the division on this step is not happening on the "ok" clusters. so the okay clusters should still exist in the temp. 
    #if(sum(abs(id.ok-id.sum)<0.1)==0){
    #newtemp<-temp[!(temp$id%in%id.ok),]
    #if(!id.I%in%id.ok){ 
    temp.points <- rbind(temp.points, c(i, ok.points));
    temp.step <- rbind(temp.step, c(i, stepI.x.el.step));
    temp.point.step <- rbind(temp.point.step, c(i, sum(ok.points & stepI.x.el.step)))
    
    
    print(ok.points)
    print(sum(ok.points & stepI.x.el.step))
    if(sum(ok.points & stepI.x.el.step) == 0){     ###meaning stepI.x.el is not in ok.points 
      
      stepII.breaks <- setdiff(stepII[,-c(4,6)], stepI[,-c(4,6)])
      if(nrow(stepII.breaks)==1){
        print("why did this happen????")
        stepII.breaks <- rbind(stepII.breaks, stepII.breaks)
      }			
      stepII.num <- stepII[, 4]
      stepII.breaks <- cbind(stepII.breaks, stepII[1:2,c(4,6)])   #just add the step number back to the matrix
      stepII.breaks.x.el <- setdiff(stepII.x.el, stepI.x.el)      
           for(ii in 1:2){  
         if(stepII.breaks[ii, 'okay'] ==1 | stepII.breaks$var[ii]<var.bound) {        #okay cluster   APRIL 15
            
          str.f <- c(str.f, stepII.breaks[ii,1])
          #print("here")
          #print(str)
          str <- rbind(str, stepII.breaks[ii,])
          print(stepII.breaks.x.el)
          cutstr2 <- rbind(cutstr2, stepII.breaks.x.el[ii,])
          #print(str)
          
          #id.ok<-c(id.ok, stepII.breaks$id[ii])
          ok.points <- ok.points + stepII.breaks.x.el[ii,]
          #print(ok.points)
        }        
      }
    }
    
    #print("LAST PRINT")
    if(sum(ok.points>1)>0){
      print("ERROR")
    }
    if(sum(ok.points) == n){
      flag<-TRUE
    } else {
      i<-(v-1)
    }  
    #if(v==i){iii<-TRUE}else{i<-(v-1)}
  }
  
  for( i in 1:nrow(str)){
    #print(i)
    xel <- mag.output$x.el[[str$step[i]]]
    #print(xel)
    for(j in 1:nrow(xel)){
      m1 <- mean(mag.output$freq.s[xel[j,]==1])    ##changed to sorted. after oc.9#### OCT 30: THis should be changed for OC11 since there is no sorting in that
      if(abs(m1-str[i,1]) < 0.001){
        cutstr<-rbind(cutstr, xel[j,])  ## what's the difference between cutstr and cutstr2????????
      }
    }
  }
  
  
  time1 <- as.numeric(Sys.time())
  cat("took ", (time1 - time0), " seconds.\n");
  #print(cutstr)
  colors<-colSums(cutstr*c(1:nrow(cutstr)))+1
  final.data<-cbind(mag.output$data.sorted,colors)
  results<-list(str=str, cutclust=cutstr, cutclust2=cutstr2, colors=colors, final.data=final.data)
  return(results)  
}














###############################
###############################
############################### Reduce                                

#BEFORE FEB 26 2018 version

#####
##### mag.reduce does single clustering on the data then do the reduction. 
##### mag.reduce.fast does the reduction without doing single sample clusterin. We use this one.
##### in the mag.reduce.1 I should change all mag.reduce, to mag.reduce.fast   ##FEB 12








#### factors to check: * optim.bound  * sizethreshold. 
#### 
mbhc.reduce.graph<- function(arg.data,x.el=data.frame() ,sampleNum=1, n=-100){
  #x.el<- reduce.out$x.el.red
  #data<- reduce.out$data
  #data.sorted<- reduce.out$data.sorted
  #IDs<- data.sorted[, dim(data.sorted)[2]]
  
  data<- arg.data
  data.sorted<- arg.data[order(arg.data[,sampleNum]),]
  IDs<- data.sorted[, dim(data.sorted)[2]]
  
  if(dim(x.el)[1]==0){
    l<- dim(data)[1]
    x.el<- data.frame(diag(l))
    names(x.el)<- IDs
  }
  
  
  mean.frqs<-t(apply( x.el, 1, function(x){mean(data.sorted[as.logical(x),sampleNum])}))
  x.el.sorted<-x.el[order(mean.frqs),]
  x.el.final<- data.frame()
  x.nb<- data.sorted[,sampleNum]
  ss.other<- data.sorted[,-c(sampleNum, dim(data.sorted)[2]), drop=F]
  
  size.threshold<-5
  
  l<-nrow(x.el)
  full<-floor(l/size.threshold)
  rem<-l%%size.threshold
  row.fold<-0
  
  if( rem < 5 ){ full<- full-1; rem<- rem+size.threshold}
  
  
  optim.value.bound<-fit.variants.weight(c(0.4,0.46,0.05))[1]
  
  
  
  
  ind.full<- c(1:size.threshold)
  combn.el<- t(combn(ind.full,2))
  
  if(full>0){
    for(i in 1:full){
      
      fold<- c((1+(i-1)*size.threshold):(i*size.threshold))  
      x.el.fold<-x.el.sorted[fold, ] 
      #print(i)
      initial.full<-apply(combn.el, 1, function(x){fit.elements.weight.reduce( x.el.fold[x[1],], 
                                                                               x.el.fold[x[2],],
                                                                               x.nb, ss.other )})
      tempMatLoop<-matrix(1000, nrow=size.threshold, ncol=size.threshold)
      tempMatLoop.2<-lower.tri(diag(0, nrow=size.threshold, ncol=size.threshold))
      
      tempMatLoop[tempMatLoop.2] <- initial.full[1, ]
      mat.loop<-t(tempMatLoop)
      mat.loop.g<- mat.loop
      mat.loop.g[mat.loop <= optim.value.bound]<- 1
      mat.loop.g[mat.loop> optim.value.bound]<- 0
      sum(mat.loop.g==1)
      g1<-graph_from_adjacency_matrix(mat.loop.g, mode = 'undirected')
      complete_points<- max_cliques(g1, min=2)
      if(length(complete_points)>0){
        x.el.graph<-c()
        el.rm.graph<-c()
        for(i in 1:length(complete_points)){
          el.rm<- as.numeric(complete_points[[i]])
          el.rm<- el.rm[!el.rm%in%el.rm.graph]
          if(length(el.rm)>0){
            x.el.graph<- rbind(x.el.graph, colSums(x.el.fold[el.rm,]))
            el.rm.graph<- c(el.rm.graph, el.rm)
          }}
        x.el.fold<- rbind(x.el.fold[-el.rm.graph, , drop=F],x.el.graph)
        colSums(x.el.fold)}
      x.el.final<- rbind(x.el.final, x.el.fold)
      {  #while(  min(mat.loop)< optim.value.bound & nrow(mat.loop)>2){#
        
        # el.rm <- which(mat.loop==min(mat.loop), arr.ind=T)
        #  el.rm <- el.rm[1,]   # april 14 
        
        #el.rm<- unique(as.numeric(el.rm))
        
        # mat.loop<- mat.loop[-el.rm, -el.rm, drop=F] 
        #if(sum(dim(mat.loop))==0){mat.loop=matrix(c(0,0,0,0),nrow = 2)}
        #x.el.fold<- rbind( x.el.fold[-el.rm, , drop=F], colSums(x.el.fold[el.rm,]))
        #colnames(x.el.fold)<-row.names(data.sorted)
        #}
      }
    }}
  
  rems<- c((l-rem+1):l)
  #print(rems)
  x.el.rem<- x.el.sorted[rems,]
  ind.rem<- c(1:rem)
  combn.rem<- t(combn(ind.rem, 2))
  initial.rems<-apply(combn.rem, 1, function(x){fit.elements.weight.reduce( x.el.rem[x[1],], 
                                                                            x.el.rem[x[2],],
                                                                            x.nb, ss.other )})
  
  ###may 6th
  tempMatLoop<-matrix(1000, nrow=rem, ncol=rem)
  tempMatLoop.2<-lower.tri(diag(0, nrow=rem, ncol=rem))
  tempMatLoop[tempMatLoop.2] <- initial.rems[1, ]
  mat.loop<-t(tempMatLoop)
  #diag(mat.loop) <- 1e5;
  mat.loop.g<- mat.loop
  mat.loop.g[mat.loop <= optim.value.bound]<- 1
  mat.loop.g[mat.loop> optim.value.bound]<- 0
  sum(mat.loop.g==1)
  g1<-graph_from_adjacency_matrix(mat.loop.g, mode = 'undirected')
  complete_points<- max_cliques(g1, min=2)
  x.el.graph<-c()
  el.rm.graph<-c()
  
  #print(complete_points)
  if(length(complete_points)>0){
    for(i in 1:length(complete_points)){
      el.rm<- as.numeric(complete_points[[i]])
      el.rm<- el.rm[!el.rm%in%el.rm.graph]
      if(length(el.rm)>0){
        x.el.graph<- rbind(x.el.graph, colSums(x.el.rem[el.rm,]))
        el.rm.graph<- c(el.rm.graph, el.rm)
      }}
    x.el.rem<- rbind(x.el.rem[-el.rm.graph, , drop=F],x.el.graph)
    colSums(x.el.rem)}
  x.el.final<- rbind(x.el.final, x.el.rem)
  {
    #while(  min(mat.loop)< -10 &nrow(mat.loop)>2){#
    #  el.rm<- which(mat.loop==min(mat.loop), arr.ind=T)
    #  el.rm<- unique(as.numeric(el.rm))
    
    # mat.loop<- mat.loop[-el.rm, -el.rm, drop=F] 
    #  if(sum(dim(mat.loop))==0){mat.loop=matrix(c(0,0,0,0),nrow = 2)}
    #  x.el.rem<- rbind( x.el.rem[-el.rm, , drop=F], colSums(x.el.rem[el.rm,]))
    #  colnames(x.el.fold)<-row.names(data.sorted)
    #}
    #x.el.final<- rbind(x.el.final, x.el.rem)
  }
  
  results<- list(x.el.red=x.el.final, data=arg.data, data.sorted=data.sorted);
  nn <- nrow(results$x.el.red)
  print("jj")
  print(jj)
  jj<-jj+1
  if(nn < size.threshold | nn == n) {
    return(results)
  } else {
    #mbhc.reduce.2(results, sampleNum=sampleNum,n=nn)
    mbhc.reduce.graph(arg.data=data, sampleNum=sampleNum,x.el=x.el.final,n=nn)
  }
  
  
  
}







#####



###############################
###############################
############################### MAG MULTIPLE                                    




#### clean up: 
mbhc.multiple.cut <- function(arg.data,x.el.cut=data.frame(), sampleNum=1) {        
  
  if(sampleNum>dim(arg.data)[2]-1){print("stop the algorithm and check the sampleNum; run mbhc.prepdata() first")}
  
  time0<-as.numeric(Sys.time())
  params<-c();  x.values<-c();  x.steps<-c();  row.loop<-c();x.step<-c()
  x.elements.list <- list()
  
  weights.list <- list()
  weights.dist <- c()
  arg.data <- as.matrix(arg.data);
  arg.data.other <- c();
  
  if(ncol(arg.data) == 2){        ##### change '1' to '2' because the 2nd column is always id
    print('executing single sample clustering....'); flush.console();
    return(mag.single(arg.data[,1])) #feb 26: the second column is the IDs
  } else {
    arg.data.other<- arg.data[, -c(sampleNum,dim(arg.data)[2]), drop=F]    #other samples in the data     FEB26: aarg.datauming last column is ID
  }
  
  #feb 27: took these out of the if statement. So every version is done based on sorted data.
  order<-order(arg.data[,sampleNum])
  arg.data.other<-arg.data.other[order,,drop=F]
  data.sorted<- arg.data[order, ]
  x.nb<- arg.data[order, sampleNum]
  IDs<- arg.data[order,dim(arg.data)[2]] #FeB26
  
  if( sum(dim(x.el.cut))==0){
    print("Mag.reduce has not been run.")
    x.elements<-data.frame(diag(length(x.nb)))
    #names(x.elements)<-c(1:length(x.nb))   feb 26
    names(x.elements)<- IDs    
  }else{
    x.elements<- x.el.cut
  } 
  
  
  
  #print(class(arg.data.other))
  
  #x.elements<- x.el.cut
  
  
  #work on the target sample--- use other samples to find weight  
  x.loop <- x.nb   
  l <- length(x.loop)
  # generate pairs of the indexes... ----------------------------------------
  indx<-c(1:length(x.nb))
  combn.idx <- t(combn(indx, 2))
  weights <- apply(combn.idx, 1, function(x) 
  { temp.elements <- rep(0, length(x.nb)); temp.elements[x] <- 1; 
  return(get.weights.within(temp.elements, arg.data.other)) } );  ### Now I Have the initial weights 
  combn.x <- t(combn(x.nb, 2))     
  initial <- cbind(combn.x, weights)
  weights.list[[1]]<-initial
  #print(initial)
  
  
  l<- nrow(x.elements)
  
  #### new matloop with reduced x.el:
  indx.cut<-c(1:nrow(x.elements))
  combn.cut <- t(combn(indx.cut, 2))
  
  ###no need to run weights...
  cat('initial fitting ... '); flush.console();
  time1 <- as.numeric(Sys.time())     
  initial.1<-apply(combn.cut, 1, function(x){fit.elements.weight( x.elements[x[1],], x.elements[x[2],],
                                                                  x.nb, arg.data.other )})
  time2 <- as.numeric(Sys.time())     
  cat('took ', time2-time1, ' seconds.\n'); flush.console();
  
  tempMatLoop<-matrix(1000, nrow=l, ncol=l)  #may 6th
  tempMatLoop.2<-lower.tri(diag(0, nrow=l, ncol=l))
  tempMatLoop[tempMatLoop.2] <- initial.1[1, ]
  mat.loop<-t(tempMatLoop)
  #diag(mat.loop) <- 1e5;
  
  
  tempPar1Loop<-lower.tri(diag(0, nrow=l, ncol=l))
  #tempPar2Loop<-lower.tri(diag(0, nrow=l, ncol=l))
  tempPar1Loop[tempPar1Loop] <- initial.1[2, ]
  #tempPar2Loop[tempPar2Loop] <- initial.1[3, ]
  
  par.1.loop<-t(tempPar1Loop)
  #par.2.loop<-t(tempPar2Loop)  
  weights.step<-c()
  var.mean <- c()
  
  
  time2<-as.numeric(Sys.time())
  s<-1
  while(nrow(x.elements) > 1){
    cat(s, '...', nrow(x.elements), '...\n'); flush.console();
    
    x.values[s]<-min(mat.loop)
    el1<-which(mat.loop==min(mat.loop), arr.ind = T)[1,1]
    el2<-which(mat.loop==min(mat.loop), arr.ind = T)[1,2]
    el.rm<-c(el1,el2)
    params<-c(params, par.1.loop[el1,el2])#, par.2.loop[el1,el2]))
    
    #what happenes at each step
    x.step <- x.elements[el.rm[1], ] + x.elements[el.rm[2],]   ###this is very clever! I dont need to update the xloop.
    ###only use x.elements rows for each cluster!
    x.steps <- rbind(x.steps, x.step)
    # print(x.steps)
    x.elements <- rbind(x.elements[-el.rm,], x.step)
    x.elements.list[[s]]<-x.elements
    
    var.mean.temp <- t(apply(x.elements, 1, get.var, x.nb=x.nb, s=s))
    weights.temp <- apply(x.elements, 1, get.weights.within, arg.data=arg.data.other)
    var.mean.temp <- cbind(var.mean.temp, weights.temp)
    var.mean <- rbind(var.mean, var.mean.temp) 
    
    #####UPDATE MATRIX HERE
    mat.loop<-mat.loop[-el.rm, -el.rm, drop=F]      
    par.1.loop<-par.1.loop[-el.rm, -el.rm, drop=F]
    
    
    mat.column.update<-c()
    par1.column.update<-c()
    #par2.column.update<-c()
    
    x1 <- x.elements[1:(nrow(x.elements)-1), ]
    x2 <- x.elements[nrow(x.elements), ]   #the ones that are put together
    #print(x1)
    #print(x2)
    
    fit.between <- t(apply(x1, 1, fit.elements.weight, element.2=x2, arg.data=x.loop, arg.data.other=arg.data.other))
    mat.column.update <- fit.between[, 1]        
    par1.column.update <- fit.between[, 2]
    #par2.column.update <- fit.between[, 3]
    
    #print(mat.loop)
    if(nrow(mat.loop) != 0){
      mat.loop<-cbind(mat.loop, mat.column.update )
      mat.loop<-rbind(mat.loop, rep(1000,dim(mat.loop)[2]))    #just to keep mat.loop square. 
      par.1.loop<-cbind(par.1.loop, par1.column.update)
      par.1.loop<-rbind(par.1.loop, rep(0, dim(par.1.loop)[2]))
      # par.2.loop<-cbind(par.2.loop, par2.column.update)
      # par.2.loop<-rbind(par.2.loop, rep(0, dim(par.2.loop)[2]))
      #print("MATLOOP UPDATE SHODE")
      #print(mat.loop)
      #print("END")
    }
    #print(mat.loop)
    
    s<-s+1 
  }  
  colnames(var.mean)<-c("mean", "var","NumberOfPoints", "step", "weights")
  
  time3 <- as.numeric(Sys.time())
  cat('loop took ', time3 - time2, ' seconds.\n')
  cat('total took ', time3 - time0, ' seconds.\n')
  
  result <- list(x.el=x.elements.list,data.sorted=data.sorted, data=arg.data, params=params, var.mean=var.mean, wights=weights.list ,freq.s=x.nb ,freq=x.nb)
  return(result)
}




###############################
###############################
############################### FIND CUT OFF MULTIPLE                                   




#MAR 9      add multiple sample cut off to find different cutoffs
mbhc.find.cut.off.single <- function(var.mean, data, x.el.list){
  #detach('package:plyr')
  library(dplyr)
  
  time0 <- as.numeric(Sys.time())
  #print(dim(mag.output))
  #var.mean <- mag.output$var.mean
  var.mean <- data.frame(var.mean, row.names = NULL)    
  #var.mean$id<-round(var.mean$mean*var.mean$NumberOfPoints,5)
  var.mean[is.na(var.mean)] <- 0
  var.mean.first <- var.mean 
  #print(var.mean)
  var.mean$var <- round(var.mean$var, 3)
  #var.mean$var <- var.mean$var + 1e-7
  #id.max<-max(var.mean$id)
  ##################### add   okay structure to the var mean: 
  #var.mean$okay<- 0
  #var.mean[var.mean$step==1,'okay']<- 1
  vars<- var.mean$var[var.mean$var>0]
  var.bound<- quantile(vars)[4]*1.5
  
  
  
  
  flag <- FALSE
  m <- c()
  str <- c()  #clone structure line from temp var.mean matrix
  str.f <- c()
  cutstr <- c()
  cutstr2 <- c()  
  
  #pointsMax<- max(rowSums(x.el.list[[1]]))
  
  
  i <- max(var.mean$step)
  # allPoints <- mag.output$x.el[[i]]
  allPoints<- x.el.list[[i]]   #feb28
  ok.points <- allPoints - allPoints  ##### To generate clean slate
  n <- length(allPoints)
  #	print(allPoints)
  id.ok <-c() #clusters with okay variance! 
  temp.mv <- c();
  temp.points <- c();
  temp.step <- c();
  temp.point.step <- c();
  while(flag == FALSE){
    #print(flag)
    print(i); flush.console();
    v <- i
    
    ############    
    ## I need to use x.el to get the breaks, because using StepI would cause problem 
    ## in identical clusters and mutations.
    ## TO FIX THIS I JUST KEPT MY ORIGINAL METHOD - 
    ## I just duplicated the rows in the case of identical mutations
    ## because there should be onlytwo rows at each step 
    ## and having only 1 row means identical clusters. 
    ############
    stepI <- var.mean[var.mean$step==i,]
    stepI.x.el <- x.el.list[[i]] 
    stepI.x.el.step <- stepI.x.el[nrow(stepI.x.el), ]  #Only look at the most recent cluster
    
    stepII <- var.mean[var.mean$step==i-1, ]  ## previous step
    #stepII$var<-stepII$var,)
    stepII.x.el <- x.el.list[[i-1]]
    
    #if( !(id.sum%in%id.ok)){    ###this means if the division on this step is not happening on the "ok" clusters. so the okay clusters should still exist in the temp. 
    #if(sum(abs(id.ok-id.sum)<0.1)==0){
    #newtemp<-temp[!(temp$id%in%id.ok),]
    #if(!id.I%in%id.ok){ 
    temp.points <- rbind(temp.points, c(i, ok.points));
    temp.step <- rbind(temp.step, c(i, stepI.x.el.step));
    temp.point.step <- rbind(temp.point.step, c(i, sum(ok.points & stepI.x.el.step)))
    
    
    print('here')
    #print(ok.points)
    #print(stepI.x.el.step)
    #print(sum( ok.points & stepI.x.el.step)==0)
    
    if(sum(ok.points & stepI.x.el.step) == 0){     ###meaning stepI.x.el is not in ok.points 
      
      #stepII.breaks<- stepII[!stepII$id%in%stepI,]      
      print("stepI ")
      #print(stepI )
      #print("stepII ")
      #print(stepII )			
      #stepII.breaks <- setdiff(stepII[,-c(4,6)], stepI[,-c(4,6)])
      stepII.breaks <- setdiff(stepII[,-4], stepI[,-4])
      if(nrow(stepII.breaks)==1){
        print("why did this happen????")
        stepII.breaks <- rbind(stepII.breaks, stepII.breaks)
      }			
      stepII.num <- unique(stepII[, 4])
      # stepII.breaks <- cbind(stepII.breaks, stepII[1:2,c(4,6)])   #just add the step number back to the matrix        #############MARCh 14 add "okay"
      stepII.breaks <- cbind(stepII.breaks, step=stepII.num)   #MAY 9th 
      stepII.breaks.x.el <- setdiff(stepII.x.el, stepI.x.el)
     
      #for( ii in 1:dim(newtemp)[1]){              ###stepII.breaks always has 2 rows
      for(ii in 1:2){  
        #sample.sim <- cluster.sim(edep=depth, efrq=stepII.breaks[ii,1], shape=shape, n=10000)
        #var.sim <- c()
        #size <- ifelse(stepII.breaks$NumberOfPoints[ii]<5, 5, stepII.breaks$NumberOfPoints[ii])
        #print(size)
        #for(j in 1:500){
         # sample.sim.2 <- sample.sim[sample(1:10000, size=size),]   #temp[1,3]
          #var.sim <- c(var.sim, sd(sample.sim.2[, 4]))
        #}
       # m <- mean(var.sim) + 1e-6
      #  vv <- sd(var.sim)
       # vmax <- max(var.sim)
      #  temp.mv <- rbind(temp.mv, c(i, ii, m, vv, vmax, stepII.breaks[ii,2]))
 
    
        
        #if(stepII.breaks[ii, 'okay'] ==1 | sqrt(stepII.breaks[ii,2])*(1+stepII.breaks[ii, "weights"]) < m+vv) {        #okay cluster
        #if(stepII.breaks[ii, 'okay'] ==1 | sqrt(stepII.breaks[ii,2])*(1+stepII.breaks[ii, "weights"]) < m+vv) {        #okay cluster    ####change april 11
        if(stepII.breaks[ii, 'okay'] ==1 | stepII.breaks$var[ii]<var.bound) {        #okay cluster    ####change 6/2/2018
            
          str.f <- c(str.f, stepII.breaks[ii,1])
          #print("here")
          #print(str)
          str <- rbind(str, stepII.breaks[ii,])
          cutstr2 <- rbind(cutstr2, stepII.breaks.x.el[ii,])
          #print(str)
          
          #id.ok<-c(id.ok, stepII.breaks$id[ii])
          ok.points <- ok.points + stepII.breaks.x.el[ii,]
          #print('OKAY LAST')
          #print(sum(ok.points))
          #print(ok.points)
          
        }        
      }
    }
    
    ###stop condition
    #if(abs(sum(str$id)-id.max)<0.000001){flag<-TRUE}else{i<-(v-1)} 
    #print(i)      ###new condition: sum of id.okay== id.max
    
    #if(abs(sum(str$id)-id.max)<0.000001){flag<-TRUE}else{i<-(v-1)} 
    
    ####new condition    oct 31
    #print("LAST PRINT")
    if(sum(ok.points>1)>0){
      print("ERROR")
    }
    if(sum(ok.points) == n){
      flag<-TRUE
    } else {
      i<-(v-1)
    }  
    #if(v==i){iii<-TRUE}else{i<-(v-1)}
  }
  
  for( i in 1:nrow(str)){
    #print(i)
    xel <- x.el.list[[str$step[i]]]     #feb 28
    #print(xel)
    for(j in 1:nrow(xel)){
      m1 <- mean(data[xel[j,]==1])   #feb 28 
      ##changed to sorted. after oc.9#### OCT 30: THis should be changed for OC11 since there is no sorting in that
      if(abs(m1-str[i,1]) < 0.001){
        cutstr<-rbind(cutstr, xel[j,])  ## what's the difference between cutstr and cutstr2????????
      }
    }
  }
  
  
  time1 <- as.numeric(Sys.time())
  cat("took ", (time1 - time0), " seconds.\n");
  #print(cutstr)
  colors<-colSums(cutstr*c(1:nrow(cutstr)))+1
  final.data<-cbind(data,colors)
  results<-list(str=str, cutclust=cutstr, cutclust2=cutstr2, colors=colors, final.data=final.data)
  return(results)  
}

mbhc.find.cut.off.multiple <- function(mbhc.var){
  library(dplyr)
  if(is.null(mbhc.var[['var.mean.all']])) return('incomplte, run mbhc.var() first');
  
  final.data.clusters<- mbhc.var$mbhc.out$data.sorted
  cuts<-list()
  colors.all<-c()
  x.el<- mbhc.var$mbhc.out$x.el
  
  for(i in 1:length(mbhc.var[['var.mean.all']])) {
    var.mean<- mbhc.var$var.mean.all[[i]]
    data<- mbhc.var$mbhc.out$data.sorted[,i,drop=F]
    
    #print(shape)
    #print(depth)
    #shape=shape
    #depth=depth
    cut.temp<-mbhc.find.cut.off.single(var.mean = var.mean, data = data,x.el.list = x.el)
    cuts[[i]]<-cut.temp
    name.temp<-paste(colnames(data),'.color',sep = "")
    
    colors.all<-cbind(colors.all, cut.temp$colors*10^(2*i-2))
    final.data.clusters<-cbind(final.data.clusters, cut.temp$colors)
    
    colnames(colors.all)[dim(colors.all)[2]]<- name.temp
    colnames(final.data.clusters)[dim(final.data.clusters)[2]]<-name.temp
  }
  
  #print('colorsall')
  #print(colors.all)
  
  colors.all<-data.frame(colors.all)
  colors.all$final<- rowSums(colors.all)
  #colors.all$final
  library(plyr)
  #print(colors.all$final)
  updated.colors<-mapvalues(colors.all$final, from = unique(colors.all$final), to=c(1:length(unique(colors.all$final))))
  colors.all$colors.final<- updated.colors
  detach('package:plyr')
  final.data.assignments<- cbind(final.data.clusters, cluster=updated.colors)
  final.data.clusters<- cbind( mbhc.var$mbhc.out$data.sorted, cluster=updated.colors)
  
  return(list( final.data.clusters=final.data.clusters,final.data.all.samples=final.data.assignments, final.colors=colors.all))
  
}


#to run this, data only contains the frequensies. 






###############################
###############################
############################### MAG VAR                                  




mbhc.var<- function(mbhc.out, sampleNum=1){
  
  #varmean is geneated based on sorted data
  x.el.list<- mbhc.out$x.el
  data<- data.frame(mbhc.out$data.sorted)
  #this also assumes that the data was ran on prep.data code so the last column has ID.
  num<- dim(data)[2]-1
  var.mean.list<-list()
  
  ss.other<- data[,-c(sampleNum,dim(data)[2]),drop=F]
  x.el.red<- x.el.list[[1]]
  
  x.el.red.okay<- apply(x.el.red, 1, function(x) paste(x, collapse = ''))  #april 15
  
  
  
  for ( i in 1:num){
    var.mean<-c()  
    x.data<- data[,i]
    names(data)
    
    #s<-7
    for( s in 1:length(x.el.list)){
      x.el<- x.el.list[[s]]
      
      temp<- t(apply(x.el, 1, get.var, x.nb=x.data, s=s))
      weights.temp <- apply(x.el, 1, get.weights.within, arg.data=ss.other)
      #okay<- apply( x.el, 1, function(x) nrow(merge(t(x),x.el.red)))
      #april 16
      okay<- apply(x.el, 1, function(x) sum(paste(x,collapse = '')%in%x.el.red.okay))
      
      
      temp <- cbind(temp, weights.temp, okay)
      var.mean<- rbind(var.mean, temp)
    }
    
    colnames(var.mean)<-c("mean", "var","NumberOfPoints", "step","weights","okay")  
    var.mean.list[[i]]<- var.mean 
    names(var.mean.list)[i]<-paste(names(data)[i],".var.mean",sep="")
    
    
  }
  
  return(list(var.mean.all=var.mean.list, mbhc.out=mbhc.out))
}











###############################
###############################
############################### MAG RUN                                  


#may 20th
### input should have ID. 
### input should only be the frequencies. Each column is one sample.
### sampleFreqs: character vector of names of the samples. (s1, s2, s3, ...)
mbhc.run<- function(data,name="",sampleFreqs= NA,sampleNum=1, depth=200){
  if(length(data$ID)==0){return(print("RUN mbhc.prepdata"))}
  if(sum(is.na(sampleFreqs))>0){ data=data }else{data=cbind(data[,sampleFreqs], ID=data$ID)}
  
  #data.prep<- mbhc.prepdata(data)
  # data should have id.
  data.prep<- data
  dim<-dim(data)[2]
  s<- sampleNum
  
  if(dim==2){
    print('Running Single Sample Clustering')
    mbhc.out.single<- mbhc.single(data.prep, sampleNum=s)
    mbhc.cut.single<- mbhc.find.cut.off(mbhc.out.single, depth=depth)
    results<- list(final.dataframe<- mbhc.cut.single$final.data, mbhc.output<- mbhc.out.single,
                   data.prep<- data.prep)
    names(results)<- c(paste(name,".final.data", sep = ""),paste(name,".mbhc.single.output", sep = ""),
                       paste(name,".data.prep",sep = ""))
    
  }
  
  if(dim>2){
    print('Running Multiple Sample Clustering')
    mbhc.red<- mbhc.reduce.graph(data.prep)
    mbhc.out<- mbhc.multiple.cut(data.prep, x.el.cut = mbhc.red$x.el.red, sampleNum=s)
    mbhc.mv <- mbhc.var(mbhc.out,sampleNum = s)
    mbhc.cut<- mbhc.find.cut.off.multiple(mbhc.mv, shape=shape, depth=depth)
    results <- list(final.dataframe<-mbhc.cut$final.data.clusters ,data.prep<- data.prep, 
                   reduce.output<- mbhc.red, mbhc.output<- mbhc.out, mbhc.cutoff<- mbhc.cut, mbhc.mvar= mbhc.mv )
    names(results)<- c(paste(name,".final.data", sep = ""), paste(name,".data.prep",sep = ""), 
                       paste(name,".reduce.output", sep = ""),
                       paste(name,".mbhc.out", sep = ""), paste(name,".mbhc.cut", sep = ""),
                       paste(name,".mbhc.vars", sep=""))
    
  }
  return(results)
  
}







#test on ss

ss2<- rbind(ss2,ss*2)
ss2<- rbind(ss2, ss[1:7,])
ss2
mag.out<- test.ss.2

test.prep<- mag.prepdata(test)

test.run<- mag.run(test.prep)

test.run$.final.data

test.var<-mag.var(test.ss.2)

cut.mult<- find.cut.off.multiple(test.var, 8, 200)

test.red<- mag.reduce.1(test.prep)
test.red.2<- mag.reduce.2(test.red)



mag.out<- mag.single(test.prep[,-2,drop=F ])
var.mean<-mag.out$var.mean
mbhc.cut.test<- mbhc.find.cut.off(mag.out)
mbhc.cut.test$str
mbhc.cut.test$colors



###
mag.out<- mbhc.multiple.cut(test.prep)
mbhc.var<-mbhc.var(mag.out)
test.cut<- mbhc.find.cut.off.multiple(mbhc.var)
test.cut$final.data.clusters

###testt mag.run































