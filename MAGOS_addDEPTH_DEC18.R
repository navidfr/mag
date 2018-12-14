#MAGOS v3

#use 
#HI

sim50x4_1

head(sim50x4_1)


test.data<- sim50x4_1[1:30,-1]
head(test.data)
class(test.data)


########################################
########################################          Prep data

mag.prepdata.v3<- function(arg.data){
  num<- dim(arg.data)[2]/2   # number of samples. 
  cat('Number of samples: ', num)
  vaf.data<- c()
  depth.data<- c()
  for(i in 1:num){
    depth<-arg.data[,i*2]+arg.data[,i*2-1]
    vaf<- arg.data[,i*2]/depth
    
    depth.data<- cbind(depth.data, depth)
    colnames(depth.data)[i]<-paste("depth.", i, sep='')
    vaf.data<- cbind(vaf.data, vaf)
    colnames(vaf.data)[i]<-paste("vaf.", i, sep='')
  }
  
  ### arg data should have the frequencies in rows and columns for samples. for 3 samples, 100 mutation: 3 columns,100rows. 
  vaf.data <- as.matrix(vaf.data)
  vaf.data <- round(vaf.data,3)
  vaf.data <- ifelse(vaf.data<1e-3, 1e-3, vaf.data)
  vaf.data <- ifelse(vaf.data> 1-(1e-3), 1-(1e-3), vaf.data)
  
  
  counts.data<- data.frame(arg.data, row.names=NULL)
  vaf.data<- data.frame(vaf.data, row.names=NULL)
  depth.data<- data.frame(depth.data, row.names=NULL)
  
  ID<- c(1:dim(counts.data)[1])
  
  counts.data$ID<- ID
  vaf.data$ID<- ID
  depth.data$ID<- ID
  
  results<- list(counts=counts.data,vafs=vaf.data, depths=depth.data)
  return(results)
}

test<- mag.prepdata.v3(test.data)
test$counts
test$vafs
test$depths



prep.data<-mag.prepdata.v3(test.data[,c(1,2)]) 


########################################
########################################          MAG SINGLE UPDATE
mag.output<- test.ms


mag.single.v3 <- function(prep.data){
  
  if(is.null(prep.data$vafs$ID)){ return(print('Run mag.prepdata on the data'))}
  time0 <- as.numeric(Sys.time())
  
  
  vaf.data<- prep.data$vafs
  vaf.data<- round(vaf.data,2)    #change 1 (aug6,2018)
  vaf.data[,-dim(vaf.data)[2]]<- ifelse(vaf.data[,-dim(vaf.data)[2]] < 1e-3, 1e-3, vaf.data[,-dim(vaf.data)[2]])
  vaf.data[,-dim(vaf.data)[2]]<- ifelse(vaf.data[,-dim(vaf.data)[2]] > 1-1e-3, 1-1e-3, vaf.data[,-dim(vaf.data)[2]])
  
  depth.data<- prep.data$depths
  
  
  x.nb <- vaf.data;
  if(!is.null(ncol(vaf.data)) && ncol(vaf.data) > 0) {
    x.nb <- vaf.data$vaf.1
  }
  params<-c();  x.values<-c();  x.steps<-c();  row.loop<-c(); x.step<-c()
  x.elements.list<-list()
  x.elements<-data.frame(diag(length(x.nb)))
  
  
  
  sort.x<-sort(vaf.data$vaf.1)
  
  #vaf.depth<- merge(prep.data$vafs,prep.data$depths, by="ID")
  
  #x.loop<-sort.x   
  sort.vaf.data <- vaf.data[order(vaf.data$vaf.1),]
  sort.depth.data<- depth.data[order(vaf.data$vaf.1),]
  
  
  names(x.elements)<-sort.vaf.data$ID
  #names(x.el.first)<- sort.vaf.data$ID
  
  
  
  x.elements.preprc<- data.frame()
  for(i in unique(sort.x)){
    temp<-colSums(x.elements[sort.x==i,])
    x.elements.preprc<- rbind(x.elements.preprc, temp)
    
  }
  names(x.elements.preprc)<- sort.vaf.data$ID
  ### names checked! correct! 
  
  
  
  #### LOOK HERE:       THIS CHANGES TO X.element.preprc
  #x.el.first<- data.frame(diag(length(x.nb)))
  #x.el.first.okay<- apply(x.el.first, 1, function(x) paste(x, collapse = '')) # april 15
  x.el.first.okay<- apply(x.elements.preprc, 1, function(x) paste(x, collapse = '')) # april 15
  x.el.first.okay.id<- c(1:dim(x.elements.preprc)[1]) 
  l<- dim(x.elements.preprc)[1]
  
  mat.loop<-matrix(0, ncol=l, nrow=l)
  par.1.loop<-matrix(0, ncol=l, nrow=l)
  
  s.pairs <- cbind(1:(l-1), 2:l)  
  #print(s.pairs)
  #s.likls <- t(apply(s.pairs, 1, fit.variants))  #calculated likelihoods on the pairs
  s.likls <- t(apply(s.pairs, 1, function(x){
    #print(x)
    #x.elem
   # x<-c(1,2)
    #fit.elements(x.elements.preprc[x[1],], x.elements.preprc[x[2],], sort.x)   
    #fit.elements.v2(as.logical(x.elements.preprc[s.pairs[1,1],]+
     #                            x.elements.preprc[s.pairs[1,1],]), sort.vaf.data)# changed to fit.elements.v2 dec 6 2018
    fit.elements.v3(as.logical(x.elements.preprc[x[1],]+
                                 x.elements.preprc[x[2],]), sort.vaf.data, sort.depth.data)
  }))
  #mat.loop<-diag(0, nrow=l, ncol=l)
  for(i in 2:l){
    mat.loop[i-1,i] <- s.likls[i-1, 1]      ### the s.likls is the likelihood between i-1)th element and the i th. 
    par.1.loop[i-1,i] <- s.likls[i-1, 2]
    # par.2.loop[i-1,i] <- s.likls[i-1, 3]
  }  
  
  time1 <- as.numeric(Sys.time())
  var.mean<-c() 
  
  ### aug12 add the step 0 to var.mean
  temp<- t(apply(x.elements.preprc, 1, get.var, x.nb=sort.x, s=0 ))
  #okay<- apply(x.elements.preprc, 1, function(x) sum(paste(x, collapse = '')%in%x.el.first.okay))
  
  okay<- apply(x.elements.preprc, 1, function(x) #dec 6 2018
    {y=paste(x, collapse = '');
    ind=x.el.first.okay%in%y;
    res=ifelse(sum(ind)==0, 0, x.el.first.okay.id[ind])
    return(res)
    })

  
  weights<- rep(0, length(okay))  ##### just to have weight  APRIL 15
  
  ###add depth to var.mean    dec 6 2018
  depths<- apply(x.elements.preprc,1, function(x) round(mean(sort.depth.data[as.logical(x),'depth.1'])))
  
  
  temp<- cbind(temp,depths,weights ,okay)
  
  var.mean<- rbind(var.mean, temp)
  s<-1
  x.elements<- x.elements.preprc
  
  m1=as.matrix(x.elements)
  m.preprc=Matrix(m1, sparse = T)   ###just to save the preproces matrix
  #x.elements.list[[1]]<-m2
  ###need to create new "n' for number of mutations. to find neighbours
  n<-length(sort.x)
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
    
    
    ########## change to sparse matrix saving:
    
    m1<- as.matrix(x.elements)
    m2<- Matrix(m1, sparse = T)
    x.elements.list[[s]]<- m2
    
    #### april 14 add "okay"
    temp<- t(apply(x.elements, 1, get.var, x.nb=sort.x, s=s ))          ####for oc.v9 I changed x.nb=x.nb to x.nb=
    #####
    #okay<- apply(x.elements, 1, function(x) nrow(merge(t(x), x.el.first)))               ################# APRIL 14 2018   takeeeesss a long time
    ####################### faster version? :
   # okay<- apply(x.elements, 1, function(x) sum(paste(x, collapse = '')%in%x.el.first.okay))
    okay<- apply(x.elements, 1, function(x) #dec 6 2018
    {y=paste(x, collapse = '');
    ind=x.el.first.okay%in%y;
    res=ifelse(sum(ind)==0, 0, x.el.first.okay.id[ind])
    return(res)
    })
    
    
    weights<- rep(0, length(okay))  ##### just to have weight  APRIL 15
    
    
    ###add depth to var.mean    dec 6 2018
    depths<- apply(x.elements,1, function(x) round(mean(sort.depth.data[as.logical(x),'depth.1'])))
    
    temp<- cbind(temp,depths,weights ,okay)
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
    if(min == 1 & max != n){
      neiIND<-which(x.elements[,max+1]==1)
      nei<-x.elements[neiIND,]
      #			print("if1")
    }else if(max == n & min != 1){
      neiIND<-which(x.elements[,min-1]==1)
      nei<-x.elements[neiIND,]
      #			print("if2")
    }else if( max!=n & min!=1){
      nei1<-which(x.elements[,min-1]==1)
      nei2<-which(x.elements[,max+1]==1)
      neiIND<-c(nei1,nei2)
      nei<-x.elements[neiIND,]
      #			print("if3")
    }
    
    if(length(nei) > 0){
      x.last <- x.elements[nrow(x.elements), ]

      temp <- apply(nei, 1, function(x){
        fit.elements.v3(as.logical(x+x.last), vaf.data=sort.vaf.data,depth.data=sort.depth.data)
      })
      
      ###try to only look at cluster neighbours 
      ###instead of dim i will use sqrt of length. Because mat.loop is always square matrix   WORKS!!! SEP12
      l2<-sqrt(length(mat.loop))
      mat.column.update<-rep(0,l2)
      mat.column.update[neiIND]<-temp[1,]
      par1.column.update<-rep(0, l2)
      par1.column.update[neiIND]<-temp[2,]
      #par2.column.update<-rep(0, l2)
      #par2.column.update[neiIND]<-temp[3,]
      
      mat.loop<-cbind(mat.loop, mat.column.update )
      mat.loop<-rbind(mat.loop, rep(0,dim(mat.loop)[2]))    #just to keep mat.loop square. 
      par.1.loop<-cbind(par.1.loop, par1.column.update)
      par.1.loop<-rbind(par.1.loop, rep(0, dim(par.1.loop)[2]))
      #par.2.loop<-cbind(par.2.loop, par2.column.update)
      #par.2.loop<-rbind(par.2.loop, rep(0, dim(par.2.loop)[2]))
    }
    
    s<-s+1 
  }
  
  #var.mean <-cbind(var.mean, 1/(var.mean[,3]))
  colnames(var.mean)<-c("mean", "var","NumberOfPoints", "step","depth", "weights", 'okay')
  
  time2 <- as.numeric(Sys.time())
  cat("while took ", (time2 - time1), " seconds.\n")
  cat("total took ", (time2 - time0), " seconds.\n")
  
  result <- list(x.el=x.elements.list,x.el.preprocess=m.preprc, 
                 params=params, var.mean=var.mean, prep.data=prep.data, 
                 vaf.sorted=sort.vaf.data,depth.sorted= sort.depth.data,
                 freq.s= sort.x, time=list(time0,time1, time2))
  return(result)
}
find.cut.off.v3 <- function(mag.output){
  time0 <- as.numeric(Sys.time())
  #print(dim(mag.output))
  var.mean <- mag.output$var.mean
  head(var.mean)
  var.mean <- data.frame(var.mean, row.names = NULL)    
  #var.mean$id<-round(var.mean$mean*var.mean$NumberOfPoints,5)
  var.mean[is.na(var.mean)] <- 0
  #var.mean.first <- var.mean 
  #print(var.mean)
  var.mean$var <- round(var.mean$var, 6)
  var.mean$var <- var.mean$var + 1e-7
  #id.max<-max(var.mean$id)
  flag <- FALSE
  m <- c()
  str <- c()  #clone structure line from temp var.mean matrix
  str.f <- c()
  cutstr <- c()
  cutstr2 <- c()  
  
  i <- max(var.mean$step)
  ##added aug 12:
  if(i==1){
    flag=T
    str<-var.mean
  }
  print(flag)
  allPoints <- mag.output$x.el[[i]]   ### added +1 with the AUG6th2018 version, because the x.element.list has an extra element.
  allPoints<- as.data.frame(as.matrix(allPoints))
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
    print('enter')
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
    stepI.x.el<- as.data.frame(as.matrix(stepI.x.el))
    
    stepI.x.el.step <- stepI.x.el[nrow(stepI.x.el), ]  #Only look at the most recent cluster
    
    stepII <- var.mean[var.mean$step==i-1, ]  ## previous step
    #stepII$var<-stepII$var,)
    if(i>1){
      stepII.x.el <- mag.output$x.el[[i-1]]
      stepII.x.el<- as.data.frame(as.matrix(stepII.x.el))
    }
    if(i==1){
      stepII.x.el <- mag.output$x.el.preprocess
      stepII.x.el<- as.data.frame(as.matrix(stepII.x.el))
      
    }
    print("II")
    #print(stepII.x.el)
    
    #if( !(id.sum%in%id.ok)){    ###this means if the division on this step is not happening on the "ok" clusters. so the okay clusters should still exist in the temp. 
    #if(sum(abs(id.ok-id.sum)<0.1)==0){
    #newtemp<-temp[!(temp$id%in%id.ok),]
    #if(!id.I%in%id.ok){ 
    temp.points <- rbind(temp.points, c(i, ok.points));
    temp.step <- rbind(temp.step, c(i, stepI.x.el.step));
    temp.point.step <- rbind(temp.point.step, c(i, sum(ok.points & stepI.x.el.step)))
    
    
    #print(ok.points)
    #print(sum(ok.points & stepI.x.el.step))
    if(sum(ok.points & stepI.x.el.step) == 0){     ###meaning stepI.x.el is not in ok.points 
      
      #stepII.breaks<- stepII[!stepII$id%in%stepI,]      
      #print("stepI ")
      #print(stepI )
      #print("stepII ")
      #print(stepII )			
      #stepII.breaks <- setdiff(stepII[,-c(4,7)], stepI[,-c(4,7)])
      stepII.breaks <- dplyr::setdiff(stepII[,-4], stepI[,-4])   #dec 6 2018
      if(nrow(stepII.breaks)==1){
        print("why did this happen????")
        stepII.breaks <- rbind(stepII.breaks, stepII.breaks)
      }			
      stepII.num <- stepII[, 4]
      stepII.breaks <- cbind(stepII.breaks, stepII[1:2,4])   #just add the step number back to the matrix
      stepII.breaks.x.el <- dplyr::setdiff(stepII.x.el, stepI.x.el)      
      
      #for( ii in 1:dim(newtemp)[1]){              ###stepII.breaks always has 2 rows
      for(ii in 1:2){  
        #sample.sim <- cluster.sim(edep=depth, efrq=stepII.breaks[ii,1], shape=shape, n=10000)
        var.sim <- c()
        size <- ifelse(stepII.breaks$NumberOfPoints[ii]<5, 5, stepII.breaks$NumberOfPoints[ii])
        #print(size)
        #for(j in 1:500){
        # sample.sim.2 <- sample.sim[sample(1:10000, size=size),]   #temp[1,3]
        #  var.sim <- c(var.sim, var(sample.sim.2[, 4]))
        #}
        var.sim<- mag.exp.var.v2(efrq=stepII.breaks[ii,1], edep=stepII.breaks[ii,'depth'], num=size)
        
        
        m <- mean(var.sim, na.rm=T) + 1e-6
        vv <- sd(var.sim, na.rm=T)
        vmax <- max(var.sim, na.rm=T)
        #temp.mv <- rbind(temp.mv, c(i, ii, m, vv, vmax, stepII.breaks[ii,2]))
        
        #if(stepII.breaks[ii, 'NumberOfPoints'] < 11 | stepII.breaks[ii,2]*(1+stepII.breaks[ii, "weights"]) < m+vv) {        #okay cluster
        #if(stepII.breaks[ii, 'okay'] ==1 | (stepII.breaks[ii,2])*(1+stepII.breaks[ii, "weights"]) < m+vv) {        #okay cluster   APRIL 15
        if(stepII.breaks[ii, 'okay'] >0 |
           stepII.breaks[ii,2] < m+3*vv ) {   #dec 6 2018
          str.f <- c(str.f, stepII.breaks[ii,1])
          print("here")
          #print(str)
          str <- rbind(str, stepII.breaks[ii,])
          #print(stepII.breaks.x.el)
          cutstr2 <- rbind(cutstr2, stepII.breaks.x.el[ii,])
          #print(str)
          
          #id.ok<-c(id.ok, stepII.breaks$id[ii])
          ok.points <- ok.points + stepII.breaks.x.el[ii,]
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
    if(str$step[i]!=0){
      xel <- mag.output$x.el[[str$step[i]]]}
    if(str$step[i]==0){
      xel <- mag.output$x.el.preprocess}
    
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
  final.data<-cbind(mag.output$vaf.sorted,colors)
  results<-list(str=str, cutclust=cutstr, cutclust2=cutstr2, colors=colors, final.data=final.data)
  return(results)  
}


hist(test.prep$vafs$vaf.1,50)
test.prep<- mag.prepdata.v3(sim50x4_1[,c(2,3)])
test.ms<- mag.single.v3(test.prep)
test.cut<- find.cut.off.v3(test.ms)
plot(test.cut$final.data$vaf.1, test.ms$depth.sorted$depth.1, col=test.cut$final.data$colors, pch=19)


### dec 2018 version: DONT NEED OPTIM ANYMORE
optim.obj.v3 <- function(elements,vaf.depth){     
  
  vafs<- vaf.depth[elements,'vaf']
  vaf<- mean(vafs)
  depth<- mean(vaf.depth[elements, 'depth'])
  
  s1<- depth*vaf
  s2<- depth-s1
  
  R <-dbeta(vafs, shape1=s1, shape2=s2, log=T) #log should be F- because if Log=T we could have positive and negative values. We should add log

  v <- sd(vafs)
  v<- var(vafs)
  range<- max(vafs)-min(vafs)

  v <- ifelse(v<0.005, 0.005, v)
  range<- ifelse(range<0.005,0.005,range)
  mad<- ifelse(mad<0.005,0.005,mad)
  
  #R <- R/v
  R <- R/(range*v)

  results <- -mean(R)#/exp(5*(weight+0.0001))
  #print(paste(results, 'end'))
  return(results)  
}  #MAY 6th version. Uses range- only one par.


mag.dist.v3<- function(vafs, depths){
  vaf<- mean(vafs)
  depth<- mean(depths)
  
  s1<- depth*vaf
  s2<- depth-s1
  
  R <-dbeta(vafs, shape1=s1, shape2=s2, log=T) #log should be F- because if Log=T we could have positive and negative values. We should add log
  
  v <- sd(vafs)
  v<- var(vafs)
  range<- max(vafs)-min(vafs)
  
  v <- ifelse(v<0.005, 0.005, v)
  range<- ifelse(range<0.005,0.005,range)
  
  #R <- R/v
  R <- R/(range*v)
  
  dist <- -mean(R)#/exp(5*(weight+0.0001))
  dist<- ifelse(dist>0, -0.001,dist)
  results<- list(dist=dist, s1=s1, s2=s2)
  #print(paste(results, 'end'))
  return(results)  
}
fit.elements.v3<- function(elements,vaf.data, depth.data){
#  fit.obj<- c()
 # fit.minval<- c()
  dist.val<-c()
  dist.shape<- c()
  for(i in 1:(dim(vaf.data)[2]-1)){
    #print(i)
    variants<- vaf.data[elements, i]
    depths<- depth.data[elements, i]
    
    dist<- mag.dist.v3(variants,depths)
    dist.val<- c(dist.val, dist$dist)
    dist.shape<- c(dist.shape, dist$s1)
    
    #fit1<-optimize(y=variants,exp.depth=exp.depth ,f=optim.obj, interval = c(1,exp.depth-1),tol=0.001)
    #fit.obj<- c(fit.obj, fit1$object)
    #fit.minval<- c(fit.minval, fit1$minimum)
  }
  
  #return(c(max(fit.obj),fit.minval[which.max(fit.obj)]))
  return(c(max(dist.val),dist.shape[which.max(dist.val)]))
}










#### MAG MULTIPLE UPDATE



test.data.mult<- mag.prepdata.v3(test.data[,c(1,2,3,4,5,6)])
prep.data<- test.data.mult

###reduce
mag.reduce.graph.v3<-function(prep.data,x.el=data.frame(), n=-100){
  
  
  #sampleNum<-1  #added oct 11 2018
  
  order.var<- apply(prep.data$vafs[,-which(names(prep.data$vafs)=='ID')],1,prod)
  order.ind<- order(order.var)
  sorted.order.var<- order.var[order.ind]
  
  vaf.data<- prep.data$vafs
  depth.data<- prep.data$depths
  
  vaf.sorted<- vaf.data[order.ind,]
  depth.sorted<- depth.data[order.ind,]
  
  IDs<- vaf.sorted[, 'ID']
  
  if(dim(x.el)[1]==0){
    l<- dim(vaf.data)[1]
    x.el<- data.frame(diag(l))
    names(x.el)<- IDs
  }
  
  
  mean.frqs<-t(apply( x.el, 1, function(x){mean(sorted.order.var[as.logical(x)])}))
  x.el.sorted<-x.el[order(mean.frqs),]
  x.el.final<- data.frame()
  
  #x.nb<- data.sorted[,sampleNum]
  #ss.other<- data.sorted[,-c(sampleNum, dim(data.sorted)[2]), drop=F]
  
  size.threshold<-10
  
  l<-nrow(x.el)
  full<-floor(l/size.threshold)
  rem<-l%%size.threshold
  row.fold<-0
  
  if( rem < 5 ){ full<- full-1; rem<- rem+size.threshold}
  
  
  # optim.value.bound<-fit.variants.weight(c(0.4,0.46,0.05))[1]
  
  
  #optim.value.bound<-fit.elements.v2(c(TRUE,TRUE,TRUE), cbind(c(0.4,0.46, 0.43),c(1,2,3)))[1]

  
  depth.bound<- min(apply(depth.sorted[,-dim(depth.sorted)[2]],2,mean))
  
  
  
  optim.value.bound<- fit.elements.v3(c(TRUE,TRUE), cbind(c(0.4,0.44),c(1,2)),
                                      cbind(c(depth.bound,depth.bound),c(1,2)))[1]
    
  
  
  
  ind.full<- c(1:size.threshold)
  combn.el<- t(combn(ind.full,2))
  
  if(full>0){
    for(i in 1:full){
      
      fold<- c((1+(i-1)*size.threshold):(i*size.threshold))  
      x.el.fold<-x.el.sorted[fold, ] 
      #print(i)
      #initial.full<-apply(combn.el, 1, function(x){fit.elements.weight.reduce( x.el.fold[x[1],], 
      #                                                                        x.el.fold[x[2],],
      #                                                                       x.nb, ss.other )})
      
      
      initial.full<- apply(combn.el, 1, function(x){
      #fit.elements.v2(as.logical(x.el.fold[x[1],]+x.el.fold[x[2],]),data.sorted)})
      fit.elements.v3(as.logical(x.el.fold[x[1],]+x.el.fold[x[2],]),vaf.sorted, depth.sorted)})

      
      tempMatLoop<-matrix(0, nrow=size.threshold, ncol=size.threshold)
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
  # initial.rems<-apply(combn.rem, 1, function(x){fit.elements.weight.reduce( x.el.rem[x[1],], 
  #                                                                           x.el.rem[x[2],],
  #                                                                           x.nb, ss.other )})
  initial.rems<- apply(combn.rem, 1, function(x){
    #fit.elements.v2(as.logical(x.el.rem[x[1],]+x.el.rem[x[2],]),arg.data)})
    fit.elements.v3(as.logical(x.el.rem[x[1],]+x.el.rem[x[2],]),vaf.sorted, depth.sorted)})
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
  
  results<- list(x.el.red=x.el.final, vaf.data=vaf.data,vaf.sorted=vaf.sorted,
                 depth.data=depth.data,depth.sorted=depth.sorted);
  nn <- nrow(results$x.el.red)
  #print("jj")
  #print(jj)
  #jj<-jj+1
  if(nn < size.threshold | nn == n) {
    return(results)
  } else {
    #mag.reduce.2(results, sampleNum=sampleNum,n=nn)
    mag.reduce.graph.v3(prep.data = prep.data,x.el=x.el.final,n=nn)
  }
  
  
  
}


test.red<- mag.reduce.graph.v3(prep.data)
test.red$x.el.red
test.red$vaf.sorted





mag.multiple.v3<- function(prep.data,x.el.cut=data.frame()) {        
  
  
  if(1>dim(prep.data$vafs)[2]-1){print("stop the algorithm and check the sampleNum; run mag.prepdata() first")}
  
  time0<-as.numeric(Sys.time())
  params<-c();  x.values<-c();  x.steps<-c();  row.loop<-c();x.step<-c()
  x.elements.list <- list()
  
  
  #weights.list <- list()
  #weights.dist <- c()
 # arg.data <- as.matrix(arg.data);
  #arg.data.other <- c();
  
  #sampleNum<-1 #used for ordering only
  

  
  if(ncol(prep.data$vafs) == 2){        ##### change '1' to '2' because the 2nd column is always id
    print('executing single sample clustering....'); flush.console();
    #return(mag.single(arg.data[,1])) #feb 26: the second column is the IDs
  } #else {
    #arg.data.other<- arg.data[, -c(sampleNum,dim(arg.data)[2]), drop=F]    #other samples in the data     FEB26: aarg.datauming last column is ID
  #}
  
  #### DO I NEED arg.data.other??  DEC 2018
  
  
  #feb 27: took these out of the if statement. So every version is done based on sorted data.
  #order<-order(arg.data[,sampleNum])
  
  #data.sorted<- arg.data[order, ]
  #IDs<- arg.data[order,dim(arg.data)[2]] #FeB26
  
  
  order.var<- apply(prep.data$vafs[,-which(names(prep.data$vafs)=='ID')],1,prod)
  order.ind<- order(order.var)
  sorted.order.var<- order.var[order.ind]
  
  vaf.data<- prep.data$vafs
  depth.data<- prep.data$depths
  
  vaf.sorted<- vaf.data[order.ind,]
  depth.sorted<- depth.data[order.ind,]
  
  IDs<- vaf.sorted[, 'ID']
  
  
  
  
  
  
  #oct 18 2018    I WILL KEEP THE SORTED VERSION. 

  n<- dim(vaf.data)[1]
  
  
  # not need? OCT 2018
  #arg.data.other<-arg.data.other[order,,drop=F]
  #x.nb<- arg.data[order, sampleNum]
  
  if( sum(dim(x.el.cut))==0){
    print("Mag.reduce has not been run.")
    x.elements<-data.frame(diag(n))
    #names(x.elements)<-c(1:length(x.nb))   feb 26
    names(x.elements)<- IDs    
  }else{
    x.elements<- x.el.cut
  } 
  
  
  

  
  l<- nrow(x.elements)
  
  #### new matloop with reduced x.el:
  indx.cut<-c(1:nrow(x.elements))
  combn.cut <- t(combn(indx.cut, 2))
  
  ###no need to run weights...
  cat('initial fitting ... '); flush.console();
  time1 <- as.numeric(Sys.time())     
  #initial.1<-apply(combn.cut, 1, function(x){fit.elements.weight( x.elements[x[1],], x.elements[x[2],],
  #                                                               x.nb, arg.data.other )})
  
  initial.1<- apply(combn.cut, 1, function(x){
    #fit.elements.v2(as.logical(c(x.elements[x[1],]+x.elements[x[2],])),data.sorted)})
    fit.elements.v3(as.logical(c(x.elements[x[1],]+x.elements[x[2],])),vaf.sorted, depth.sorted)})

  
  
  time2 <- as.numeric(Sys.time())     
  cat('took ', time2-time1, ' seconds.\n'); flush.console();
  
  tempMatLoop<-Matrix(0, nrow=l, ncol=l, sparse = T)  #may 6th   DEC 2018: added sparse
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
  #weights.step<-c()
  #var.mean <- c()
  
  
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
    
    ### OCT 2018
    #var.mean.temp <- t(apply(x.elements, 1, get.var, x.nb=x.nb, s=s))
    #weights.temp <- apply(x.elements, 1, get.weights.within, arg.data=arg.data.other)
    #var.mean.temp <- cbind(var.mean.temp, weights.temp)
    #var.mean <- rbind(var.mean, var.mean.temp) 
    
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
    
    #fit.between <- t(apply(x1, 1, fit.elements.weight, element.2=x2, 
    #                      arg.data=x.loop, arg.data.other=arg.data.other))
    
    #oct 2018
    fit.between<- apply(x1, 1, function(x){
      #fit.elements.v2(as.logical(x+x2), data.sorted)} )
      fit.elements.v3(as.logical(x+x2), vaf.sorted, depth.sorted)} )
    
    mat.column.update <- fit.between[1,]        
    par1.column.update <- fit.between[2,]
    #par2.column.update <- fit.between[, 3]
    
    #print(mat.loop)
    if(nrow(mat.loop) != 0){
      mat.loop<-cbind(mat.loop, mat.column.update )
      mat.loop<-rbind(mat.loop, rep(0,dim(mat.loop)[2]))    #just to keep mat.loop square. 
      par.1.loop<-cbind(par.1.loop, par1.column.update)
      par.1.loop<-rbind(par.1.loop, rep(0, dim(par.1.loop)[2]))
      # par.2.loop<-cbind(par.2.loop, par2.column.update)
      # par.2.loop<-rbind(par.2.loop, rep(0, dim(par.2.loop)[2]))
      #print("MATLOOP UPDATE SHODE")
      #print(mat.loop)
      #print("END")
    }
    # print(mat.loop)
    
    s<-s+1 
  }  
  # colnames(var.mean)<-c("mean", "var","NumberOfPoints", "step", "weights")
  
  time3 <- as.numeric(Sys.time())
  cat('loop took ', time3 - time2, ' seconds.\n')
  cat('total took ', time3 - time0, ' seconds.\n')
  
  result <- list(x.el=x.elements.list,vaf.sorted=vaf.sorted,depth.sorted=depth.sorted,
                 prep.data=prep.data)#,var.mean=var.mean, wights=weights.list ,freq.s=x.nb ,freq=x.nb)
  return(result)
}


test.multiple<- mag.multiple.v3(prep.data)

mag.out<- test.multiple


mag.var.v3<- function(mag.out){
  
  #varmean is geneated based on sorted data
  x.el.list<- mag.out$x.el
  data<- data.frame(mag.out$vaf.sorted)
  #this also assumes that the data was ran on prep.data code so the last column has ID.
  num<- dim(data)[2]-1
  var.mean.list<-list()
  
  #ss.other<- data[,-c(sampleNum,dim(data)[2]),drop=F]
  x.el.red<- x.el.list[[1]]
  x.el.red.okay<- apply(x.el.red, 1, function(x) paste(x, collapse = ''))  #april 15
  x.el.red.IDs<- as.numeric(names(x.el.list[[1]]))
  
  
  for ( i in 1:num){
    var.mean<-c()  
    x.data<- data[,i]
    names(data)
    depth.sorted.each<- mag.out$depth.sorted[,i]
    #s<-7
    for( s in 1:length(x.el.list)){
      x.el<- x.el.list[[s]]
      
      temp<- t(apply(x.el, 1, get.var, x.nb=x.data, s=s))
      #weights.temp <- apply(x.el, 1, get.weights.within, arg.data=ss.other)
      #okay<- apply( x.el, 1, function(x) nrow(merge(t(x),x.el.red)))
      #april 16
      #okay<- apply(x.el, 1, function(x) sum(paste(x,collapse = '')%in%x.el.red.okay))
      okay<- apply(x.el, 1, function(x){y=paste(x,collapse = ''); ind=x.el.red.okay%in%y;
      res=ifelse(sum(ind)==0,0,x.el.red.IDs[ind])
      return(res)})
      
      depths<- apply(x.el,1, function(x) round(mean(depth.sorted.each[as.logical(x)])))
      
      
      
      #temp <- cbind(temp, weights.temp, okay)
      temp <- cbind(temp, depths ,okay)
      var.mean<- rbind(var.mean, temp)
    }
    
    #colnames(var.mean)<-c("mean", "var","NumberOfPoints", "step","weights","okay")  
    colnames(var.mean)<-c("mean", "var","NumberOfPoints", "step","depth","okay")  
    var.mean.list[[i]]<- var.mean 
    names(var.mean.list)[i]<-paste(names(data)[i],".var.mean",sep="")
    
    
  }
  
  return(list(var.mean.all=var.mean.list, mag.out=mag.out))
}


test.var<- mag.var.v3(mag.out)
test.var$var.mean.all$vaf.1.var.mean

test.multiple$x.el
plot(prep.data$vafs[,c(3,2)], col=mag.step.check(test.multiple$x.el[[28]]))


mag.var<- test.var
mag.var$var.mean.all

find.cut.off.multiple.v3<-function(mag.var, depth){
  
  if(is.null(mag.var[['var.mean.all']])) return('incomplte, run mag.var() first');
  
  final.data.clusters<- mag.var$mag.out$vaf.sorted
  cuts<-list()
  colors.all<-c()
  x.el<- mag.var$mag.out$x.el
  
  for(i in 1:length(mag.var[['var.mean.all']])) {
    var.mean<- mag.var$var.mean.all[[i]]
    data<- mag.var$mag.out$data.sorted[,i,drop=F]
    
    #print(shape)
    #print(depth)
    #shape=shape
    depth=depth
    #print('zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz')
    cut.temp<-find.cut.off.single.v2(var.mean = var.mean, data = data,x.el.list = x.el, depth )
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
  
  #print(colors.all$final)
  updated.colors<-plyr::mapvalues(colors.all$final, from = unique(colors.all$final), to=c(1:length(unique(colors.all$final))))
  colors.all$colors.final<- updated.colors
  final.data.assignments<- cbind(final.data.clusters, cluster=updated.colors)
  final.data.clusters<- cbind( mag.var$mag.out$data.sorted, cluster=updated.colors)
  
  return(list( final.data.clusters=final.data.clusters,final.data.all.samples=final.data.assignments, final.colors=colors.all))
  
}

x.el.list<- mag.out$x.el
var.mean<- test.var$var.mean.all$vaf.1.var.mean
unique(var.mean[,4])
depth.sorted<- mag.out$depth.sorted


cut.off.multiple.v3  <-function(mag.var){
  
  time0 <- as.numeric(Sys.time())
  #print(dim(mag.output))
  #var.mean <- mag.output$var.mean
  mag.var.list<-mag.var$var.mean.all
  for( i in 1:length(mag.var.list)){
    temp<- mag.var.list[[i]]
    temp<- data.frame(temp, row.names = NULL)  
    temp[is.na(temp)]<-0
    temp$var<- round(temp$var, 6)
    temp$var<- temp$var+1e-7
    mag.var.list[[i]]<-temp
  }
  names(mag.var.list)<- names(mag.var$var.mean.all)
  
  x.el.list<- mag.var$mag.out$x.el
  
  #var.mean <- data.frame(var.mean, row.names = NULL)  
  #var.mean[is.na(var.mean)] <- 0
  #var.mean.first <- var.mean 
  #var.mean$var <- round(var.mean$var, 6)
  #var.mean$var <- var.mean$var + 1e-7
  
  flag <- FALSE
  m <- c()
  str <- c()  #clone structure line from temp var.mean matrix
  str.f <- c();  cutstr <- c();cutstr2 <- c()  

  i <- max(mag.var.list[[1]]$step)
  # allPoints <- mag.output$x.el[[i]]
  allPoints<- x.el.list[[i]]   #feb28
  ok.points <- allPoints - allPoints  ##### To generate clean slate
  n <- length(allPoints)
  #	print(allPoints)
  id.ok <-c() #clusters with okay variance! 
  temp.mv <- c(); temp.points <- c(); temp.step <- c();temp.point.step <- c();
  
  while(flag == FALSE){
    print(i); flush.console();
    ############    
    ## I need to use x.el to get the breaks, because using StepI would cause problem 
    ## in identical clusters and mutations.
    ## TO FIX THIS I JUST KEPT MY ORIGINAL METHOD - 
    ## I just duplicated the rows in the case of identical mutations
    ## because there should be onlytwo rows at each step 
    ## and having only 1 row means identical clusters. 
    ############
    v <- i
    stepI.x.el <- x.el.list[[i]] 
    stepI.x.el.step <- stepI.x.el[nrow(stepI.x.el), ]  #Only look at the most recent cluster
    stepII.x.el <- x.el.list[[i-1]]
    stepII.breaks.x.el <- dplyr::setdiff(stepII.x.el, stepI.x.el)
    
    temp.points <- rbind(temp.points, c(i, ok.points));
    temp.step <- rbind(temp.step, c(i, stepI.x.el.step));
    temp.point.step <- rbind(temp.point.step, c(i, sum(ok.points & stepI.x.el.step)))
    
    print(sum( ok.points & stepI.x.el.step)==0)
    
    keep<- 0
    for( j in 1:length(mag.var.list)){
      var.mean<- mag.var.list[[j]]
      stepI <- var.mean[var.mean$step==i,]
      stepII <- var.mean[var.mean$step==i-1, ]  ## previous step
      
      if(sum(ok.points & stepI.x.el.step) == 0){     ###meaning stepI.x.el is not in ok.points 
        
        #stepII.breaks <- setdiff(stepII[,-c(4,6)], stepI[,-c(4,6)])
        stepII.breaks <- dplyr::setdiff(stepII[,-4], stepI[,-4])
        if(nrow(stepII.breaks)==1){
          print("why did this happen????")
          stepII.breaks <- rbind(stepII.breaks, stepII.breaks)
        }			
        stepII.num <- unique(stepII[, 4])
        stepII.breaks <- cbind(stepII.breaks, step=stepII.num)   #MAY 9th 
        # stepII.breaks <- cbind(stepII.breaks, stepII[1:2,c(4,6)])   #just add the step number back to the matrix        #############MARCh 14 add "okay"
        #for( ii in 1:dim(newtemp)[1]){              ###stepII.breaks always has 2 rows
        for(ii in 1:2){  
          #sample.sim <- cluster.sim(edep=depth, efrq=stepII.breaks[ii,1], shape=shape, n=10000)
          var.sim <- c()
          size <- ifelse(stepII.breaks$NumberOfPoints[ii]<5, 5, stepII.breaks$NumberOfPoints[ii])
          # for(j in 1:500){
          #   sample.sim.2 <- sample.sim[sample(1:10000, size=size),]   #temp[1,3]
          #   var.sim <- c(var.sim, sd(sample.sim.2[, 4]))
          # }
          var.sim<- mag.exp.var.v2(efrq = stepII.breaks[ii,1], edep = stepII.breaks[ii,'depth'], num=size)
          
          m <- mean(var.sim) + 1e-6
          vv <- sd(var.sim)
          vmax <- max(var.sim)
          #       if(stepII.breaks[ii,2]>=vmax+vv){i<-(i-1)} 
          #       if(stepII.breaks[ii, 'NumberOfPoints'] == 1 | stepII.breaks[ii,2] < m+vv) {        #okay cluster
          #       if(stepII.breaks[ii, 'okay'] ==1 | sqrt(stepII.breaks[ii,2])*(1+stepII.breaks[ii, "weights"]) < m+vv) {        #okay cluster
          if(stepII.breaks[ii, 'okay'] >0 |
             stepII.breaks[ii,2] < m+3*vv ) {
            keep<-keep+1
          }#changed to vmax on oct 14 2018
          
          cutstr2 <- rbind(cutstr2, stepII.breaks.x.el[ii,])
          #id.ok<-c(id.ok, stepII.breaks$id[ii])
          ok.points <- ok.points + stepII.breaks.x.el[ii,]
        }        
      }
    }
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





