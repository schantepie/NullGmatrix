#ped : list of pedigree. Each pedigree should have 3 columns 'id', 'dam', 'sire'
#G : list of G matrix. Each G matrix correspond to animal part of MCMCglmm animal models. 
#animalToExtractByPop : list of vector of animal to extract in each pop. By default all individuals with an entry in the pedigree are returned.
##Example below the function




simul_Null_Gmatrix<-function(ped,G,animalToExtractByPop=NULL){
  
  require(pedigree)
  require(mvtnorm)
  #initialize containers and parameters   
  if (length(G)!=length(ped)) stop("pb number pedigree vs G")
  nbpop=length(G)
  # for univariate case, transforms G in matrix
  if (is.null(dim(G[[1]]))) {
    for(popNum in 1:nbpop) {
      G[[popNum]]=as.matrix(G[[popNum]])
    }
  }
  
  nbTraits=sqrt(dim(G[[1]])[2])
  nbsample=dim(G[[1]])[1]

  pb <- txtProgressBar(min = 0,      
                       max = nbsample, 
                       style = 3,   
                       width = 50,  
                       char = "=") 
  fullpedigree=NULL
  sizepedigree=NULL
  
  for(popNum in 1:nbpop) {
    # if(!all(colnames(ped[[i]][1:3])==c("id","dam","sire"))) stop("pb in pedigree name")
    pedi= ped[[popNum]]
    
    #unsure that same id in different pop keep different
    pedi$id=paste("pop",popNum,pedi$id,sep="_")
    pedi$sire=paste("pop",popNum,pedi$sire,sep="_")
    pedi$dam=paste("pop",popNum,pedi$dam,sep="_")
    
    ###### Transform NA to zero in pedigree 
    pedi[pedi==paste("pop",popNum,NA,sep="_")] <-0
    
    ###### Add dummyparents as founders
    pedi=rbind(cbind(paste("pop",popNum,"dummyDam",1:length(which(pedi[,"dam"]==0 & pedi[,"sire"]!=0)),sep="_"),0,0),
               cbind(paste("pop",popNum,"dummySire",1:length(which(pedi[,"dam"]!=0 & pedi[,"sire"]==0)),sep="_"),0,0),
               as.matrix(pedi))
    pedi[which(pedi[,"dam"]==0 & pedi[,"sire"]!=0),"dam"]=paste("pop",popNum,"dummyDam",1:length(which(pedi[,"dam"]==0 & pedi[,"sire"]!=0)),sep="_")
    pedi[which(pedi[,"dam"]!=0 & pedi[,"sire"]==0),"sire"]=paste("pop",popNum,"dummySire",1:length(which(pedi[,"dam"]!=0 & pedi[,"sire"]==0)),sep="_")
    pedi=as.data.frame(pedi)
    
    ###### estimate inbreeding coef
    F <- calcInbreeding(pedi[,1:3])
    
    ###### bind all pedigree ,add columns for location of parents in pedigree, inbreeding coef and breeding values
    parent1=NA
    parent2=NA
    BreedingValue=matrix(NA,ncol=nbTraits,nrow=dim(pedi)[1])
    colnames(BreedingValue)=paste("BreedingValue",(1:nbTraits),sep="")
    generation=countGen(pedi)
    pedi=cbind(pedi[,1:3],generation,popNum,parent1,parent2,F,BreedingValue)
    fullpedigree=rbind(fullpedigree,pedi)
    sizepedigree= c(sizepedigree,dim(pedi)[1])
  }
  
  
  #inialise the countainer of breeding values for all pop and all sample
  sizepedigree=sum(sizepedigree)
  fullpedigreeSample=array(NA,c(sizepedigree,nbTraits,nbsample))
  
  cpt_progressBar=0
  nsample=1
  for (nsample in 1:nbsample){
    cpt_progressBar=cpt_progressBar+1
    #reset BreedingValue to avoid undetected error
    fullpedigree[,grep("BreedingValue",colnames(fullpedigree))]=NA
    
    for (popNum in 1:nbpop){
      Founders=which(fullpedigree$generation==0 & fullpedigree$popNum==popNum)
      ## predict breeding values of founders
      fullpedigree[Founders,grep("BreedingValue",colnames(fullpedigree))]=rmvnorm(n = length(Founders), mean=rep(0,nbTraits), sigma=matrix(G[[popNum]][nsample,],nbTraits,nbTraits))
    }
    ##### locate parent 
    fullpedigree$parent1=match(fullpedigree$sire,fullpedigree$id)
    fullpedigree$parent2=match(fullpedigree$dam,fullpedigree$id)
    ##### store the rows corresponding to founders and predict their breeding values
    Founders_pos=which(fullpedigree$generation==0)
    ##### store the rows corresponding to non-founders
    nonFounders_pos=c(1:dim(fullpedigree)[1])[-Founders_pos]
    ##### Randomize founders among populations
    fullpedigree[Founders_pos,grep("BreedingValue",colnames(fullpedigree))]=fullpedigree[sample(Founders_pos,size = length(Founders_pos),replace = FALSE),grep("BreedingValue",colnames(fullpedigree))]
    ##### Estimate G based on the breeding values of randomized founders
    GAfterRand=matrix(NA,nrow = length(G),ncol=dim(G[[1]])[2])
    for (i in 1:nbpop) {
      if (dim(G[[i]])[2]==1){
        GAfterRand[i,]=var(fullpedigree[fullpedigree$popNum==i & fullpedigree$sire==0 & fullpedigree$dam==0 ,grep("BreedingValue",colnames(fullpedigree))])
      }else{
        GAfterRand[i,]=c(cov(fullpedigree[fullpedigree$popNum==i & fullpedigree$sire==0 & fullpedigree$dam==0 ,grep("BreedingValue",colnames(fullpedigree))]))
      }
    }
    
    ##### Predict breeding value of non-founders according to standard rules of polygenic inheritance 
    # Estimate scaling factor of G using the inbreeding coeff of parents
    varscale=matrix(NA,dim(fullpedigree)[1],1)
    varscale[nonFounders_pos]=as.matrix(1-((fullpedigree[fullpedigree[nonFounders_pos,"parent1"],"F"]+fullpedigree[fullpedigree[nonFounders_pos,"parent2"],"F"])/2))

    for (gen in (unique(fullpedigree$generation)[-1])){ #go through generations 
      
      ##### Estimate mean BV for all indiv at generation gen 
      fullpedigree[which(fullpedigree$generation==gen),grep("BreedingValue",colnames(fullpedigree))]=
                (fullpedigree[fullpedigree[which(fullpedigree$generation==gen),"parent1"],grep("BreedingValue",colnames(fullpedigree))]+
                 fullpedigree[fullpedigree[which(fullpedigree$generation==gen),"parent2"],grep("BreedingValue",colnames(fullpedigree))])/2
      
      # initialize pop number having the generation equals to gen
      popgen=unique(fullpedigree[which(fullpedigree$generation==gen),"popNum"])
             
      ##### Add within family deviation that depends on the G matrices of randomized founders and the scaling factor
      for (pop in popgen) {
          if(length(which(fullpedigree$popNum==pop & fullpedigree$generation==gen))>0){ #all population do not have the same number of generation
          fullpedigree[which(fullpedigree$popNum==pop & fullpedigree$generation==gen),grep("BreedingValue",colnames(fullpedigree))]=
            fullpedigree[which(fullpedigree$popNum==pop & fullpedigree$generation==gen),grep("BreedingValue",colnames(fullpedigree))]+ #retrieve mean BV
            rmvnorm(n = length(which(fullpedigree$popNum==pop & fullpedigree$generation==gen)), mean=rep(0,nbTraits), sigma=matrix(GAfterRand[pop,]/2,nbTraits,nbTraits))* #draw MVT using GAfterRand to add to mean BV
            matrix(rep(t(sqrt(as.matrix(varscale[which(fullpedigree$popNum==pop & fullpedigree$generation==gen)]))),nbTraits),ncol=nbTraits) # previously draw MVT deviation is scaled using mean F of parents (varscale)
          }
        }
    }
    fullpedigreeSample[,,nsample]=as.matrix(fullpedigree[,grep("BreedingValue",colnames(fullpedigree))])
    setTxtProgressBar(pb,cpt_progressBar) 
  }
  ### Shape output results
  nam= do.call(rbind,lapply(strsplit(as.character(fullpedigree[,1]),"_",fixed=TRUE), function(x){x[3]}))
  dimnames(fullpedigreeSample)=list(nam,colnames(fullpedigree)[grep("BreedingValue",colnames(fullpedigree))])
  random_bv=list()
  if(!is.null(animalToExtractByPop)){
    for (pop in 1:nbpop) {
      random_bv[[pop]]=fullpedigreeSample[match(animalToExtractByPop[[pop]],rownames(fullpedigreeSample)),,,drop=FALSE]
    }
  }else{
    for (pop in 1:nbpop) {
      random_bv[[pop]]=fullpedigreeSample[which(fullpedigree$popNum==pop&(match(nam,ped[[pop]][,1]))),,,drop=FALSE]   
    }
  }
  return(random_bv)
}


###########################################################################################################################
##################################### EXAMPLE based on Morrissey et al. (2019)#############################################
###########################################################################################################################

#Code to simulate data from Morrissey et al. (2019) : Morrissey, M. B., Hangartner, S., & Monro, K. (2019). A note on simulating null distributions for G matrix comparisons. Evolution, 73(12), 2512-2517.
makeSimDat<-function(nsires = 100, noff =5, Va=0.5, Ve=0.5){
  # set up parental information in a data frame called d
  # along the lines of a dams-within-sires design
  d<-as.data.frame(expand.grid(1:noff,1:nsires))
  names(d)<-c("id","sire")
  d$dam<-nsires+1:(nsires*noff)
  d$id<-(nsires+nsires*noff)+1:(nsires*noff)
  # the full pedigree needs the mums and dads bolted
  # on to the top:
  si<-data.frame(id=1:(nsires),sire=NA,dam=NA)
  da<-data.frame(id=d$dam,sire=NA,dam=NA)
  ped<-rbind(si,da,d[,c("id","sire","dam")])
  # clone the experimental design as represented in d
  # into two structures,to contain phenotypes
  # from each of two simulated datasets from separate
  # populations
  d1<-d2<-d
  # generate breeding values using the rbf() function
  # from MCMCglmm, then compose phenotypes
  a1<-rbv(ped,Va)
  a2<-rbv(ped,Va)
  d1$a<-a1[d1$id]
  d2$a<-a2[d2$id]
  d1$z<-d1$a+rnorm(nsires*noff,0,sqrt(Ve))
  d2$z<-d2$a+rnorm(nsires*noff,0,sqrt(Ve))
  # return phenotypic data frames from the dams-within-
  # sires study in each of the two popualtions, and a
  # pedigree object that serves for both datasets
  return(list(d1=d1,d2=d2,ped=ped))
}

sim<-makeSimDat()

p<-list(G=list(G1=list(V=1,nu=0.002)),R=list(V=1,nu=0.002))
m1<-MCMCglmm(z~1,random=~sire,data=sim$d1,
             nitt=130000,thin=100,burnin=30000,prior=p,verbose=FALSE)
m2<-MCMCglmm(z~1,random=~sire,data=sim$d2,
             nitt=130000,thin=100,burnin=30000,prior=p,verbose=FALSE)

#Va=4*Vsire
Va_pop1=4*m1$VCV[,1]
Va_pop2=4*m2$VCV[,1]
#pedigree are similar for both populations
ped1=sim$ped
ped2=sim$ped


NullmodelBV= simul_Null_Gmatrix(ped=list(ped1,ped2),G=list(Va_pop1,Va_pop2),animalToExtractByPop = list(sim$d1$id,sim$d1$id))

#then phenotype can be computed by adding random and fixed part of the model animal model 
# For the morrissey design Ve=Vp-Va
# Ve_m1 = (m1$VCV[,2]+m1$VCV[,1])-(4*m1$VCV[,1])
