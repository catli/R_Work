library(MASS)
#####################################################################################
####################################################################################
#DEFINE FUNCTIONS HERE

betasimul<-function(betamastS,betarefS,betanamesS,wkday,monthlistS,monthnoS,stationlistS,stationnoS,nsimulS) {

    betadata<-betamastS[betarefS[wkday]:(betarefS[wkday+1]-1),1:(betarefS[wkday+1]-betarefS[wkday]-1)] 
    colnames(betadata)<-betanamesS[wkday,1:(betarefS[wkday+1]-betarefS[wkday]-1)]

    betamu<-as.numeric(betadata[1,])
    betacov<-as.numeric(betadata[2,])
        for (i in 3:dim(betadata)[1]) {
            betacov<-rbind(betacov,as.numeric(betadata[i,])) }
    colnames(betacov)<-colnames(betadata)
    rownames(betacov)<-colnames(betadata)
  
    montebeta<-mvrnorm(nsimulS,betamu,betacov)

    for (i in 1:length(monthnoS))                                                 #change the names of the month and station variables to make it easier to reference 
        { colnames(montebeta)[which(colnames(montebeta)==monthlistS[i])]=monthnoS[i] }

    stationmiss<-rep(0,66-50)                                                       
    for (i in 1:length(stationnoS)) {                                              #find the reference position for missing stations
        if(length(which(colnames(montebeta)==stationlistS[i]))<1)  {stationmiss[i]=1}
        colnames(montebeta)[which(colnames(montebeta)==stationlistS[i])]=stationno[i] 
        }

    refmiss<-which(stationmiss==1)

    for (i in 1:length(refmiss)) {                                                #add the missing stations back in as fixed effects
        stationfei<-stationfe[wkday,1+refmiss[i]]
        stationadd<-rep(stationfei,dim(montebeta)[1])
        if(refmiss[i]==1) {refadd<-which(colnames(montebeta)==12)}
        else {refadd<-which(colnames(montebeta)==49+refmiss[i])}
        montebeta<-cbind(montebeta[,1:refadd],stationadd,montebeta[,(refadd+1):dim(montebeta)[2]])
        colnames(montebeta)[refadd+1]<-50+refmiss[i] 
        }
    montebeta
    }   

forecastyrvar<-function(forecastbaseS,montesdcovS,matrixcovS,growthbaseS,yrS,inflationrefS)  {                            #function to select out forecast year data

    forecastsim<-forecastbaseS
    for (i in 1:dim(forecastsim)[1]) {                                             #function to simulate forecast data
        forecastsim[i,3:dim(forecastsim)[2]]<-forecastsim[i,3:dim(forecastsim)[2]]+montesdcovS[1,]} #INCLUDE *sqrt(diag(matrixcovS))?
 
    curmos<-which(forecastsim[,which(colnames(forecastsim)=="year")]==yrS)     #reference for the months selected for montecarlo 
    prevmos<-which(growthbaseS$year==yrS-1)
    prevmos<-prevmos[(1+length(prevmos)-length(curmos)):length(prevmos)]
    
    forecastnew<-matrix(0,ncol=dim(forecastsim)[2]-2,nrow=length(curmos))      #create the forecast data for the selected year
    forecastnew<-growthbaseS[prevmos,3:dim(growthbaseS)[2]]*(1+forecastsim[curmos,3:dim(forecastsim)[2]])
    forecastnew<-cbind(forecastsim[curmos,1:2],forecastnew)
    
    cangasr<-forecastnew$cangas*inflationrefS/forecastnew$inflation
    forecastnew$cangas<-cangasr
    
    bcemplag<-growthbaseS[prevmos-1,3]*(1+forecastsim[curmos-1,3])              #create lag variables for bcemp and cangas
    inflationlag<-growthbaseS[prevmos-1,8]*(1+forecastsim[curmos-1,8])
    cangaslag<-growthbaseS[prevmos-1,7]*(1+forecastsim[curmos-1,7])*inflationrefS/inflationlag
    
    while (length(bcemplag)<dim(forecastnew)[1]) {                              #delete first observation if there is a mismatch between lag and unlagged observations
        forecastnew<-forecastnew[(1+dim(forecastnew)[1]-length(bcemplag)):dim(forecastnew)[1],] }     
        
    refmonth<-forecastnew[1,2]                                                  #add fixed variables - service and nominal fare
    fixdatasm<-fixdata[which(fixdata[,1]==yrS),]
    fixdatasm<-fixdatasm[which(fixdatasm[,2]==refmonth):dim(fixdatasm)[1],]

    forecastnew<-cbind(forecastnew,bcemplag,cangaslag,fixdatasm[,3:4])
   
    fareadj<-forecastnew$farenom*inflationref/forecastnew$inflation             #add adjusted fare
    
    forecastnew<-cbind(forecastnew,fareadj)
    forecastnew[which(forecastnew<0)]<-.001                                     #make sure there are no negative numbers
    
    forecastnew[,3:dim(forecastnew)[2]]<-log(forecastnew[,3:dim(forecastnew)[2]]) #log variables   
    
    forecastnew }

matchvar<-function(forecastnewS,montebetaS,varkeyS,stationk) {          #function to select out forecast year data - input forecast for relevant year and the simulated model coefficients
        
    forecastkey<-colnames(montebetaS)[1:(which(colnames(montebetaS)==1)-1)]     #forecastkey: variable names of the relevant model
    covmatrix<-matrix(0,nrow=dim(forecastnewS)[1],ncol=length(forecastkey))
    colnames(covmatrix)<-forecastkey                                            #create covariate matrix that matches the model for relevant year
    for (i in 1:length(forecastkey)) {
        ref<-which(varkeyS==forecastkey[i])
        if(length(ref)<1) next
        covmatrix[,i]<-forecastnewS[,ref] }

    monthvar<-matrix(0,ncol=11,nrow=dim(forecastnewS)[1])                        # add month dummies
    for (i in 1:10) {
      monthvar[,i]=as.numeric(forecastnewS$month==i) }
      monthvar[,11]=as.numeric(forecastnewS$month==12)
    
    stationpck<-rep(0,66-50)                                                     # add station dummies, will change with loop
    stationpck[stationk-50]=1                
    stationvar<-t(matrix(rep(stationpck,dim(covmatrix)[1]),ncol=dim(covmatrix)[1]))
    
    covmatrix<-cbind(covmatrix,monthvar,stationvar,1)
    
    covmatrix  }

###################################################################
#####READ DATA AND CREATE GROWTH COVARIANCE MATRIX FOR SIMULATION

#COEFFICIENTS
betamon<-as.matrix(read.delim("matrixmon.csv",sep=",",head=TRUE))                          #read estimated coeffiicent data
betatues<-as.matrix(read.delim("matrixtues.csv",sep=",",head=TRUE))
betawed<-as.matrix(read.delim("matrixwed.csv",sep=",",head=TRUE))
betathurs<-as.matrix(read.delim("matrixthurs.csv",sep=",",head=TRUE))
betafri<-as.matrix(read.delim("matrixfri.csv",sep=",",head=TRUE))
betasat<-as.matrix(read.delim("matrixsat.csv",sep=",",head=TRUE))
betasun<-as.matrix(read.delim("matrixsun.csv",sep=",",head=TRUE))

betaref<-c(1,dim(betamon)[1],dim(betatues)[1],dim(betawed)[1],dim(betathurs)[1],dim(betafri)[1],dim(betasat)[1],dim(betasun)[1])
dimmax<-max(betaref)-1
    for (i in 2:length(betaref)) { betaref[i]<-sum(betaref[(i-1):i]) }

betamast<-matrix(0,nrow=max(betaref)-1,ncol=dimmax)                               #collect all coefficients into one matrix 
	⁃	betanames<-matrix(0,nrow=7,ncol=dimmax)                                         #collect all col names into matrix
betamast[betaref[1]:(betaref[2]-1),1:dim(betamon)[2]]<-betamon[,]
betanames[1,1:dim(betamon)[2]]<-colnames(betamon)
betamast[betaref[2]:(betaref[3]-1),1:dim(betatues)[2]]<-betatues[,]
betanames[2,1:dim(betatues)[2]]<-colnames(betatues)
betamast[betaref[3]:(betaref[4]-1),1:dim(betawed)[2]]<-betawed[,]
betanames[3,1:dim(betawed)[2]]<-colnames(betawed)
betamast[betaref[4]:(betaref[5]-1),1:dim(betathurs)[2]]<-betathurs[,]
betanames[4,1:dim(betathurs)[2]]<-colnames(betathurs)
betamast[betaref[5]:(betaref[6]-1),1:dim(betafri)[2]]<-betafri[,]
betanames[5,1:dim(betafri)[2]]<-colnames(betafri)
betamast[betaref[6]:(betaref[7]-1),1:dim(betasat)[2]]<-betasat[,]
betanames[6,1:dim(betasat)[2]]<-colnames(betasat)
betamast[betaref[7]:(betaref[8]-1),1:dim(betasun)[2]]<-betasun[,]
betanames[7,1:dim(betasun)[2]]<-colnames(betasun)


stationfe<-read.delim("stationfe.csv",sep=",",head=TRUE)


#GROWTH COVARIANCE
covdata<-read.delim("Covdistribution.csv",header=TRUE, sep=",")                 #read covariance matrix

buildperm<-covdata$buildpermits
for (i in 36:length(buildperm)) {
    buildperm[i]<-mean(covdata$buildpermits[(i-35):i]) }

basec<-c(1,2,5:9,11)
growthdata<-covdata[,basec]
growthdata$buildpermits<-buildperm
ref<-which(is.na(rowSums(growthdata))==FALSE)
growthdata<-growthdata[ref,]


bcempperc<-rep(0,which(growthdata$year==2012)[1]-12)                            #filter for historical data
gasperc<-rep(0,which(growthdata$year==2012)[1]-12)
tourperc<-rep(0,which(growthdata$year==2012)[1]-12)
accomoperc<-rep(0,which(growthdata$year==2012)[1]-12)
buildperc<-rep(0,which(growthdata$year==2012)[1]-12)
infperc<-rep(0,which(growthdata$year==2012)[1]-12)

for (i in 13:which(growthdata$year==2012)[1]) {                                 #find the base percentage growth
  bcempperc[i]<-growthdata$bcemp[i]/growthdata$bcemp[i-12]-1  
  gasperc[i]<-growthdata$usgasprice_n[i]/growthdata$usgasprice_n[i-12]-1  
  tourperc[i]<-growthdata$cantourist[i]/growthdata$cantourist[i-12]-1  
  accomoperc[i]<-growthdata$vancaccomo[i]/growthdata$vancaccomo[i-12]-1  
  buildperc[i]<-growthdata$buildpermits[i]/growthdata$buildpermits[i-12]-1  
  infperc[i]<-growthdata$inflation[i]/growthdata$inflation[i-12]-1 
   }
growthmatrix<-cbind(bcempperc,tourperc,accomoperc,buildperc,gasperc,infperc)
matrixcov<-cov(growthmatrix)                                                    #find covariance of growth


#BASE GROWTH
basec<-c(1,2,5:8,10,11)
growthbase<-covdata[,basec]                       
ref<-which(is.na(rowSums(growthbase))==FALSE)
growthbase<-growthbase[ref,]

forecastbase<-matrix(0,ncol=dim(growthbase)[2],nrow=dim(growthbase)[1]-which(growthbase$year==2012)[1]+1)   
rstart<-which(growthbase$year==2012)[1]
forecastbase[,1]<-growthbase[(rstart):dim(growthbase)[1],1]
forecastbase[,2]<-growthbase[(rstart):dim(growthbase)[1],2]

for (i in 1:dim(forecastbase)[1]) {                                             #generate base forecast
    for (j in 3:dim(forecastbase)[2]) {
        forecastbase[i,j]<-growthbase[rstart+i,j]/growthbase[rstart+i-12,j]-1 
        }
    }
colnames(forecastbase)<-colnames(growthbase)

# FIXED VECTOR
fixdata<-covdata[,c(1:4)]                                                       #generate fixed vectors of service and nominal fare
ref<-which(is.na(rowSums(fixdata))==FALSE)
fixdata<-fixdata[ref,]
rstart<-which(fixdata$year==2012)[1]
fixdata<-fixdata[(rstart):dim(fixdata)[1],]

#STANDARD ERROR
stderdata<-as.numeric(read.delim("standarderror.csv",sep=",",head=TRUE))           #mc for standard error


#DEFINE REFERENCE OBJECTS
monthlist<-c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X12")           
monthno<-c(1:10,12)
stationlist<-c("X51","X52","X53","X54","X55","X56","X57","X58","X59","X60","X61","X62","X63","X64","X65","X66")
stationno<-c(51:66)

inflationref<-covdata$inflation[which(covdata$year==2011)[12]]       
yr<-2012 
varkey<-c("year","month","bcemplog","touristlog","accomolog","buildlog","cangaslog","inflation","bcemplaglog","cangaslaglog","farenom","servhat","fare")  


    
################################################################################
##MONTECARLO SIMULATION#########################################################
nsimul=10000

yrlist<-c(2012,2013,2014,2015,2020,2025,2030,2040)
yrmolist<-c(11,rep(12,7))

output2012<-matrix(0,ncol=16*11*7,nrow=nsimul)
colnames(output2012)<-rep(rep(51:66,each=dim(output2012)[2]/7/16),7)
output2013<-matrix(0,ncol=16*12*7,nrow=nsimul)
colnames(output2013)<-rep(rep(51:66,each=dim(output2013)[2]/7/16),7)
output2040=output2030=output2025=output2020=output2015=output2014=output2013

outputmast<-cbind(output2012,output2013,output2014,output2015,output2020,output2025,output2030,output2040)
refout<-c(1,dim(output2012)[2],dim(output2013)[2],dim(output2014)[2],dim(output2015)[2],dim(output2020)[2],dim(output2025)[2],dim(output2030)[2],dim(output2040)[2])
    for (i in 2:length(refout)) {refout[i]<-sum(refout[(i-1):i]) }
        
for (d in 1:dim(betanames)[1]) {   
    montebeta<-betasimul(betamast,betaref,betanames,d,monthlist,monthno,stationlist,stationno,nsimul)        #mc for model coefficient 
    montesdcov<-mvrnorm(nsimul,rep(0,dim(matrixcov)[1]),matrixcov)                     #mc for covariate growth standard error
    
    for (y in 1:length(yrlist)) {      
        yr<-yrlist[y] 
        montestder<-rnorm(nsimul,0,1-exp(-mean(stderdata)*sqrt(yr-2011)))
        montestder[which((montestder+1)<0)]=-1
        forecastnew<-forecastyrvar(forecastbase,montesdcov,matrixcov,growthbase,yr,inflationref)     #select relevant forecast year for growth
        refyr<-which(yrlist==yr)
        refmo<-yrmolist[refyr]
        outputmat<-matrix(0,ncol=16*refmo,nrow=nsimul)

        for (stk in 51:66) {
           covmatrix<-matchvar(forecastnew,montebeta,varkey,stk)                        
           yhat<-montebeta%*%t(covmatrix)                                 #only for one year    
           forecast<-log(exp(yhat)*(1+montestder))
           outputmat[,(refmo*(stk-51)+1):(refmo*(stk-50))]<-forecast   
           }
        refcol<-refout[refyr]+(refmo*16)*(d-1)
        outputmast[,refcol:(refcol+(refmo*16)-1)]<-outputmat
        }
    }
                       
#output2012<-outputmast[,refout[1]:(refout[2]-1)]
#output2013<-outputmast[,refout[2]:(refout[3]-1)]
#output2014<-outputmast[,refout[3]:(refout[4]-1)]
#output2015<-outputmast[,refout[4]:(refout[5]-1)]
#output2020<-outputmast[,refout[5]:(refout[6]-1)]
#output2025<-outputmast[,refout[6]:(refout[7]-1)]
#output2030<-outputmast[,refout[7]:(refout[8]-1)]
#output2040<-outputmast[,refout[8]:(refout[9]-1)]

#write.table(output2012,file="output2012.csv",col.names=TRUE,row.names=TRUE,sep=",")
#write.table(output2013,file="output2013.csv",col.names=TRUE,row.names=TRUE,sep=",")
#write.table(output2014,file="output2014.csv",col.names=TRUE,row.names=TRUE,sep=",")
#write.table(output2015,file="output2015.csv",col.names=TRUE,row.names=TRUE,sep=",")
#write.table(output2020,file="output2020.csv",col.names=TRUE,row.names=TRUE,sep=",")
#write.table(output2025,file="output2025.csv",col.names=TRUE,row.names=TRUE,sep=",")
#write.table(output2030,file="output2030.csv",col.names=TRUE,row.names=TRUE,sep=",")
#write.table(output2040,file="output2040.csv",col.names=TRUE,row.names=TRUE,sep=",")

##############OUTPUTS###########################################################
###############################################################################
#DISTRIBUTION TABLES
futdays<-read.delim("futdaystot.csv",header=TRUE, sep=",")       
annualoutput<-matrix(0,ncol=length(yrlist),nrow=nsimul)   
colnames(annualoutput)<-yrlist

for (i in 1:length(yrlist)) { 
    yr=yrlist[i]
    reffut<-which(futdays$year==yr)
    daysmat=0
        for (j in 1:7) { daysmat<-c(daysmat,futdays[reffut,2+j]) }
    daysmat<-daysmat[2:length(daysmat)]
    daysmat<-rep(daysmat,each=16)
    outputi<-exp(outputmast[,refout[i]:(refout[i+1]-1)]) 
    outputtot<-t(t(outputi)*daysmat)
    annualtot<-rowSums(outputtot)   
    annualoutput[,i]<-annualtot }
 
par(mfrow=c(3,3))
for (i in 1:8) {  
    graphdata<-annualoutput[,i] 
    hist(graphdata,breaks=100,main=paste("Histogram of Annual Counts for Year",yrlist[i]))
    }

write.table(annualoutput,file="outputannual.csv",col.names=TRUE,row.names=TRUE,sep=",")


int<-16*11
par(mfrow=c(3,3))
for (i in 1:7) {  
    graphdata<-output2012[,(int*(i-1)+1):(int*i)] 
    graphdata<-rowSums(exp(graphdata))
    hist(graphdata,breaks=100,main=paste("Year 2012 Histogram of Day",i))
    }
    
int<-16*12    
par(mfrow=c(3,3))
for (i in 1:7) {  
    graphdata<-output2013[,(int*(i-1)+1):(int*i)] 
    graphdata<-rowSums(exp(graphdata))
    hist(graphdata,breaks=100,main=paste("Year 2013 Histogram of Day",i))
    }

par(mfrow=c(3,3))
for (i in 1:7) {  
    graphdata<-output2014[,(int*(i-1)+1):(int*i)] 
    graphdata<-rowSums(exp(graphdata))
    hist(graphdata,breaks=100,main=paste("Year 2014 Histogram of Day",i))
    }

par(mfrow=c(3,3))
for (i in 1:7) {  
    graphdata<-output2015[,(int*(i-1)+1):(int*i)] 
    graphdata<-rowSums(exp(graphdata))
    hist(graphdata,breaks=100,main=paste("Year 2015 Histogram of Day",i))
    }

par(mfrow=c(3,3))
for (i in 1:7) {  
    graphdata<-output2020[,(int*(i-1)+1):(int*i)] 
    graphdata<-rowSums(exp(graphdata))
    hist(graphdata,breaks=100,main=paste("Year 2020 Histogram of Day",i))
    }

par(mfrow=c(3,3))
for (i in 1:7) {  
    graphdata<-output2030[,(int*(i-1)+1):(int*i)] 
    graphdata<-rowSums(exp(graphdata))
    hist(graphdata,breaks=100,main=paste("Year 2030 Histogram of Day",i))
    }
    
par(mfrow=c(3,3))
for (i in 1:7) {  
    graphdata<-output2040[,(int*(i-1)+1):(int*i)] 
    graphdata<-rowSums(exp(graphdata))
    hist(graphdata,breaks=100,main=paste("Year 2040 Histogram of Day",i))
    }
    