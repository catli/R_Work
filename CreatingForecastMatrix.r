data<-read.csv("ForecastData.csv", header = TRUE,sep=",")
data<-as.matrix(data[,5:dim(data)[2]])
data2<-matrix(as.numeric(data),ncol=dim(data)[2])

nanum<-rep(0,dim(data2)[2])
for (i in 1:dim(data2)[2]) {
  nanum[i]<-sum(is.na(data2[,i])==TRUE) }
nanum
#data2[which(is.na(data2)==TRUE)]<-0

model<-read.csv("ForecastModel.csv", header = TRUE,sep=",")
model2<-as.matrix(model)
model2[which(is.na(model2)==TRUE )]<-0
 
forecast<-data2%*%model2

write(t(forecast),file="forecast.csv", sep=",",ncol=dim(forecast)[2])