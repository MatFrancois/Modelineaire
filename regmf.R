# don<- readRDS("starting_k/data.rds")
# 
# gs <- c("BRAF",          # Tumor supressors
#   "FASTKD1", "ND2", "ND3", "ND4"   # Mitochondiral activity
# )
# colnames(don)
# reg.mf<-lm(BRAF~age+sex+tissue_status+project+tissue,data=don)
# step(reg.mf)
# reg.mf<-lm(formula = BRAF ~ age + tissue_status + project + tissue, data = don)
# summary(reg.mf)
# table(df.nona$tissue_status)
# df.nona<-don[is.na(don$t)==FALSE &
#                is.na(don$n)==FALSE &
#                is.na(don$m)==FALSE &
#                is.na(don$tnm_grade)==FALSE &
#                is.na(don$tnm_stage)==FALSE,]
# reg.mf<-lm(BRAF~t+n+m+project+tissue,data=df.nona)
# summary(reg.mf)
# step(reg.mf)
# reg.step<-lm(formula = BRAF ~ m + tissue, data = df.nona)
# summary(reg.step)

######
######Prédiction des gênes en fonctions des autres gênes
######

###Copie du jeu de données
don<- readRDS("starting_k/data.rds")
don.imput<-don##jeu de données sans nans
gs <- c(
  "ATAD2", "SYCP3",                # Chromatin 
  "H19", "IGF2", "NNAT", "BLCAP",  # inprinted genes
  "BRD4", "BRDT", "NUTM1",         # Testis restricted
  "MAGEB6", "TUBA3C",              # Testis specific
  "SMYD3", "MAP3K2", "KDR",        # K methyl transferase
  "TP53", "KRAS", "BRAF",          # Tumor supressors
  "FASTKD1", "ND2", "ND3", "ND4"   # Mitochondiral activity
)

##prediction des valeurs manquantes par un simple lm en fonction des variables biologiques
pred <- sapply(gs, function(g, don.imput){
  # m = lm(don.imput[[g]] ~ d$tissue)
  m<-lm(don.imput[[g]] ~ don.imput$age+
           don.imput$sex+
           don.imput$tissue+
           don.imput$tissue_status+
           don.imput$project)
  m<-step(m)
  predict(m, don.imput)
}, don.imput)

summary(pred)
r2<-NULL
for (i in gs){
  m<-lm(don.imput[[i]] ~ don.imput$age+
          don.imput$sex+
          don.imput$tissue+
          don.imput$tissue_status+
          don.imput$project)
  r2<-c(r2,summary(m)$r.squared)
}
##imputation des na par la valeur prédite
for (i in 1:length(gs)){
  don.imput[is.na(don.imput[,gs[i]]),gs[i]]<-pred[is.na(don.imput[,gs[i]]),gs[i]]
}

##suppression des valeurs de merdes
don<-don[,!colnames(don)%in%c("t","n","m","tnm_stage","tnm_grade")]
# don.imput[,gs]<-pred##remplacement des valeurs prédites
don.imput<-don.imput[,!colnames(don.imput)%in%c("t","n","m","tnm_stage","tnm_grade")]
# res.pred <- sapply(gs, function(g){
#   m<-step(lm(don.imput[[g]] ~ .,data=don.imput))
#   predict(m)
# })

##modele pour chaque gene
require(lmtest)##pour hypothese d'homogénéité des variances
pred1<-step(lm(don.imput$ATAD2~.,data=don.imput))
shapiro.test(sample(pred1$residuals,5000)) #p-value < 0.05
gqtest(pred1) #p-value = 1
hist(pred1$residuals)

pred2<-step(lm(don.imput$SYCP3~.,data=don.imput))
shapiro.test(sample(pred2$residuals,5000)) #p-value < 0.05
gqtest(pred2) #p-value = 1
hist(pred2$residuals)  

pred3<-step(lm(don.imput$BRDT~.,data=don.imput))
shapiro.test(sample(pred3$residuals,5000)) #p-value < 0.05
gqtest(pred3) #p-value = 1
hist(pred3$residuals)  

pred4<-step(lm(don.imput$BRD4~.,data=don.imput))
shapiro.test(sample(pred4$residuals,5000)) #p-value < 0.05
gqtest(pred4) #p-value = 0.69
hist(pred4$residuals)  

pred5<-step(lm(don.imput$NUTM1~.,data=don.imput))
shapiro.test(sample(pred5$residuals,5000)) #p-value < 0.05
gqtest(pred5) #p-value = 1
hist(pred5$residuals)  

pred6<-step(lm(don.imput$MAGEB6~.,data=don.imput))
shapiro.test(sample(pred6$residuals,5000)) #p-value < 0.05
gqtest(pred6) #p-value = 1
hist(pred6$residuals)  

pred7<-step(lm(don.imput$TUBA3C~.,data=don.imput))
shapiro.test(sample(pred7$residuals,5000)) #p-value < 0.05
gqtest(pred7) #p-value = 1
hist(pred7$residuals)  

pred8<-step(lm(don.imput$H19~.,data=don.imput))
shapiro.test(sample(pred8$residuals,5000)) #p-value < 0.05
gqtest(pred8) #p-value = 0.00014
hist(pred8$residuals)  

pred9<-step(lm(don.imput$IGF2~.,data=don.imput))
shapiro.test(sample(pred9$residuals,5000)) #p-value < 0.05
gqtest(pred9) #p-value = 0.38
hist(pred9$residuals)  

pred10<-step(lm(don.imput$NNAT~.,data=don.imput))
shapiro.test(sample(pred10$residuals,5000)) #p-value < 0.05
gqtest(pred10) #p-value = 0.21
hist(pred10$residuals)  

pred11<-step(lm(don.imput$BLCAP~.,data=don.imput))
shapiro.test(sample(pred11$residuals,5000)) #p-value < 0.05
gqtest(pred11) #p-value = 1
hist(pred11$residuals)  

pred12<-step(lm(don.imput$SMYD3~.,data=don.imput))
shapiro.test(sample(pred12$residuals,5000)) #p-value < 0.05
gqtest(pred12) #p-value = 1
hist(pred12$residuals)  

pred13<-step(lm(don.imput$MAP3K2~.,data=don.imput))
shapiro.test(sample(pred13$residuals,5000)) #p-value < 0.05
gqtest(pred13) #p-value = 1
hist(pred13$residuals)  

pred14<-step(lm(don.imput$KDR~.,data=don.imput))
shapiro.test(sample(pred14$residuals,5000)) #p-value < 0.05
gqtest(pred14) #p-value < 0.05
hist(pred14$residuals)  

pred15<-step(lm(don.imput$TP53~.,data=don.imput))
shapiro.test(sample(pred15$residuals,5000)) #p-value < 0.05
gqtest(pred15) #p-value = 1
hist(pred15$residuals)  

pred16<-step(lm(don.imput$KRAS~.,data=don.imput))
shapiro.test(sample(pred16$residuals,5000)) #p-value < 0.05
gqtest(pred16) #p-value = 0.001469
hist(pred16$residuals)  

pred17<-step(lm(don.imput$BRAF~.,data=don.imput))
shapiro.test(sample(pred17$residuals,5000)) #p-value < 0.05
gqtest(pred17) #p-value = 1
hist(pred17$residuals)  

pred18<-step(lm(don.imput$FASTKD1~.,data=don.imput))
shapiro.test(sample(pred18$residuals,5000)) #p-value < 0.05
gqtest(pred18) #p-value = 0.5255
hist(pred18$residuals)  

pred19<-step(lm(don.imput$ND2~.,data=don.imput))
shapiro.test(sample(pred19$residuals,5000)) #p-value < 0.05
gqtest(pred19) #p-value = 1
hist(pred19$residuals)  

pred20<-step(lm(don.imput$ND3~.,data=don.imput))
shapiro.test(sample(pred20$residuals,5000)) #p-value < 0.05
gqtest(pred20) #p-value = <0.05
hist(pred20$residuals)  

pred21<-step(lm(don.imput$ND4~.,data=don.imput))
shapiro.test(sample(pred21$residuals,5000)) #p-value < 0.05
gqtest(pred21) #p-value = 0.2
hist(pred21$residuals)  

##prediction
res.predM<-cbind(predict(pred1),predict(pred2),predict(pred3),predict(pred4)
                ,predict(pred5),predict(pred6),predict(pred7)
                ,predict(pred8),predict(pred9),predict(pred10),predict(pred11)
                ,predict(pred12),predict(pred13),predict(pred14),predict(pred15)
                ,predict(pred16),predict(pred17),predict(pred18),predict(pred19)
                ,predict(pred20),predict(pred21))

##rendu
colnames(res.predM)<-colnames(don)[6:26]
gs<-colnames(don)[6:26]
res.predM[!is.na(as.matrix(don[,gs]))] = NA # Sparse, only submit NA values!
saveRDS(res.predM, "results.rds")
zip_filename = paste(sep="",  "results_", format(Sys.time(), format="%m_%d_%Y_%s"), ".zip")
zip(zip_filename, "results.rds")
print(zip_filename)

#####
#####evolution R2

mod2<-c(summary(pred1)$r.squared,summary(pred2)$r.squared,summary(pred3)$r.squared,summary(pred4)$r.squared
      ,summary(pred5)$r.squared,summary(pred6)$r.squared,summary(pred7)$r.squared
      ,summary(pred8)$r.squared,summary(pred9)$r.squared,summary(pred10)$r.squared,summary(pred11)$r.squared
      ,summary(pred12)$r.squared,summary(pred13)$r.squared,summary(pred14)$r.squared,summary(pred15)$r.squared
      ,summary(pred16)$r.squared,summary(pred17)$r.squared,summary(pred18)$r.squared,summary(pred19)$r.squared
      ,summary(pred20)$r.squared,summary(pred21)$r.squared)

df.r2<-(rbind(r2,mod2))
colnames(df.r2)<-gs
gene<-gs
df.r2<-data.frame(t(df.r2),gene)
save(df.r2,file = "evor2.RData")
plot(df.r2$r2[order(df.r2$mod2)],col="blue",ylim = c(min(c(df.r2$r2,df.r2$mod2)),max(c(df.r2$r2,df.r2$mod2))+0.2),pch=16,cex=df.r2$r2[order(df.r2$mod2)],ylab="R²")
text(df.r2$r2[order(df.r2$mod2)],labels = df.r2$gene[order(df.r2$mod2)],col="blue",cex=df.r2$r2[order(df.r2$mod2)],adj = c(0.5,-1))
points(sort(df.r2$mod2),col="red",pch=16,cex=sort(df.r2$mod2))
text(sort(df.r2$mod2),labels = df.r2$gene[order(df.r2$mod2)],col="red",cex=sort(df.r2$mod2),adj = c(0.5,-1))


####
####matrice des parametres de chaque modele
param<-function(obj){
  noms<-NULL
  for (i in 3:length(attr(obj$terms,"variables"))){
    noms<-c(noms,paste(attr(obj$terms,"variables")[[i]]))
  }
  return(noms)
}
###matrice
gs<-c("sex","tissue_status","project","tissue",gs)
param.mat<-matrix(nrow = length(gs),ncol = length(gs))


param.mat<-rbind(ifelse(gs%in%param(pred1),"x",""),
      ifelse(gs%in%param(pred2),"x",""),
      ifelse(gs%in%param(pred3),"x",""),
      ifelse(gs%in%param(pred4),"x",""),
      ifelse(gs%in%param(pred5),"x",""),
      ifelse(gs%in%param(pred6),"x",""),
      ifelse(gs%in%param(pred7),"x",""),
      ifelse(gs%in%param(pred8),"x",""),
      ifelse(gs%in%param(pred9),"x",""),
      ifelse(gs%in%param(pred10),"x",""),
      ifelse(gs%in%param(pred11),"x",""),
      ifelse(gs%in%param(pred12),"x",""),
      ifelse(gs%in%param(pred13),"x",""),
      ifelse(gs%in%param(pred14),"x",""),
      ifelse(gs%in%param(pred15),"x",""),
      ifelse(gs%in%param(pred16),"x",""),
      ifelse(gs%in%param(pred17),"x",""),
      ifelse(gs%in%param(pred18),"x",""),
      ifelse(gs%in%param(pred19),"x",""),
      ifelse(gs%in%param(pred20),"x",""),
      ifelse(gs%in%param(pred21),"x",""))
row.names(param.mat)<-gs[-c(1:4)]
colnames(param.mat)<-gs
