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
shapiro.test(sample(pred1$residuals,5000))
gqtest(pred1)
hist(pred1$residuals)
pred2<-step(lm(don.imput$SYCP3~.,data=don.imput))
pred3<-step(lm(don.imput$BRDT~.,data=don.imput))
pred4<-step(lm(don.imput$BRD4~.,data=don.imput))
pred5<-step(lm(don.imput$NUTM1~.,data=don.imput))
pred6<-step(lm(don.imput$MAGEB6~.,data=don.imput))
pred7<-step(lm(don.imput$TUBA3C~.,data=don.imput))
pred8<-step(lm(don.imput$H19~.,data=don.imput))
pred9<-step(lm(don.imput$IGF2~.,data=don.imput))
pred10<-step(lm(don.imput$NNAT~.,data=don.imput))
pred11<-step(lm(don.imput$BLCAP~.,data=don.imput))
pred12<-step(lm(don.imput$SMYD3~.,data=don.imput))
pred13<-step(lm(don.imput$MAP3K2~.,data=don.imput))
pred14<-step(lm(don.imput$KDR~.,data=don.imput))
pred15<-step(lm(don.imput$TP53~.,data=don.imput))
pred16<-step(lm(don.imput$KRAS~.,data=don.imput))
pred17<-step(lm(don.imput$BRAF~.,data=don.imput))
pred18<-step(lm(don.imput$FASTKD1~.,data=don.imput))
pred19<-step(lm(don.imput$ND2~.,data=don.imput))
pred20<-step(lm(don.imput$ND3~.,data=don.imput))
pred21<-step(lm(don.imput$ND4~.,data=don.imput))

##prediction
res.pred<-cbind(predict(pred1),predict(pred2),predict(pred3),predict(pred4)
                ,predict(pred5),predict(pred6),predict(pred7)
                ,predict(pred8),predict(pred9),predict(pred10),predict(pred11)
                ,predict(pred12),predict(pred13),predict(pred14),predict(pred15)
                ,predict(pred16),predict(pred17),predict(pred18),predict(pred19)
                ,predict(pred20),predict(pred21))

##rendu
colnames(res.pred)<-colnames(don)[6:26]
gs<-colnames(don)[6:26]
res.pred[!is.na(as.matrix(don[,gs]))] = NA # Sparse, only submit NA values!
saveRDS(res.pred, "results.rds")
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
df.r2<-data.frame(t(df.r2))
plot(df.r2)
df.r2$gene<-gs
df.new<-c(df.r2$r2,df.r2$mod2)
df.new<-data.frame(r2=df.new,gene=c(df.r2$gene,df.r2$gene))
