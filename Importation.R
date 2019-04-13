####Importation

don<- readRDS("starting_k/data.rds")

gs <- c(
  "ATAD2", "SYCP3",                # Chromatin 
  "H19", "IGF2", "NNAT", "BLCAP",  # inprinted genes
  "BRD4", "BRDT", "NUTM1",         # Testis restricted
  "MAGEB6", "TUBA3C",              # Testis specific
  "SMYD3", "MAP3K2", "KDR",        # K methyl transferase
  "TP53", "KRAS", "BRAF",          # Tumor supressors
  "FASTKD1", "ND2", "ND3", "ND4"   # Mitochondiral activity
)

pred <- sapply(gs, function(g, don){
  # m<-lm(don[[g]] ~ don$tissue)
  m<-lm(don[[g]] ~ don$age+d$sex+don$tissue+don$tissue_status+don$project)
  predict(m, don)
}, don)

res<-matrix(nrow = length(gs),ncol = 1)
colnames(res)<-c("age","sex","tissue","status","project")
rownames(res)<-gs
lescall<-NULL
for (i in gs){
  res<-lm(new.don[[i]] ~ age+sex+tissue+tissue_status+project,data = new.don)
  lescall<-c(lescall,step(res)$call)
}

don.pred<-cbind(predict(lm(formula = new.don[[gs[1]]] ~ age + sex + tissue + tissue_status + 
             project, data = new.don),don),
      predict(lm(formula = new.don[[gs[2]]] ~ tissue + tissue_status + project, 
                 data = new.don),don),
      predict(lm(formula = new.don[[gs[3]]] ~ age + tissue + tissue_status + project, 
                 data = new.don),don),
      predict(lm(formula = new.don[[gs[4]]] ~ age + tissue + tissue_status + project, 
                 data = new.don),don),
      predict(lm(formula = new.don[[gs[5]]] ~ sex + tissue + tissue_status + project, 
                 data = new.don),don),
      predict(lm(formula = new.don[[gs[6]]] ~ age + sex + tissue + tissue_status + 
                   project, data = new.don),don),
      predict(lm(formula = new.don[[gs[7]]] ~ age + tissue + tissue_status + project, 
                 data = new.don),don),
      predict(lm(formula = new.don[[gs[8]]] ~ tissue + tissue_status + project, 
                 data = new.don),don),
      predict(lm(formula = new.don[[gs[9]]] ~ tissue + project, data = new.don),don),
      predict(lm(formula = new.don[[gs[10]]] ~ age + sex + tissue + tissue_status + 
                   project, data = new.don),don),
      predict(lm(formula = new.don[[gs[11]]] ~ age + tissue_status + project, data = new.don),don),
      predict(lm(formula = new.don[[gs[12]]] ~ age + sex + tissue + tissue_status + 
                   project, data = new.don),don),
      predict(lm(formula = new.don[[gs[13]]] ~ tissue_status + project, data = new.don),don),
      predict(lm(formula = new.don[[gs[14]]] ~ sex + tissue + tissue_status + project, 
                 data = new.don),don),
      predict(lm(formula = new.don[[gs[15]]] ~ age + tissue + tissue_status + project, 
                 data = new.don),don),
      predict(lm(formula = new.don[[gs[16]]] ~ tissue + tissue_status + project, 
                 data = new.don),don),
      predict(lm(formula = new.don[[gs[17]]] ~ age + tissue + tissue_status + project, 
                 data = new.don),don),
      predict(lm(formula = new.don[[gs[18]]] ~ tissue + tissue_status + project, 
                 data = new.don),don),
      predict(lm(formula = new.don[[gs[19]]] ~ age + sex + tissue_status + project, 
                 data = new.don),don),
      predict(lm(formula = new.don[[gs[20]]] ~ age + sex + tissue + tissue_status + 
                   project, data = new.don),don),
      predict(lm(formula = new.don[[gs[21]]] ~ age + tissue_status + project, data = new.don),don))
colnames(don.pred)<-gs      
      
don.pred[!is.na(as.matrix(don[,gs]))] = NA # Sparse, only submit NA values!
saveRDS(don.pred, "results.rds")
zip_filename = paste(sep="",  "results_", format(Sys.time(), format="%m_%d_%Y_%s"), ".zip")
zip(zip_filename, "results.rds")
print(zip_filename)
      
      
      
      
      
      
