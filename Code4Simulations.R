library(CLEAN)
library(doParallel)
library(foreach)
library(Matrix)
library(data.table)
library(ciftiTools)
library(Rcpp)
library(RcppArmadillo)
args=(commandArgs(TRUE))
job_num=as.numeric(args[1])
path=args[2] #path to files
cluster=as.logical(args[3]) #set cluster enhancement or not
max.radius=as.numeric(args[4]) #set maximum radius
spatial.cov=as.logical(args[5]) #set spatial autocorrelation modelling or not
ntwins=as.numeric(args[6]) #set number of true twin pairs

Rcpp::sourceCpp(paste0(path,"/Codes/Heritability/Clean_support.cpp"))
source(paste0(path,"/Codes/Heritability/CLEANV.R"))


total <- 120. #N=120 for heritability

ciftiTools.setOption("wb_path", paste0(path,"/workbench"))

#Read the independent subjects behavior information
#Read the twins behavior information
#Read the independent subjects imaging data
#Read the twins imaging data
dif_info <- readRDS(paste0(path,"/Codes/Heritability/Data/dif_info.rds")) 
twins_info <- readRDS(paste0(path,"/Codes/Heritability/Data/Etwins_info.rds"))
twins_data <- readRDS(paste0(path,"/Codes/Heritability/Data/twins_10000.rds"))
dif_data <- readRDS(paste0(path,"/Codes/Heritability/Data/dif_10000.rds"))

#Read the left distance matrix 
Ldismat <- fread(file=paste0(path,"/Data/Dismat/l_exinflated_dismat10000.txt"),sep=" ", data.table =F,header=F)
l_cortex <-which(twins_data$meta$cortex$medial_wall_mask$left)
ldismat <-  unname(as.matrix(Ldismat[l_cortex,l_cortex]))

#Read the right distance matrix 
Rdismat <- fread(file=paste0(path,"/Data/Dismat/r_exinflated_dismat10000.txt"),sep=" ", data.table =F,header=F)
r_cortex <-which(twins_data$meta$cortex$medial_wall_mask$right)
rdismat <- unname(as.matrix(Rdismat[r_cortex, r_cortex]))

##Sample TWINS
set.seed(job_num)
twins_id <- sample(nrow(twins_info)/2, ntwins)
twins_row <- sort(c(twins_id*2-1,twins_id*2))
ldata_twins <-  twins_data$data$cortex_left[,twins_row]
rdata_twins <- twins_data$data$cortex_right[,twins_row]
binfo <- twins_info[twins_row,4:5]

##Construct K matrix based on the images' relationship
KMat <- Matrix(kronecker(diag(1, total/2), matrix(c(1,1,1,1),nrow=2,ncol=2)),sparse=T)
if (ntwins != 0 ) {
  for (i in c(1:ntwins)) {
    if(twins_info[twins_id[i],'ZygosityGT'] == 'DZ') {
      KMat[2*i, 2*i-1] = 0.5
      KMat[2*i-1, 2*i] = 0.5
    } else {
      KMat[2*i, 2*i-1] = 1
      KMat[2*i-1, 2*i] = 1
    }
  }
}

#Sample NON-Twins
#NON-Twins
difsize <- nrow(dif_info)
dif_id <- sample(difsize, total-ntwins*2)
ldata_dif <- dif_data$data$cortex_left[,dif_id]
rdata_dif <- dif_data$data$cortex_right[,dif_id]
binfodif <- dif_info[dif_id,2:3]

#Combine twins sample and non-twins sample
Binfo <- rbind(binfo, binfodif)
ldata <- cbind(ldata_twins, ldata_dif)
rdata <-  cbind(rdata_twins, rdata_dif)
mod0 <- model.matrix(~Binfo$Age_in_Yrs+Binfo$Gender)

#Modelling the data based on settings
lres <- CleanV(ldata, ldistmat, K=KMat, max.radius =max.radius, spatial=spatial.cov,
               mod=mod0,
               nperm = 2000,
               alpha = 0.05,
               ncores = 1)

rres <-CleanV(rdata, rdistmat, K=KMat, max.radius = max.radius, spatial=spatial.cov,
              mod=mod0,
              nperm = 2000,
              alpha = 0.05,
              ncores = 1)

res <- list(left = lres,
            right = rres)

saveRDS(res,file=paste0(path,"/Codes/Heritability/Result/job_num=",job_num,
                        "setting",
                        ntwins,
                        cluster,
                        max.radius,
                        spatial.cov,"power.rds"))


