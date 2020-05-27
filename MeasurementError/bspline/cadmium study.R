
# Data handling -----------------------------------------------------------

library(dplyr)
cadmium = read.csv("./cadmium.csv", header=TRUE)
study = cadmium$studycode   # study code
w = cbind(ifelse(cadmium$ethnicity == 1, 1, 0),
          ifelse(cadmium$age >= 50, 1, 0),
          ifelse(cadmium$gender == 1, 1, 0),
          ifelse(cadmium$gender == 0, 1, 0))
colnames(w)=c('asian','age>=50','female','male')
y = log(cadmium$b2_GM)
x = log(cadmium$Ucd_GM)

# sort the data w.r.t. studycode, to make correspondence with MCMC
cadmium_bind = cbind.data.frame(study,w,y,x) %>% arrange(study)
study = cadmium_bind$study
y = cadmium_bind$y
x = cadmium_bind$x; ord = order(x)
w = cadmium_bind %>% dplyr::select(-study,-y,-x) %>% as.matrix()
label = unique(study)
Z = matrix(0, nrow(w), length(label))
for(idx in 1:length(label)) Z[which(study==label[idx]),idx] = 1


# Random Intercept MER ----------------------------------------------------

source("./randintmer.R")
start = as.numeric(Sys.time())
vb_result = randintmer(y,w,x,Z,n_intknot=2)
print(paste("VB", round(as.numeric(Sys.time()) - start, 3), "seconds elapsed"))
plot(vb_result$lb, main="ELBO", type="l", xlab="Iterations", ylab="")
residual = y-cbind(1,w)%*%vb_result$mubeta.q-Z%*%vb_result$mub.q
plot(x,residual,main="Estimated mean curve",xlab=expression(logUcd[GM]))
ord = order(vb_result$xgrid)
lines(vb_result$xgrid[ord],vb_result$post_curve[ord],lwd=2,col=2)
lines(vb_result$xgrid[ord],vb_result$post_upper[ord],lwd=2,col=3)
lines(vb_result$xgrid[ord],vb_result$post_lower[ord],lwd=2,col=3)
