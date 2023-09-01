#Input Data Excel
library(readxl)
idx30 <- read_excel("Febrina/Skripsi/skripsi.xlsx", sheet = "Data2021-2022")
SBI<- read_excel("Febrina/Skripsi/skripsi.xlsx", sheet = "SBI2021", range = "C1:C25")

#Return Saham
## Return ADRO
adro<-idx30$ADRO
n<-length(adro)
ADRO<-matrix(adro,nrow=n,ncol=1)
rADRO<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rADRO[i]<-log(ADRO[i+1,]/ADRO[i,])}
rADRO 
## Return ANTM
antm<-idx30$ANTM
n<-length(antm)
ANTM<-matrix(antm,nrow=n,ncol=1)
rANTM<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rANTM[i]<-log(ANTM[i+1,]/ANTM[i,])}
rANTM 
## Return ASII
asii<-idx30$ASII
n<-length(asii)
ASII<-matrix(asii,nrow=n,ncol=1)
rASII<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rASII[i]<-log(ASII[i+1,]/ASII[i,])}
rASII
## Return BBCA
bbca<-idx30$BBCA
n<-length(bbca)
BBCA<-matrix(bbca,nrow=n,ncol=1)
rBBCA<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rBBCA[i]<-log(BBCA[i+1,]/BBCA[i,])}
rBBCA
## Return BBNI
bbni<-idx30$BBNI
n<-length(bbni)
BBNI<-matrix(bbni,nrow=n,ncol=1)
rBBNI<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rBBNI[i]<-log(BBNI[i+1,]/BBNI[i,])}
rBBNI
## Return BBRI
bbri<-idx30$BBRI
n<-length(bbri)
BBRI<-matrix(bbri,nrow=n,ncol=1)
rBBRI<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rBBRI[i]<-log(BBRI[i+1,]/BBRI[i,])}
rBBRI
## Return BMRI
bmri<-idx30$BMRI
n<-length(bmri)
BMRI<-matrix(bmri,nrow=n,ncol=1)
rBMRI<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rBMRI[i]<-log(BMRI[i+1,]/BMRI[i,])}
rBMRI
## Return CPIN
cpin<-idx30$CPIN
n<-length(cpin)
CPIN<-matrix(cpin,nrow=n,ncol=1)
rCPIN<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rCPIN[i]<-log(CPIN[i+1,]/CPIN[i,])}
rCPIN
## Return ICBP
icbp<-idx30$ICBP
n<-length(icbp)
ICBP<-matrix(icbp,nrow=n,ncol=1)
rICBP<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rICBP[i]<-log(ICBP[i+1,]/ICBP[i,])}
rICBP
## Return INDF
indf<-idx30$INDF
n<-length(indf)
INDF<-matrix(indf,nrow=n,ncol=1)
rINDF<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rINDF[i]<-log(INDF[i+1,]/INDF[i,])}
rINDF
## Return INKP
inkp<-idx30$INKP
n<-length(inkp)
INKP<-matrix(inkp,nrow=n,ncol=1)
rINKP<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rINKP[i]<-log(INKP[i+1,]/INKP[i,])}
rINKP
## Return KLBF
klbf<-idx30$KLBF
n<-length(klbf)
KLBF<-matrix(klbf,nrow=n,ncol=1)
rKLBF<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rKLBF[i]<-log(KLBF[i+1,]/KLBF[i,])}
rKLBF
## Return MDKA
mdka<-idx30$MDKA
n<-length(mdka)
MDKA<-matrix(mdka,nrow=n,ncol=1)
rMDKA<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rMDKA[i]<-log(MDKA[i+1,]/MDKA[i,])}
rMDKA
## Return PGAS
pgas<-idx30$PGAS
n<-length(pgas)
PGAS<-matrix(pgas,nrow=n,ncol=1)
rPGAS<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rPGAS[i]<-log(PGAS[i+1,]/PGAS[i,])}
rPGAS
## Return PTBA
ptba<-idx30$PTBA
n<-length(ptba)
PTBA<-matrix(ptba,nrow=n,ncol=1)
rPTBA<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rPTBA[i]<-log(PTBA[i+1,]/PTBA[i,])}
rPTBA
## Return SMGR
smgr<-idx30$SMGR
n<-length(smgr)
SMGR<-matrix(smgr,nrow=n,ncol=1)
rSMGR<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rSMGR[i]<-log(SMGR[i+1,]/SMGR[i,])}
rSMGR
## Return TBIG
tbig<-idx30$TBIG
n<-length(tbig)
TBIG<-matrix(tbig,nrow=n,ncol=1)
rTBIG<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rTBIG[i]<-log(TBIG[i+1,]/TBIG[i,])}
rTBIG
## Return TLKM
tlkm<-idx30$TLKM
n<-length(tlkm)
TLKM<-matrix(tlkm,nrow=n,ncol=1)
rTLKM<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rTLKM[i]<-log(TLKM[i+1,]/TLKM[i,])}
rTLKM
## Return TOWR
towr<-idx30$TOWR
n<-length(towr)
TOWR<-matrix(towr,nrow=n,ncol=1)
rTOWR<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rTOWR[i]<-log(TOWR[i+1,]/TOWR[i,])}
rTOWR
## Return UNTR
untr<-idx30$UNTR
n<-length(untr)
UNTR<-matrix(untr,nrow=n,ncol=1)
rUNTR<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rUNTR[i]<-log(UNTR[i+1,]/UNTR[i,])}
rUNTR
## Return UNVR
unvr<-idx30$UNVR
n<-length(unvr)
UNVR<-matrix(unvr,nrow=n,ncol=1)
rUNVR<-matrix(nrow=n-1,ncol=1)
for (i in 1:n-1)
{rUNVR[i]<-log(UNVR[i+1,]/UNVR[i,])}
rUNVR

#Expected Return dan Standar deviasi Saham
mean(rADRO)
mean(rANTM)
mean(rASII)
mean(rBBCA)
mean(rBBNI)
mean(rBBRI)
mean(rBMRI)
mean(rCPIN)
mean(rICBP)
mean(rINDF)
mean(rINKP)
mean(rKLBF)
mean(rMDKA)
mean(rPGAS)
mean(rPTBA)
mean(rSMGR)
mean(rTBIG)
mean(rTLKM)
mean(rTOWR)
mean(rUNTR)
mean(rUNVR)

sd(rADRO)
sd(rANTM)
sd(rASII)
sd(rBBCA)
sd(rBBNI)
sd(rBBRI)
sd(rBMRI)
sd(rCPIN)
sd(rICBP)
sd(rINDF)
sd(rINKP)
sd(rKLBF)
sd(rMDKA)
sd(rPGAS)
sd(rPTBA)
sd(rSMGR)
sd(rTBIG)
sd(rTLKM)
sd(rTOWR)
sd(rUNTR)
sd(rUNVR)

#Normalitas Univariat
library(nortest)
ks.test(rADRO,"pnorm",mean(rADRO),sd(rADRO))
ks.test(rBBCA,"pnorm",mean(rBBCA),sd(rBBCA))
ks.test(rBBNI,"pnorm",mean(rBBNI),sd(rBBNI))
ks.test(rBBRI,"pnorm",mean(rBBRI),sd(rBBRI))
ks.test(rBMRI,"pnorm",mean(rBMRI),sd(rBMRI))
ks.test(rICBP,"pnorm",mean(rICBP),sd(rICBP))
ks.test(rKLBF,"pnorm",mean(rKLBF),sd(rKLBF))
ks.test(rMDKA,"pnorm",mean(rMDKA),sd(rMDKA))
ks.test(rPGAS,"pnorm",mean(rPGAS),sd(rPGAS))
ks.test(rPTBA,"pnorm",mean(rPTBA),sd(rPTBA))
ks.test(rTBIG,"pnorm",mean(rTBIG),sd(rTBIG))
ks.test(rTLKM,"pnorm",mean(rTLKM),sd(rTLKM))
ks.test(rTOWR,"pnorm",mean(rTOWR),sd(rTOWR))

#Korelasi antar return saham
Q=matrix(c(rADRO,rBBCA,rBBNI,rBBRI,rBMRI,rICBP,rKLBF,rMDKA,rPGAS,rPTBA,rTLKM,rTOWR),ncol=12)
cor(Q)

#Hitung Rf
Rf=(mean(SBI$SBI))/(52)
Rf

#Metode CCM
##Hitung ERS
#BMRI
ERS_rBMRI<-(mean(rBMRI)-Rf)/sd(rBMRI)
ERS_rBMRI
#ICBP
ERS_rICBP<-(mean(rICBP)-Rf)/sd(rICBP)
ERS_rICBP
#KLBF
ERS_rKLBF<-(mean(rKLBF)-Rf)/sd(rKLBF)
ERS_rKLBF
#MDKA
ERS_rMDKA<-(mean(rMDKA)-Rf)/sd(rMDKA)
ERS_rMDKA
#PGAS
ERS_rPGAS<-(mean(rPGAS)-Rf)/sd(rPGAS)
ERS_rPGAS
#TLKM
ERS_rTLKM<-(mean(rTLKM)-Rf)/sd(rTLKM)
ERS_rTLKM
#TOWR
ERS_rTOWR<-(mean(rTOWR)-Rf)/sd(rTOWR)
ERS_rTOWR

##Menghitung Korelasi rata-rata (p*)
P=matrix(c(rBMRI,rKLBF,rMDKA,rPGAS,rTLKM,rTOWR),ncol=6)
cor(P)
kor=sum(c(0.2871907,0.2500558,0.3566650,0.3181231,-0.1013199,0.1436972,0.1965193,0.1097404,-0.1116666,0.1542619,0.1039276,0.1052117,0.2250230,-0.1140887,0.1219574))
kor
r=kor/((6*(6-1)/2))
r

##Menghitung Ci
C_rBMRI<-(r/(1-r+(1*r)))*(ERS_rBMRI)
C_rBMRI
C_rKLBF<-(r/(1-r+(2*r)))*(ERS_rBMRI+ERS_rKLBF)
C_rKLBF
C_rMDKA<-(r/(1-r+(3*r)))*(ERS_rBMRI+ERS_rKLBF+ERS_rMDKA)
C_rMDKA
C_rTOWR<-(r/(1-r+(4*r)))*(ERS_rBMRI+ERS_rKLBF+ERS_rMDKA+ERS_rTOWR)
C_rTOWR
C_rTLKM<-(r/(1-r+(5*r)))*(ERS_rBMRI+ERS_rKLBF+ERS_rMDKA+ERS_rTOWR+ERS_rTLKM)
C_rTLKM
C_rPGAS<-(r/(1-r+(6*r)))*(ERS_rBMRI+ERS_rKLBF+ERS_rMDKA+ERS_rTOWR+ERS_rTLKM+ERS_rPGAS)
C_rPGAS

##Menentukan Cut-off Rate(C*)
Cstar<-max(C_rBMRI,C_rKLBF,C_rMDKA,C_rTOWR,C_rTLKM,C_rPGAS)
Cstar

##Menghitung nilai Z
print(Z_rBMRI<-(1/((1-r)*sd(rBMRI))*(ERS_rBMRI-Cstar)))
print(Z_rKLBF<-(1/((1-r)*sd(rKLBF))*(ERS_rKLBF-Cstar)))
print(Z_rMDKA<-(1/((1-r)*sd(rMDKA))*(ERS_rMDKA-Cstar)))
print(Z<-sum(Z_rBMRI,Z_rKLBF,Z_rMDKA))

##Menghitung Bobot (w)
print(w_rBMRI<-Z_rBMRI/Z)
print(w_rKLBF<-Z_rKLBF/Z)
print(w_rMDKA<-Z_rMDKA/Z)

##Uji Normalitas Multivariat
mnorm.test<-function(x) 
{ 
  rata2 <- apply(x, 2, mean) 
  mcov <- var(x) 
  ds<-sort(mahalanobis(x, center = rata2,cov = mcov)) 
  n<-length(ds) 
  p <- (1:n-0.5)/n 
  chi <- qchisq(p, df = ncol(x)) 
  win.graph() 
  plot(ds, chi, type = "p") 
  return(ks.test(ds,chi,df=ncol(x))) } 
A=matrix(c(rBMRI,rKLBF,rMDKA),ncol=3)
mnorm.test(A)

##Perhitungan Banyak Pengulangan
myu<-mean(c(sum(rBMRI),sum(rKLBF),sum(rMDKA)))
myu
strBMRI=(sum(rBMRI)-myu)^2
strKLBF=(sum(rKLBF)-myu)^2
strMDKA=(sum(rMDKA)-myu)^2
std<-sqrt((sum(strBMRI,strKLBF,strMDKA))/3) #standardeviasi
std 
E=myu/(1/0.05) #Nilai Kesalahan
E
K=round(((3*std)/E)^2) # Banyak Pengulangan
K

#Perhitungan ES dengan Simulasi Monte Carlo 1 kali
library(MASS)
A=cbind(rBMRI,rKLBF,rMDKA)
A
w=matrix(c(w_rBMRI,w_rKLBF,w_rMDKA))
w
n=104
V0=1
alpha=0.05
t=1
MeanR=colMeans(A)
SigmaR=var(A)
return.ccm=mvrnorm(n,MeanR,SigmaR)
Rp=return.ccm%*%w
MeanRp=mean(Rp)
MeanRp
sdRp=sd(Rp)
sdRp
Rstar=quantile(Rp,1-alpha)
Rstar
VaR=V0*Rstar*sqrt(t)
VaR
ES=V0*(MeanRp+(sdRp*dnorm(qnorm(VaR))/alpha))*sqrt(t)
ES

#Perhitungan ES dengan Simulasi Monte Carlo
A=cbind(rBMRI,rKLBF,rMDKA)
A
w=matrix(c(w_rBMRI,w_rKLBF,w_rMDKA))
w
n=104
V0=1
alpha=0.05
t=1
VaR=matrix(nrow=K,ncol=1)
ES=matrix(nrow=K,ncol=1)
for(i in 1:K)
{
  MeanR=colMeans(A)
  SigmaR=var(A)
  return.ccm=mvrnorm(n,MeanR,SigmaR)
  Rp=return.ccm%*%w
  MeanRp=mean(Rp)
  sdRp=sd(Rp)
  Rstar=quantile(Rp,1-alpha)
  VaR[i]=V0*Rstar*sqrt(t)
  ES[i]=V0*(MeanRp+(sdRp*dnorm(qnorm(VaR[i]))/alpha))*sqrt(t)
}
MeanRp
MeanVaR=print(mean(VaR))
MeanES=print(mean(ES))


