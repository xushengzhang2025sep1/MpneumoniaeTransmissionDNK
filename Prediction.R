####################################### Date 26 December 2025 #####################################
# Base on the calibration model to predict the epidemic up to 2032
# we consider two potential behaviors in future from March 2025 on
# 1. completely return to normal as before COVID-19 pandemic
# 2. the behvaior of April 2024 to March 2025 will be repeated in future


## To predic the future epidemics, input your assumption about the change in births post March 2025
#  0   	births in future keep the same as that in 2024
# -1   	births in future decrease at annual rate of 7.95%
#  1    	births in future increase at annual rate of 7.60%
ID_births =0

## to draw epdiemics only (1 -- part) or epidmcsi plus relative transmisison (2 --- whole)
WoP =2


######## Estimating model parameters by calibrating to observations from Jan 2011 to just before NPIs (March 2020) ##########
tsiRparameter<-function(Data) {
 require(kernlab)

DenMark_preNPI = Data[(Data$time<2020.24),]  

## model calibration to the data prior to NPIs
parms <- estpars(
              DenMark_preNPI,
              xreg="cumcases",
              seasonality="standard",
              IP=3.0588,
              regtype="gaussian",
              sigmamax=3,
              family="gaussian",	
              link ="identity",
              userYhat= numeric(),
              alpha=NULL,
              sbar =NULL,
 	        printon= F);       

sim <- simulatetsir(
              DenMark_preNPI,
              IP=3.0588,   
              parms=parms,
              method = "deterministic",
              epidemics = "cont",  
              pred ="forward",
              threshold=1,  		
              inits.fit=FALSE,
              add.noise.sd =0,
              mul.noise.sd =0)

 return(list(1/sim$rho, parms$contact$beta, parms$alpha,sim$inits))
}
#####-----------------------------------------------------------------------------------------------------------#######


##--------------------------------------------------------------------------------------------------------##
##== Model predictions to asssess the effect of NPIs from its lifting to March 2025 and further to 2032 ==##
##--------------------------------------------------------------------------------------------------------##
NPI_Effect<-function(inputparameter) {    

 reduce_NPI = inputparameter[1:Npoint];					# reduction in contact rate due to NPIs among Npoint tri-weeks
 epsilon    = inputparameter[1+Npoint];					# increase in reporting rate from NPIs 
 KEY        = inputparameter[2+Npoint];					# choice for contact pattern post March 2025: 1 -- return to normal, 2 -- repeat last year(4/2024-3/2025) 

 times <- seq(1,55+45+8+1-86.940515-7, by = 1/ (52/IP))[1:TWeek]   	
 controlStart= (2020-2011)*(52/IP) + controlWeekStart
 controlEnd  =  controlStart       + controlWeekLength
 tEnd       = (2020-2011)*(52/IP) + (2025.176-2020)*17;            # Fixing the time point after which NPI effect disappears
 Dnpi=ceiling(controlWeekLength)

if(KEY==1) {                                                                                                       #1 resume the beta proir to NPIs (1)
  pred <- NPIeffectModel(times = times, births = births, beta = BETA, alpha = ALPHA,
                        S0 =floor(S_0), I0 = floor(I_0), nsim = 2, stochastic = F, 
                        controlStart = controlStart, controlEnd =controlEnd, betachangeNPI = 1-reduce_NPI, tEnd=tEnd)
} else if(KEY==2) {
  reduce_NPIbeyond =c(reduce_NPI,rep(reduce_NPI[(Npoint-16):Npoint],length(times)))[1:length(times)];    #2 repeat last year pattern 
  tEnd = (2020-2011)*(52/IP) + TWeek;   
  pred <- NPIeffectModel(times = times, births = births, beta = BETA, alpha = ALPHA,
                        S0 =floor(S_0), I0 = floor(I_0), nsim = 2, stochastic = F, 
                        controlStart = controlStart, controlEnd =controlEnd, betachangeNPI = 1-reduce_NPIbeyond, tEnd=tEnd)
}

 subsetI0 <- pred$I$mean
 subsetS  <- pred$S$mean

 mRrate =mean(ReportR)
 Rrate=c(rep(mRrate,4),ReportR,rep(epsilon*ReportR,12))[1:length(subsetI0)]/mRrate;    
 subsetI <-Rrate*subsetI0;						# reported cases  (detections)

 return(list(subsetS,subsetI))
}


################################# plot model result ##########################################
Drawepidemics<-function(ITW.quant,Ymax1,STW.quant,Ymax2,T3,Labelling) {
 tM=length(ITW.quant[3,])
 timeuse <- seq(2011,2100,1/(52/IP))[1:tM] 

par(mar=c(3,4,1,4))  			
Ymax=Ymax1; 
plot(timeuse,ITW.quant[2,]/1000,type="n",col="orange",lwd =1,ylim= c(0,Ymax/1000), xlab="",ylab="",xaxs="i",xlim=c(2011,2032))  
polygon(c(timeuse[ceiling((52/IP)*9+controlWeekStart)],timeuse[ceiling((52/IP)*9+controlWeekStart)],timeuse[ceiling((52/IP)*9 + controlWeekStart + controlWeekLength)],
          timeuse[ceiling((52/IP)*9 + controlWeekStart + controlWeekLength)]), c(-90000,Ymax/1000,Ymax/1000,-90000), border = NA,col="gray90")  
polygon(c(timeuse[1:(tM-T3)],rev(timeuse[1:(tM-T3)])),c(ITW.quant[1,1:(tM-T3)]/1000,rev(ITW.quant[3,1:(tM-T3)]/1000)),col=mycol3,border=NA)			

 lines(timeuse[1:(tM-T3)],ITW.quant[2,1:(tM-T3)]/1000,col="orange",lwd =1,xlab="",ylab=        ""); 

polygon(c(timeuse[(tM+1-T3):tM],rev(timeuse[(tM+1-T3):tM])),c(ITW.quant[1,(tM+1-T3):tM]/1000,rev(ITW.quant[3,(tM+1-T3):tM]/1000)),col=mycol4,border=NA)	
 lines(timeuse[(tM+1-T3):tM],ITW.quant[2,(tM+1-T3):tM]/1000,col="purple",lwd =1,xlab="",ylab=        ""); 

points(timeuse[1:242],dnkcases[1:242]/1000,cex=0.5,col="grey32",pch=18);  #cex=0.75,  242-- March 2025
points(timeuse[243:length(dnkcases)],dnkcases[243:length(dnkcases)]/1000,cex=0.5,col="grey32",pch=19);  #cex=0.75,
abline(v = 2025.176,lty=2,col="black",,lwd =1.5)

abline(v = seq(2011,2100,1),lty=3,col="gray")
par(new = TRUE)
plot(timeuse, STW.quant[2,]/1000, col="blue",type="l",xlab="",ylab="",axes =F,xlim=c(2011,2032),ylim=c(0,Ymax2/1000), lwd =2, lty =2)
polygon(c(timeuse,rev(timeuse)),c(STW.quant[1,]/1000,rev(STW.quant[3,]/1000)),col=mycol2,border=NA)			
lines(timeuse, STW.quant[2,]/1000, col="blue", lwd =1.1, lty =1)  #lty =2

axis(side = 4)
title(xlab="Year", line = 2)
title(ylab="Cases (n in 1000)", line = 2,col.lab="orange")
mtext(side = 4, line = 2, 'Susceptibles (n in 1000)',col="blue",cex = 0.8)
 text(2010,max(STW.quant[3,]/1000)*0.9,Labelling,pos = 4,cex = 0.8)   

}
##########------------------------------------------------------------######


################## draw the relative transmission rate due to NPIs and reporting rate ########################
DrawNPIeffect<-function(median.NPI,low.NPI,high.NPI,T3,Labelling) {    
 tM= length(median.NPI)
 timeuse <- seq(2011,2100,1/(52/IP))[1:tM] 

par(mar=c(3,4,1,4))  			
plot(timeuse,median.NPI,type="n",col="orange",lwd =1.2,ylim=c(0,2.0),xlab="",ylab="",xaxs="i",xlim=c(2011,2032))  
polygon(c(timeuse[ceiling((52/IP)*9+controlWeekStart)],timeuse[ceiling((52/IP)*9+controlWeekStart)],timeuse[ceiling((52/IP)*9 + controlWeekStart + controlWeekLength)],
          timeuse[ceiling((52/IP)*9 + controlWeekStart + controlWeekLength)]), c(-1,2,2,-1), border = NA,col="gray90")    
polygon(c(timeuse[1:(tM-T3)],rev(timeuse[1:(tM-T3)])),c(low.NPI[1:(tM-T3)],rev(high.NPI[1:(tM-T3)])),col=mycol3,border=NA)			
lines(timeuse[1:(tM-T3)], median.NPI[1:(tM-T3)], col="orange", lwd =1.2, lty =1)

polygon(c(timeuse[(tM+1-T3):tM],rev(timeuse[(tM+1-T3):tM])),c(low.NPI[(tM+1-T3):tM],rev(high.NPI[(tM+1-T3):tM])),col="grey",border=NA)	
lines(timeuse[(tM+1-T3):tM], median.NPI[(tM+1-T3):tM], col="dark grey", lwd =1.2, lty =1)   

abline(v = seq(2017,2100,1),lty=3,col="gray")
abline(v = 2025.176,lty=2,col="black",lwd =1.5)

par(new = TRUE)   # 2nd Y axis to draw the reporting rate
 mRrate =mean(ReportR);  mepsilon <- mean(epsilon);
  Rrate=100*c(rep(mRrate,4),ReportR,rep(mepsilon*ReportR,12))[1:TWeek]; 
plot(timeuse, Rrate, col="navy",type="l",xlab="",ylab="", axes =F,xlim=c(2011,2032),ylim=c(0.7*min(Rrate),max(Rrate)), lwd =0.8, lty =1);  
abline(h = 100*mRrate,lty=1,col="navy")

axis(side = 4)
title(xlab="Year", line = 2)
title(ylab="relative transmission", line = 2,col.lab="orange") #"grey32")
mtext(side = 4, line = 2, 'Reporting rate (%)',col="navy",cex = 0.8);  
 text(2010,max(Rrate)*0.9,Labelling,pos = 4,cex = 0.8)  
}


###-------------------------------------------------------------------------------------------###
############################################# main ##############################################
##--------------------------------------------------------------------------------------------###
 #setwd("C:/Users/xu-sheng.zhang/Documents/PCR_OCT_2025/manuscript/PNASLetter/PNAS_code/data_code_repository") 
 IP =3.058824;   							# infectious period of MP=52/17=3.058824 weeks
 TWeek=375;                          			# 22*17+1=375  tri-weeks from 2011 to 2032 with 51 weeks per year approximatedly
 controlWeekStart =14/IP;      			      # from April 2020
 controlWeekLength =97/IP;  					# up to 31/1/2022;  

#### load the estimates of model parameters  ####
 load("NPI_effect_Denmark.RData");   #obtained via MP_NPIs_MCMCsampling.R

 enpi    = Sampleenpi		# estimate of NPI-effect on Npoint=83 -- time point of triweeks from the start of NPI to March 0f 2025 (2025.176)
 epsilon = Sampleepsilon		# increase in the testing effort ( reporting rate of cases) from NPI lifitings
 Npoint = dim(enpi)[1];		#total number of points at which reduction due to NPI is estimated 
 NSample= dim(enpi)[2];   	#total number of samples


# load the packages for tsiR analysis
library(tsiR)
library(kernlab)
library(ggplot2)
library(fields)
library("pals")

 require("fields")
 source('NPIeffectModel.R', encoding = 'UTF-8');  

#####-----------------------------------------------------------------------------####
#### Transparent colors    ## Mark Gardener 2015  ## www.dataanalytics.org.uk-----####
t_col <- function(color, percent = 50, name = NULL) { # color = color name;  percent = % transparency;  name = an optional name for the color
## Get RGB values for named color
rgb.val <- col2rgb(color)

## Make new color using input color as base and alpha set by transparency
t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
             max = 255,
             alpha = (100 - percent) * 255 / 100,
             names = name)

## Save the color
invisible(t.col)
}
## END
 mycol2 <- t_col("blue",  perc = 80, name = "lt.blue")     
 mycol3 <- t_col("orange",perc = 60, name = "lt.orange")    
 mycol4 <- t_col("purple",perc = 70, name = "lt.purple")
##---------------------------------------------------------------------------------##


### Denmark national survey data from January 2011 to December 2025: the bith, population and observed data ###
 DenMark_tWK<- read.csv(file=paste("DenmarkMPdetectionsJan2011_Dec2025.csv",sep=""),head=TRUE)  
 pop 	  <- DenMark_tWK$pop[1];				    						  	# total population Data2011
 L_data = dim(DenMark_tWK)[1]
 births0 <- c(DenMark_tWK$births,rep(DenMark_tWK$births[(L_data-2):L_data],100))[1:TWeek];	# total tri-weekly births  births in future are assumed to be the same as 2024-2025

 aP = 0.076                      # largest increase rate per year
 aN = 0.0795                     # largest decrease rate per year during the period from 2017 to 2024

c_factorP=c(rep(1,17*15),rep((1+aP),17),rep((1+aP)^2,17),rep((1+aP)^3,17),rep((1+aP)^4,17),rep((1+aP)^5,17),rep((1+aP)^6,17),rep((1+aP)^7,18))
c_factorN=c(rep(1,17*15),rep((1-aN),17),rep((1-aN)^2,17),rep((1-aN)^3,17),rep((1-aN)^4,17),rep((1-aN)^5,17),rep((1-aN)^6,17),rep((1-aN)^7,18))

if(ID_births==0) {  births =births0*1; }   			# births in future keep the same as that in 2024
if(ID_births==-1) { births =births0*c_factorN;  }      	# births in future decreases at annual rate of 7.95%
if(ID_births==1) {  births =births0*c_factorP;  }      	# births in future increases at annual rate of 7.60%

Parameters=tsiRparameter(DenMark_tWK)		#use tsiR to estimate model parameters

 ReportR  = Parameters[[1]];				# reporting rate
 BETA     = Parameters[[2]];				# 17 triweek transmission rate
 ALPHA    = Parameters[[3]];        		# contact index alpha
 S_0     = Parameters[[4]][1];      		# initial number of the susceptible 
 I_0     = Parameters[[4]][2];    			# initial number of infections

# load denmark national surveillance data
dnk <- DenMark_tWK
dnk <- dnk[order(dnk$year, dnk$week),]

dfcorrect0 <- data.frame(week = c(rep(1+3*0:16, 14),1+3*0:16), year = c(rep(seq(2011,2024,1), each = ceiling(52/IP)),rep(2025,4+13)))        
dfcorrect <- dfcorrect0[1:dim(dfcorrect0)[1],];										#start from Jan 2011
dffin <- merge(dfcorrect,dnk, by=c("week","year"), all.x=T)
dffin <- dffin[order(dffin$year, dffin$week),]
dnkts <- dffin$cases;  		#dffin$percent_specimen_positive

 Npoint= ceiling((2025.176-2020)*17-as.numeric(14)/IP); 	# time point of triweeks from the start of NPI to March 0f 2025 (2025.176)
 NPITW =floor(as.numeric(97)/IP);                   		# duration of NPIs in triweek

################################### Variable for storing the epidemic curves #######################################
 ITWSample <- matrix(0,nrow=TWeek,ncol=NSample);     		# mean nos of infections of TWeek triweekly  
 STWSample <- matrix(0,nrow=TWeek,ncol=NSample);     		# mean nos of the susceptible of TWeek triweekly  

 ENPI1.quant <-apply(enpi[,1:NSample],1,quantile,probs=c(0.025,0.5,0.975));  			#the 95%CI & median for NPI effect on MP 

   t1 = (6+3)*17+ceiling(controlWeekStart);  			# from Jan 2011 to the start of NPIs
   t2=Npoint                                    		# NPI period
   t3=max(0,TWeek-t1-t2)                                	# NPI lifted

for(Key in 1:2) {  
  for(ii in 1:NSample) {		
    NewP=c(enpi[,ii],epsilon[ii],Key)
    ReEmergenceSeries<-  NPI_Effect(NewP); 			# tsiR model to generate the time series of re-emerged MP
    ITWSample[,ii] <- ReEmergenceSeries[[2]]; 			# the time series of reported cases by model		 
    STWSample[,ii] <- ReEmergenceSeries[[1]]; 			# the time series of the susceptible by model 	  	
  }  ##NSample

   Itweek<- ReEmergenceSeries[[2]]; 				# the time series of reported cases by model
   dnkcases <- (dnkts/(mean(dnkts[1:220], na.rm=T)))*mean(Itweek[1:floor(3*52/IP)])  # scale positives to match simulated incidence  1:220 from 2011 to March 2020

######----------------------------------------------------------------------------------#######

  if(Key==1) {
    ITW1.quant <-apply(ITWSample[,1:NSample],1,quantile,probs=c(0.025,0.5,0.975));  	# the 95%CI & median of infections per tri-week
    STW1.quant <-apply(STWSample[,1:NSample],1,quantile,probs=c(0.025,0.5,0.975));  	# the 95%CI & median of the susceptible 
### NPI-effect  --- return to normal (1) from April 2025

    median.NPI1=c(rep(1,t1),1-ENPI1.quant[2,],rep(1,t3))
    low.NPI1   =c(rep(1,t1),1-ENPI1.quant[1,],rep(1,t3))
    high.NPI1  =c(rep(1,t1),1-ENPI1.quant[3,],rep(1,t3))
  }
  if(Key==2) {
    ITW2.quant <-apply(ITWSample[,1:NSample],1,quantile,probs=c(0.025,0.5,0.975));  	# the 95%CI & median of infections per tri-week
    STW2.quant <-apply(STWSample[,1:NSample],1,quantile,probs=c(0.025,0.5,0.975));  	# the 95%CI & median of the susceptible 
### NPI-effect  --- repeat what were during 4/2024-3/2025
    median.NPI2=c(rep(1,t1),1-ENPI1.quant[2,],rep(1-ENPI1.quant[2,(Npoint-16):Npoint],t3)[1:t3])
    low.NPI2   =c(rep(1,t1),1-ENPI1.quant[1,],rep(1-ENPI1.quant[1,(Npoint-16):Npoint],t3)[1:t3])
    high.NPI2  =c(rep(1,t1),1-ENPI1.quant[3,],rep(1-ENPI1.quant[3,(Npoint-16):Npoint],t3)[1:t3])
   }
}

# Align to the same size for the two situation for the sake of visual comparison
YM1= max(ITW1.quant[2,],ITW2.quant[2,],max(dnkcases,na.rm=T));   	
YM2= max(STW1.quant[3,],STW2.quant[3,]);					

if(WoP==1) { 
  if(ID_births==0) {  pdf(paste0("Fig1_Prediction2032_birthsSame.pdf"),width=6,height=4);  }   	 # births in future keep the same as that in 2024
  if(ID_births==-1) { pdf(paste0("Fig1_Prediction2032_birthsDecrease.pdf"),width=6,height=4);  }   # births in future decreases at annual rate of 7.95%
  if(ID_births==1) {  pdf(paste0("Fig1_Prediction2032_birthsIncrease.pdf"),width=6,height=4);  }   # births in future increases at annual rate of 7.60%

  par(mfrow=c(2,1))
  Drawepidemics(ITW1.quant,YM1,STW1.quant,YM2,t3,"contact resumes that before NPIs");       		
  Drawepidemics(ITW2.quant,YM1,STW2.quant,YM2,t3,"contact repeats last year")  

  dev.off()
}

if(WoP==2) { 
  if(ID_births==0) {  pdf(paste0("Fig1_Prediction2032_birthsSameW.pdf"),width=6,height=4);  }   	  # births in future keep the same as that in 2024
  if(ID_births==-1) { pdf(paste0("Fig1_Prediction2032_birthsDecreaseW.pdf"),width=6,height=4);  }   # births in future decreases at annual rate of 7.95%
  if(ID_births==1) {  pdf(paste0("Fig1_Prediction2032_birthsIncreaseW.pdf"),width=6,height=4);  }   # births in future increases at annual rate of 7.60%
 
par(mfrow=c(2,2))
  Drawepidemics(ITW1.quant,YM1,STW1.quant,YM2,t3,"a)");     Drawepidemics(ITW2.quant,YM1,STW2.quant,YM2,t3,"c)")  
  DrawNPIeffect(median.NPI1,low.NPI1,high.NPI1,t3,"b)"); 	DrawNPIeffect(median.NPI2,low.NPI2,high.NPI2,t3,"d)"); 

dev.off()
}
