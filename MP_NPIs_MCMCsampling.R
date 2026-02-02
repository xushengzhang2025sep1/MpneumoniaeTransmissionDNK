
#### Estimating model parameters by calibrating to observations from January 2011 to just before NPIs (March 2020) #####
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


###################################################################################################################
##== Model predictions to asssess the effect of NPIs from the liting of NPIs to March 2025 and further to 2032 ==##
###################################################################################################################
NPI_Effect<-function(inputparameter) {    

 reduce_NPI = inputparameter[1:Npoint];					# reduction in contact due to NPIs control among Npoint tri-weeks
 epsilon    = inputparameter[1+Npoint];					# increase in reporting rate from NPIs 

 times <- seq(1,55+45+8+1-86.940515-7, by = 1/ (52/IP))[1:TWeek]   
 controlStart= (2020-2011)*(52/IP) + controlWeekStart
 controlEnd  =  controlStart       + controlWeekLength
 tEnd       = (2020-2011)*(52/IP) +  (2025.176-2020)*17;           # Fixing the time point after which NPI effect is assumed to disappear

pred <- NPIeffectModel(times = times, births = births, beta = BETA, alpha = ALPHA,
                        S0 =floor(S_0), I0 = floor(I_0), nsim = 2, stochastic = F, 
                        controlStart = controlStart, controlEnd =controlEnd, betachangeNPI = 1-reduce_NPI, tEnd=tEnd)
 subsetI0<- pred$I$mean
 subsetS <- pred$S$mean

 mRrate =mean(ReportR)
 Rrate=c(rep(mRrate,4),ReportR, rep(epsilon*ReportR,12))[1:length(subsetI0)]/mRrate;   
 subsetI <-Rrate*subsetI0;						# cases reported in accordance with the reporting pattern over the 2011-2020 period

 return(list(subsetS,subsetI))
}
###########################################################################################################


###########################################################################################################
################## calculate negative log NB likelihood for outbreak data:Mcases ##########################
LikelihoodNB<-function(Mcases,Ocases,eta) {  	# eta dispersions of NB distribution 
  aa1<-log(eta);      			ab1 <- log(1.-1./eta);	
  rr1<- Mcases[Tc:242]/(eta-1.0);	

  LLike<- -(LogMcasesplus1) 
  LLike<- LLike +ab1*sum(Ocases[Tc:242])+ sum(lgamma(Ocases[Tc:242]+rr1)-lgamma(rr1)-rr1*aa1)  

	return (-LLike);								
}


############# Sample new parameters with Normal walk around the old parameters for HM algoritham of MCMC ################
SampleParameter<-function(OldPara,stDev,Ich) { 	# Npoint enpi (NPI effect),epsilon, and eta (dispension)

 Para<-rep(-1,ndim);
 for(i in 1:ndim) {
  if(i==Ich) {  Para[i] = rnorm(1)*stDev[Ich] + OldPara[Ich]; 
  } else { Para[i] =OldPara[i]; }
 }

return(Para);
}


############## Sample new parameters from priori distribution for Bayesian framework ###################
SampleParameterSeeds<-function() {
 NewPara<-rep(-1,ndim);
 for(i in 1:ndim) NewPara[i] <- runif(1,min=PI_shape[i],max=PI_rate[i]);  
 return (NewPara);
}


Draw_Sample_Density<-function(Name_Parameter,MCMC_Sample,Prior_L_U_MLE) {
##----------------------- MCMC samples ---------------------------##
  Qant_MCMC<-quantile(MCMC_Sample,probs=c(0.025,0.5,0.975));
  plot(1:length(MCMC_Sample),MCMC_Sample,col="pink",xlab="#Sample",ylab=Name_Parameter,lty=5,type="b",cex=.4,pch=20,
  main=paste(round(Qant_MCMC[2],3)," [",round(Qant_MCMC[1],3),",",round(Qant_MCMC[3],3),"]",sep="")); 

## ------------ hist of posterior distribution --------------##
 hist(MCMC_Sample, xlab=Name_Parameter,ylab="Number of simulations",main=paste("MLE:",round(Prior_L_U_MLE[3],3)))
 abline(v=Prior_L_U_MLE[1],lty=2,col="red")
 abline(v=Prior_L_U_MLE[2],lty=2,col="red"); abline(v=Prior_L_U_MLE[3],lty=2,col="green")
}


##---------------------------------------------------------------------------------------------##
################################################## main #########################################
##---------------------------------------------------------------------------------------------##
 #setwd("C:/Users/xu-sheng.zhang/Documents/PCR_OCT_2025/manuscript/PNASLetter/PNAS_code/data_code_repository") 
 today=Sys.Date() 						# because reports are done one day after later
 IP =3.058824;   							# infectious period of MP=52/17=3.058824 weeks
 TWeek=242;                                           # (2025-2011)*17+4   for model fit up to March 2025
 controlWeekStart =14/IP; 	   				# from April 2020
 controlWeekLength =97/IP;  					# up to 31/1/2022;  

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
t_col <- function(color, percent = 50, name = NULL) { # color = color name; percent = % transparency;  name = an optional name for the color
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
 mycol2 <- t_col("green", perc = 90, name = "lt.green")
 mycol3 <- t_col("Magenta",   perc = 90, name = "lt.Magenta")
##---------------------------------------------------------------------------------##


## Denmark national survey data from 2011 to March 2025: 
 Data_tWK<- read.csv(file=paste("DenmarkMPdetectionsJan2011_Dec2025.csv"),head=TRUE) 
 DenMark_tWK  = Data_tWK[(Data_tWK$time<2025.235),]  

 pop 	  <- DenMark_tWK$pop[1];				    						  	# total population Data2011
 L_data = dim(DenMark_tWK)[1]
 births <- c(DenMark_tWK$births,rep(DenMark_tWK$births[(L_data-2):L_data],100))[1:TWeek];	 	# total tri-weekly births  births in future are assumed to be the same as in 2025

 Tc = 0+(6+3)*17+floor(controlWeekStart);  		# from Jan 2011 to the start of NPIs    	# point to start to calibrate with data

Parameters=tsiRparameter(DenMark_tWK)		#use tsiR to eatimate model parameters

 ReportR  = Parameters[[1]];				# reporting rate
 BETA     = Parameters[[2]];				# 17 triweek transmissiob rate
 ALPHA    = Parameters[[3]];  		      # contact index alpha
 S_0     = Parameters[[4]][1];      		# initial number of the susceptible 
 I_0     = max(10,Parameters[[4]][2]);    	# initial number of infections

# load denmark national surveillance data
dnk <- DenMark_tWK
dnk <- dnk[order(dnk$year, dnk$week),]

dfcorrect0 <- data.frame(week = c(rep(1+3*0:16, 14),1+3*0:3), year = c(rep(seq(2011,2024,1), each = ceiling(52/IP)),rep(2025,4)))        
dfcorrect <- dfcorrect0[1:dim(dfcorrect0)[1],];										#start from Jan 2011
dffin <- merge(dfcorrect,dnk, by=c("week","year"), all.x=T)
dffin <- dffin[order(dffin$year, dffin$week),]
dnkts <- dffin$cases;  		


####################### setup for MCMC ##############################
#_____________________________adaptive parameter for MCMC sampling___________________________________#
 mixingthresholdUB<-0.40; 					# The upper boundary  for the accepted ratio
 mixingthresholdLB<-0.12; 					# The lower boundary  for the accepted ratio
 Increasingfactor <-1.2;  					# increase the stDev if accepted ratio is lower than thresholdUB
 Decreasingfactor <-0.8;  					# decrease the stDev if accepted ratio is higher than thresholdLB

#____________________________parameters controling the running of the MCMC___________________________#
 IOP  <-  200;       					# steps for output  "thinning of sampling"
 Iadp <-  150;   						# steps for adaptive change of stdeviation of ndim parameters
 Loop <-  2000; 	      	      		# size of samples
 Lburnin<-1000000;  						# burn-in steps

 Npoint= ceiling((2025.176-2020)*17-as.numeric(14)/IP); #time point of triweeks from the start of NPI to March 0f 2025 (2025.176)
 NPITW =floor(as.numeric(99+1-3)/IP);                 #duration of NPIs in triweek


 ndim<-Npoint+1+1;   	 					#Npoint enpi (NPI effect),epsilon, and eta (dispension of NB likelihood function)
##_______________________hyperpareameter for prior distribution of model parameters____________________
#            enpi              and             epsilon eta--- uniformly distributed [PI_shape,PI_rate]
 PI_shape<-c(rep(0.25,NPITW),rep(0.1,Npoint-NPITW),2.5,1.5);       
 PI_rate <-c(rep(0.75,Npoint),                     4.5,100000);   	

 para_a  <-c(rep(0.00,NPITW),rep(-0.5,Npoint-NPITW),1.0, 1.01);   		#lower boundary 
 para_b  <-c(rep(1.0,Npoint),                      10., 1000000);   	#upper boundary  
 stDev   <-c(rep(0.3,Npoint),                      2.5, 50000); 		#SD for sampling new parameters in the MCMC process (HM algoritham) 

################################### Variables for storing the epidemic curves #######################################
ITWSample    <- matrix(0,nrow=TWeek,ncol=Loop);     	# nos of infections of TWeek triweekly  
STWSample    <- matrix(0,nrow=TWeek,ncol=Loop);     	# nos of the susceptible of TWeek triweekly  

###---  ---  	---   	&&&   	****** 	&&&  		--- 	--- 	----    ##
 Itweek   <- rep(0,TWeek);     	 	  	 	# tri-weekly number of reported cases predicted by model
 Stweek   <- rep(0,TWeek);     	 	   		# tri-weekly number of susceptible individuals
    
logAlph<-0;
OldP <-rep(0,ndim);					NewP <-rep(0,ndim);       			# sample of parameters 
OldQ <-rep(0,ndim);					NewQ <-rep(0,ndim);       			# the inverse of the proposal function
PriorOld<-rep(0,ndim);					PriorNew<-rep(0,ndim);    			# prior density 

Sampleenpi <-matrix(0,nrow=Npoint,ncol=Loop);     	Sampleepsilon<-rep(0,ncol=Loop);	
Sampleeta  <-rep(0,ncol=Loop);				SampleLL     <-rep(0,ncol=Loop);

 OldP<-SampleParameterSeeds();									
 ReEmergenceSeries<- NPI_Effect(OldP[1:(Npoint+1)]); 							# tsiR model to generate the time series of outbreak due to NPIs
 Itweek[1:TWeek]  <- ReEmergenceSeries[[2]]; 								# the time series of reported cases from model

 dnkcases <- (dnkts/(mean(dnkts[1:220], na.rm=T)))*mean(Itweek[1:floor(3*52/IP)]) 		# scale % pos to match simulated incidence  1:220 from Jan 2011 to Mar 2020
 LogMcasesplus1  <-sum(lgamma(dnkcases[Tc:242]+1)); 							# for Negative binomial likelihood function

 prevLogLike<- LikelihoodNB(Itweek,dnkcases,OldP[ndim]) 						#calculate neg log Negative Binomial likelihood of tweek with outbreak data  

for(i in 1:length(stDev)) {			   								# for the use in calculation of acceptance rate
  OldQ[i] = pnorm((para_b[i]-OldP[i])/stDev[i])-pnorm((para_a[i]-OldP[i])/stDev[i]);  	# the inverse of the proposal function (normal walk)
  PriorOld[i]<- dunif(OldP[i],para_a[i],para_b[i]); 							# print(paste("PriorOld[",i,"]=",PriorOld[i]));
}

##  ----------------------------------- Markov Chains -------------------------------------------
print("enpi1[i], i=1...Npoint   epsilon  eta  LL");
SI <-0;                                   					#no of accepted samples
j  <-0;    											#total accepted steps
while(SI<Loop) {
 numberAccepted<-rep(0,ndim);							    	#empty the accepted number box to monitor the accepted ratio
for(jj in 1:Iadp) {    			  	    	     	    			#adaptive
   if(SI==Loop) break;
   for(ii in 1:IOP) {						    			#output "thinning"

      for(i1 in 1:ndim) {
         	NewP<-SampleParameter(OldP,stDev,i1);   	    			#generate new proposal NewP from the current point OldP
  if(NewP[i1]>=para_a[i1]&&NewP[i1]<=para_b[i1]) {
       PriorNew[i1]<- dunif(NewP[i1],para_a[i1],para_b[i1]);

       ReEmergenceSeries<-  NPI_Effect(NewP[1:(Npoint+1)]); 		# tsiR model to generate the time series of MP cases
       Stweek[1:TWeek]<- ReEmergenceSeries[[1]]; 				# the time series of the susceptible
       Itweek[1:TWeek]<- ReEmergenceSeries[[2]]; 				# the time series of reported cases 
       dnkcases <- (dnkts/(mean(dnkts[1:220], na.rm=T)))*mean(Itweek[1:floor(3*52/IP)]) 		# scale positives to match simulated incidence  1:220 from 2011 to March 2020

       postLogLike<-LikelihoodNB(Itweek,dnkcases,NewP[ndim]) 			
	 NewQ[i1] = pnorm((para_b[i1]-NewP[i1])/stDev[i1])-pnorm((para_a[i1]-NewP[i1])/stDev[i1]);      
           
		logAlph = -(postLogLike - prevLogLike);
              if(runif(1)<=(PriorNew[i1]/PriorOld[i1])*exp(logAlph)*(OldQ[i1]/NewQ[i1])) { 
                  prevLogLike = postLogLike;
			Itweek_select<-Itweek;		Stweek_select<-Stweek;			
                  OldP[i1]    =NewP[i1];        OldQ[i1]    =NewQ[i1];
			PriorOld[i1]=PriorNew[i1];
                  numberAccepted[i1]<- numberAccepted[i1]+1;
			j<-j+1; 		
              }
  }  ### if(Newp[i1]....)
       }   ## for(i1 in 1:ndim) 
   }       #for(ii in 1:IOP) 

  if(j>Lburnin) {    
    SI <- SI+1;
    ITWSample[,SI] <-Itweek_select[1:TWeek];		 STWSample[,SI] <-Stweek_select[1:TWeek]; 	  	

    Sampleenpi[,SI]<- round(OldP[1:Npoint],6);	 	Sampleepsilon[SI]<- round(OldP[ndim-1],4); 		
    Sampleeta[SI]  <- round(OldP[ndim],4); 		SampleLL[SI]     <- round(prevLogLike,3);

   print(paste(SI,Sampleenpi[1,SI],Sampleenpi[10,SI],Sampleenpi[Npoint,SI],Sampleeta[SI],Sampleepsilon[SI],SampleLL[SI],sep=" "));
  }   ### j

} ## Iadp  

  for(i0 in 1:ndim) {   		## control the walk steps by judging the step sizes
     	ratio = numberAccepted[i0]/(IOP*Iadp);   						  
	if (ratio>mixingthresholdUB) { stDev[i0] = stDev[i0]*Increasingfactor; }
     	else if (ratio<mixingthresholdLB) { stDev[i0]=stDev[i0]*Decreasingfactor; }
  }

}  ##Loop
########################################################################################################


ITW.quant <-apply(ITWSample[,1:SI],1,quantile,probs=c(0.025,0.5,0.975));  	# the 95%CI & median for infections per tri-week
STW.quant <-apply(STWSample[,1:SI],1,quantile,probs=c(0.025,0.5,0.975));  	# the 95%CI & median for the susceptible 
ENPI.quant <-apply(Sampleenpi[,1:SI],1,quantile,probs=c(0.025,0.5,0.975));  	# the 95%CI & median for NPI effect on MP transmission

####################    save the MCMC samples of model parameters    ########################
save(Sampleenpi,Sampleepsilon,Sampleeta,SampleLL, file = paste0("NPI_effect_Denmark.RData"))


################################# plot model result ##########################################
 timeuse <- seq(2011,2100,1/(52/IP))[1:length(ITW.quant[3,])] 

pdf(paste0("Figure_Modelfit.pdf"),width=6,height=4)
par(mar=c(3,4,1,4))  			
Ymax=max(ITW.quant[3,],max(dnkcases,na.rm=T))
plot(timeuse,ITW.quant[2,],type="n",col="#F64740",lwd =1,ylim= c(0,Ymax), xlab="",ylab="",xaxs="i",bty = "n",xlim=c(2011,2032))   

polygon(c(timeuse[ceiling((52/IP)*9+controlWeekStart)],timeuse[ceiling((52/IP)*9+controlWeekStart)],timeuse[ceiling((52/IP)*9 + controlWeekStart + controlWeekLength)],
          timeuse[ceiling((52/IP)*9 + controlWeekStart + controlWeekLength)]), c(-90000,Ymax,Ymax,-90000), border = NA,col="#D1D3DF")
polygon(c(timeuse,rev(timeuse)),c(ITW.quant[1,],rev(ITW.quant[3,])),col=mycol3,border=NA)			

 lines(timeuse,ITW.quant[2,],col="#F64740",lwd =1,xlab="",ylab=        ""); 

points(timeuse[1:length(dnkcases)],dnkcases,cex=0.75,col="grey32",pch=18)
abline(v = seq(2011,2100,1),lty=3,col="gray")
par(new = TRUE)
plot(timeuse, STW.quant[2,], col="green",type="l",xlab="",ylab="",axes =F,xlim=c(2011,2032),ylim=c(0.7*min(STW.quant[1,]),max(STW.quant[3,])), lwd =2, lty =2)
polygon(c(timeuse,rev(timeuse)),c(STW.quant[1,],rev(STW.quant[3,])),col=mycol2,border=NA)			
lines(timeuse, STW.quant[2,], col="green", lwd =2, lty =2)

 axis(side = 4)
 title(xlab="Year", line = 2)
 title(ylab="Cases (n)", line = 2,col.lab="#F64740")
 mtext(side = 4, line = 2, 'Susceptibles (n)',col="green",cex = 0.65)
 text(2011,max(STW.quant[3,])*0.9,"Denmark",pos = 4,cex = 1.0)


# display NPI effects
par(mar=c(3,3,1,3))
t1 = (6+3)*17+ceiling(controlWeekStart);  		# from Jan 2011 to the start of NPIs
t2=Npoint                                    		# NPI effect period
t3=max(0,TWeek-t1-t2)                                	# NPI effect disappear
median.NPI=c(rep(1,t1),1-ENPI.quant[2,],rep(1,t3))
low.NPI=c(rep(1,t1),1-ENPI.quant[1,],rep(1,t3))
high.NPI=c(rep(1,t1),1-ENPI.quant[3,],rep(1,t3))

plot(timeuse,median.NPI,type="n",col="#F64740",lwd =2,ylim=c(0,2.0),xlab="",ylab="",bty = "n",xaxs="i",xlim=c(2011,2032))
polygon(c(timeuse[ceiling((52/IP)*9+controlWeekStart)],timeuse[ceiling((52/IP)*9+controlWeekStart)],timeuse[ceiling((52/IP)*9 + controlWeekStart + controlWeekLength)],
          timeuse[ceiling((52/IP)*9 + controlWeekStart + controlWeekLength)]), c(-1,2,2,-1), border = NA,col="#D1D3DF")
polygon(c(timeuse,rev(timeuse)),c(low.NPI,rev(high.NPI)),col=mycol3,border=NA)			
lines(timeuse, median.NPI, col="#F64740", lwd =2, lty =1)

abline(v = seq(2011,2100,1),lty=3,col="gray")
title(xlab="Year", line = 2)
title(ylab="relative transmission rate", line = 2,col.lab="grey32") #"#F64740")
text(2011,1.35,paste0("Denmark_NPI effect on MP"),pos = 4,cex = 0.8)

dev.off()

