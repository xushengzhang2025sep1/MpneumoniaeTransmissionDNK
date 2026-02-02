# MpneumoniaeTransmissionDNK
Data and Code for Improving prediction of future Mycoplasma pneumoniae epidemics using transmission-based models

This readme file provides an explanation for data, model information and code used in the Letter

Improving prediction of future Mycoplasma pneumoniae epidemics using transmission-based models

by Xu-Sheng Zhang, Semjon Sidorov, Ursula Dalrymple, Hanne-Dorthe Emborg, Søren Anker Uldum, Michael L. Beeton, Patrick M. Meyer Sauteur.

Data

National Mycoplasma pneumoniae detection data, along with birth and population data for Denmark spanning January 2011 to March 2025, are provided in the file:
	PNAS_Letter_ESGMAC-MAPS_Repository_Data-Denmark-Jan-2011-Dec-2025.csv

Model

To investigate the effect of NPIs on the transmission dynamics of M. pneumoniae, we assume that transmission rate changed as:

beta(t) = betahat(t)            when t<t_sNPI;

beta(t) = delta(t)*betahat(t)   when t_sNPI <= t < March 2025

beta(t) = betahat(t)            when t>= tMarch 2025   	(1)

and reporting rate changed as:

rho(t) = rhohat(t)                when t <= t_sNPI

rho(t) = epsilon*rhohat(t)        when t > t_sNPI        (2)

Here betahat(t) represents the seasonal transmission rate and rhohat(t) the reporting rate from January 2011 to the start of NPIs obtained by applying TSIR model. The parameter tsNPI represents the start time of NPIs (and teNPI the end date of NPIs). The parameter delta(t) represents the relative transmission rate reflecting the effects of NPIs on the transmission of M. pneumoniae and epsilon the change in reporting rate from the start of NPIs, which is simply assumed to be a constant. March 2025 is the date of most recent data collection of the ESGMAC M. pneumoniae Surveillance (MAPS) study. 

In this study, we use the generation time of 3 weeks of M. pneumoniae as the time step and one year is approximated as 17 triweeks. Consequently, the relative transmission during the period [tsNPI, March 2025] is approximated by a vector of Npoint= (March2025-TsNPI)/3weeks point: delta(t), t=1,2,…Npoint.  The Npoint +1 parameters (delta(t) and epsilon) are estimated using Bayesian inference via MCMC sampling. As the purpose of NPIs suggested (reducing the contact rate between people), we assume the prior distribution of delta(t), t = tsNPI,…,teNPI, during the NPI period as U[0,1] (i.e., unform distribution from 0 to 1), and their prior distributions from NPI lifting to March 2025 (i.e., t = teNPI+1, …,March 2025) as U[0,1.5] in view of the potential large fluctuation in population contact after NPI lifting. The prior distribution of epsilon is set as U[1,4].

Codes

	NPIeffectModel.R
	generates the transmission dynamics incorporating the NPI effect on transmission rate, which will be called in the following two R codes:

	MP_NPIs_MCMCsampling.R
	generates the posterior distribution of NPI effects and reporting change in equations (1) and (2). The output (estimates of model parameters) is saved in NPI_effect_Denmark.RData.

	prediction.R
	generates the predictions of M pneumoniae detections from 2011 to 2032 using estimates of model parameters saved in NPI_effect_Denmark.RData.
