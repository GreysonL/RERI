UML, CML and EB Wald Tests of G-E interaction in multiplicative and additive scale

We provide an R function which returns point estimates, 95% CIs and p-values for UML,CML and EB estimator of relative excess of risk due tointeraction (RERI). You can simply call GE.wald.test(G0,G1,E0,E1,data,response.var,snp.var,main.vars,int.vars,strata.var,modelnum) to obtain those estimates. The function requires loading package "CGEN". The input is as same as "snp.logistic" in CGEN package but additionally needs G0,G1,E0,E1, which means the relative excess risk due to interaction when environmental risk factor changes from E0 to E1 and genetic risk factor changes from G0 to G1 but other covariates are held constant. Only one genetic risk factor is allowed and it must be binary or trinary. There can be multiple environmental risk factors to interact with G and any type of E, namely, binary, categorical and continuous can be handled in this function. Because we take advantage of snp.logistic() in the function, it may be convenient to return the summary result in multiplicative scale togenther with the result in additive scale, so that you don't have to run it twice.

For the input of snp.logistic, the arguments are copied here:
data: Data frame containing all the data. No default.

response.var: Name of the binary response variable coded as 0 (controls) and 1 (cases). No default.

snp.var: Name of the SNP variable, which must be coded 0-1-2 (or 0-1). The SNP will be included as a main effect in the model. No default.

main.vars: Character vector of variable names or a formula for all covariates of interest which need to be included in the model as main effects. The default is NULL, so that only the SNP variable will be included as a main effect in the model.

int.vars: Character vector of variable names or a formula for all covariates of interest that will interact with the SNP variable. The default is NULL, so that no interactions will be in the model.

strata.var: Name of the stratification variable or a formula (see details for more info). If strata.var="SVAR", where "SVAR" is a factor or character variable in data, then "SVAR" will be treated as categorical. Otherwise, "SVAR" is treated as a continuous variable. The default is NULL (1 stratum).

modelnum: 0-3: The genetic model for the SNP. 0=additive, 1=dominant, 2=recessive, 3=general (co-dominant).
