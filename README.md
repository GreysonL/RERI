# RERI: UML, CML and EB Estimators and Wald Tests of Relative Excess Risk due to Interaction.

We provide an R function which returns point estimates, 95% CIs and p-values for UML,CML and EB estimator of RERI. You can simply call RERI.test(G0,G1,E0,E1,data,response.var,snp.var,main.vars,int.vars,strata.var) to obtain those estimates. The function requires loading package "CGEN". The input is as same as "snp.logistic" in CGEN package but additionally needs G0,G1,E0,E1, which means the relative excess risk due to interaction when environmental risk factor changes from E0 to E1 and genetic risk factor changes from G0 to G1 but other covariates z are held constant. In this version of code, only one genetic risk factor is allowed and it must be binary or trinary. There can be multiple environmental risk factors to interact with G and the function can deal with all types of E variable: binary, categorical and continuous. 

For the input of snp.logistic, the arguments are copied here:
data: Data frame containing all the data. No default.
response.var: Name of the binary response variable coded as 0 (controls) and 1 (cases). No default.
snp.var: Name of the SNP variable, which must be coded 0-1-2 (or 0-1). The SNP will be included as a main effect in the model. No default.
main.vars: Character vector of variable names or a formula for all covariates of interest which need to be included in the model as main effects. The default is NULL, so that only the SNP variable will be included as a main effect in the model.
int.vars: Character vector of variable names or a formula for all covariates of interest that will interact with the SNP variable. The default is NULL, so that no interactions will be in the model.
strata.var: Name of the stratification variable or a formula (see details for more info). If strata.var="SVAR", where "SVAR" is a factor or character variable in data, then "SVAR" will be treated as categorical. Otherwise, "SVAR" is treated as a continuous variable. The default is NULL (1 stratum).
