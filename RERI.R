library(CGEN)
GE.wald.test=function(G0,G1,E0,E1,data,response.var,snp.var,main.vars=NULL,int.vars=NULL,strata.var=NULL,modelnum){
  elevel=c()
  for(i in 1:length(int.vars)){
    if(is.factor(data[,int.vars[i]])){
      elevel=c(elevel,length(unique(data[,int.vars[i]]))-sum(is.na(data[,int.vars[i]])))
    }
    if(is.factor(data[,int.vars[i]])==F){
      elevel=c(elevel,1)
    }
  }
  if(modelnum==0|modelnum==3){glevel=length(unique(data[,snp.var]))-sum(is.na(data[,snp.var]))}
  if(modelnum==1|modelnum==2){glevel=2}
  main_e_name=c()
  for(i in 1:length(int.vars)){
    if(elevel[i]<=2){main_e_name=c(main_e_name,int.vars[i])}
    if(elevel[i]>2){
      for(j in 1:(elevel[i]-1)){
        main_e_name=c(main_e_name,paste0(int.vars[i],"_",j))
      }
    }
  }
  if(glevel==2){main_g_name=snp.var}
  if(glevel==3 & modelnum==0){main_g_name=snp.var}
  if(glevel==3 & modelnum==3){main_g_name=c(paste0(snp.var,"1"),paste0(snp.var,"2"))}
  inter_name=c()
  for(i in 1:length(main_g_name)){
    for(j in 1:length(main_e_name)){
      inter_name=c(inter_name,paste0(main_g_name[i],":",main_e_name[j]))
    }
  }
  fit=snp.logistic(data,response.var,snp.var,main.vars,int.vars,strata.var,op=list(genetic.model=modelnum))
  index=match(c(main_g_name,main_e_name,inter_name),names(fit$UML$parms))
  beta_uml=fit$UML$parms[index]
  beta_cml=fit$CML$parms[index]
  
  ############################
  if(is.null(G0)==F & is.null(G1)==F & is.null(E0)==F & is.null(E1)==F){
    E0_new=c()
    E1_new=c()
    for(i in 1:length(int.vars)){
      if(elevel[i]<=2){E0_new=c(E0_new,E0[i]);E1_new=c(E1_new,E1[i])}
      if(elevel[i]>2){
        temp0=rep(0,elevel[i]-1)
        temp1=rep(0,elevel[i]-1)
        temp0[E0[i]]=1
        temp1[E1[i]]=1
        E0_new=c(E0_new,temp0)
        E1_new=c(E1_new,temp1)
      }
    }
    G1_E1=c();G1_E0=c();G0_E1=c();G0_E0=c()
    if(modelnum==3){
      G0_new=rep(0,glevel-1)
      G0_new[G0]=1
      G1_new=rep(0,glevel-1)
      G1_new[G1]=1
      for(i in 1:(glevel-1)){
        G1_E1=c(G1_E1,G1_new[i]*E1_new)
        G1_E0=c(G1_E0,G1_new[i]*E0_new)
        G0_E1=c(G0_E1,G0_new[i]*E1_new)
        G0_E0=c(G0_E0,G0_new[i]*E0_new)
      }
    }
    if(modelnum==0){
      G0_new=G0
      G1_new=G1
      G1_E1=E1_new*G1_new
      G1_E0=E0_new*G1_new
      G0_E1=E1_new*G0_new
      G0_E0=E0_new*G0_new
    }
    if(modelnum==1){ ##dominant AA v.s. Aa+aa ##
      G0_new=min(G0,1)
      G1_new=min(G1,1)
      G1_E1=E1_new*(G1_new!=0)
      G1_E0=E0_new*(G1_new!=0)
      G0_E1=E1_new*(G0_new!=0)
      G0_E0=E0_new*(G0_new!=0)
    }
    if(modelnum==2){ ##recessive AA+Aa v.s. aa ##
      G0_new=max(G0,1)-1
      G1_new=max(G1,1)-1
      G1_E1=E1_new*(G1_new!=0)
      G1_E0=E0_new*(G1_new!=0)
      G0_E1=E1_new*(G0_new!=0)
      G0_E0=E0_new*(G0_new!=0)
    }
    if(sum((G0_new-G1_new)^2)>0 & sum((E1_new-E0_new)^2)>0){
      reri_uml=exp(beta_uml%*%c(G1_new-G0_new,E1_new-E0_new,G1_E1-G0_E0))-exp(beta_uml%*%c(G1_new-G0_new,E0_new-E0_new,G1_E0-G0_E0))-exp(beta_uml%*%c(G0_new-G0_new,E1_new-E0_new,G0_E1-G0_E0))+1
      reri_cml=exp(beta_cml%*%c(G1_new-G0_new,E1_new-E0_new,G1_E1-G0_E0))-exp(beta_cml%*%c(G1_new-G0_new,E0_new-E0_new,G1_E0-G0_E0))-exp(beta_cml%*%c(G0_new-G0_new,E1_new-E0_new,G0_E1-G0_E0))+1
      deriv_G_uml=(exp(beta_uml%*%c(G1_new-G0_new,E1_new-E0_new,G1_E1-G0_E0))-exp(beta_uml%*%c(G1_new-G0_new,E0_new-E0_new,G1_E0-G0_E0)))*(G1_new-G0_new)
      deriv_E_uml=(exp(beta_uml%*%c(G1_new-G0_new,E1_new-E0_new,G1_E1-G0_E0))-exp(beta_uml%*%c(G0_new-G0_new,E1_new-E0_new,G0_E1-G0_E0)))*(E1_new-E0_new)
      deriv_GE_uml=exp(beta_uml%*%c(G1_new-G0_new,E1_new-E0_new,G1_E1-G0_E0))*(G1_E1-G0_E0)-exp(beta_uml%*%c(G0_new-G0_new,E1_new-E0_new,G0_E1-G0_E0))*(G0_E1-G0_E0)-exp(beta_uml%*%c(G1_new-G0_new,E0_new-E0_new,G1_E0-G0_E0))*(G1_E0-G0_E0)
      B1=c(deriv_G_uml,deriv_E_uml,deriv_GE_uml)
      deriv_G_cml=(exp(beta_cml%*%c(G1_new-G0_new,E1_new-E0_new,G1_E1-G0_E0))-exp(beta_cml%*%c(G1_new-G0_new,E0_new-E0_new,G1_E0-G0_E0)))*(G1_new-G0_new)
      deriv_E_cml=(exp(beta_cml%*%c(G1_new-G0_new,E1_new-E0_new,G1_E1-G0_E0))-exp(beta_cml%*%c(G0_new-G0_new,E1_new-E0_new,G0_E1-G0_E0)))*(E1_new-E0_new)
      deriv_GE_cml=exp(beta_cml%*%c(G1_new-G0_new,E1_new-E0_new,G1_E1-G0_E0))*(G1_E1-G0_E0)-exp(beta_cml%*%c(G0_new-G0_new,E1_new-E0_new,G0_E1-G0_E0))*(G0_E1-G0_E0)-exp(beta_cml%*%c(G1_new-G0_new,E0_new-E0_new,G1_E0-G0_E0))*(G1_E0-G0_E0)
      B0=c(deriv_G_cml,deriv_E_cml,deriv_GE_cml)
      V1=fit$UML$cov[index,index]
      V0=fit$CML$cov[index,index]
      var1=t(B1)%*%V1%*%B1
      var0=t(B0)%*%V0%*%B0
      reri_eb=(reri_uml-reri_cml)^2/((reri_uml-reri_cml)^2+var1)*reri_uml+var1/((reri_uml-reri_cml)^2+var1)*reri_cml
      R1=reri_uml;R0=reri_cml
      deriv_uml=(3*R1^2-4*R0*R1+R0^2)/(var1+(R1-R0)^2)-(2*var1*R0*(R1-R0)+2*R1*(R1-R0)^3)/(var1+(R1-R0)^2)^2
      deriv_cml=(var1^2-var1*(R1-R0)^2)/(var1+(R1-R0)^2)^2
      A=c(deriv_uml,deriv_cml)
      sigma=fit$EB$UML.CML.cov[index,index]
      var=t(A)%*%matrix(c(var1,t(B1)%*%sigma%*%B0,t(B1)%*%sigma%*%B0,var0),nrow=2,ncol=2)%*%A
      EB=c(reri_eb,reri_eb-1.96*sqrt(var),reri_eb+1.96*sqrt(var),2*(1-pnorm(abs(reri_eb/sqrt(var)))))
      UML=c(reri_uml,reri_uml-1.96*sqrt(var1),reri_uml+1.96*sqrt(var1),2*(1-pnorm(abs(reri_uml/sqrt(var1)))))
      CML=c(reri_cml,reri_cml-1.96*sqrt(var0),reri_cml+1.96*sqrt(var0),2*(1-pnorm(abs(reri_cml/sqrt(var0)))))
      names(EB)=c("stat","lower","upper","p-value")
      names(UML)=c("stat","lower","upper","p-value")
      names(CML)=c("stat","lower","upper","p-value")
      add=list(UML,CML,EB)
      names(add)=c("UML","CML","EB")
      result=list(add,summary(fit))
      names(result)=c("add","mult")
    }
    if(sum((G0_new-G1_new)^2)==0|sum((E1_new-E0_new)^2)==0){
      print("Please make sure there are changes in both G and E. Be careful that if you choose a dominant genetic model,then G=1 and G=2 are equivalent; if you choose a recessive model, then G=0 and G=1 are equivalent. Only multiplicative interactions are displayed below.")
      result=summary(fit)
    }
  }
  if(is.null(G0)|is.null(G1)|is.null(E0)|is.null(E1)){
    print("Please specify G0,G1,E0,E1, otherwise only multiplicative interactions will be displayed")
    result=summary(fit)
  }
  return(result)
}
