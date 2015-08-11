

#Function to return sensitivity of lambda to age-specific mortality from a
#stage-classified model
#Code adapted  by Rob Salguero-Gomez from Hal Caswell's matlab scripts
#Last modification: Aug 11th 2015

library(Matrix)
library(Deducer)
 
matU=matrix(c(.1,0,0,.3,0.2,0,0,.3,.4),nrow=3,byrow=T)
matF=matrix(c(0,0,30,0,0,0,0,0,0),nrow=3,byrow=T)



vecperm <-function(n,m){

P = matrix(rep(0,(m*n)^2),nrow=m*n)
I = diag(m*n)

k = 1
for (i in 1:n){
    for (j in 1:m*n){
        P[k,] = I[j,]
        k = k+1
    }
}

return(P)
}



ageStage <- function(matU,matF){
 
    #matU=transition matrix
    #matF=fertility matrix
    #matA=projection matrix

        matA = matU + matF 
        s = dim(matA)[1]
        sigma = colSums(matU)
        mort = 1 - sigma
 
    #Survival matrix
        S=diag(sigma)
    #Conditional transition matrix
        G=matU*diag(1/sigma)
 
    #Fundamental matrix
        N=solve(diag(s)-matU)

    #Mean life expectancy
        meaneta=colSums(N)
 
    #Maximum age= omega
        om = 2
    #Proportion of stable population in age omega
        prop = 1
 
    while (prop>.01){
    #while (om<101){
        om=om+1
        
        print(om)
        
        #Create some necessary matrices
        #Sparse identity matrix
        #Iom=speye(om);
            Iom = diag(om)
            Iom = Matrix(c(Iom),byrow=T,nrow=om,sparse=T)
        #Is=speye(s)
            Is = diag(s)
            Is = Matrix(c(Is),byrow=T,nrow=s,sparse=T)
        #Isom=speye(s*om);
            Isom = diag(s*om)
            Isom = Matrix(c(Isom),byrow=T,nrow=s*om,sparse=T)
        
        #Age transition matrix (no mortality)
            du = Matrix(rep(0,s*s),nrow=s,sparse=T)
            du[row(du)-1 == col(du)] = 1
            du[dim(du)[1],dim(du)[2]]=1
        
        #Age transition matrix for reproduction
            #df=sparse(zeros(om));
            df = Matrix(rep(0,om*om),nrow=om,sparse=T)
            df[1,]=1
        
        #block diagonal matrices
        #U=kron(Iom,u);
            U = kronecker(Iom, matU, make.dimnames = TRUE)
        #F=kron(Iom,f);
            F = kronecker(Iom, matF, make.dimnames = TRUE)
        #Du=kron(Is,du);
            Du = kronecker(Is,du)
        #Df=kron(Is,df);
            Df = kronecker(Is,df)
        #A=kron(Iom,a);
            A = kronecker(Iom,matA)
        
        #Vec-permutation matrix
        #K=spvecperm(s,om);
        #K = vecperm(s,om)
            K = vecperm(s,om)
        
        #Big transient matrix
            Utilde = t(K)*Du*K*U
        
        #Big fertility matrix
            Ftilde = t(K)*Df*K*F
        
        #Big projection matrix
            Atilde = Utilde + Ftilde
        
        #clear some un-needed matrices
        rm(Utilde,Ftilde)
        
        #eigenvalues and eigenvectors
        
        lambda = max(eigen(Atilde)$values)
        pick = which(eigen(Atilde)$values==max(eigen(Atilde)$values))
        w = eigen(Atilde)$vectors[,pick]
        pick = which(eigen(t(Atilde))$values==max(eigen(t(Atilde))$values))
        v = eigen(t(Atilde))$vectors[,pick]
            
        #Rescale eigenvectors
        w=w/(t(w)*v)
        w[is.nan(w)]=0
        
        #Proportion of population in last age class
        xx=matrix(w/sum(w),s,om)
        prop=sum(xx[,dim(xx)[2]])
        
        #Continue loop until that proportion gets small
        
    }


    #Sensitivity of lambda to A
        dlam_dvecA=kronecker(t(w),t(v))
 
    #Sensitivity of U to survival
        ones=matrix(1,s,1)
        dvecu_dtheta=kronecker(Is,G)*diag(as.numeric(Is))*kronecker(ones,Is)
 
    #some constant matrices needed
        Qu = t(K)*Du*K
        aa = kronecker(Isom,Qu)
        bb = kronecker(as.numeric(Iom),speye(s^2))
        cc = kronecker(K,Is)
 
    #Generate sensitivity of lambda to survival each age
        for (i in 1:om) {
            dlam_dtheta[i,] = dlam_dvecA*aa*kronecker((Iom(:,i)*Iom[i,]),cc)*bb*dvecu_dtheta
            elam_etheta[i,] = (1/lambdamax)*dlam_dtheta[i,]*diag(sigma)
        }
 
    out$u=u
    out$f=f
    out$G=G
    out$S=S
    out$dlam_dtheta=dlam_dtheta
    out$elam_etheta=elam_etheta
    out$v=v
    out$w=w
    out$prop=prop
    return(out)
}



set(0,'DefaultAxesColorOrder',[0 0 0],...
     'DefaultAxesLineStyleOrder','-o|-v|-s|-^',...
     'DefaultAxesFontSize',10)
colormap gray
 
%Sensitivity of lambda to survival in each stage
 figure
 bar(sum(dlam_dtheta(:,:,1)))
 xlabel('Stage')
 ylabel('Sensitivity')
 colormap gray
 nameSave=strcat(species, '_fig1')
 print('-depsc','-tiff','-r300',nameSave)
 
%Elasticity of lambda to survival in each stage
 figure
 bar(sum(elam_etheta(:,:,1)))
 xlabel('Stage')
 ylabel('Elasticity')
 colormap gray
 nameSave=strcat(species, '_fig2')
 print('-depsc','-tiff','-r300',nameSave)
 
%Sensitivity of lambda to survival in stage
 figure
 bar(dlam_dtheta(:,:,1)')
 xlabel('Stage')
 ylabel('Sensitivity of \lambda to survival')
 nameSave=strcat(species, '_fig3')
 print('-depsc','-tiff','-r300',nameSave)
 
%Elasticity of lambda to survival in stage
 figure
 bar(elam_etheta(:,:,1)')
 xlabel('Stage')
 ylabel('Elasticity of \lambda to survival')
 nameSave=strcat(species, '_fig4')
 print('-depsc','-tiff','-r300',nameSave)
 
%Sensitivity of lambda to survival in each stage and age
 figure
 bar3(fliplr(dlam_dtheta(:,:,1)))
 set(gca,'xticklabel',(7:-1:1))
 set(gca,'ylim',[0  om+1])
 ylabel('Age')
 xlabel('Stage')
 zlabel('Sensitivity')
 colormap(flipud(jet))
 nameSave=strcat(species, '_fig5')
 print('-depsc','-tiff','-r300',nameSave)
     
 %Elasticity of lambda to survival in each stage and age
 figure
 bar3(fliplr(elam_etheta(:,:)))
 set(gca,'xticklabel',(7:-1:1))
 set(gca,'ylim',[0  om+1])
 colormap(flipud(jet))
 ylabel('Age')
 xlabel('Stage')
 zlabel('Elasticity')
 nameSave=strcat(species, '_fig6')
 print('-depsc','-tiff','-r300',nameSave)
 
nameSave=strcat(species, '_sens2.csv')
csvwrite(nameSave,dlam_dtheta)
nameSave=strcat(species, '_elas')
csvwrite(nameSave,elam_etheta)
nameSave=strcat(species, '_maneta')
csvwrite(nameSave,maneta)
 
 


