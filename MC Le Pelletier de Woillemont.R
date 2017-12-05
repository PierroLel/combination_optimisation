## Ripley, 2003, souligne que l'AIC a été introduit pour retenir des
## variables pertinentes lors de prévisions, et que le critère BIC
## vise la sélection de variables statistiquement significative dans
## le modèle.
############################ FONCTIONS ##################################
############## Création de la base ###########
tab=function(n,p,p1) {
  p2=p-p1;
  xt=matrix(runif(n = n*p1,min=-1000,max=1000), nrow = n ,ncol =p1)
  beta=runif(n = p1,min=-100,max=100);length(beta)
  brui=rnorm(n = n,mean = 0,sd = 300)
  yt=xt%*%beta+brui
  xf1=matrix(rnorm(n = n*round(p2/4)+20,mean=300,sd=500), nrow = n ,ncol =round(p2/4)+20)
  xf2=matrix(rgamma(n = n*round(p2/4),shape = 2,scale = 100), nrow = n ,ncol =round(p2/4))
  xf3=matrix(1000*rbeta(n = n*round(p2/4),shape1 = 3.5,shape2 = 3.5,ncp = 5), nrow = n ,ncol =round(p2/4))
  xf4=matrix(100*(rt(n = n*round(p2/4),df=10)-10*rt(n = n*round(p2/4),df=100)), nrow = n ,ncol =round(p2/4))
  x=cbind(yt,xt,xf1,xf2,xf3,xf4);x=x[,1:(p+1)]
  return(x)
}

################# Noyau GIBBS ################
algo_recuit_gibbs=function(x,y) {
  p=length(x[1,]);n=length(x[,1]) # nombre de regresseurs et d'observation
  b=rbind(matrix(data=1,nrow =p1 ,ncol =1 ),matrix(data =0 ,nrow = p2,ncol = 1)) #la véritable valeur que la fonction doit trouver
  s0=rbinom(n =p,size = 1,prob = ifelse(n>p,0.5,13/p)) # état initial 
  s_0=x[,s0!=0];modele_s0=lm(yt~0+.,data = as.data.frame(s_0)) #modele associé à l'état initial
  e=BIC(modele_s0) # energie associé au modele inital
  Tp0=e/100 # temprature initiale, proportionnelle à l'energie initiale 
  Tp=Tp0
  energ=c(e) # marqueur de progression de l'energie (donc du BIC)
  ecart=c()  # marqueur de progression du nombre de composante en commun, ecart signifie en réalité nombre de composante en commun avec la solution que la fonction doit trouver
  mark_stop=c(1:10) # marqueur d'arret, si on a le meme BIC 10 itérations de suite, on peut penser qu'on est arriver et qu'il ne sert à rien de continuer à fair tourner l'algortihme
  k=1
  kmax=500 # nombre d'iteration maximale
  while ((k<=kmax)& (length(unique(mark_stop))!=1) ) { 
    sn=s0;
    for(i in 1:p){
      # composante i:=0
      s1=sn;
      s1[i]=0;
      s_1=x[,s1!=0];
      modele_s1=lm(yt~0+.,data = as.data.frame(s_1))
      e1=BIC(modele_s1);
      
      # composante i:=1
      s2=sn;
      s2[i]=1;
      s_2=x[,s2!=0];
      modele_s2=lm(yt~0+.,data = as.data.frame(s_2));
      e2=BIC(modele_s2); 
      
      prob_condi=1/(exp((e2-e1)/Tp)+1); # calcul de la proba conditionnelle
      sn[i]=rbinom(n = 1,size = 1,prob = prob_condi); # remplacement de la composante i en fonction de la proba conditionnelle associé à cette composante
    }
    s_n=x[,sn!=0];modele_sn=lm(yt~0+.,data = as.data.frame(s_n)) # nouvel état associé
    en=BIC(modele_sn) # energie associé au nouvel état
    if ((en<e) | (runif(n = 1,min = 0,max = 1)<(exp(-(en-e)/Tp)) )) { # recuit: on prend le nouvel état si il est meilleur. Si il est moins bon on le prend avec une certaine proba decroissante en fnction du temps et de l'ecart entre l'etat actuel et le nouvel etat
      e=en;s0=sn;
    }
    k=k+1
    mark_stop=c(e,mark_stop)[1:10] # mise à jour du marqueur d'arret
    energ=c(energ,e);ecart=c(ecart,nb_comm(s0,b))
    Tp=Tp/log(1+k); # mise à jour de la temperature
    plot(ecart,type="l")
  }
  return(list(s0,ecart))
}

############ Noyau Metropolis ##########
algo_recuit_metro=function(x,y) {
  p=length(x[1,]);n=length(x[,1]) # nombre de regresseurs et d'observation
  b=rbind(matrix(data=1,nrow =p1 ,ncol =1 ),matrix(data =0 ,nrow = p2,ncol = 1)) #la véritable valeur que la fonction doit trouver
  s0=rbinom(n =p,size = 1,prob = ifelse(n>p,0.5,13/p)) # état initial 
  s_0=x[,s0!=0];modele_s0=lm(yt~0+.,data = as.data.frame(s_0)) #modele associé à l'état initial
  e=BIC(modele_s0) # energie associé au modele inital
  Tp0=e/100 # temprature initiale, proportionnelle à l'energie initiale 
  Tp=Tp0
  energ=c(e) # marqueur de progression de l'energie (donc du BIC)
  ecart=c()  # marqueur de progression du nombre de composante en commun, ecart signifie en réalité nombre de composante en commun avec la solution que la fonction doit trouver
  mark_stop=c(1:(p*3)) # marqueur d'arret; si on a le meme BIC 3*p itérations de suite, on peut penser qu'on est arriver et qu'il ne sert à rien de continuer à fair tourner l'algortihme
  k=0
  if (p<40) {N=2000} else {N=0.4755*p^2-25*p+9000} #adaptation du nombre d'iteration maximal en fonction de la taille de la base
  while ( (Tp>Tp0*(0.99**N)) & (length(unique(mark_stop))!=1) ) { 
    sn=s0
    for ( j in 1:round(max(round(p/50),1))) {
      dumm=sample(1:p,1)
      sn[dumm]=abs(sn[dumm]-1) # on change une (ou pls) composante au hasard. une ou plusieurs si p est très grand
    }
    s_n=x[,sn!=0];modele_sn=lm(y~0+.,data = as.data.frame(s_n)) # nouvel etat voisin
    en=BIC(modele_sn) # energie associee au nouvel état, si nombre de regresseurs > nombre d'observation on peut avoir des BIC=-infini, la boucle while ci dessou permet d'eviter de consider ces cas
    while (en==-Inf) {
      sn=s0
      for ( j in 1:round(max(round(p/50),1))) {
        sn[sample(1:p,1)]=abs(sn[sample(1:p,1)]-1) ## on change plusieurs composante au hasard, on determine le nombre de composante à change en fonction de la taille de la base. Si p=50, on change une composante par une composante; mais si p=5000 changer une composante par une composante est trop long, donc on en change 20 par 20 (par exemple) dans une premier temps
      }
      s_n=x[,sn!=0];modele_sn=lm(y~0+.,data = as.data.frame(s_n)) # nouvel etat voisin
      en=BIC(modele_sn) # energie associee au nouvel état
    }
    if ((en<e) | (runif(n = 1,min = 0,max = 1)<exp(-(en-e)/Tp))) { # recuit: on prend le nouvel état si il est meilleur. Si il est moins bon on le prend avec une certaine proba decroissante en fnction du temps et de l'ecart entre l'etat actuel et le nouvel etat
      e=en;s0=sn;
    }
    Tp=Tp*0.99      #decroissance log de la temperature
    k=k+1
    energ=c(energ,e);ecart=c(ecart,nb_comm(s0,b))
    if (k%in%(c(1:(N/100))*100)) {
      plot(ecart,type="l")
      print(paste(round(100*k/N),"%"))
    }
    mark_stop=c(e,mark_stop)[1:(3*p)] # historique des 3p dernier e choisit (pour eviter de trouner dans le vide vers la fin de l'algo)
  }
  mark_stop=c(1:(p*10),mark_stop) # on est arrivé  un stade ou on est proche de la solution et donc il faut affiner en ne changeant qu'1 seule composante, et ce quel que soit p.
  # le boucle while ci dessous est construite exactement comme celle d'au dessus, à la seule difference que dans celle ci on change uniquement une seule composante par une composante
  while ( (Tp>Tp0*(0.99**N)) & (length(unique(mark_stop))!=1) ) { 
    sn=s0
    sn[sample(1:p,1)]=abs(sn[sample(1:p,1)]-1)
    s_n=x[,sn!=0];modele_sn=lm(y~0+.,data = as.data.frame(s_n)) # voisin
    en=BIC(modele_sn)
    while (en==-Inf) {
      sn=s0
      sn[sample(1:p,1)]=abs(sn[sample(1:p,1)]-1)
      s_n=x[,sn!=0];modele_sn=lm(y~0+.,data = as.data.frame(s_n)) # voisin
      en=BIC(modele_sn)
    }
    if ((en<e) | (runif(n = 1,min = 0,max = 1)<exp(-(en-e)/Tp))) {
      e=en;s0=sn;
    }
    Tp=Tp*0.99      #decroissance log de la temperature
    k=k+1
    energ=c(energ,e);ecart=c(ecart,nb_comm(s0,b))
    if (k%in%(c(1:(N/100))*100)) {
      par(mfrow = c(1,1))
      plot(ecart,type="l")
      print(paste(round(100*k/N),"%"))
    }
    mark_stop=c(e,mark_stop)[1:(13*p)]
  }
  return(list(s0,ecart))
}

############ algo CE ##########
algo_ce=function(x,y) {
  p=length(x[1,]);n=length(x[,1]) # nombre de regresseurs et d'observation
  b=rbind(matrix(data=1,nrow =p1 ,ncol =1 ),matrix(data =0 ,nrow = p2,ncol = 1))
  theta0=matrix(data = 0.5,nrow = p,ncol = 1) # parametre initale de la bernoulli
  N=p*7 # definition de la taille de l'échantillon simulé
  tmax=15;t=1; # nombre max d'iteration
  distance=c() # marqueur de progression
  while ((t<tmax)&(sum(theta0)!=sum(theta0!=0))) { # la boucle s'arretera quand le parametre de la Bernoulli sera un vecteur binaire
    theta1=theta0
    ## parametres initiaux de la boucle d'echantillonnage
    BIC1=c();X=data.frame(X_=1:p)
    ## Création de l'échantillon de taille N de vecteur binaire suivant une bernoulli de parametre theta1
    for (i in 1:N) {
      X_=rbinom(n =p,size = 1,prob = theta1)
      X=data.frame(X,X_)
      x_i=x[,X_!=0]
      modele_i=lm(y~0+.,data = as.data.frame(x_i))
      while (BIC(modele_i)==-Inf) { #cete boucle while permet d'éviter d'avoir des BIC=-infini (ce qui arrive lorsque nombre regresseurs > nombre observations)
        X_=rbinom(n =p,size = 1,prob = theta1)
        X=data.frame(X,X_)
        x_i=x[,X_!=0]
        modele_i=lm(y~0+.,data = as.data.frame(x_i))
        }
      BIC1=c(BIC1,BIC(modele_i))
      }
    X=X[,-1]
    ## trie en fonction du BIC
    score=data.frame(id=1:N,bic=BIC1)
    score=score[order(score[,"bic"],decreasing=F), ] 
    if (p<100) { # dertermination de la taille de l'échantillon "preferentiel", en fonction de la taille de l'échantillon générée au dessus
      rho=round((log(p)/p)*N)
      } else {
        rho=round(N-(1-0.001)*N)
        }
    for (j in 1:p) { # estimation du nouveau paramtre de la loi bernoulli (theta0)
      pj=0
      for (i in 1:rho) {
        pj=pj+X[j,paste("X_.",score[i,1],sep="")]
        }
      theta0[j]=(1/rho)*pj
      }
    #print(mean(BIC1));print(sum(theta0));print(head(theta0,6));print(sum(theta0!=0)) #si on souhaite suivre en détails l'avancé de l'algo
    t=t+1
    norme_p(b,theta0)
    distance=c(distance,norme_p(b,theta0)) # marqueur de progression
    plot(distance,type="l")
  }
  return(list(theta0,distance))
}

########## nombre en commun ########
nb_comm=function(x,y) {
  return(length(x)-sum(abs(x-y)))
}

########## norme P ###########
norme_p=function(a,b) {
  (sum(abs(a-b)**length(a)))**(1/length(a))
}

############################ COMPARAISON ##############################
N=1;n=2000;p=50;p1=5;p2=p-p1
Ecart_ce=data.frame(id=1:5)
Temps_ce=c()
Resultats_ce=data.frame(matrix(data = 0,nrow = p,ncol = N))

Ecart_recuit_metro=data.frame(id=1:5)
Temps_recuit_metro=c()
Resultats_recuit_metro=data.frame(matrix(data = 0,nrow = p,ncol = N))

Ecart_recuit_gibbs=data.frame(id=1:5)
Temps_recuit_gibbs=c()
Resultats_recuit_gibbs=data.frame(matrix(data = 0,nrow = p,ncol = N))

for (i in 1:N) {
  X=tab(n = n,p = p,p1 = p1)
  yt=X[,1];x=X[,2:(p+1)]
  
  # Algo CE
  T1_ce=Sys.time() 
  resul_ce=algo_ce(x = x,y = yt)
  T2_ce=Sys.time()
  Temps_ce=c(Temps_ce,difftime(T2_ce,T1_ce,units = "secs"))
  ecart_ce_i=data.frame(id=1:length(resul_ce[[2]]),distance=resul_ce[[2]])
  Ecart_ce=merge(Ecart_ce,ecart_ce_i,by="id",all.x=T,all.y=T)
  Resultats_ce[,i]=resul_ce[[1]]
  
  # Algo Metro
  T1_recuit_metro=Sys.time() 
  resul_recuit_metro=algo_recuit_metro(x = x,y = yt)
  T2_recuit_metro=Sys.time()
  Temps_recuit_metro=c(Temps_recuit_metro,difftime(T2_recuit_metro,T1_recuit_metro,units = "secs"))
  ecart_recuit_metro_i=data.frame(id=1:length(resul_recuit_metro[[2]]),distance=resul_recuit_metro[[2]])
  Ecart_recuit_metro=merge(Ecart_recuit_metro,ecart_recuit_metro_i,by="id",all.x=T,all.y=T)
  Resultats_recuit_metro[,i]=resul_recuit_metro[[1]]
  
  # Algo Gibbs
  T1_recuit_gibbs=Sys.time() 
  resul_recuit_gibbs=algo_recuit_gibbs(x = x,y = yt)
  T2_recuit_gibbs=Sys.time()
  Temps_recuit_gibbs=c(Temps_recuit_gibbs,difftime(T2_recuit_gibbs,T1_recuit_gibbs,units = "secs"))
  ecart_recuit_gibbs_i=data.frame(id=1:length(resul_recuit_gibbs[[2]]),distance=resul_recuit_gibbs[[2]])
  Ecart_recuit_gibbs=merge(Ecart_recuit_gibbs,ecart_recuit_gibbs_i,by="id",all.x=T,all.y=T)
  Resultats_recuit_gibbs[,i]=resul_recuit_gibbs[[1]]
}

Resultat_p_50_p1_5=list(CE=list(Ecart=Ecart_ce,Temps=Temps_ce,Resultats=Resultats_ce),
                        Metro=list(Ecart=Ecart_recuit_metro,Temps=Temps_recuit_metro,Resultats=Resultats_recuit_metro),
                        Gibbs=list(Ecart=Ecart_recuit_gibbs,Temps=Temps_recuit_gibbs,Resultats=Resultats_recuit_gibbs))
save(Resultat_p_50_p1_5,file="Resultat_p_50_p1_5")
save(Resultat_p_100_p1_10,file="Resultat_p_100_p1_10")
save(Resultat_p_150_p1_15,file="Resultat_p_150_p1_15")
## la boucle d'au dessus peut etre très longue à faire tourner selon les paramtres, N et p choisit
## c'est pour cette raison qu'on a sauvegarder les resultats dans des listes.

load("Resultat_p_50_p1_5")
load("Resultat_p_100_p1_10")
load("Resultat_p_150_p1_15")

## Boxplots liÃ©s au temps

#p50, p1 = 5

boxplot(Resultat_p_50_p1_5$CE$Temps,Resultat_p_50_p1_5$Gibbs$Temps,Resultat_p_50_p1_5$Metro$Temps,names=c("CE","Gibbs","Metropolis"),border=c("blue","red","green"),col="wheat",ylab="Temps de RÃ©solution (en sec)",main="p=50 | p1=5", las=1)
abline(h=mean(Resultat_p_50_p1_5$CE$Temps), col="blue")
abline(h=mean(Resultat_p_50_p1_5$Metro$Temps), col="green")
abline(h=mean(Resultat_p_50_p1_5$Gibbs$Temps), col="red")

#p100, p1 = 10

boxplot(Resultat_p_100_p1_10$CE$Temps,Resultat_p_100_p1_10$Gibbs$Temps,Resultat_p_100_p1_10$Metro$Temps,names=c("CE","Gibbs","Metropolis"),border=c("blue","red","green"),col="wheat",ylab="Temps de RÃ©solution (en sec)",main="p=100 | p1=10",las=1)
abline(h=mean(Resultat_p_100_p1_10$CE$Temps), col="blue")
abline(h=mean(Resultat_p_100_p1_10$Metro$Temps), col="green")
abline(h=mean(Resultat_p_100_p1_10$Gibbs$Temps), col="red")

#p150, p1 = 15

boxplot(Resultat_p_150_p1_15$CE$Temps,Resultat_p_150_p1_15$Gibbs$Temps,Resultat_p_150_p1_15$Metro$Temps,names=c("CE","Gibbs","Metropolis"),border=c("blue","red","green"),col="wheat",ylab="Temps de RÃ©solution (en sec)",main="p=150 | p1=15",las=1)
abline(h=mean(Resultat_p_150_p1_15$CE$Temps), col="blue")
abline(h=mean(Resultat_p_150_p1_15$Metro$Temps), col="green")
abline(h=mean(Resultat_p_150_p1_15$Gibbs$Temps), col="red")

## Score d'exactitude 

#p1=5

true_vect1=rep(1,5)
true_vect2=rep(0,45)
true_vect=c(true_vect1,true_vect2)

##CE
score_CE_5 <- NULL
for (i in 1:30) {
  score_CE_5[i] <- nb_comm(true_vect,Resultat_p_50_p1_5$CE$Resultats[i])
}

##Metropolis
score_met_5 <- NULL
for (i in 1:30) {
  score_met_5[i] <- nb_comm(true_vect,Resultat_p_50_p1_5$Metro$Resultats[i])
}

##Gibbs
score_gibbs_5 <- NULL
for (i in 1:30) {
  score_gibbs_5[i] <- nb_comm(true_vect,Resultat_p_50_p1_5$Gibbs$Resultats[i])
}

#p1=15

true_vect1=rep(1,15)
true_vect2=rep(0,35)
true_vect=c(true_vect1,true_vect2)

##CE
score_CE_15 <- NULL
for (i in 1:30) {
  score_CE_15[i] <- nb_comm(true_vect,Resultat_p_50_p1_15$CE$Resultats[i])
}

##Metropolis
score_met_15 <- NULL
for (i in 1:30) {
  score_met_15[i] <- nb_comm(true_vect,Resultat_p_50_p1_15$Metro$Resultats[i])
}

##Gibbs
score_gibbs_15 <- NULL
for (i in 1:30) {
  score_gibbs_15[i] <- nb_comm(true_vect,Resultat_p_50_p1_15$Gibbs$Resultats[i])
}

#p1=25

true_vect1=rep(1,25)
true_vect2=rep(0,25)
true_vect=c(true_vect1,true_vect2)

##CE
score_CE_25 <- NULL
for (i in 1:30) {
  score_CE_25[i] <- nb_comm(true_vect,Resultat_p_50_p1_25$CE$Resultats[i])
}

##Metropolis
score_met_25 <- NULL
for (i in 1:30) {
  score_met_25[i] <- nb_comm(true_vect,Resultat_p_50_p1_25$Metro$Resultats[i])
}

##Gibbs
score_gibbs_25 <- NULL
for (i in 1:30) {
  score_gibbs_25[i] <- nb_comm(true_vect,Resultat_p_50_p1_25$Gibbs$Resultats[i])
}

#p1=40

true_vect1=rep(1,40)
true_vect2=rep(0,10)
true_vect=c(true_vect1,true_vect2)

##CE
score_CE_40 <- NULL
for (i in 1:30) {
  score_CE_40[i] <- nb_comm(true_vect,Resultat_p_50_p1_40$CE$Resultats[i])
}

##Metropolis
score_met_40 <- NULL
for (i in 1:30) {
  score_met_40[i] <- nb_comm(true_vect,Resultat_p_50_p1_40$Metro$Resultats[i])
}

##Gibbs
score_gibbs_40 <- NULL
for (i in 1:30) {
  score_gibbs_40[i] <- nb_comm(true_vect,Resultat_p_50_p1_40$Gibbs$Resultats[i])
}

## Calculs des moyennes 

mean(score_CE_5)
mean(score_CE_15)
mean(score_CE_25)
mean(score_CE_40)

mean(score_gibbs_5)
mean(score_gibbs_15)
mean(score_gibbs_25)
mean(score_gibbs_40)

mean(score_met_5)
mean(score_met_15)
mean(score_met_25)
mean(score_met_40)

## Calcul des Ã©carts-types

sd(score_CE_5)
sd(score_CE_15)
sd(score_CE_25)
sd(score_CE_40)

sd(score_gibbs_5)
sd(score_gibbs_15)
sd(score_gibbs_25)
sd(score_gibbs_40)

sd(score_met_5)
sd(score_met_15)
sd(score_met_25)
sd(score_met_40)

############################ Application Vraie base ##################################
## on test l'aglo recuit metropolis
data<- read.csv("OnlineNewsPopularity.csv");x = data[,3:60];y = data[,61]
p=length(x[1,])
s0=rbinom(n =p,size = 1,prob = ifelse(n>p,0.5,13/p)) ## 13/p pour eviter bic=-inf dès que p>100
s_0=x[,s0!=0];modele_s0=lm(y~.,data = as.data.frame(s_0)) #etat initial
e=BIC(modele_s0)  #energie initiale
Tp=e/1000          #temperature initiale, elle est proportionnel au bic en gros
Tp0=e/1000
energ=c(e)        #marqueur de progression
ecart=c()         #marqueur de progression
mark_stop=c(1:(p*3))
k=0
if (p<40) {N=2000} else {N=0.4755*p^2-25*p+9000}
while ( (Tp>Tp0*(0.99**N)) & (length(unique(mark_stop))!=1) ) { 
  sn=s0
  for ( j in 1:round(max(round(p/50),1))) {
    sn[sample(1:p,1)]=abs(sn[sample(1:p,1)]-1) # on change une (ou pls) composante au hasard. une ou plusieurs si p est très grand
  }
  s_n=x[,sn!=0];modele_sn=lm(y~.,data = as.data.frame(s_n)) # voisin
  en=BIC(modele_sn)
  if ((en<e) | (runif(n = 1,min = 0,max = 1)<exp(-(en-e)/Tp))) {
    e=en;s0=sn;
  }
  Tp=Tp*0.99      #decroissance log de la temperature
  k=k+1
  energ=c(energ,e)
  mark_stop=c(e,mark_stop)[1:(3*p)]
}
mark_stop=c(1:(p*10),mark_stop) # on est arrivé  un stade ou on est proche de la solution et donc il faut affiner en ne changeant qu'1 seule composante, et ce quel que soit p
while ( (Tp>Tp0*(0.99**N)) & (length(unique(mark_stop))!=1) ) { 
  sn=s0
  sn[sample(1:p,1)]=abs(sn[sample(1:p,1)]-1)
  s_n=x[,sn!=0];modele_sn=lm(y~.,data = as.data.frame(s_n)) # voisin
  en=BIC(modele_sn)
  while (en==-Inf) {
    sn=s0
    sn[sample(1:p,1)]=abs(sn[sample(1:p,1)]-1)
    s_n=x[,sn!=0];modele_sn=lm(y~.,data = as.data.frame(s_n)) # voisin
    en=BIC(modele_sn)
  }
  if ((en<e) | (runif(n = 1,min = 0,max = 1)<exp(-(en-e)/Tp))) {
    e=en;s0=sn;
  }
  Tp=Tp*0.99      #decroissance log de la temperature
  k=k+1
  energ=c(energ,e);
  if (k%in%(c(1:(N/100))*100)) {
    par(mfrow = c(1,1))
    plot(energ,type="l")
    print(paste(round(100*k/N),"%"))
  }
  mark_stop=c(e,mark_stop)[1:(13*p)]
}

s_0=x[,s0!=0];modele_s0=lm(y~.,data = as.data.frame(s_0)) #etat initial
summary(modele_s0) 
## le probleme d'une application de nos algorithme à de vraies "grosses" bases est le fait que
## dans de telle base il y a des variables binaire qui ne se traitent pas de la meme façon que des variables continues, ou meme discretes
