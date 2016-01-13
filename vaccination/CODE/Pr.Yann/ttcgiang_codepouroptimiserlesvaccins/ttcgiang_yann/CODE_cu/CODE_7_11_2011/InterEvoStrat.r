############################## Interface de VACCINATION utilisant l'algorithme (1+1)-ES  ###################################################################
# Le but de cette fonction est de fonctionner le programmes "../EvolStrat.cpp"                                                                            ##
# qui evalue les politiques de vaccinations en utilisant l'algorithme (1+1)-ES                                                                            ##
# Par      TRAN Thi Cam Giang - P15 - IFI - Hanoi - Vietname                                                                                              ##
############################################################################################################################################################
InterEvoStrat <- function(
			tmax = 2000,initN=10000,initI=100,nu=0.01,rmu=1/(10*365),beta=1000/365,sigma2 = 0, sigma=1/7,
			gamma=1/7, epsilon =0.10, nbVilles = 1,topology=1,graine=10, villeDessiner=1,
			periode =90, tstart = 0, tfinal = 1000, tauxVaccToutesVilles = 80, nombreMaxVaccins = 2000,nStrategies = 5, nPolitiques =10,
			EvoStratInput = "EvoStrat_Input.txt",EvoStratOutput = "Output/EvoStratOutput.txt",...
			) {

tmax0 = tmax
initN0=initN
initI0=initI
nu0=nu
rmu0=rmu
beta0=beta
sigma0=sigma
gamma0=gamma
epsilon0=epsilon
nbVilles0=nbVilles
topology0=topology
graine0=graine
villeDessiner0=villeDessiner
tstart0=tstart
tfinal0=tfinal
periode0=periode
tauxVaccToutesVilles0=tauxVaccToutesVilles
nombreMaxVaccins0 =  nombreMaxVaccins;
nStrategies0 = nStrategies;
nPolitiques0 = nPolitiques;
EvoStratInput0 = EvoStratInput
EvoStratOutput0 = EvoStratOutput

lAr=23
lNouVal = 25
lValDef = 25
lSig = 50

parametresDeFaut<- c("Argument            |   Valeur défaut          |     Signification ",
        	paste(stringFixe("tmax",lAr),stringFixe(tmax0,lValDef),stringFixe("Temps pour fonctionner une simulation",lSig),sep=""),
		paste(stringFixe("initN",lAr),stringFixe(initN0,lValDef),stringFixe("Nombre de population initial",lSig),sep=""),
		paste(stringFixe("initI",lAr),stringFixe(initI0,lValDef),stringFixe("Nombre de personnes infectées initial",lSig),sep=""),
		paste(stringFixe("rmu",lAr),stringFixe(rmu0,lValDef),stringFixe("Taux de mortalite et aussi le taux de connaissance",lSig),sep=""),
		paste(stringFixe("beta",lAr),stringFixe(beta0,lValDef),stringFixe("Taux de contact S->E initial",lSig),sep=""),
		paste(stringFixe("sigma",lAr),stringFixe(sigma0,lValDef),stringFixe("Taux de transmission E->I",lSig),sep=""),
		paste(stringFixe("gamma",lAr),stringFixe(gamma0,lValDef),stringFixe("Taux de transmission I->R",lSig),sep=""),
		paste(stringFixe("nbVilles",lAr),stringFixe(nbVilles0,lValDef),stringFixe("Nombre de villes",lSig),sep=""),
		paste(stringFixe("epsilon",lAr),stringFixe(epsilon0,lValDef),stringFixe("Taux de transmission entre deux villes",lSig),sep=""),
		paste(stringFixe("topology",lAr),stringFixe(topology0,lValDef),stringFixe("Type de topology entre des villes",lSig),sep=""),
		paste(stringFixe("graine",lAr),stringFixe(graine0,lValDef),stringFixe("Chiffre de graine",lSig),sep=""),
		paste(stringFixe("villeDessiner",lAr),stringFixe(villeDessiner0,lNouVal),stringFixe("Quelle ville que vous voulez saver ses informations",lSig),sep=""),
		paste(stringFixe("nu",lAr),stringFixe(nu0,lValDef),stringFixe("Taux de l'infection a l'exterieur",lSig),sep=""),
                paste(stringFixe("tstart",lAr),stringFixe(tstart0,lValDef),stringFixe("Temps ou on commence a faire la vaccination",lSig),sep=""),
		paste(stringFixe("tfinal",lAr),stringFixe(tfinal0,lValDef),stringFixe("Temps ou on finit la vaccination",lSig),sep=""),
		paste(stringFixe("periode",lAr),stringFixe(periode0,lValDef),stringFixe("Temps entre deux fois de vaccinations",lSig),sep=""),          
		paste(stringFixe("tauxVaccToutesVilles",lAr),stringFixe(tauxVaccToutesVilles0,lValDef),stringFixe("(%S) - Chiffre de graine du pourcentage de susceptibles vaccinees",lSig),sep=""),
		paste(stringFixe("nombreMaxVaccins",lAr),stringFixe(nombreMaxVaccins0,lValDef),stringFixe("Nombre de vaccins maximal",lSig),sep=""),
		paste(stringFixe("nStrategies",lAr),stringFixe(nStrategies0,lValDef),stringFixe("Nombre de stratefies",lSig),sep=""),
		paste(stringFixe("nPolitiques",lAr),stringFixe(nPolitiques0,lValDef),stringFixe("Nombre de politiques dans une strategie",lSig),sep=""),
		paste(stringFixe("EvoStratInput",lAr),stringFixe(EvoStratInput0,lValDef),stringFixe("Fichier qui contient 5 valeurs de 5 elements du vecteur politique(w0,w1,w2,w3,w4)",lSig),sep=""),
		paste(stringFixe("EvoStratOutput",lAr),stringFixe(EvoStratOutput0,lValDef),stringFixe("Fichier qui contient le resulat apres avoir execute le programme",lSig),sep="")
	       )

menu(parametresDeFaut,graphics=FALSE,title="\n table 1: VALEURS DES PARAMETRES PAR ANNEE, par defaut");
print("*********************************************************************************")	
switch(menu(c("S   Selectionner les paramatres par defaut", "M   Modifier parametres"),graphics=FALSE,title="\n table 2: DONNER votre choix!") + 1,  cat("Rien de fait\n"), "S", "M")->choix
if(choix == "S")
{
print("Vous avez choisi: S - Selectionner les paramatres par defaut!")
print("Le programme est en train de marcher....")
}
else if (choix == "M")
{
print("Tapez quelle parametre que vous modifiez sa valeur!")
print("Note : - Parametres de Simulation SEIR: tmax, initN, initI, rmu, beta0, sigma2, sigma, gamma, nbVilles, epsilon, topology, villeDessiner")
print("       - Type de parametres est un chiffre ou une fraction comme sigma=0.13343 ou gamma=23/100.")

repeat
{
print("\n Est-ce que vous voullez continuer a changer la valer defaut de parametres (y(Y)/n(N))?")
scan(what="character") -> lettrechoisie
k<- length(lettrechoisie)
lettre<-lettrechoisie[k]

if((lettre=='y')||(lettre=='Y'))
{
nouvellesparametres<- c("Argument             | Nouvelle valeur      |   Valeur défaut          |     Signification ",
        	paste(stringFixe("tmax",lAr),stringFixe(tmax,lNouVal),stringFixe(tmax0,lValDef),stringFixe("Temps pour fonctionner une simulation",lSig),sep=""),
		paste(stringFixe("initN",lAr),stringFixe(initN,lNouVal),stringFixe(initN0,lValDef),stringFixe("Nombre de population initial",lSig),sep=""),
		paste(stringFixe("initI",lAr),stringFixe(initI,lNouVal),stringFixe(initI0,lValDef),stringFixe("Nombre de personnes infectées initial",lSig),sep=""),
		paste(stringFixe("rmu",lAr),stringFixe(rmu,lNouVal),stringFixe(rmu0,lValDef),stringFixe("Taux de mortalite et aussi le taux de connaissance",lSig),sep=""),
		paste(stringFixe("beta",lAr),stringFixe(beta,lNouVal),stringFixe(beta0,lValDef),stringFixe("Taux de contact S->E initial",lSig),sep=""),
		paste(stringFixe("sigma",lAr),stringFixe(sigma,lNouVal),stringFixe(sigma0,lValDef),stringFixe("Taux de transmission E->I",lSig),sep=""),
		paste(stringFixe("gamma",lAr),stringFixe(gamma,lNouVal),stringFixe(gamma0,lValDef),stringFixe("Taux de transmission I->R",lSig),sep=""),
		paste(stringFixe("nbVilles",lAr),stringFixe(nbVilles,lNouVal),stringFixe(nbVilles0,lValDef),stringFixe("Nombre de villes",lSig),sep=""),
		paste(stringFixe("epsilon",lAr),stringFixe(epsilon,lNouVal),stringFixe(epsilon0,lValDef),stringFixe("Taux de transmission entre deux villes",lSig),sep=""),
		paste(stringFixe("topology",lAr),stringFixe(topology,lNouVal),stringFixe(topology0,lValDef),stringFixe("Type de topology entre des villes",lSig),sep=""),
		paste(stringFixe("graine",lAr),stringFixe(graine,lNouVal),stringFixe(graine0,lValDef),stringFixe("Chiffre de graine",lSig),sep=""),
		paste(stringFixe("villeDessiner",lAr),stringFixe(villeDessiner,lNouVal),stringFixe(villeDessiner0,lNouVal),stringFixe("Quelle ville que vous voulez saver ses informations",lSig),sep=""),
		paste(stringFixe("nu",lAr),stringFixe(nu,lNouVal),stringFixe(nu0,lValDef),stringFixe("Taux de l'infection a l'exterieur",lSig),sep=""),
                paste(stringFixe("tstart",lAr),stringFixe(tstart,lNouVal),stringFixe(tstart0,lValDef),stringFixe("Temps ou on commence a faire la vaccination",lSig),sep=""),
		paste(stringFixe("tfinal",lAr),stringFixe(tfinal,lNouVal),stringFixe(tfinal0,lValDef),stringFixe("Temps ou on finit la vaccination",lSig),sep=""),
		paste(stringFixe("periode",lAr),stringFixe(periode,lNouVal),stringFixe(periode0,lValDef),stringFixe("Temps entre deux fois de vaccinations",lSig),sep=""),          
		paste(stringFixe("tauxVaccToutesVilles",lAr),stringFixe(tauxVaccToutesVilles,lNouVal),stringFixe(tauxVaccToutesVilles0,lValDef),stringFixe("(%S) - Chiffre de graine du pourcentage de susceptibles vaccinees",lSig),sep=""),
		paste(stringFixe("nombreMaxVaccins",lAr),stringFixe(nombreMaxVaccins,lNouVal),stringFixe(nombreMaxVaccins0,lValDef),stringFixe("Nombre de vaccins maximal",lSig),sep=""),
		paste(stringFixe("nStrategies",lAr),stringFixe(nStrategies,lNouVal),stringFixe(nStrategies0,lValDef),stringFixe("Nombre de vaccins maximal",lSig),sep=""),
		paste(stringFixe("nPolitiques",lAr),stringFixe(nPolitiques,lNouVal),stringFixe(nPolitiques0,lValDef),stringFixe("Nombre de politiques dans une strategie",lSig),sep=""),
		paste(stringFixe("EvoStratInput",lAr),stringFixe(EvoStratInput,lNouVal),stringFixe(EvoStratInput0,lValDef),stringFixe("Fichier qui contient 5 valeurs de 5 elements du vecteur politique(w0,w1,w2,w3,w4)",lSig),sep=""),
		paste(stringFixe("EvoStratOutput",lAr),stringFixe(EvoStratOutput,lNouVal),stringFixe(EvoStratOutput0,lValDef),stringFixe("Fichier qui contient le resulat apres avoir execute le programme",lSig),sep="")
	       )
switch(menu(nouvellesparametres,graphics=FALSE,title="\n table 3: VALEUR DE PARAMETRES"),
				"Rien de fait!","tmax","initN","initI","rmu","beta","sigma","gamma",
				"nbVilles","epsilon","topology","graine","villeDessiner","nu",
				"tstart","tfinal","periode","tauxVaccToutesVilles","nombreMaxVaccins","nStrategies","nPolitiques","EvoStratInput","EvoStratOutput") -> parametrechoisie;
print(paste(parametrechoisie, " = "))
if(parametrechoisie == "tmax")
{
parametrechoisie=scan(what="integer", flush="TRUE")
n=length(parametrechoisie)
tmax = sapply(parametrechoisie[n], as.integer)
}
else if(parametrechoisie == "initN")
{
parametrechoisie=scan(what="integer", flush="TRUE")
n=length(parametrechoisie)
initN = sapply(parametrechoisie[n], as.integer)
}
else if(parametrechoisie == "initI")
{
parametrechoisie=scan(what="integer", flush="TRUE")
n=length(parametrechoisie)
initI = sapply(parametrechoisie[n], as.integer)
}
else if(parametrechoisie == "rmu")
{
parametrechoisie=scan(what="numeric", flush="TRUE")
n=length(parametrechoisie)
rmu <- ConStrToNum(parametrechoisie[n])
}
else if(parametrechoisie == "beta")
{
parametrechoisie=scan(what="numeric", flush="TRUE")
n=length(parametrechoisie)
beta <- ConStrToNum(parametrechoisie[n])
}
else if(parametrechoisie == "sigma")
{
parametrechoisie=scan(what="integer", flush="TRUE")
n=length(parametrechoisie)
sigma <- ConStrToNum(parametrechoisie[n])
}
else if(parametrechoisie == "gamma")
{
parametrechoisie=scan(what="integer", flush="TRUE")
n=length(parametrechoisie)
gamma <- ConStrToNum(parametrechoisie[n])
} 
else if(parametrechoisie == "nu")
{
parametrechoisie=scan(what="integer", flush="TRUE")
n=length(parametrechoisie)
nu <- ConStrToNum(parametrechoisie[n])
} 
else if(parametrechoisie == "nbVilles")
{
parametrechoisie=scan(what="integer", flush="TRUE")
n=length(parametrechoisie)
nbVilles = sapply(parametrechoisie[n], as.integer)
}
else if(parametrechoisie == "graine")
{
parametrechoisie=scan(what="integer", flush="TRUE")
n=length(parametrechoisie)
graine = sapply(parametrechoisie[n], as.integer)
}
else if(parametrechoisie == "topology")
{
parametrechoisie=scan(what="integer", flush="TRUE")
n=length(parametrechoisie)
topology = sapply(parametrechoisie[n], as.integer)
}
else if(parametrechoisie == "epsilon")
{
parametrechoisie=scan(what="integer", flush="TRUE")
n=length(parametrechoisie)
epsilon <- ConStrToNum(parametrechoisie[n])
}
else if(parametrechoisie == "tstart")
{
parametrechoisie=scan(what="integer", flush="TRUE")
n=length(parametrechoisie)
tstart = sapply(parametrechoisie[n], as.integer)
}
else if(parametrechoisie == "tfinal")
{
parametrechoisie=scan(what="integer", flush="TRUE")
n=length(parametrechoisie)
tfinal = sapply(parametrechoisie[n], as.integer)
}
else if(parametrechoisie == "periode")
{
parametrechoisie=scan(what="integer", flush="TRUE")
n=length(parametrechoisie)
periode = sapply(parametrechoisie[n], as.integer)
}
else if(parametrechoisie == "tauxVaccToutesVilles")
{
parametrechoisie=scan(what="integer", flush="TRUE")
n=length(parametrechoisie)
tauxVaccToutesVilles = sapply(parametrechoisie[n], as.integer)
}
else if(parametrechoisie == "nombreMaxVaccins")
{
parametrechoisie=scan(what="integer", flush="TRUE")
n=length(parametrechoisie)
nombreMaxVaccins = sapply(parametrechoisie[n], as.integer)
}
else if(parametrechoisie == "nStrategies")
{
parametrechoisie=scan(what="integer", flush="TRUE")
n=length(parametrechoisie)
nStrategies = sapply(parametrechoisie[n], as.integer)
}
else if(parametrechoisie == "nPolitiques")
{
parametrechoisie=scan(what="integer", flush="TRUE")
n=length(parametrechoisie)
nPolitiques = sapply(parametrechoisie[n], as.integer)
}
else if(parametrechoisie == "EvoStratInput")
{
parametrechoisie=scan(what="integer", flush="TRUE")
n=length(parametrechoisie)
EvoStratInput = parametrechoisie[n]
}
else if(parametrechoisie == "EvoStratOutput")
{
parametrechoisie=scan(what="integer", flush="TRUE")
n=length(parametrechoisie)
EvoStratOutput = parametrechoisie[n]
}
else
{
print("il faut revoir votre parametre!")
}
#The end of if loop
}
else
{
menu(nouvellesparametres,graphics=FALSE,title="\n table 3: VALEUR DE PARAMETRES");
break
}
#The end repeat
}
#The end of M
}
else
{
print("Il n'y a pas de lettre comme votre choix!")
print("Le programme marche par defaut.....")
}
	
#Recuperer les valeurs des parametres
nbVilles <- formatC(nbVilles,mode="integer")
text1<-paste(" -nbVilles ",nbVilles)
text2<-paste(" -tmax ",tmax)
xI <- formatC(initI,mode="integer")
text3<-paste(" -initI ",xI)
xN <- formatC(initN,mode="integer")
text4<-paste(" -initN ",xN)
text5<-paste(" -beta ",beta)
text6<-paste(" -epsilon ",epsilon)
text7<-paste(" -rmu ",rmu)
text8<-paste(" -sigma ",sigma)
text9<-paste(" -gamma ",gamma)
text10<-paste(" -topology ",topology)
text11<-paste(" -tstart ",tstart)
text12<-paste(" -tfinal ",tfinal)
text13<-paste(" -periode ", periode)
#text15<-paste(" -MaxCycle ", MaxCycle) 
text16<-paste(" -graine ",format(graine,mode="integer"))
text17<-paste(" -nu ", nu)
text18<-paste(" -tauxVaccToutesVilles ", tauxVaccToutesVilles) 
text19<-paste(" -nombreMaxVaccins ", nombreMaxVaccins) 
text20<-paste(" -nStrategies ", nStrategies) 
text21<-paste(" -nPolitiques ", nPolitiques) 
text22<-paste(" -EvoStratInput ", EvoStratInput) 
text23<-paste(" -EvoStratOutput ", EvoStratOutput) 


#Fonctionner le programme avec vaccination en utilisant (1+1)-ES sous C++
#Creer une ligne de texte des parametres
textTotal=paste(text1, text2, text3,text4,text5,text6,text7,text8,text9,text10, text11, text12, text13, text16,text17,text18,text19, text20, text21, text22,text23)

#Creer un chemin pour fonctionner le programme 
EvoStrat<-paste("./EvoStrat ",textTotal)
system(EvoStrat)



# The end
}


