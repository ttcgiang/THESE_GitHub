############################## Interface de VACCINATION utilisant l'algorithme (1+1)-ES  ###################################################################
# Le but de cette fonction est de fonctionner le programmes "../EvolStrat.cpp"                                                                            ##
# qui evalue les politiques de vaccinations en utilisant l'algorithme (1+1)-ES                                                                            ##
# Par      TRAN Thi Cam Giang - P15 - IFI - Hanoi - Vietname                                                                                              ##
############################################################################################################################################################
Simulator = "./BatchVacc"
source("stringFixe.r")
source("ConStrToNum.r")
strPDFLocal = paste("ResulBatchVaccLocal",Sys.Date(),".pdf",seq="")
strPDFGeneral = paste("ResulBatchVaccGeneral",Sys.Date(),".pdf",seq="")
fBatchParms <- paste("./BatchParms",Sys.Date(),".txt",seq="")
couleurs<-c("red","blue","green","yellow","violet","black")
nCol = length(couleurs)


InterBatchVacc <- function(
			tmax = 2000,initN=10000,initI=100,nu=0.01,rmu=1/(10*365),beta=1000/365,sigma2 = 0, sigma=1/7,
			gamma=1/7, epsilon =0.10, nbVilles = 1,topology=1,graine=10, villeDessiner=1,
			...
			) {
#Creer un fichier en PDF pour sauvegarder l'image de resultats
pdf(strPDFLocal)

#Sauvegarder les valeurs par defaut des arguments
tmax0                = tmax
initN0               = initN
initI0               = initI
nu0                  = nu
rmu0                 = rmu
beta0                = beta
sigma0               = sigma
gamma0               = gamma
epsilon0             = epsilon
nbVilles0            = nbVilles
topology0            = topology
graine0              = graine
villeDessiner0       = villeDessiner

#Creer les valeurs des largeurs pour afficher le menu plus beau
lAr                  = 23
lNouVal              = 25
lValDef              = 25
lSig                 = 50

#Executer BatchVacc.cpp jusqu'à ne pas vouloir fonctionner
repeat
{
print("\n Voull-vous executer le programme BatchVacc.cpp(y(Y)/n(N))?")
scan(what="character") -> lettrechoisie
k<- length(lettrechoisie)
lettre<-lettrechoisie[k]

if((lettre=='y')||(lettre=='Y'))
{
#Creer le munu des arguments avec les valeurs par defaut
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
		paste(stringFixe("nu",lAr),stringFixe(nu0,lValDef),stringFixe("Taux de l'infection a l'exterieur",lSig),sep="")
                
	       )

menu(parametresDeFaut,graphics=FALSE,title="\n table 1: VALEURS DES PARAMETRES PAR ANNEE, par defaut\n Pouvez choisir n'importe quel chiffre");
print("*********************************************************************************")	
switch(menu(c("Selectionner les paramatres par defaut", "Modifier parametres"),graphics=FALSE,title="\n table 2: DONNER votre choix!") + 1,  cat("Rien de fait\n"), "S", "M")->choix
if(choix == "S")
{
print("Vous avez choisi: Selectionner les paramatres par defaut!")
print("Le programme est en train de marcher....")
}
else if (choix == "M")
{
print("Tapez quelle parametre que vous modifiez sa valeur!")
print("Note : - Parametres de Simulation SEIR: tmax, initN, initI, rmu, beta0, sigma2, sigma, gamma, nbVilles, epsilon, topology, villeDessiner")
print("       - Type de parametres est un chiffre ou une fraction comme sigma=0.13343 ou gamma=23/100.")

repeat
{
print("Voullez-vous changer la valer defaut de parametres (y(Y)/n(N))?")
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
		paste(stringFixe("nu",lAr),stringFixe(nu,lNouVal),stringFixe(nu0,lValDef),stringFixe("Taux de l'infection a l'exterieur",lSig),sep="")
	       )
switch(menu(nouvellesparametres,graphics=FALSE,title="\n table 3: VALEUR DE PARAMETRES"),
				"Rien de fait!","tmax","initN","initI","rmu","beta","sigma","gamma",
				"nbVilles","epsilon","topology","graine","villeDessiner","nu"
				) -> parametrechoisie;
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

else
{
print("il faut revoir votre parametre!")
}
#The end of if loop
}
else
{
#menu(nouvellesparametres,graphics=FALSE,title="\n table 3: VALEUR DE PARAMETRES");
print("Executant le programme!")
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
text1<-paste(" -tmax ",tmax)
xI <- formatC(initI,mode="integer")
text2<-paste(" -initI ",xI)
xN <- formatC(initN,mode="integer")
text3<-paste(" -initN ",xN)
text4<-paste(" -beta ",beta)
text5<-paste(" -rmu ",rmu)
text6<-paste(" -sigma ",sigma)
text7<-paste(" -gamma ",gamma)
text8<-paste(" -nu ", nu)

nbVilles <- formatC(nbVilles,mode="integer")
text9<-paste(" -nbVilles ",nbVilles)
text10<-paste(" -epsilon ",epsilon)
text11<-paste(" -topology ",topology)
text12<-paste(" -graine ",format(graine,mode="integer"))



#Fonctionner le programme avec vaccination en utilisant (1+1)-ES sous C++
#Creer une ligne de texte des parametres
BatchParm=paste(text1, text2, text3,text4,text5,text6,text7,text8,text9,text10, text11, text12)
StartBatchProcess<-paste(Simulator,BatchParm)
system(StartBatchProcess)
#Afficher les resultats sur une image qui est sauvegardee dans le fichier ResulBatchVaccLocal en pdf
nameFile<-paste("SimVil",villeDessiner,"Eps",sprintf("%.3f", epsilon),"Top",topology,".csv",sep="")
#Lire le resultat
dataSansVac<- read.table(nameFile)
plot(dataSansVac[,1],dataSansVac[,4],type="n",xlab="Temps (jours)",ylab="nombre d'Infectes")
maxSansVac<-max(dataSansVac[,4])
#Dessiner une courbe qui a le nom "textOut"
textOut <- paste(" Simulation sans Vaccination","\n","nb Villes:",nbVilles,",beta:",sprintf("%.3f",beta),",gamma:",sprintf("%.3f",gamma),",sigma:",sprintf("%.3f",sigma),"\nI0:",initI,",N0:",initN,",epsilon",sprintf("%.3f",epsilon),",rmu:",sprintf("%.5f",rmu),",nu:",sprintf("%.3f",nu),",tmax:",tmax,",topology:",topology,sep="  ")
title(textOut,cex.main=0.6)
for(ville in 1:nbVilles)
{
#Creer le nom de fichier qui contient le resultat
nameFileSansVac<-paste("SimVil",ville,"Eps",sprintf("%.3f", epsilon),"Top",topology,".csv",sep="")
#Lire le resultat
dataSansVac<- read.table(nameFileSansVac)
x<-dataSansVac[,1]
y<-dataSansVac[,4]
lines(x,y,type="l",col=couleurs[ville%%nCol])
ymin<- ville*15;
text(tmax-0.5*tmax,maxSansVac-ymin,paste("ville ",ville,"(",couleurs[ville],")"," : somme de I = ", sum(dataSansVac[,4])),cex=0.6)
}
dev.off()

#Sauvegarder les valeurs des arguments dans fichier.txt
print("Est-ce que vous voulez sauvegarder les valeurs des arguments que vous vien d'utiliser pour simulateur, au fichier general BatchParms.txt (y(Y)/n(N)?")
lettres <-scan(what="complex")
nlettres <- length(lettres)
lettre <- lettres[nlettres]
if((lettre =='y')||(lettre=='Y'))
{
write(BatchParm, file = fBatchParms,append=TRUE, sep = "\n")
}
}
else 
{
print("avoir sauvegarde !")
break
}
}

#pdf general
print("Est-ce que vous voulez executer toutes les valeurs des arguments des simulateurs que vous avez sauvegardes dans le fichier general BatchParms.txt (y(Y)/n(N)?")
lettres <-scan(what="complex")
nlettres <- length(lettres)
lettre <- lettres[nlettres]
if((lettre =='y')||(lettre=='Y'))
{
pdf(strPDFGeneral)
arrBatchParms <- readLines(fBatchParms, n=-1)
nBatchParms <- length(arrBatchParms)
for( iBat in 1:nBatchParms)
{
StartBatchProcess<-paste(Simulator,arrBatchParms[iBat])
system(StartBatchProcess)
#Pour nbVilles villes
#Fonctionner le programme sans vaccination 
nameFile<-paste("SimVil",villeDessiner,"Eps",sprintf("%.3f", epsilon),"Top",topology,".csv",sep="")

#Lire le resultat
dataSansVac<- read.table(nameFile)
plot(dataSansVac[,1],dataSansVac[,4],type="n",xlab="Temps (jours)",ylab="nombre d'Infectes")
maxSansVac<-max(dataSansVac[,4])

#Dessiner une courbe qui a le nom "textOut"
textmainTP <- unlist(strsplit(arrBatchParms[iBat],split="-nu"))
textmain <- paste(textmainTP[1],"\n","-nu",textmainTP[2])

textOut <- paste(" Simulation sans Vaccination","\n",textmain)
title(textOut,cex.main=0.6)

#Prendre la valeur de nbVilles
arrCharacters<- scan(fBatchParms,character(0), sep="")
#arrCharacters <- strsplit(arrBatchParms[iBat], split="\\  ", fixed = FALSE, perl = FALSE, useBytes = FALSE) 
arrIndexNbVils <- which(arrCharacters=="-nbVilles")
nbVilles <- sapply(arrCharacters[sapply(arrIndexNbVils[iBat],as.integer) + 1], as.integer)
print(nbVilles)

for(ville in 1:nbVilles)
{
#Creer le nom de fichier qui contient le resultat
nameFileSansVac<-paste("SimVil",ville,"Eps",sprintf("%.3f", epsilon),"Top",topology,".csv",sep="")
#Lire le resultat
dataSansVac<- read.table(nameFileSansVac)
x<-dataSansVac[,1]
y<-dataSansVac[,4]
lines(x,y,type="l",col=couleurs[ville%%nCol])
ymin<- ville*15;
text(tmax-0.5*tmax,maxSansVac-ymin,paste("ville ",ville,"(",couleurs[ville],")"," : somme de I = ", sum(dataSansVac[,4])),cex=0.6)
}
}
dev.off()
#The end if general
}


print("Finir le programme InetrBatchVacc.r!")
#print(nameFile)
nameFile

# The end
}


