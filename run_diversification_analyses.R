# If using this analytical pipeline, please cite:
# Condamine F.L., Antonelli A., Lagomarsino L.P., Hoorn C., Liow L.H. 2017. Teasing apart mountain uplift, climate change and biotic drivers of diversification. In: Mountains, Climate, and Biodiversity (eds. Hoorn C., Antonelli A.). Wiley Blackwell, pp. xx-xx.
# Author: Fabien L. Condamine (email: fabien.condamine@gmail.com). Please don't hesitate to contact me if you have any issue or question with this code.

#Required R packages (will install them automatically if not yet installed)
if ( ! ("DDD" %in% installed.packages())) {install.packages("DDD", dependencies=T)}
if ( ! ("picante" %in% installed.packages())) {install.packages("picante", dependencies=T)}
if ( ! ("pspline" %in% installed.packages())) {install.packages("pspline", dependencies=T)}
if ( ! ("TreePar" %in% installed.packages())) {install.packages("TreePar", dependencies=T)}
library("DDD")
library("picante")
library("pspline")
library("TreePar")

#R codes and functions
source("diversification_library/fit_bd.R")
source("diversification_library/fit_env_bd.R")
source("diversification_library/integrate.R")
source("diversification_library/likelihood_bd.R")
source("diversification_library/Phi.R")
source("diversification_library/Psi.R")
source("diversification_library/tables.summary.R")

no.extension <- function(filename)
{
	if (substr(filename, nchar(filename), nchar(filename))==".") {return(substr(filename, 1, nchar(filename)-1))} 
	else {no.extension(substr(filename, 1, nchar(filename)-1))} 
}

################################################################################
# LIST of FUNCTIONS
# 1. run_Morlon_models (tree_file, sampling_fraction=1, number_of_trees=1)
# 2. run_PaleoEnv (tree_file, env_data_file,sampling_fraction=1, number_of_trees=1)
# 3. run_DDD(tree_file, total_richness=Ntip(tree_file), number_of_trees=1)
# 4. run_TreePar (tree_file, sampling_fraction=1, number_of_trees=1)
################################################################################


## Morlon's models: time-dependent diversification models with continuous rates through time
# If using this appraoch, please cite:
# Morlon H., Parsons T.L., Plotkin J. 2011. Reconciling molecular phylogenies with the fossil record. Proc. Natl. Acad. Sci. USA 108:16327–16332.

run_Morlon_models <- function (tree_file, sampling_fraction=1, number_of_trees=1)
{
	tree_file_name <- tree_file
	tree_file<-read.nexus(tree_file)
	posteriors<-sample(tree_file, number_of_trees)

	finaltree_file<-list()

	for (i in 1:length(posteriors))
{
	print(i)

	if (length(posteriors)==1){phyloi<-tree_file} else {phyloi<-posteriors[[i]]}

	tot_time<-max(node.age(phyloi)$ages)
	f<-sampling_fraction
	cond="crown"
		

# BCST (Pure birth)
print("BCST")
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){0}
lamb_par<-c(0.1)
mu_par<-c()
cst.lamb=T; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T

	treei_BCST<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BCST)


# BCST DCST (constant Birth-death)
print("BCST DCST")
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){y[1]}
lamb_par<-c(treei_BCST$lamb_par[1])
mu_par<-c(0.01)
cst.lamb=T; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

	treei_BCSTDCST<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BCSTDCST)

###############################################
###### Time Dependence (exponential variation) ######
###############################################

	print(i)
	
# BTimeVar EXPO
print("BTimeVar EXPO")
f.lamb<-function(x,y){y[1]*exp(y[2]*x)}
f.mu<-function(x,y){0}
lamb_par<-c(treei_BCSTDCST$lamb_par[1],0.01)
mu_par<-c()
cst.lamb=F; cst.mu=T; expo.lamb=T; expo.mu=F; fix.mu=T

	treei_BTimeVar_EXPO<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTimeVar_EXPO)


# BTimeVar DCST EXPO
print("BTimeVar DCST EXPO")
f.lamb<-function(x,y){y[1]*exp(y[2]*x)}
f.mu<-function(x,y){y[1]}
lamb_par<-c(treei_BTimeVar_EXPO$lamb_par[1],treei_BTimeVar_EXPO$lamb_par[2])
#mu_par<-c(treei_BCSTDCST$mu_par[1])
mu_par<-c(0.01)
cst.lamb=F; cst.mu=T; expo.lamb=T; expo.mu=F; fix.mu=F

	treei_BTimeVarDCST_EXPO<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTimeVarDCST_EXPO)


# BCST DTimeVar EXPO
print("BCST DTimeVar EXPO")
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(treei_BCSTDCST$lamb_par[1])
mu_par<-c(0.01,0.001)
cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=T; fix.mu=F

	treei_BCSTDTimeVar_EXPO<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BCSTDTimeVar_EXPO)


# BTimeVar DTimeVar EXPO
print("BTimeVar DTimeVar EXPO")
f.lamb<-function(x,y){y[1]*exp(y[2]*x)}
f.mu<-function(x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(treei_BTimeVarDCST_EXPO$lamb_par[1],0.001)
mu_par<-c(0.05,0.001)
cst.lamb=F; cst.mu=F; expo.lamb=T; expo.mu=T; fix.mu=F

	treei_BTimeVarDTimeVar_EXPO<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTimeVarDTimeVar_EXPO)


##########################################
###### Time Dependence (linear variation) ######
##########################################

	print(i)
	
# BTimeVar LIN
print("BTimeVar LIN")
f.lamb<-function(x,y){y[1]+(y[2]*x)}
f.mu<-function(x,y){0}
lamb_par<-c(0.01,0.001)
mu_par<-c()
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T

	treei_BTimeVar_LIN<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTimeVar_LIN)


# BTimeVar DCST LIN
print("BTimeVar DCST LIN")
f.lamb<-function(x,y){y[1]+(y[2]*x)}
f.mu<-function(x,y){y[1]}
lamb_par<-c(abs(treei_BTimeVar_LIN$lamb_par[1]),treei_BTimeVar_LIN$lamb_par[2])
mu_par<-c(0.01)
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

	treei_BTimeVarDCST_LIN<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTimeVarDCST_LIN)


# BCST DTimeVar LIN
print("BCST DTimeVar LIN")
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){y[1]+(y[2]*x)}
lamb_par<-c(treei_BCSTDCST$lamb_par[1])
mu_par<-c(0.001,0.001)
cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

	treei_BCSTDTimeVar_LIN<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BCSTDTimeVar_LIN)


# BTimeVar DTimeVar LIN
print("BTimeVar DTimeVar LIN")
f.lamb<-function(x,y){y[1]+(y[2]*x)}
f.mu<-function(x,y){y[1]+(y[2]*x)}
lamb_par<-c(abs(treei_BTimeVarDCST_LIN$lamb_par[1]),0.001)
mu_par<-c(0.05,-0.001)
cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

	treei_BTimeVarDTimeVar_LIN<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTimeVarDTimeVar_LIN)

	############# RESULTS ###########################################
	
	results<-matrix(NA,10,8)
	colnames(results)<-c("Models","NP","logL","AICc","Lambda","AlphaTime","Mu","BetaTime")

#Models
	results[,1]<-c("BCST","BCSTDCST","BTimeVar_EXPO","BTimeVarDCST_EXPO","BCSTDTimeVar_EXPO","BTimeVarDTimeVar_EXPO","BTimeVar_LIN","BTimeVarDCST_LIN","BCSTDTimeVar_LIN","BTimeVarDTimeVar_LIN")

#NP
	results[1,2]<-1
	results[2,2]<-2
	results[3,2]<-2
	results[4,2]<-3
	results[5,2]<-3
	results[6,2]<-4
	results[7,2]<-2
	results[8,2]<-3
	results[9,2]<-3
	results[10,2]<-4

#logL
	results[1,3]<-round(treei_BCST$LH,3)
	results[2,3]<-round(treei_BCSTDCST$LH,3)
	results[3,3]<-round(treei_BTimeVar_EXPO$LH,3)
	results[4,3]<-round(treei_BTimeVarDCST_EXPO$LH,3)
	results[5,3]<-round(treei_BCSTDTimeVar_EXPO$LH,3)
	results[6,3]<-round(treei_BTimeVarDTimeVar_EXPO$LH,3)
	results[7,3]<-round(treei_BTimeVar_LIN$LH,3)
	results[8,3]<-round(treei_BTimeVarDCST_LIN$LH,3)
	results[9,3]<-round(treei_BCSTDTimeVar_LIN$LH,3)
	results[10,3]<-round(treei_BTimeVarDTimeVar_LIN$LH,3)

#AICc
	results[1,4]<-round(treei_BCST$aicc,3)
	results[2,4]<-round(treei_BCSTDCST$aicc,3)
	results[3,4]<-round(treei_BTimeVar_EXPO$aicc,3)
	results[4,4]<-round(treei_BTimeVarDCST_EXPO$aicc,3)
	results[5,4]<-round(treei_BCSTDTimeVar_EXPO$aicc,3)
	results[6,4]<-round(treei_BTimeVarDTimeVar_EXPO$aicc,3)
	results[7,4]<-round(treei_BTimeVar_LIN$aicc,3)
	results[8,4]<-round(treei_BTimeVarDCST_LIN$aicc,3)
	results[9,4]<-round(treei_BCSTDTimeVar_LIN$aicc,3)
	results[10,4]<-round(treei_BTimeVarDTimeVar_LIN$aicc,3)
	
#Lambda0
	results[1,5]<-round(abs(treei_BCST$lamb_par[1]),4)
	results[2,5]<-round(abs(treei_BCSTDCST$lamb_par[1]),4)
	results[3,5]<-round(abs(treei_BTimeVar_EXPO$lamb_par[1]),4)
	results[4,5]<-round(abs(treei_BTimeVarDCST_EXPO$lamb_par[1]),4)
	results[5,5]<-round(abs(treei_BCSTDTimeVar_EXPO$lamb_par[1]),4)
	results[6,5]<-round(abs(treei_BTimeVarDTimeVar_EXPO$lamb_par[1]),4)
	results[7,5]<-round(abs(treei_BTimeVar_LIN$lamb_par[1]),4)
	results[8,5]<-round(abs(treei_BTimeVarDCST_LIN$lamb_par[1]),4)
	results[9,5]<-round(abs(treei_BCSTDTimeVar_LIN$lamb_par[1]),4)
	results[10,5]<-round(abs(treei_BTimeVarDTimeVar_LIN$lamb_par[1]),4)

#Alpha Time
	results[3,6]<-round(treei_BTimeVar_EXPO$lamb_par[2],5)
	results[4,6]<-round(treei_BTimeVarDCST_EXPO$lamb_par[2],5)
	results[6,6]<-round(treei_BTimeVarDTimeVar_EXPO$lamb_par[2],5)
	results[7,6]<-round(treei_BTimeVar_LIN$lamb_par[2],5)
	results[8,6]<-round(treei_BTimeVarDCST_LIN$lamb_par[2],5)
	results[10,6]<-round(treei_BTimeVarDTimeVar_LIN$lamb_par[2],5)

#Mu0
	results[2,7]<-round(abs(treei_BCSTDCST$mu_par[1]),5)
	results[4,7]<-round(abs(treei_BTimeVarDCST_EXPO$mu_par[1]),5)
	results[5,7]<-round(abs(treei_BCSTDTimeVar_EXPO$mu_par[1]),5)
	results[6,7]<-round(abs(treei_BTimeVarDTimeVar_EXPO$mu_par[1]),5)
	results[8,7]<-round(abs(treei_BTimeVarDCST_LIN$mu_par[1]),5)
	results[9,7]<-round(abs(treei_BCSTDTimeVar_LIN$mu_par[1]),5)
	results[10,7]<-round(abs(treei_BTimeVarDTimeVar_LIN$mu_par[1]),5)

#Beta Time
	results[5,8]<-round(treei_BCSTDTimeVar_EXPO$mu_par[2],5)
	results[6,8]<-round(treei_BTimeVarDTimeVar_EXPO$mu_par[2],5)
	results[9,8]<-round(treei_BCSTDTimeVar_LIN$mu_par[2],5)
	results[10,8]<-round(treei_BTimeVarDTimeVar_LIN$mu_par[2],5)
	
	finaltree_file[[i]]<-results
	#print(results)
}

	final_table_tree_file<-tables.summary(finaltree_file)

	fname <- no.extension(basename(tree_file_name))
	outfile <- paste(dirname(tree_file_name), "/", fname, "_results_Morlon.txt", sep="")
	out_R <- paste(dirname(tree_file_name), "/", fname, "_complete_results_Morlon.Rdata", sep="")

	write.table(final_table_tree_file,file=outfile,quote=FALSE,sep="\t",row.names=FALSE)
	save(final_table_tree_file,file=out_R)
	
	print(final_table_tree_file)
	print("Analyses with the time-dependent (Morlon) models completed ! Check your results in the working directory.")
}



## Condamine's models: paleoenvironment-dependent diversification models with continuous rates through time
# If using this appraoch, please cite:
# Condamine F.L., Rolland J., Morlon H. 2013. Macroevolutionary perspectives to environmental change. Ecol. Lett. 16:72–85.

run_PaleoEnv <- function (tree_file, env_data_file, sampling_fraction=1, number_of_trees=1)
{
	env_data<-read.table(env_data_file,header=T)
	tree_file_name <- tree_file
	tree_file<-read.nexus(tree_file)
	posteriors<-sample(tree_file, number_of_trees)

	finaltree_file<-list()

	for (i in 1:length(posteriors))
	{
	print(i)

	if (length(posteriors)==1){phyloi<-tree_file} else {phyloi<-posteriors[[i]]}
		
	tot_time<-max(node.age(phyloi)$ages)
	f<-sampling_fraction
	cond="crown"
	
	# BCST DCST (constant Birth-death)
	print("BCST DCST")
	f.lamb<-function(x,y){y[1]}
	f.mu<-function(x,y){y[1]}
	lamb_par<-c(0.1)
	mu_par<-c(0.01)
	cst.lamb=T; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

		treei_BCSTDCST<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
		print(treei_BCSTDCST)

	
	###############################################
	###### Temp Dependence (exponential variation) ######
	###############################################

		print(i)
	
	# BEnv.Var EXPO
	print("BEnv.Var EXPO")
	f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
	f.mu<-function(t,x,y){0}
	lamb_par<-c(0.1, 0.0)
	mu_par<-c()
	cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T

		treei_BEnv.Var_EXPO<-fit_env_bd(phyloi,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
		print(treei_BEnv.Var_EXPO)


	# BEnv.Var DCST EXPO
	print("BEnv.Var DCST EXPO")
	f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
	f.mu<-function(t,x,y){y[1]}
	lamb_par<-c(abs(treei_BEnv.Var_EXPO$lamb_par[1]), 0.0)
	mu_par<-c(0.01)
	cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

		treei_BEnv.VarDCST_EXPO<-fit_env_bd(phyloi,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
		print(treei_BEnv.VarDCST_EXPO)


	# BCST DEnv.Var EXPO
	print("BCST DEnv.Var EXPO")
	f.lamb<-function(t,x,y){y[1]}
	f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
	lamb_par<-c(treei_BCSTDCST$lamb_par[1])
	mu_par<-c(0.01,0.001)
	cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

		treei_BCSTDEnv.Var_EXPO<-fit_env_bd(phyloi,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
		print(treei_BCSTDEnv.Var_EXPO)


	# BEnv.Var DEnv.Var EXPO
	print("BEnv.Var DEnv.Var EXPO")
	f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
	f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
	lamb_par<-c(abs(treei_BEnv.VarDCST_EXPO$lamb_par[1]), 0.0)
	mu_par<-c(0.01, 0.0)
	cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

		treei_BEnv.VarDEnv.Var_EXPO<-fit_env_bd(phyloi,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
		print(treei_BEnv.VarDEnv.Var_EXPO)


	##########################################
	###### Temp Dependence (linear variation) ######
	##########################################

		print(i)
	
	# BEnv.Var LIN
	print("BEnv.Var LIN")
	f.lamb<-function(t,x,y){y[1]+y[2]*x}
	f.mu<-function(t,x,y){0}
	lamb_par<-c(abs(treei_BEnv.Var_EXPO$lamb_par[1]), 0.0)
	mu_par<-c()
	cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T

		treei_BEnv.Var_LIN<-fit_env_bd(phyloi,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
		print(treei_BEnv.Var_LIN)


	# BEnv.Var DCST LIN
	print("BEnv.Var DCST LIN")
	f.lamb<-function(t,x,y){y[1]+y[2]*x}
	f.mu<-function(t,x,y){y[1]}
	lamb_par<-c(abs(treei_BEnv.Var_LIN$lamb_par[1]), 0.0)
	mu_par<-c(0.01)
	cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

		treei_BEnv.VarDCST_LIN<-fit_env_bd(phyloi,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
		print(treei_BEnv.VarDCST_LIN)


	# BCST DEnv.Var LIN
	print("BCST DEnv.Var LIN")
	f.lamb<-function(t,x,y){y[1]}
	f.mu<-function(t,x,y){y[1]+y[2]*x}
	lamb_par<-c(treei_BCSTDCST$lamb_par[1])
	mu_par<-c(0.02, 0.0)
	cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

		treei_BCSTDEnv.Var_LIN<-fit_env_bd(phyloi,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
		print(treei_BCSTDEnv.Var_LIN)


	# BEnv.Var DEnv.Var LIN
	print("BEnv.Var DEnv.Var LIN")
	f.lamb<-function(t,x,y){y[1]+y[2]*x}
	f.mu<-function(t,x,y){y[1]+y[2]*x}
	lamb_par<-c(abs(treei_BEnv.VarDCST_LIN$lamb_par[1]), 0.0)
	mu_par<-c(0.02, 0.0)
	cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

		treei_BEnv.VarDEnv.Var_LIN<-fit_env_bd(phyloi,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
		print(treei_BEnv.VarDEnv.Var_LIN)


	############# RESULTS ###########################################

		results<-matrix(NA,8,8)
		colnames(results)<-c("Models","NP","logL","AICc","Lambda","AlphaEnv.","Mu","BetaEnv.")

	#Models
		results[,1]<-c("BEnv.Var_EXPO","BEnv.VarDCST_EXPO","BCSTDEnv.Var_EXPO","BEnv.VarDEnv.Var_EXPO",
		"BEnv.Var_LIN","BEnv.VarDCST_LIN","BCSTDEnv.Var_LIN","BEnv.VarDEnv.Var_LIN")

	#NP
		results[1,2]<-2
		results[2,2]<-3
		results[3,2]<-3
		results[4,2]<-4
		results[5,2]<-2
		results[6,2]<-3
		results[7,2]<-3
		results[8,2]<-4

	#logL
		results[1,3]<-treei_BEnv.Var_EXPO$LH
		results[2,3]<-treei_BEnv.VarDCST_EXPO$LH
		results[3,3]<-treei_BCSTDEnv.Var_EXPO$LH
		results[4,3]<-treei_BEnv.VarDEnv.Var_EXPO$LH
		results[5,3]<-treei_BEnv.Var_LIN$LH
		results[6,3]<-treei_BEnv.VarDCST_LIN$LH
		results[7,3]<-treei_BCSTDEnv.Var_LIN$LH
		results[8,3]<-treei_BEnv.VarDEnv.Var_LIN$LH

	#AICc
		results[1,4]<-treei_BEnv.Var_EXPO$aicc
		results[2,4]<-treei_BEnv.VarDCST_EXPO$aicc
		results[3,4]<-treei_BCSTDEnv.Var_EXPO$aicc
		results[4,4]<-treei_BEnv.VarDEnv.Var_EXPO$aicc
		results[5,4]<-treei_BEnv.Var_LIN$aicc
		results[6,4]<-treei_BEnv.VarDCST_LIN$aicc
		results[7,4]<-treei_BCSTDEnv.Var_LIN$aicc
		results[8,4]<-treei_BEnv.VarDEnv.Var_LIN$aicc

	#Lambda0
		results[1,5]<-abs(treei_BEnv.Var_EXPO$lamb_par[1])
		results[2,5]<-abs(treei_BEnv.VarDCST_EXPO$lamb_par[1])
		results[3,5]<-abs(treei_BCSTDEnv.Var_EXPO$lamb_par[1])
		results[4,5]<-abs(treei_BEnv.VarDEnv.Var_EXPO$lamb_par[1])
		results[5,5]<-abs(treei_BEnv.Var_LIN$lamb_par[1])
		results[6,5]<-abs(treei_BEnv.VarDCST_LIN$lamb_par[1])
		results[7,5]<-abs(treei_BCSTDEnv.Var_LIN$lamb_par[1])
		results[8,5]<-abs(treei_BEnv.VarDEnv.Var_LIN$lamb_par[1])

	#Alpha Env.
		results[1,6]<-treei_BEnv.Var_EXPO$lamb_par[2]
		results[2,6]<-treei_BEnv.VarDCST_EXPO$lamb_par[2]
		results[4,6]<-treei_BEnv.VarDEnv.Var_EXPO$lamb_par[2]
		results[5,6]<-treei_BEnv.Var_LIN$lamb_par[2]
		results[6,6]<-treei_BEnv.VarDCST_LIN$lamb_par[2]
		results[8,6]<-treei_BEnv.VarDEnv.Var_LIN$lamb_par[2]	

	#Mu0
		results[2,7]<-abs(treei_BEnv.VarDCST_EXPO$mu_par[1])
		results[3,7]<-abs(treei_BCSTDEnv.Var_EXPO$mu_par[1])
		results[4,7]<-abs(treei_BEnv.VarDEnv.Var_EXPO$mu_par[1])
		results[6,7]<-abs(treei_BEnv.VarDCST_LIN$mu_par[1])
		results[7,7]<-abs(treei_BCSTDEnv.Var_LIN$mu_par[1])
		results[8,7]<-abs(treei_BEnv.VarDEnv.Var_LIN$mu_par[1])

	#Beta Env.
		results[3,8]<-treei_BCSTDEnv.Var_EXPO$mu_par[2]
		results[4,8]<-treei_BEnv.VarDEnv.Var_EXPO$mu_par[2]
		results[7,8]<-treei_BCSTDEnv.Var_LIN$mu_par[2]
		results[8,8]<-treei_BEnv.VarDEnv.Var_LIN$mu_par[2]

		finaltree_file[[i]]<-results
		#print(results)
}

	final_table_tree_file<-tables.summary(finaltree_file)

	fname <- no.extension(basename(tree_file_name))
	if (env_data_file=="./PaleoEnv/PastTemperature.txt") 
	{
	outfile <- paste(dirname(tree_file_name), "/", fname, "_results_PastTemperature.txt", sep="")
	out_R <- paste(dirname(tree_file_name), "/", fname, "_complete_results_PastTemperature.Rdata", sep="")
	write.table(final_table_tree_file,file=outfile,quote=FALSE,sep="\t",row.names=FALSE)
	save(final_table_tree_file,file=out_R)
	}
	if (env_data_file=="./PaleoEnv/PastSeaLevel.txt") 
	{
	outfile <- paste(dirname(tree_file_name), "/", fname, "_results_PastSeaLevel.txt", sep="")
	out_R <- paste(dirname(tree_file_name), "/", fname, "_complete_results_PastSeaLevel.Rdata", sep="")
	write.table(final_table_tree_file,file=outfile,quote=FALSE,sep="\t",row.names=FALSE)
	save(final_table_tree_file,file=out_R)
	}
	if (env_data_file=="./PaleoEnv/PastAndeanAltitude.txt") 
	{
	outfile <- paste(dirname(tree_file_name), "/", fname, "_results_PastAndeanAltitude.txt", sep="")
	out_R <- paste(dirname(tree_file_name), "/", fname, "_complete_results_PastAndeanAltitude.Rdata", sep="")
	write.table(final_table_tree_file,file=outfile,quote=FALSE,sep="\t",row.names=FALSE)
	save(final_table_tree_file,file=out_R)
	}

	print(final_table_tree_file)
	print("Analyses with the paleoenvironment-dependent models completed ! Check your results in the working directory.")
}



## DDD models: diversity-dependent diversification models
# If using this appraoch, please cite:
# Etienne R.S., Haegeman B., Stadler T., Aze T., Pearson P.N., Purvis A., Phillimore A.B. 2012. Diversity-dependence brings molecular phylogenies closer to agreement with the fossil record. Proc. Roy. Soc. Lond. B 279:1300–1309.

run_DDD<-function(tree_file, total_richness=Ntip(tree_file), number_of_trees=1)
{
	tree_file_name<-tree_file
	tree_file<-read.nexus(tree_file)
	posteriors<-sample(tree_file, number_of_trees)

	final<-list()

	for (i in 1:number_of_trees)
	{
		print(i)
		
		if (length(posteriors)==1){phyloi<-tree_file} else {phyloi<-posteriors[[i]]}
		brtsi<-getx(phyloi)
		missing.lineages <- total_richness - Ntip(phyloi)
		resi<-10*(total_richness)
	
print("Linear dependence of speciation rate without extinction")
	DDD_1<-dd_ML(brtsi, ddmodel=1, initparsopt=c(0.3, total_richness), idparsopt=c(1,3), idparsfix=c(2), parsfix=c(0), res=resi, missnumspec=missing.lineages, cond=1, btorph=1, soc=2)
	print(DDD_1)

print("Linear dependence of speciation rate with extinction")
	DDD_2<-dd_ML(brtsi, ddmodel=1, initparsopt=c(0.3,0.1, total_richness), res=resi, missnumspec=missing.lineages, cond=1, btorph=1, soc=2)
	print(DDD_2)

print("Exponential dependence of speciation rate with extinction")
	DDD_3<-dd_ML(brtsi, ddmodel=2, initparsopt=c(0.1,0.01, total_richness), res=resi, missnumspec=missing.lineages, cond=1, btorph=1, soc=2)
	print(DDD_3)

print("Linear dependence of extinction rate")
	DDD_4<-dd_ML(brtsi, ddmodel=3, initparsopt=c(0.3,0.1, total_richness), res=resi, missnumspec=missing.lineages, cond=1, btorph=1, soc=2)
print(DDD_4)

print("Exponential dependence of extinction rate")
	DDD_5<-dd_ML(brtsi, ddmodel=4, initparsopt=c(0.1,0.01, total_richness), res=resi, missnumspec=missing.lineages, cond=1, btorph=1, soc=2)
print(DDD_5)

print("Linear dependence of speciation and extinction rates")
	DDD_6<-dd_ML(brtsi, ddmodel=5, initparsopt=c(0.5,0.1, total_richness,0.001), res=resi, missnumspec=missing.lineages, cond=1, btorph=1, soc=2)
print(DDD_6)

	############# RESULTS ###########################################

	results<-matrix(NA,6,8)
		
	colnames(results)<-c("Model","NP","logL","AICc","Lambda","Mu","K","r")
	results[,1]<-c("DDL","DDL+E","DDX+E","DD+EL","DD+EX","DDL+EL")
		
	#NP
	results[1,2]<-round(as.numeric(DDD_1[5]))
	results[2,2]<-round(as.numeric(DDD_2[5]))
	results[3,2]<-round(as.numeric(DDD_3[5]))
	results[4,2]<-round(as.numeric(DDD_4[5]))
	results[5,2]<-round(as.numeric(DDD_5[5]))
	results[6,2]<-round(as.numeric(DDD_6[6]))
	
	#logL
	results[1,3]<-round(as.numeric(DDD_1[4]),4)
	results[2,3]<-round(as.numeric(DDD_2[4]),4)
	results[3,3]<-round(as.numeric(DDD_3[4]),4)
	results[4,3]<-round(as.numeric(DDD_4[4]),4)
	results[5,3]<-round(as.numeric(DDD_5[4]),4)
	results[6,3]<-round(as.numeric(DDD_6[5]),4)

	#AICc
	results[1,4]<-round((2*(-round(as.numeric(DDD_1[4]),4))+2*round(as.numeric(DDD_1[5]))+(2*round(as.numeric(DDD_1[5]))*(round(as.numeric(DDD_1[5]))+1))/(Ntip(phyloi)-round(as.numeric(DDD_1[5]))-1)),3)
	results[2,4]<-round((2*(-round(as.numeric(DDD_2[4]),4))+2*round(as.numeric(DDD_2[5]))+(2*round(as.numeric(DDD_2[5]))*(round(as.numeric(DDD_2[5]))+1))/(Ntip(phyloi)-round(as.numeric(DDD_2[5]))-1)),3)
	results[3,4]<-round((2*(-round(as.numeric(DDD_3[4]),4))+2*round(as.numeric(DDD_3[5]))+(2*round(as.numeric(DDD_3[5]))*(round(as.numeric(DDD_3[5]))+1))/(Ntip(phyloi)-round(as.numeric(DDD_3[5]))-1)),3)
	results[4,4]<-round((2*(-round(as.numeric(DDD_4[4]),4))+2*round(as.numeric(DDD_4[5]))+(2*round(as.numeric(DDD_4[5]))*(round(as.numeric(DDD_4[5]))+1))/(Ntip(phyloi)-round(as.numeric(DDD_4[5]))-1)),3)
	results[5,4]<-round((2*(-round(as.numeric(DDD_5[4]),4))+2*round(as.numeric(DDD_5[5]))+(2*round(as.numeric(DDD_5[5]))*(round(as.numeric(DDD_5[5]))+1))/(Ntip(phyloi)-round(as.numeric(DDD_5[5]))-1)),3)
	results[6,4]<-round((2*(-round(as.numeric(DDD_6[5]),4))+2*round(as.numeric(DDD_6[6]))+(2*round(as.numeric(DDD_6[6]))*(round(as.numeric(DDD_6[6]))+1))/(Ntip(phyloi)-round(as.numeric(DDD_6[6]))-1)),3)
	
	#Lambda
	results[1,5]<-round(as.numeric(DDD_1[1]),4)
	results[2,5]<-round(as.numeric(DDD_2[1]),4)
	results[3,5]<-round(as.numeric(DDD_3[1]),4)
	results[4,5]<-round(as.numeric(DDD_4[1]),4)
	results[5,5]<-round(as.numeric(DDD_5[1]),4)
	results[6,5]<-round(as.numeric(DDD_6[1]),4)

	#Mu
	results[2,6]<-round(as.numeric(DDD_2[2]),5)
	results[3,6]<-round(as.numeric(DDD_3[2]),5)
	results[4,6]<-round(as.numeric(DDD_4[2]),5)
	results[5,6]<-round(as.numeric(DDD_5[2]),5)
	results[6,6]<-round(as.numeric(DDD_6[2]),5)

	#K
	results[1,7]<-round(as.numeric(DDD_1[3]),2)
	results[2,7]<-round(as.numeric(DDD_2[3]),2)
	results[3,7]<-round(as.numeric(DDD_3[3]),2)
	results[4,7]<-round(as.numeric(DDD_4[3]),2)
	results[5,7]<-round(as.numeric(DDD_5[3]),2)
	results[6,7]<-round(as.numeric(DDD_6[3]),2)
	
	#r
	results[6,8]<-round(as.numeric(DDD_6[4]),4)
	
		final[[i]]<-results
		print(results)
}
	
	final_table_tree_file<-tables.summary(final)

	fname <- no.extension(basename(tree_file_name))
	outfile <- paste(dirname(tree_file_name), "/", fname, "_results_DDD.txt", sep="")
	out_R <- paste(dirname(tree_file_name), "/", fname, "_complete_results_DDD.Rdata", sep="")

	write.table(final_table_tree_file,file=outfile,quote=FALSE,sep="\t",row.names=FALSE)
	save(final_table_tree_file,file=out_R)
	
	print(final_table_tree_file)
	print("Analyses with the diversity-dependent (DDD) models completed ! Check your results in the working directory.")
}



## TreePar models: time-dependent diversification models with discrete rates through time
# If using this appraoch, please cite:
# Stadler T. 2011. Mammalian phylogeny reveals recent diversification rate shifts. Proc. Natl. Acad. Sci. USA 108:6187–6192.

run_TreePar <- function (tree_file, sampling_fraction=1, grid=0.1, number_of_trees=1)
{
	tree_file_name<-tree_file
	tree_file<-read.nexus(tree_file)
	posteriors<-sample(tree_file, number_of_trees)

	final<-list()

	for (i in 1:number_of_trees)
	{
		print(i)
		
		if (length(posteriors)==1){phyloi<-tree_file} else {phyloi<-posteriors[[i]]}
		brtsi<-getx(phyloi)

		BD_shifts<-bd.shifts.optim(brtsi,c(sampling_fraction,1,1,1,1),grid=grid,start=0,end=ceiling(max(brtsi)),yule=FALSE,ME=FALSE,all=FALSE,posdiv=FALSE)
	
		res<-BD_shifts[[2]]

	############# RESULTS ###########################################
	
		results<-matrix(NA,5,18)
	
		colnames(results)<-c("Model","NP","logL","AICc","DivRate1","Turnover1","ShiftTime1","DivRate2","Turnover2","ShiftTime2","DivRate3","Turnover3","ShiftTime3","DivRate4","Turnover4","ShiftTime4","DivRate5","Turnover5")
		results[,1]<-c("NoShiftTime","1ShiftTime","2ShiftTimes","3ShiftTimes","4ShiftTimes")
	
	#NP
		results[1,2]<-2
		results[2,2]<-5
		results[3,2]<-8
		results[4,2]<-11
		results[5,2]<-14

	#logL
		results[1,3]<-round(-res[[1]][1],4)
		results[2,3]<-round(-res[[2]][1],4)
		results[3,3]<-round(-res[[3]][1],4)
		results[4,3]<-round(-res[[4]][1],4)
		results[5,3]<-round(-res[[5]][1],4)

	#AICc
		results[1,4]<-round((2*(-round(-res[[1]][1],4))+2*2+(2*2*(2+1))/((length(brtsi)+1)-2-1)),3)
		results[2,4]<-round((2*(-round(-res[[2]][1],4))+2*5+(2*5*(5+1))/((length(brtsi)+1)-5-1)),3)
		results[3,4]<-round((2*(-round(-res[[3]][1],4))+2*8+(2*8*(8+1))/((length(brtsi)+1)-8-1)),3)
		results[4,4]<-round((2*(-round(-res[[4]][1],4))+2*11+(2*11*(11+1))/((length(brtsi)+1)-11-1)),3)
		results[5,4]<-round((2*(-round(-res[[2]][1],4))+2*14+(2*14*(14+1))/((length(brtsi)+1)-14-1)),3)

	#DivRate1
		results[1,5]<-round(res[[1]][3],4)
		results[2,5]<-round(res[[2]][4],4)
		results[3,5]<-round(res[[3]][5],4)
		results[4,5]<-round(res[[4]][6],4)
		results[5,5]<-round(res[[5]][7],4)
	#Turnover1
		results[1,6]<-round(res[[1]][2],4)
		results[2,6]<-round(res[[2]][2],4)
		results[3,6]<-round(res[[3]][2],4)
		results[4,6]<-round(res[[4]][2],4)
		results[5,6]<-round(res[[5]][2],4)
	#ShitTime1
		results[2,7]<-round(res[[2]][6],4)
		results[3,7]<-round(res[[3]][8],4)
		results[4,7]<-round(res[[4]][10],4)
		results[5,7]<-round(res[[5]][12],4)

	#DivRate2
		results[2,8]<-round(res[[2]][5],4)
		results[3,8]<-round(res[[3]][6],4)
		results[4,8]<-round(res[[4]][7],4)
		results[5,8]<-round(res[[5]][8],4)
	#Turnover2
		results[2,9]<-round(res[[2]][3],4)
		results[3,9]<-round(res[[3]][3],4)
		results[4,9]<-round(res[[4]][3],4)
		results[5,9]<-round(res[[5]][3],4)
	#ShitTime2
		results[3,10]<-round(res[[3]][9],4)
		results[4,10]<-round(res[[4]][11],4)
		results[5,10]<-round(res[[5]][13],4)

	#DivRate3
		results[3,11]<-round(res[[3]][7],4)	
		results[4,11]<-round(res[[4]][8],4)
		results[5,11]<-round(res[[5]][9],4)
	#Turnover3
		results[3,12]<-round(res[[3]][4],4)
		results[4,12]<-round(res[[4]][4],4)
		results[5,12]<-round(res[[5]][4],4)
	#ShitTime3
		results[4,13]<-round(res[[4]][12],4)
		results[5,13]<-round(res[[5]][14],4)

	#DivRate4
		results[4,14]<-round(res[[4]][9],4)
		results[5,14]<-round(res[[5]][10],4)
	#Turnover4
		results[4,15]<-round(res[[4]][5],4)
		results[5,15]<-round(res[[5]][5],4)
	#ShitTime4
		results[5,16]<-round(res[[5]][15],4)

	#DivRate5
		results[5,17]<-round(res[[5]][11],4)
	#Turnover5
		results[5,18]<-round(res[[5]][6],4)
	
		final[[i]]<-results
		#print(results)
	}

	final_table_tree_file<-tables.summary(final)

	fname <- no.extension(basename(tree_file_name))
	outfile <- paste(dirname(tree_file_name), "/", fname, "_results_TreePar.txt", sep="")
	out_R <- paste(dirname(tree_file_name), "/", fname, "_complete_results_TreePar.Rdata", sep="")

	write.table(final_table_tree_file,file=outfile,quote=FALSE,sep="\t",row.names=FALSE)
	save(final_table_tree_file,file=out_R)
	
	print(final_table_tree_file)
	print("Analyses with the time-dependent (TreePar) models completed ! Check your results in the working directory.")
}
