# This script is used to set up an XML file to run a latent liability analysis of discrete and combined trait data.
# The XML file will run an analysis of ONLY the latent liability model, if you want a combined analysis with, e.g., phylogeny estimation, 
#     you will need to copy and paste the elements of this file into the corresponding components of an XML set up for the phylogeny estimation
# This script can only set up an analysis for BINARY discrete characters
# This script is set up to in principle allow missing discrete or continuous data. It is advisable, however, to remove species that are missing data 

# arguments: 
#			file: file name to which to output
#			latent.liability.info: the output from prepForLatentLiability()
#			tree: the (fixed) tree with branch lengths you want to analyze, as an object of class phy
#			ngen: the desired number of MCMC cycles for which to run
#			log.every: the frequency, in number of MCMC cycles, with which to sample the MCMC
#			walk.or.scale: "walk"|"scale", should the MCMC propose changes to the latent liability model using a randomWalk operator or a scaleOperator
#			name.for.traits: the ID reference for the traits in the latent liability block, not particularly important, default "Traits"
#			jitter: single value to be used in the XML for jittering, default 0.01
#			precision: single value to be used for the diagonal of the precision matrix prior, default 0.01
#			wishart: single value to be used for the diagonal of the wishart prior, default 1
#			autoOptimize: boolean, should autoOptimize be on or off, default TRUE

# WARNING: this script has minimal checking for user errors, xml files should be tested locally before they are assumed to work properly

library(ape)

printLatentLiability = function(file,latent.liability.info,tree,log.name,ngen,log.every,walk.or.scale,name.for.traits="Traits",jitter=0.01,precision=0.01,wishart=1,autoOptimize=TRUE,...){
	# recover()
  table = latent.liability.info
	taxa = row.names(table)
	types = names(table)
	d.or.c <- names(table)
	
	ntraits= length(types)
	
	if (length(taxa) != length(tree$tip.label)) {
		stop("ERROR: tree size not equal to number of species in table")
	}
	if (any(grepl(".",tree$tip.label,fixed=T)) || any(grepl(".",taxa,fixed=T))) {
	  stop("ERROR: taxon names in tree file and data frame must not contain \".\" in")
	}
	
	tree.string = write.tree(tree)
	tree.string = gsub("(","",tree.string,fixed=TRUE)
	tree.string = gsub(")","",tree.string,fixed=TRUE)
	tree.string = gsub("[0-9]","",tree.string,fixed=FALSE)
	tree.string = gsub(".","",tree.string,fixed=TRUE)
	tree.string = gsub(":","",tree.string,fixed=TRUE)
	tree.string = gsub(";","",tree.string,fixed=TRUE)
	ordered.taxa = strsplit(tree.string,split=",")[[1]]
  
	alignment = c()
	# missing data coded as ?
	for (i in 1:length(taxa)) {
		tmp = table[i,]
		tmp[which(tmp >= 0)] = 1
		tmp[which(tmp < 0)] = 0
		tmp[which(table[i,] == "?")] = "?"
		alignment = rbind(alignment,tmp)
	}
	
	types[which(tolower(types) == "c")] = "1.0"
	types[which(tolower(types) == "d")] = "2.0"
	# recover()
	if (hasArg(is.multistate)) {
		is.multistate <- list(...)$is.multistate
  	if (hasArg(is.unordered)) {
  	  is.unordered <- list(...)$is.unordered
  	} else {is.unordered <- 0}
		given.states <- as.numeric(unique(table[,is.multistate]))
		if (min(given.states) != 0 || max(given.states) != length(given.states)-1) {
		  stop("Please ensure that the multistate trait has values ranging from 0 to n (the number of states)")
		}
		alignment[,is.multistate] <- as.numeric(table[,is.multistate])
		# alignment[which(table[,is.multistate] == "0.0"),is.multistate] <- 0
		# alignment[which(table[,is.multistate] == "1.0"),is.multistate] <- 1
		# alignment[which(table[,is.multistate] == "2.0"),is.multistate] <- 2
	  types[is.multistate] = as.character(length(given.states))
		if (is.unordered != 0) {
		  if (is.unordered != is.multistate && length(is.unordered != 1)) {
		    stop("Sorry, cannot handle multiple unordered traits yet")
		  }
		  warning("Setting up this model is not fully understood.")
		  d.or.c.left <- names(table)[1:(is.unordered-1)]
		  d.or.c.right <- names(table)[(is.unordered+1):ncol(table)]
		  d.or.c.tmp <- rep("d",length(given.states)-1)
		  tmp.table.left <- table[,1:(is.unordered-1)]
		  tmp.table.right <- table[,(is.unordered+1):ncol(table)]
		  tmp.insert <- matrix(ncol=length(given.states)-1,nrow=nrow(table),data=0) 
		  multi <- table[,is.multistate]
		  for (x in 1:max(given.states)) {
		    tmp.insert[which(multi == x),x] <- 1
		  }
# 		  tmp.insert[which(multi == "1.0"),1] <- 1
# 		  tmp.insert[which(multi == "2.0"),2] <- 1
# 		  # tmp.insert[which(multi == "0.0"),1] <- 1
# # 		  tmp.insert[which(multi == "1.0"),2] <- 1
# # 		  tmp.insert[which(multi == "2.0"),2] <- 2
		  if (is.unordered == 1) {
		    table <- data.frame(tmp.insert,tmp.table.right,stringsAsFactors=F)
		    d.or.c <- c(d.or.c.tmp,d.or.c.right)
		  } else {
		    table <- data.frame(tmp.table.left,tmp.insert,tmp.table.right,stringsAsFactors=F)
		    d.or.c <- c(d.or.c.left,d.or.c.tmp,d.or.c.right)
		  }
		}
	  ntraits <- ncol(table)
	} else {
	  d.or.c <- names(table)
	}
	
	mask = c() # mask is a vector in the XML that tells BEAST what traits for which to estimate tip values, this will include all discrete traits and any traits missing data
	for (masked.taxon in ordered.taxa) {
	  m = rep(0,ntraits)
	  m[which(d.or.c == "d")] = 1
	  tmp = table[(which(taxa == masked.taxon)),]
	  m[which(tmp == "?")] <- 1
	  mask = c(mask,m)
	}

	precisionMatrix = matrix(nrow=ntraits,ncol=ntraits,data=0)
	diag(precisionMatrix) = precision
	row.names(precisionMatrix) = rep(paste('\t\t\t\t<parameter value="',sep=""),ntraits)
	
	wishartMatrix = matrix(nrow=ntraits,ncol=ntraits,data=0)
	diag(wishartMatrix) = wishart
	row.names(wishartMatrix) = rep(paste('\t\t\t\t<parameter value="',sep=""),ntraits)
	
	degreesFreedom = ntraits + 2 # this is following Moore et al who were following recommendations from Cybis et al
	
	cat('<?xml version="1.0"?>',sep="\n",file=file)
	cat('<beast>',append=TRUE,sep="\n",file=file)
	cat(paste('<!-- ntax=',length(taxa),'\t\t\t\t\t\t\t\t\t\t\t\t-->',sep=""),append=TRUE,sep="\n",file=file)
	# cat('\t<!-- Trait columns are ordered and coding\t\t\t\t\t-->',append=TRUE,sep="\n",file=file)

	cat('<taxa id="taxa">',append=TRUE,sep="\n",file=file)
	for(i in 1:length(taxa)) {
		cat(paste('\t<taxon id="',taxa[i],'"> <attr name="latent"> ',paste(table[i,],collapse=" "),' </attr></taxon>',sep=""),append=TRUE,sep="\n",file=file)
	}
	cat('</taxa>',append=TRUE,sep="\n",file=file)

	cat('<generalDataType id="multinomial">',append=TRUE,sep="\n",file=file)
	cat('\t<state code="0"/>',append=TRUE,sep="\n",file=file)
	cat('\t<state code="1"/>',append=TRUE,sep="\n",file=file)
	if(hasArg(is.multistate)) {
	  for (i in 2:max(given.states)) {
	    cat(paste0('\t<state code="',i,'"/>'),append=TRUE,sep="\n",file=file)
	  }
	}
	cat('</generalDataType>',append=TRUE,sep="\n",file=file)
	
	cat('<alignment id="alignment">',append=TRUE,sep="\n",file=file)
	cat('<dataType idref="multinomial"/>',append=TRUE,sep="\n",file=file)
	for(i in 1:length(taxa)) {
		cat(paste('\t<sequence>\n\t\t<taxon idref="',taxa[i],'"/>\n\t\t',paste(alignment[i,],collapse=""),'\n\t</sequence>',sep=""),append=TRUE,sep="\n",file=file)
	}
	cat('</alignment>',append=TRUE,sep="\n",file=file)
	
	cat('\n\t<!-- The unique patterns\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t-->',append=TRUE,sep="\n",file=file)
	cat('\t<patterns id="patterns" from="1" unique="false">',append=TRUE,sep="\n",file=file)
	cat('\t\t<alignment idref="alignment"/>',append=TRUE,sep="\n",file=file)
	cat('\t</patterns>',append=TRUE,sep="\n",file=file)

	cat('\n\t<!-- The is the fixed tree topology and branch lengths:\t\t\t\t\t\t\t\t\t\t\t\t-->',append=TRUE,sep="\n",file=file)
	cat('\t<newick id="startingTree">',append=TRUE,sep="\n",file=file)
	cat('\t\t',write.tree(tree),'\n',append=TRUE,sep="",file=file)
	cat('\t</newick>',append=TRUE,sep="\n",file=file)
	cat('\t<treeModel id="treeModel">',append=TRUE,sep="\n",file=file)
	cat('\t\t<tree idref="startingTree"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t<rootHeight>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<parameter id="treeModel.rootHeight"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t</rootHeight>',append=TRUE,sep="\n",file=file)
	cat('\t\t<nodeHeights internalNodes="true">',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<parameter id="treeModel.internalNodeHeights"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t</nodeHeights>',append=TRUE,sep="\n",file=file)
	cat('\t\t<nodeHeights internalNodes="true" rootNode="true">',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<parameter id="treeModel.allInternalNodeHeights"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t</nodeHeights>',append=TRUE,sep="\n",file=file)
	cat(paste('\t\t<nodeTraits rootNode="false" internalNodes="false" leafNodes="true" traitDimension="',ntraits,'" name="latent">',sep=""),append=TRUE,sep="\n",file=file)
	cat(paste('\t\t\t<parameter id="',name.for.traits,'"/>',sep=""),append=TRUE,sep="\n",file=file)
	cat('\t\t</nodeTraits>',append=TRUE,sep="\n",file=file)
	cat('\t</treeModel>',append=TRUE,sep="\n",file=file)
	cat('\t<report>',append=TRUE,sep="\n",file=file)
	cat('\t\tNewick Tree:',append=TRUE,sep="\n",file=file)
	cat('\t\t<tree idref="startingTree"/>',append=TRUE,sep="\n",file=file)
	cat(paste('\t\t<parameter idref="',name.for.traits,'"/>',sep=""),append=TRUE,sep="\n",file=file)
	cat('\t</report>',append=TRUE,sep="\n",file=file)

	cat('\t<!-- Begin latent-liability model\t\t\t\t\t\t\t\t\t\t\t\t\t\t-->',append=TRUE,sep="\n",file=file)
	cat('\t<multivariateDiffusionModel id="diffusionModel">',append=TRUE,sep="\n",file=file)
	cat('\t\t<precisionMatrix>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<matrixParameter id="precisionMatrix">',append=TRUE,sep="\n",file=file)
	write.table(precisionMatrix,file=file,sep=" ",eol='"/>\n',row.names=TRUE,col.names=FALSE,quote=FALSE,append=TRUE)
	cat('\t \t\t</matrixParameter>',append=TRUE,sep="\n",file=file)
	cat('\t\t</precisionMatrix>',append=TRUE,sep="\n",file=file)
	cat('\t</multivariateDiffusionModel>',append=TRUE,sep="\n",file=file)
	cat('',append=TRUE,sep="\n",file=file)
	
	cat('\t<!-- We have used the recommended value for df of k traits +2\t\t\t\t\t\t-->',append=TRUE,sep="\n",file=file)
	cat(paste('\t<multivariateWishartPrior id="precisionPrior" df="',degreesFreedom,'">'),append=TRUE,sep="\n",file=file)
	cat('\t\t<scaleMatrix>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<matrixParameter>',append=TRUE,sep="\n",file=file)
	write.table(wishartMatrix,file=file,sep=" ",eol='"/>\n',row.names=TRUE,col.names=FALSE,quote=FALSE,append=TRUE)
	cat('\t \t\t</matrixParameter>',append=TRUE,sep="\n",file=file)
	cat('\t\t</scaleMatrix>',append=TRUE,sep="\n",file=file)
	cat('\t\t<data>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<matrixParameter idref="precisionMatrix"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t</data>',append=TRUE,sep="\n",file=file)
	cat('\t</multivariateWishartPrior>',append=TRUE,sep="\n",file=file)
	cat('',append=TRUE,sep="\n",file=file)
	
	cat('\t<multivariateTraitLikelihood id="traitLikelihood" traitName="latent"',append=TRUE,sep="\n",file=file)
	cat('\t \t\t\t\t\t\t\t useTreeLength="true" scaleByTime="true"',append=TRUE,sep="\n",file=file)
	cat('\t\t\t\t\t\t\t\t cacheBranches="true"',append=TRUE,sep="\n",file=file)
	cat('\t\t\t\t\t\t\t\t reportAsMultivariate="true" reciprocalRates="true" ',append=TRUE,sep="\n",file=file)
	cat('\t\t\t\t\t\t\t\t integrateInternalTraits="true">',append=TRUE,sep="\n",file=file)
	cat('\t\t<multivariateDiffusionModel idref="diffusionModel"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t<treeModel idref="treeModel"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t<traitParameter>',append=TRUE,sep="\n",file=file)
	cat(paste('\t\t\t<parameter idref="',name.for.traits,'"/>',sep=""),append=TRUE,sep="\n",file=file)
	cat('\t\t</traitParameter>',append=TRUE,sep="\n",file=file)
	cat('\t\t',append=TRUE,sep="\n",file=file)
	cat(paste('\t\t<jitter window="',paste(rep(jitter,ntraits),collapse=" "),'" duplicatesOnly="true">',sep=""),append=TRUE,sep="\n",file=file)
	cat(paste('\t\t\t<parameter idref="',name.for.traits,'"/>',sep=""),append=TRUE,sep="\n",file=file)
	cat('\t\t</jitter>',append=TRUE,sep="\n",file=file)
	cat('\t\t',append=TRUE,sep="\n",file=file)
	cat('\t\t<conjugateRootPrior>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<meanParameter>',append=TRUE,sep="\n",file=file)
	cat(paste('\t\t\t\t<parameter value="',paste(rep("0.0",ntraits),collapse=" "),'"/>',sep=""),append=TRUE,sep="\n",file=file)
	cat('\t\t\t</meanParameter>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<priorSampleSize>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t\t<parameter value="0.01"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t</priorSampleSize>',append=TRUE,sep="\n",file=file)
	cat('\t\t</conjugateRootPrior>',append=TRUE,sep="\n",file=file)
	cat('\t</multivariateTraitLikelihood>',append=TRUE,sep="\n",file=file)
	
	cat('',append=TRUE,sep="\n",file=file)
	cat('\t<report>',append=TRUE,sep="\n",file=file)
	cat('\t\t<maskedParameter id= "latentParameters" >',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<!-- Define mask to be 1 for: a) every discrete trait b) every continuous trait missing data-->',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<mask>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t\t<parameter id="mask" value="',paste(mask,collapse=" "),'"  />\n',append=TRUE,sep="",file=file)
	cat('\t\t\t</mask>',append=TRUE,sep="\n",file=file)
	cat(paste('\t\t\t<parameter idref="',name.for.traits,'"/>',sep=""),append=TRUE,sep="\n",file=file)
	cat('\t\t</maskedParameter>',append=TRUE,sep="\n",file=file)
	cat('\t</report>',append=TRUE,sep="\n",file=file)
	cat('',append=TRUE,sep="\n",file=file)
	cat('\t<CompoundParameter id="threshold">',append=TRUE,sep="\n",file=file)
	cat('\t\t<parameter value="1"/>',append=TRUE,sep="\n",file=file)
	cat('\t</CompoundParameter>',append=TRUE,sep="\n",file=file)
	if (hasArg(is.multistate) && is.unordered != 0) {
	  cat('\t<orderedLatentLiabilityLikelihood id="liabilityLikelihood" isUnordered="true">',append=TRUE,sep="\n",file=file)
	} else {
	  cat('\t<orderedLatentLiabilityLikelihood id="liabilityLikelihood">',append=TRUE,sep="\n",file=file)
	}
	cat('\t\t<patterns idref="patterns"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t<treeModel idref="treeModel"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t<multivariateTraitLikelihood idref="traitLikelihood"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t<tipTrait>',append=TRUE,sep="\n",file=file)
	cat(paste('\t\t\t<parameter idref="',name.for.traits,'"/>',sep=""),append=TRUE,sep="\n",file=file)
	cat('\t\t</tipTrait>',append=TRUE,sep="\n",file=file)
	cat('\t\t<threshold>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<parameter idref="threshold"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t</threshold>',append=TRUE,sep="\n",file=file)
	cat('\t\t<numClasses>',append=TRUE,sep="\n",file=file)
	cat(paste('\t\t\t<parameter value = "',paste(types,collapse=" "),'"/>',sep=""),append=TRUE,sep="\n",file=file)
	cat('\t\t</numClasses>',append=TRUE,sep="\n",file=file)
	cat('\t</orderedLatentLiabilityLikelihood>',append=TRUE,sep="\n",file=file)
	cat('',append=TRUE,sep="\n",file=file)
	cat('',append=TRUE,sep="\n",file=file)
	
	cat('\t<!-- We are only proposing updates to parameters of the latent liability model\t\t  -->',append=TRUE,sep="\n",file=file)
	cat('\t<operators id="operators" optimizationSchedule="log">',append=TRUE,sep="\n",file=file)
	if ( length(grep("w",walk.or.scale)) > 0) {
		cat('\t\t<randomWalkOperator windowSize="1.0" weight="100" >',append=TRUE,sep="\n",file=file)
		cat('\t\t\t<maskedParameter idref="latentParameters"/>',append=TRUE,sep="\n",file=file)
		cat('\t\t</randomWalkOperator>',append=TRUE,sep="\n",file=file)
	} else if ( length(grep("s",walk.or.scale)) > 0) {
		cat('\t\t<scaleOperator scaleFactor="0.2" weight="100" >',append=TRUE,sep="\n",file=file)
		cat('\t\t\t<maskedParameter idref="latentParameters"/>',append=TRUE,sep="\n",file=file)
		cat('\t\t</scaleOperator>',append=TRUE,sep="\n",file=file)
	}
	else {
		return('Please use input option scale.or.walk to specify scaleOperator or randomWalkOperator')
		cat('',append=TRUE,sep="\n",file=file)
	}
	if ( hasArg(is.multistate) ) {
	  cat('\t\t<randomWalkOperator windowSize="0.01" weight="1" >',append=TRUE,sep="\n",file=file)
	  cat('\t\t\t<compountParameter idref="threshold"/>',append=TRUE,sep="\n",file=file)
	  cat('\t\t</randomWalkOperator>',append=TRUE,sep="\n",file=file)
	}
	cat('\t\t<precisionGibbsOperator weight="1">',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<multivariateTraitLikelihood idref="traitLikelihood"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<multivariateWishartPrior idref="precisionPrior"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t</precisionGibbsOperator>',append=TRUE,sep="\n",file=file)
	cat('\t</operators>',append=TRUE,sep="\n",file=file)
	cat('',append=TRUE,sep="\n",file=file)
	if (autoOptimize == T) {
		cat(paste('\t<mcmc id="mcmc" chainLength="',ngen,'" autoOptimize="true" operatorAnalysis="',log.name,'_autoOptimize.ops">',sep=""),append=TRUE,sep="\n",file=file)
	}

	else if (autoOptimize == F) {
		cat(paste('\t<mcmc id="mcmc" chainLength="',ngen,'" autoOptimize="false">',sep=""),append=TRUE,sep="\n",file=file)
	}
	cat('\t\t<posterior id="posterior">',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<prior id="prior">',append=TRUE,sep="\n",file=file)
	cat('\t\t\t\t<multivariateWishartPrior idref="precisionPrior"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t\t<gammaPrior shape="1" scale="0.5" offset="0.0">',append=TRUE,sep="\n",file=file)
	cat('\t\t\t\t\t<parameter idref="threshold"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t\t</gammaPrior>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t</prior>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<likelihood id="likelihood">',append=TRUE,sep="\n",file=file)
	cat('\t\t\t\t<multivariateTraitLikelihood idref="traitLikelihood"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t\t<MultinomialLatentLiabilityLikelihood idref="liabilityLikelihood"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t</likelihood>',append=TRUE,sep="\n",file=file)
	cat('\t\t</posterior>',append=TRUE,sep="\n",file=file)
	cat('',append=TRUE,sep="\n",file=file)
	
	cat('\t\t<operators idref="operators"/>',append=TRUE,sep="\n",file=file)
	cat('',append=TRUE,sep="\n",file=file)
	
	cat(paste('\t\t<log id="screenLog" logEvery="',log.every,'">',sep=""),append=TRUE,sep="\n",file=file)
	cat('\t\t\t<column label="Posterior" dp="4" width="12">',append=TRUE,sep="\n",file=file)
	cat('\t\t\t\t<posterior idref="posterior"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t</column>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<column label="Prior" dp="4" width="12">',append=TRUE,sep="\n",file=file)
	cat('\t\t\t\t<prior idref="prior"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t</column>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<column label="Likelihood" dp="4" width="12">',append=TRUE,sep="\n",file=file)
	cat('\t\t\t\t<likelihood idref="likelihood"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t</column>',append=TRUE,sep="\n",file=file)
# 	cat(paste('\t\t\t<column label="',name.for.traits,'" sf="6" width="12">',sep=""),append=TRUE,sep="\n",file=file)
# 	cat('\t\t\t <maskedParameter idref="latentParameters"/>',append=TRUE,sep="\n",file=file)
# 	cat(paste('\t\t\t <!--<parameter idref="',name.for.traits,'"/>\t\t  -->',sep=""),append=TRUE,sep="\n",file=file)
# 	cat('\t\t\t</column>\t  ',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<column label="threshold" sf="6" width="12">',append=TRUE,sep="\n",file=file)
	cat('\t\t\t <CompoundParameter idref="threshold"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t</column>',append=TRUE,sep="\n",file=file)
	cat('\t\t</log>',append=TRUE,sep="\n",file=file)
	cat('',append=TRUE,sep="\n",file=file)
	
	cat(paste('\t\t<log id="fileLog1" logEvery="',log.every,'" fileName="',log.name,'.log">',sep=""),append=TRUE,sep="\n",file=file)
	cat('\t\t\t<posterior idref="posterior"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<prior idref="prior"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<likelihood idref="likelihood"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<CompoundParameter idref="threshold"/>',append=TRUE,sep="\n",file=file)
	cat(paste('\t\t\t<parameter idref="',name.for.traits,'"/>',sep=""),append=TRUE,sep="\n",file=file)
	cat('\t\t\t<matrixParameter idref="precisionMatrix"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t</log>',append=TRUE,sep="\n",file=file)
	cat('\t\t',append=TRUE,sep="\n",file=file)

	cat(paste('\t\t<log id="fileLog2" logEvery="',log.every,'" fileName="',log.name,'_PRECISION.log">',sep=""),append=TRUE,sep="\n",file=file)
	cat('\t\t\t<posterior idref="posterior"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<prior idref="prior"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<likelihood idref="likelihood"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t\t<matrixParameter idref="precisionMatrix"/>',append=TRUE,sep="\n",file=file)
	cat('\t\t</log>',append=TRUE,sep="\n",file=file)
	cat('\t\t',append=TRUE,sep="\n",file=file)

	cat('\t\t',append=TRUE,sep="\n",file=file)
	cat('\t</mcmc>',append=TRUE,sep="\n",file=file)
	cat('</beast>',append=TRUE,sep="\n",file=file)
	
}


# trait.csv: the file path to a CSV containing the species names (column 1), and the traits, with a header giving the names of the traits
# transformation: a function for transforming the continuous traits, default is taking the natural log
# species: a character vecotr specifying which species in the trait table are to be included in the XML, default "all" (ie all species in trait.csv are to be analyzed)
# sep: argument to be passed to read.csv, the separator in trait.csv, default ","

# prepForLatentLiability <- function(trait.csv,exclude.columns=NA,transformation=log,species="all",sep=",") {
#   traits =  read.csv(trait.csv,sep=sep,row.names=1)
#   if (species[1] != "all") {
#     traits = traits[which(traits[,1] %in% species),]
#   }
#   if (any(is.na(traits))) {
#     warning("There are NAs in trait.csv. Missing data is usually undesirable in a latent liability analysis")
#   }
#  # recover()
#   as.data.frame(trait.csv,stringsAsFactors=FALSE)
# }
# 
# setwd("~/Desktop/beast/")
# 
# prepForLatentLiability <- function(trait.csv,exclude.columns=NA,transformation=log,species="all",sep=",") {
#   traits =  read.csv(trait.csv,sep=sep,row.names=1)
#   if (species[1] != "all") {
#     traits = traits[which(traits[,1] %in% species),]
#   }
#   if (any(is.na(traits))) {
#     warning("There are NAs in trait.csv. Missing data is usually undesirable in a latent liability analysis")
#   }
#   # recover()
#   as.data.frame(,stringsAsFactors=FALSE)
# }
#  
# prep <- prepForLatentLiability("discrete_traits.csv", species = read.nexus("~/Desktop/beast/pomacentridae_POSTERIOR_ucld_MCC.tre")$tip.label)
# 
