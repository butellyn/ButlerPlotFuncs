## R Plot Functions


## Load library(s)
#psych failed
#FOR 3.5.1: install_load('scales', 'ggplot2', 'mgcv', 'stringi', 'plyr', 'reshape2', 'gridExtra')
source("/home/ebutler/adroseHelperScripts/R/afgrHelpFunc.R")
install_load('ggplot2', 'mgcv', 'stringi', 'plyr', 'reshape2', 'gridExtra')
#Oct 15, 2018: stringi failed, need 1.1.7


## Declare a function to create mean-centered age variables (((THIS WORKS)))
ageVars <- function(dataframe, ageColNum) {
	dataframe$age_MC <- dataframe[,ageColNum] - mean(dataframe[,ageColNum])
	dataframe$ageSq_MC <- dataframe[,ageColNum]^2 - mean(dataframe[,ageColNum]^2)
	dataframe$ageCu_MC <- dataframe[,ageColNum]^3 - mean(dataframe[,ageColNum]^3)
	return(dataframe)
}


## Declare a function to remove relationship with other covariates and create z-scores of residuals (((THIS WORKS)))
# Input data types: dataframe, string, vector (put in empty vector if you just want z-scores)
regressOutVars <- function(dataframe, DV, IVColumnNums) {
	if (length(IVColumnNums) >= 1) {
		dataframe <- cbind(dataframe, NA)
		numcols = ncol(dataframe)
		names(dataframe)[numcols] <- paste0(DV, "_regressed")
		IVstring <- paste0("~",colnames(dataframe[IVColumnNums[[1]]]))
		for (i in 2:length(IVColumnNums)) {
			IVstring <- paste0(IVstring, "+", colnames(dataframe[IVColumnNums[[i]]])) 
		}
		dataframe[,numcols] <- residuals(lm(as.formula(paste(DV, IVstring)), data=dataframe, na.action = na.exclude))
		dataframe$newcol <- NA
		names(dataframe)[numcols+1] <- paste0(DV, "_regressed_zscore")
		dataframe[,numcols+1] <- dataframe[,numcols]/sd(dataframe[,numcols], na.rm=TRUE)
		return(dataframe)
	} else if (length(IVColumnNums) == 0) {
		numcols = ncol(dataframe)
		dataframe$newcol <- NA
		names(dataframe)[numcols+1] <- paste0(DV, "_zscore")
		dataframe[,numcols+1] <- (dataframe[,DV] - mean(dataframe[,DV]))/sd(dataframe[,DV], na.rm=TRUE)
		return(dataframe)
	}
}


## Declare a function that calls regressOutVars to perform this for multiple DVs (((THIS WORKS)))
# Input data types: dataframe, vector, vector
regressMultDVs <- function(dataframe, DVColumnNums, IVColumnNums) {
	for (i in 1:length(DVColumnNums)) {
		DV <- colnames(dataframe[DVColumnNums[[i]]])
		dataframe <- regressOutVars(dataframe, DV, IVColumnNums)
	}
	return(dataframe)
}


## Declare a function to create lobular volume variables
# ROI_ListofLists should have sub-elements that match column names in the dataframe
# ROI_ListofLists should contain raw values (NO z-scores)
# The -1 should only be present if the last sublist in ROI_ListofLists is not being used
lobeVolumes <- function(dataframe, ROI_ListofLists, lastList=FALSE) {
	if (lastList == TRUE) {
		max = length(ROI_ListofLists)
	} else {
		max = length(ROI_ListofLists) - 1
	}
	for (i in 1:max) {
		dataframe$newcol <- 0 
		for (j in 1:length(ROI_ListofLists[[i]])) {
			# Get the column with the name ROI_ListofLists[[i]][[j]]
			num <- which(colnames(dataframe) == ROI_ListofLists[[i]][j])
			# Add this column to newcol
			dataframe$newcol = dataframe$newcol + dataframe[ , ROI_ListofLists[[i]][j]]
		}
		# name column in dataframe
		colnames(dataframe)[colnames(dataframe)=="newcol"] <- paste0(names(ROI_ListofLists[i]), "_Vol")
	}
	return(dataframe)
}


## Declare a function to average left and right ROIs on volume (((THIS WORKS)))
averageLeftAndRight_Vol <- function(dataframe, rightString, leftString, averageString) {
	dataframe.right <- dataframe[ , grep(rightString, names(dataframe))]
	dataframe.left <- dataframe[ , grep(leftString, names(dataframe))]
	numnewcols <- ncol(dataframe.right)
	numcolsdataframe <- ncol(dataframe)
	# check if new dataframes are in the same order; if not, make it so
	right_vec <- colnames(dataframe.right)
	left_vec <- colnames(dataframe.left)
	right_neutral_vec <- gsub(pattern=rightString, replacement="", x=right_vec)
	left_neutral_vec <- gsub(pattern=leftString, replacement="", x=left_vec)
	if (identical(right_neutral_vec, left_neutral_vec)){
		# make numnewcols new columns in dataframe
		blankcol <- matrix(nrow=nrow(dataframe), ncol=1)
		for (i in 1:numnewcols) {
			dataframe <- cbind(dataframe, blankcol)
			# name the new columns by replacing rightString with averageString
			names(dataframe)[numcolsdataframe+i] <- gsub(pattern=rightString, replacement=averageString, x=names(dataframe.right)[i])
			# compute the averages
			dataframe[ , numcolsdataframe+i] <- (dataframe.right[,i] + dataframe.left[,i])/2
		}
	return(dataframe) 
	# If the vectors contain the same elements, but in different orders	
	} else if (setequal(right_neutral_vec, left_neutral_vec) && length(right_neutral_vec) == length(left_neutral_vec)) {
		# sort the dataframes
		left_ordered_vec <- gsub(pattern=rightString, replacement=leftString, x=right_vec)
		dataframe.left <- dataframe.left[left_ordered_vec]
		# subset dataframe to only include columns that aren't left or right
		dataframe <- dataframe[ , !(names(dataframe) %in% right_vec)] 
		dataframe <- dataframe[ , !(names(dataframe) %in% left_ordered_vec)]
		# merge in dataframe.right and dataframe.left
		dataframe <- cbind(dataframe, dataframe.right)
		dataframe <- cbind(dataframe, dataframe.left)
		# call the function again
		averageLeftAndRight_Vol(dataframe, rightString, leftString)
	} else {
		print("Your dataframe does not contain the same ROIs for the left and right hemispheres. Please modify your dataframe so that this is true to use this function.")
	}
}

## Declare a function to average left and right ROIs, weighted by volume
averageLeftAndRight_WeightByVol <- function(vol_dataframe, other_dataframe, volString="vol_", otherString="gmd_", rightString="_R_", leftString="_L_", averageString="_ave_") {
	vol_dataframe.right <- vol_dataframe[ , grep(rightString, names(vol_dataframe))]
	vol_dataframe.left <- vol_dataframe[ , grep(leftString, names(vol_dataframe))]
	other_dataframe.right <- other_dataframe[ , grep(rightString, names(other_dataframe))]
	other_dataframe.left <- other_dataframe[ , grep(leftString, names(other_dataframe))]
	numnewcols <- ncol(other_dataframe.right)
	numcolsdataframe <- ncol(other_dataframe)
	# check if new dataframes are in the same order; if not, make it so
	vol_right_vec <- colnames(vol_dataframe.right)
	vol_left_vec <- colnames(vol_dataframe.left)
	vol_right_neutral_vec <- gsub(pattern=rightString, replacement="", x=vol_right_vec)
	vol_right_neutral_vec <- gsub(pattern=volString, replacement="", x=vol_right_neutral_vec)
	vol_left_neutral_vec <- gsub(pattern=leftString, replacement="", x=vol_left_vec)
	vol_left_neutral_vec <- gsub(pattern=volString, replacement="", x=vol_left_neutral_vec)
	other_right_vec <- colnames(other_dataframe.right)
	other_left_vec <- colnames(other_dataframe.left)
	other_right_neutral_vec <- gsub(pattern=rightString, replacement="", x=other_right_vec)
	other_right_neutral_vec <- gsub(pattern=otherString, replacement="", x=other_right_neutral_vec)
	other_left_neutral_vec <- gsub(pattern=leftString, replacement="", x=other_left_vec)
	other_left_neutral_vec <- gsub(pattern=otherString, replacement="", x=other_left_neutral_vec)
	if ((identical(vol_right_neutral_vec, vol_left_neutral_vec)) && (identical(other_right_neutral_vec, other_left_neutral_vec)) && (identical(vol_right_neutral_vec, other_left_neutral_vec))){
		# make numnewcols new columns in dataframe
		blankcol <- matrix(nrow=nrow(other_dataframe), ncol=1)
		for (i in 1:numnewcols) {
			other_dataframe <- cbind(other_dataframe, blankcol)
			# name the new columns by replacing rightString with averageString
			names(other_dataframe)[numcolsdataframe+i] <- gsub(pattern=rightString, replacement=averageString, x=names(other_dataframe.right)[i])
			# compute the averages
			other_dataframe[ , numcolsdataframe+i] <- (other_dataframe.right[,i]*vol_dataframe.right[,i] + other_dataframe.left[,i]*vol_dataframe.left[,i])/(vol_dataframe.right[,i] + vol_dataframe.left[,i])
		}
	return(other_dataframe) 
	# If the vectors contain the same elements, but in different orders	
	} else if (setequal(vol_right_neutral_vec, vol_left_neutral_vec) && length(vol_right_neutral_vec) == length(vol_left_neutral_vec) && setequal(vol_right_neutral_vec, other_left_neutral_vec) && length(vol_right_neutral_vec) == length(other_left_neutral_vec) && setequal(other_right_neutral_vec, other_left_neutral_vec) && length(other_right_neutral_vec) == length(other_left_neutral_vec)) {
		### other
		# sort the dataframes
		other_left_ordered_vec <- gsub(pattern=rightString, replacement=leftString, x=other_right_vec)
		other_dataframe.left <- other_dataframe.left[other_left_ordered_vec]
		# subset dataframe to only include columns that aren't left or right
		other_dataframe <- other_dataframe[ , !(names(other_dataframe) %in% other_right_vec)] 
		other_dataframe <- other_dataframe[ , !(names(other_dataframe) %in% other_left_ordered_vec)]
		# merge in dataframe.right and dataframe.left
		other_dataframe <- cbind(other_dataframe, other_dataframe.right)
		other_dataframe <- cbind(other_dataframe, other_dataframe.left)
		### vol
		# sort the dataframes
		vol_left_ordered_vec <- gsub(pattern=rightString, replacement=leftString, x=other_right_vec)
		vol_left_ordered_vec <- gsub(pattern=otherString, replacement=volString, x=vol_left_ordered_vec)
		vol_dataframe.left <- vol_dataframe.left[vol_left_ordered_vec]
		# subset dataframe to only include columns that aren't left or right
		vol_dataframe <- vol_dataframe[ , !(names(vol_dataframe) %in% vol_right_vec)] 
		vol_dataframe <- other_dataframe[ , !(names(vol_dataframe) %in% vol_left_ordered_vec)]
		# merge in dataframe.right and dataframe.left
		vol_dataframe <- cbind(vol_dataframe, vol_dataframe.right)
		vol_dataframe <- cbind(vol_dataframe, vol_dataframe.left)
		# call the function again
		averageLeftAndRight_WeightByVol(vol_dataframe, other_dataframe, volString, otherString, rightString, leftString, averageString)
	} else {
		print("Your dataframes do not contain the same ROIs for the left and right hemispheres. Please modify your dataframes so that this is true to use this function.")
	}
}

## Declare a function to add a "_" at the end of each element of a character vector
addUnderScore <- function(ROIlist) {
	temp_ROIlist <- c()
	for (i in 1:length(ROIlist)) {
		roi <- paste0(ROIlist[[i]], "_")
		temp_ROIlist <- c(temp_ROIlist, roi)
	}
	ROIlist <- temp_ROIlist
	return(ROIlist)
}

## Declare a function to remove "_" from end of ROI names in output of roiLobes() when lobeDef = FALSE (((THIS WORKS)))
removeUnderScore <- function(ROI_ListofLists) {
	num_lobes <- length(ROI_ListofLists)
	temp_ROI_ListofLists <- ROI_ListofLists
	for (i in 1:num_lobes) {
		lobe_vec <- c()
		for (j in 1:length(ROI_ListofLists[[i]])) {
			withoutUnderScore <- substr(ROI_ListofLists[[i]][j], 1, nchar(ROI_ListofLists[[i]][j])-1)
			lobe_vec <- c(lobe_vec, withoutUnderScore)
		}
		temp_ROI_ListofLists[[i]] <- lobe_vec
	}
	ROI_ListofLists <- temp_ROI_ListofLists
	return(ROI_ListofLists)
}


######### Declare a function to average left and right ROIs, weighted by volume
# vol_dataframe has to have total lobe volumes for this to work
# lastList should be true if you want to use the last list in ROI_ListofLists (Other/Vol)
# only configured for type to be equal to "gmd" or "cort" at the moment
averageROIsInLobes_WeightByVol <- function(vol_dataframe, other_dataframe, ROI_ListofLists_Vol, ROI_ListofLists_Other, lastList=FALSE, type) {
	if (lastList == TRUE) {
		max = length(ROI_ListofLists_Other)
	} else {
		max = length(ROI_ListofLists_Other)-1
	}
	for (i in 1:max) {
		other_dataframe$newcol <- 0 
		for (j in 1:length(ROI_ListofLists_Other[[i]])) {
			# Get the column with the name ROI_ListofLists_Other[[i]][[j]] & Vol
			num_other <- which(colnames(other_dataframe) == ROI_ListofLists_Other[[i]][j])
			num_vol <- which(colnames(vol_dataframe) == ROI_ListofLists_Vol[[i]][j])
			# Add this column to newcol
			other_dataframe$newcol = other_dataframe$newcol + other_dataframe[ , num_other]*vol_dataframe[ , num_vol]
		}
		# divide the column by the total volume for the region
		volcol <- grep(names(ROI_ListofLists_Vol[i]), colnames(vol_dataframe), value=TRUE)
		other_dataframe$newcol = other_dataframe$newcol/vol_dataframe[,volcol]
		# name column in dataframe
		if (type == "gmd") {
			string = "_GMD"
		} else if (type == "cort") {
			string = "_CT"
		}
		colnames(other_dataframe)[colnames(other_dataframe)=="newcol"] <- paste0(names(ROI_ListofLists_Other[i]), string)
	}
	return(other_dataframe)
}


## Declare a function to remove "_" from end of ROI names in output of roiLobes() when lobeDef = FALSE (((THIS WORKS)))
removeUnderScore <- function(ROI_ListofLists) {
	num_lobes <- length(ROI_ListofLists)
	temp_ROI_ListofLists <- ROI_ListofLists
	for (i in 1:num_lobes) {
		lobe_vec <- c()
		for (j in 1:length(ROI_ListofLists[[i]])) {
			withoutUnderScore <- substr(ROI_ListofLists[[i]][j], 1, nchar(ROI_ListofLists[[i]][j])-1)
			lobe_vec <- c(lobe_vec, withoutUnderScore)
		}
		temp_ROI_ListofLists[[i]] <- lobe_vec
	}
	ROI_ListofLists <- temp_ROI_ListofLists
	return(ROI_ListofLists)
}


## Declare a function to create a dataframe with lobular definitions for each ROI (((THIS WORKS)))
# structured for one modality at a time (?)
roiLobes <- function(ROIlist, lobeDef=FALSE) {
	if (lobeDef == FALSE) { # notes: "_Thal" or "_Thalamus_Proper_"
		BasGang_vec <- c(grep("_Thal", ROIlist, value=TRUE), grep("_Putamen_", ROIlist, value=TRUE), grep("_Caudate_", ROIlist, value=TRUE), grep("_Pallidum_", ROIlist, value=TRUE), grep("_Accumbens_", ROIlist, value=TRUE), grep("_BasFor_", ROIlist, value=TRUE)) 
		Limbic_vec <- c(grep("_PHG_", ROIlist, value=TRUE), grep("_Hip_", ROIlist, value=TRUE), grep("_Hippocampus_", ROIlist, value=TRUE), grep("_PIns_", ROIlist, value=TRUE), grep("_SCA_", ROIlist, value=TRUE), grep("_AIns_", ROIlist, value=TRUE), grep("_ACgG_", ROIlist, value=TRUE), grep("_PCgG_", ROIlist, value=TRUE), grep("_Ent_", ROIlist, value=TRUE), grep("_Amygdala_", ROIlist, value=TRUE), grep("_MCgG_", ROIlist, value=TRUE))
		FrontOrb_vec <- c(grep("_FO_", ROIlist, value=TRUE), grep("_MFC_", ROIlist, value=TRUE), grep("_MOrG_", ROIlist, value=TRUE), grep("_POrG_", ROIlist, value=TRUE), grep("_OrIFG_", ROIlist, value=TRUE), grep("_TrIFG_", ROIlist, value=TRUE), grep("_AOrG_", ROIlist, value=TRUE), grep("_OpIFG_", ROIlist, value=TRUE), grep("_GRe_", ROIlist, value=TRUE), grep("_FRP_", ROIlist, value=TRUE), grep("_LOrG_", ROIlist, value=TRUE))
		FrontDors_vec <- c(grep("_PrG_", ROIlist, value=TRUE), grep("_MSFG_", ROIlist, value=TRUE), grep("_SMC_", ROIlist, value=TRUE), grep("_MFG_", ROIlist, value=TRUE), grep("_SFG_", ROIlist, value=TRUE))
		Temporal_vec <- c(grep("_FuG_", ROIlist, value=TRUE), grep("_PT_", ROIlist, value=TRUE), grep("_PP_", ROIlist, value=TRUE), grep("_ITG_", ROIlist, value=TRUE), grep("_CO_", ROIlist, value=TRUE), grep("_MTG_", ROIlist, value=TRUE), grep("_TMP_", ROIlist, value=TRUE), grep("_STG_", ROIlist, value=TRUE), grep("_TTG_", ROIlist, value=TRUE))
		Parietal_vec <- c(grep("_PCu_", ROIlist, value=TRUE), grep("_PoG_", ROIlist, value=TRUE), grep("_AnG_", ROIlist, value=TRUE), grep("_PO_", ROIlist, value=TRUE), grep("_SPL_", ROIlist, value=TRUE), grep("_MPrG_", ROIlist, value=TRUE), grep("_SMG_", ROIlist, value=TRUE), grep("_MPoG_", ROIlist, value=TRUE))
		Occipital_vec <- c(grep("_IOG_", ROIlist, value=TRUE), grep("_Cun_", ROIlist, value=TRUE), grep("_LiG_", ROIlist, value=TRUE), grep("_OFuG_", ROIlist, value=TRUE), grep("_MOG_", ROIlist, value=TRUE), grep("_Calc_", ROIlist, value=TRUE), grep("_OCP_", ROIlist, value=TRUE), grep("_SOG_", ROIlist, value=TRUE))
		Cerebellum_vec <- c(grep("_Cerebellar_Exterior_", ROIlist, value=TRUE), grep("_Cerebellar_Vermal_Lobules_I.V_", ROIlist, value=TRUE), grep("_Cerebellar_Vermal_Lobules_VI.VII_", ROIlist, value=TRUE), grep("_Cerebellar_Vermal_Lobules_VIII.X_", ROIlist, value=TRUE))
		ROILobes_list <- list(BasGang_vec, Limbic_vec, FrontOrb_vec, FrontDors_vec, Temporal_vec, Parietal_vec, Occipital_vec, Cerebellum_vec)
		names(ROILobes_list) <- list("BasGang", "Limbic", "FrontOrb", "FrontDors", "Temporal", "Parietal", "Occipital", "Cerebellum")
		return(ROILobes_list) 
	} else {
		Limbic_vec <- c(grep("Limbic", ROIlist, value=TRUE))
		Insular_vec <- c(grep("Insular", ROIlist, value=TRUE))
		Frontal_vec <- c(grep("Frontal", ROIlist, value=TRUE))
		Parietal_vec <- c(grep("Parietal", ROIlist, value=TRUE))
		Occipital_vec <- c(grep("Occipital", ROIlist, value=TRUE))
		Temporal_vec <- c(grep("Temporal", ROIlist, value=TRUE))
		ROILobes_list <- list(Limbic_vec, Insular_vec, Frontal_vec, Parietal_vec, Occipital_vec, Temporal_vec)
		names(ROILobes_list) <- list("Limbic", "Insular", "Frontal", "Parietal", "Occipital", "Temporal")
		if (length(ROILobes_list) < 6) {
			print("WARNING! You provided fewer than 6 lobes. Please check if all of the lobes you provided are present in the output of this function.")
		}
		return(ROILobes_list)
	}
}


## Declare a function to create a dataframe that can be used to plot ROI parameters
# factorsList should be names of columns from dataframe that you want in summaryDF, as a vector
# ROI_ListofLists takes output of roiLobes()
createSummaryDF <- function(dataframe, factorsList, ROI_ListofLists, stats=TRUE, simpleNames=TRUE, pattern1="vol_miccai_ave_", pattern2="_zscore", replacement="", ROIlist=c()) { #ERB: ROIlist might be a problem... unused argument
	if (length(factorsList) > 0) {
		# find the number of levels for each factor variable
		nameLevels <- list()
		for (i in 1:length(factorsList)) {
			dataframe[,factorsList[[i]]] <- as.factor(dataframe[,factorsList[[i]]])
			nameLevels <- c(nameLevels, list(levels(dataframe[,factorsList[[i]]])))
		}
		product_numLevels = 1
		for (i in 1:length(nameLevels)) {
			product_numLevels = product_numLevels*length(nameLevels[[i]])
		}
	} else {
		product_numLevels = 1
	}
	# find the number of ROI names and create a list of ROI names
	if (length(ROI_ListofLists) != 0) {
		namesROI_vec <- c()
		for (i in 1:length(ROI_ListofLists)) {
			if ((length(ROI_ListofLists[[i]]) != 0) && !(is.na(ROI_ListofLists[[i]]))) {
				namesROI_vec <- c(namesROI_vec, ROI_ListofLists[[i]])
			}
		}
	} else {
		namesROI_vec <- ROIlist
	}
	numROI= length(namesROI_vec)
	# name the columns of the matrix
	if (stats == TRUE) {
		# create a blank dataframe with the appropriate dimensions
		summaryDF <- data.frame(matrix(nrow = product_numLevels*numROI, ncol = length(factorsList) + 4))
		colnames(summaryDF) <- c("ROI_name", factorsList, "Lobe", "ROI_mean", "ROI_se")
	} else if (stats == FALSE) { #####################
		summaryDF <- data.frame(matrix(nrow = product_numLevels*numROI, ncol = length(factorsList) + 3))
		colnames(summaryDF) <- c("ROI_name", factorsList, "Lobe", "Subj_ZScore")
	}
	# ----repeat the factor levels to have a fully fleshed-out dataframe----
	ROIcolumn_vec <- c()
	for (i in 1:product_numLevels) {
		ROIcolumn_vec <- c(ROIcolumn_vec, namesROI_vec)
	}
	summaryDF$ROI_name <- ROIcolumn_vec
	if (length(factorsList) > 0) {
		# put levels into summaryDF
		num_repLevels = 1
		num_rows = nrow(summaryDF)
		# loop through all of the factors that will become columns in summaryDF
		for (i in 1:length(factorsList)) {
			if (i > 1) {
				num_repLevels = num_repLevels*length(nameLevels[[i-1]]) # number of times a complete chunk of levels will get printed
			}
			num_levelPerIter = num_rows/(length(nameLevels[[i]]))
			if (i == 1) {
				length = length(nameLevels[[i]])
			}
			if (i > 1) {
				length = length*length(nameLevels[[i]])
				num_levelPerIter = num_rows/length
			}
			# loop through the chunks
			for (j in 1:num_repLevels) {
				newcol_vec <- c()
				for (k in 1:length(nameLevels[[i]])) {
					for (p in 1:num_levelPerIter) {
						newcol_vec <- c(newcol_vec, nameLevels[[i]][k])
					}
				}
				summaryDF[ , i+1] <- c(newcol_vec)
			}
		}
	} else {
		num_rows = nrow(summaryDF)
	}
	# put in lobes
	if (length(ROI_ListofLists) != 0) {
		for (i in 1:num_rows) {
			roiName = summaryDF[i,1] #name of an roi 
			for (j in 1:length(ROI_ListofLists)) { #retrieve lobe for roiName
				if (roiName %in% ROI_ListofLists[[j]]) {
					lobe = names(ROI_ListofLists[j])
				}
			}
			summaryDF[i,"Lobe"] <- lobe
		}
	} else {
		summaryDF$Lobe <- NULL
	}
	# calculate the means and SEs #ERB: This works, but isn't running in the context of the function...
	if (stats == TRUE) {
		for (i in 1:length(namesROI_vec)) {
			summary_stats <- summarySE(dataframe, measurevar=namesROI_vec[[i]], groupvars=factorsList, na.rm=T)
			# find row numbers for this ROI
			ROI_rownums <- c()
			for (j in 1:num_rows) {
				if (summaryDF[j,1] == namesROI_vec[[i]]) {
					ROI_rownums <- c(ROI_rownums, j)
				}
			}
			# put the stats from summarySE in (ERB: will the order always match?)
			for (k in 1:length(ROI_rownums)) {
				summaryDF[ROI_rownums[[k]], "ROI_mean"] <- summary_stats[k, namesROI_vec[[i]]]
				summaryDF[ROI_rownums[[k]], "ROI_se"] <- summary_stats[k, "se"]
			}
		}	
		if (simpleNames == TRUE) {
			ROIcolumn_vec <- gsub(pattern1, replacement, ROIcolumn_vec)
			ROIcolumn_vec <- gsub(pattern2, replacement, ROIcolumn_vec)
		}
		summaryDF$ROI_name <- ROIcolumn_vec
	} else {
		if (simpleNames == TRUE) {
			ROIcolumn_vec <- gsub(pattern1, replacement, ROIcolumn_vec)
			ROIcolumn_vec <- gsub(pattern2, replacement, ROIcolumn_vec)
		}
		summaryDF$ROI_name <- ROIcolumn_vec
	}
	if (length(ROI_ListofLists) != 0) {
		# change the data type for Lobe to factor
		summaryDF$Lobe <- as.factor(summaryDF$Lobe)
		# create lobeLevels_vec
		lobeLevels_vec <- c()
		for (i in 1:length(summaryDF$Lobe)) {
			if (!(as.character(summaryDF[i, "Lobe"]) %in% lobeLevels_vec)) {
				lobeLevels_vec <- c(lobeLevels_vec, as.character(summaryDF[i, "Lobe"]))
			}
		}
		summaryDF$Lobe <- factor(summaryDF$Lobe, levels=lobeLevels_vec) ###########
	}
	# change the data type for ROI_name to factor
	summaryDF$ROI_name <- as.factor(summaryDF$ROI_name)
	# create roiLevels_vec
	roiLevels_vec <- c()
	i=1
	while (!(as.character(summaryDF[i,"ROI_name"]) %in% roiLevels_vec)) {
		roiLevels_vec <- c(roiLevels_vec, as.character(summaryDF[i,"ROI_name"]))
		i = i + 1
	}
	summaryDF$ROI_name <- factor(summaryDF$ROI_name, levels=roiLevels_vec)
	return(summaryDF)
}


## Declare a function to plot ROIs, split by some factor
# "dataframe" should be the output of createSummaryDF, but filtered so that only one factor with multiple levels remains
# "factor" should be the factor for which you want different colored lines (factor from above), data type = character
# "plotTitle" should describe the levels that you kept of the other factors from dataframe and how the z-scores were calculated (e.g., what was regressed out)
createGGPlotImage <- function(dataframe, factor = "", plotTitle, lower_order=-1, upper_order=1, increment=.2, rois=TRUE) {
	# retrieve the levels of factor
	if (factor != "") {
		levels_vec <- c()
		nrows = nrow(dataframe)
		for (i in 1:nrows) {
			if (!(dataframe[i, factor] %in% levels_vec)) {
				levels_vec <- c(levels_vec, dataframe[i, factor])
			}
		}
	} else {
		levels_vec <- c("no levels!")
	}
	# create the vector for scale_colour_manual
	if (length(levels_vec) == 1) {
		color_vec <- c("brown1")
	} else if (length(levels_vec) == 2) {
		color_vec <- c("brown1", "dodgerblue")
	} else if (length(levels_vec) == 3) {
		color_vec <- c("brown1", "dodgerblue", "darkseagreen1")
	} else if (length(levels_vec) == 4) {
		color_vec <- c("brown1", "dodgerblue", "darkseagreen1", "darkcyan")
	} else if (length(levels_vec) == 5) {
		color_vec <- c("brown1", "dodgerblue", "darkseagreen1", "darkcyan", "lightskyblue1")
	}
	# create the vector for scale_linetype_manual
	if (length(levels_vec) == 1) {
		linetype_vec <- c("solid")
	} else if (length(levels_vec) == 2) {
		linetype_vec <- c("solid", "solid")
	} else if (length(levels_vec) == 3) {
		linetype_vec <- c("solid", "solid", "solid")
	} else if (length(levels_vec) == 4) {
		linetype_vec <- c("solid", "solid", "solid", "solid")
	} else if (length(levels_vec) == 5) {
		linetype_vec <- c("solid", "solid", "solid", "solid", "solid")
	}
	# create the plot
	if (factor !=  "") {
		plotToReturn <- ggplot(dataframe, aes_string(y="ROI_mean", x="ROI_name", group=factor)) + 
			geom_line(aes_string(linetype=factor, color=factor), size=5) +
			geom_point(aes_string(shape=factor, color=factor), size=1.5) +
			scale_y_continuous(limits=c(lower_order, upper_order), breaks=round(seq(lower_order,upper_order,increment), digits=2)) +
			xlab("ROI") +
			ylab("Z-Score") +
			geom_hline(aes(yintercept=0), linetype="longdash", colour="black", size=0.5) +
			scale_colour_manual(name = factor, values = color_vec) +
			scale_linetype_manual(name = factor, values = linetype_vec) +
			theme_bw() +
			theme(legend.position="top") +
			#facet_grid(cols = vars(Lobe), scales="free", space="free_x") +
			ggtitle(plotTitle)
			if (rois == TRUE) {
				plotToReturn = plotToReturn + facet_grid(. ~ Lobe, scales="free", space="free_x") + 
				theme(text=element_text(size=20), axis.text.x = element_text(angle = 45, hjust = 1, face="bold"), axis.text.y = element_text(face="bold", size=12), axis.title.x = element_text(face="bold", size=25),
			axis.title.y = element_text(face="bold", size=25),
			plot.title = element_text(face="bold", size=25), strip.text.x = element_text(size=12))
			} else {
				plotToReturn = plotToReturn + theme(text=element_text(size=40), axis.text.x = element_text(angle = 45, hjust = 1, face="bold"), axis.text.y = element_text(face="bold", size=12), axis.title.x = element_text(face="bold", size=25),
			axis.title.y = element_text(face="bold", size=25),
			plot.title = element_text(face="bold", size=25), strip.text.x = element_text(size=12))
			}
	} else {
		plotToReturn <- ggplot(dataframe, aes(y=Subj_ZScore, x=ROI_name, group=Lobe)) + 
			geom_point(size=1.5) +
			geom_line() +
			scale_y_continuous(limits=c(lower_order, upper_order), breaks=round(seq(lower_order,upper_order,increment), digits=2)) +
			xlab("ROI") +
			ylab("Z-Score") +
			geom_hline(aes(yintercept=0), linetype="longdash", colour="black", size=0.5) +
			theme_bw() +
			theme(legend.position="top") +
			ggtitle(plotTitle)
			if (rois == TRUE) {
				plotToReturn = plotToReturn + facet_grid(. ~ Lobe, scales="free", space="free_x") + 
				theme(text=element_text(size=20), axis.text.x = element_text(angle = 45, hjust = 1, face="bold"), axis.text.y = element_text(face="bold", size=24), axis.title.x = element_text(face="bold", size=28),
			axis.title.y = element_text(face="bold", size=28),
			plot.title = element_text(face="bold", size=28), strip.text.x = element_text(size=15))
			} else {
				plotToReturn = plotToReturn + theme(text=element_text(size=40), axis.text.x = element_text(angle = 45, hjust = 1, face="bold"), axis.text.y = element_text(face="bold", size=24), axis.title.x = element_text(face="bold", size=28),
			axis.title.y = element_text(face="bold", size=28),
			plot.title = element_text(face="bold", size=28), strip.text.x = element_text(size=15))
			}
	}
	return(plotToReturn)
}


## Declare a function to create a summaryDF for subject against a normative group
# "subjectID" should be the string that represents the subject you want to compare to the rest of the sample
# "timePoint" should be 
# "columnSubjIDs" should be the column name with the subject IDs in it
# "columnTimePoint"
# "dataframe" should only contain column that specify the subject and the DVs
subjectCompare <- function(subjectID, timePoint, columnSubjID, columnTimePoint, DVColumnNums, factorsList, ROI_ListofLists, dataframe, pattern1, pattern2, simpleNames=TRUE) {
	# create a column that identifies the subject
	dataframe$subjYes <- NA
	num_rows <- nrow(dataframe)
	subj_rownum <- which(dataframe[,columnSubjID]==subjectID & dataframe[,columnTimePoint]==timePoint)
	for (i in 1:num_rows) {
		if (i != subj_rownum) {
			dataframe[i, "subjYes"] <- "nope"
		} else {
			dataframe[i, "subjYes"] <- "YESSS"
		}
	}
	# create summaryDF
	summaryDF <- createSummaryDF(data.frame(), factorsList, ROI_ListofLists, stats=FALSE, simpleNames, pattern1, pattern2)
	# create z-scores for this subject
	zscore_norm_data <- dataframe[dataframe$subjYes == "nope", ]
	# note: might want to use variables that end in "_regressed"
	subj_zscores <- data.frame(nrow=length(DVColumnNums), ncol=2)
	colnames(subj_zscores) <- c("ROI_name", "Subj_ZScore")
	for (i in 1:length(DVColumnNums)) {
		mean = mean(zscore_norm_data[ , DVColumnNums[[i]]], na.rm = TRUE)
		sd = sd(zscore_norm_data[ , DVColumnNums[[i]]], na.rm = TRUE)
		zscore = (dataframe[dataframe$subjYes == "YESSS", DVColumnNums[[i]]] - mean)/sd
		subj_zscores[i, "ROI_name"] <- colnames(zscore_norm_data[DVColumnNums[[i]]])
		subj_zscores[i, "Subj_ZScore"] <- zscore
	}
	# give subj_zscores simple names
	if (simpleNames == TRUE) {
		ROIcolumn_vec <- gsub(pattern1, replacement, ROIcolumn_vec)
		ROIcolumn_vec <- gsub(pattern2, replacement, ROIcolumn_vec)
	}
	subj_zscores$ROI_name <- ROIcolumn_vec
	# sort subj_zscores such that it's rows are in the same order as summaryDF
	target <- summaryDF[,"ROI_name"]
	subj_zscores <- subj_zscores[match(target, subj_zscores$ROI_name),]
	# put these z-scores in a summaryDF
	summaryDF[ , "Subj_ZScore"] <- subj_zscores[,"Subj_ZScore"]
	return(summaryDF)
}


## Declare a function that filters for subjects with all timepoints
# subjColumn should be a string
# timeColumn should be a string
filterAllTimepoints <- function(dataframe, subjColumn, timeColumn, numTimepoints) {
	numrows <- nrow(dataframe)
	timepointsPerSubject <- list()
	# iterate through all of the rows of dataframe, creating a vector of two-element vectors (subj, numTimepoints = 0)
	names_vec <- c()
	for (i in 1:numrows) {
		subj = as.character(dataframe[i, subjColumn])
		if (!(subj %in% names_vec)) {
			names_vec <- c(names_vec, subj)
		}
	}
	timepointsPerSubject <- as.list(matrix(0,nrow=length(names_vec)))
	names(timepointsPerSubject) <- names_vec
	for (i in 1:numrows) {
		subj = as.character(dataframe[i, subjColumn])
		num_element <- which(names(timepointsPerSubject) == subj)
		timepointsPerSubject[[num_element]] = timepointsPerSubject[[num_element]] + 1
	}
	# find the subjects that have all of the timepoints
	allTimepoints <- c()
	for (j in 1:length(timepointsPerSubject)) {
		if (timepointsPerSubject[[j]] ==  numTimepoints) {
			allTimepoints <- c(allTimepoints, names_vec[[j]])
		}
	}
	copies_allTimepoints <- c()
	for (k in 1:numTimepoints) {
		copies_allTimepoints = c(copies_allTimepoints, allTimepoints)
	}
	dataframe <- subset(dataframe, dataframe[ , subjColumn] %in% allTimepoints)
	return(dataframe)
}

## Declare a function that filters for subjects with specific timepoints
# subjColumn should be a string
# timeColumn should be a string
# timePointVec should be a vector of strings that are the timepoints you would like to select for
#filterForTimepoints <- function(dataframe, subjColumn, timeColumn, timePointVec) {
#	numrows <- nrow(dataframe)
#	timepointsPerSubject <- list()
#	# iterate through all of the rows of dataframe, creating a vector of subject names
#	names_vec <- c()
#	for (i in 1:numrows) {
#		subj = as.character(dataframe[i, subjColumn])
#		if (!(subj %in% names_vec)) {
#			names_vec <- c(names_vec, subj)
#		}
#	}
#	timepointDF <- data.frame(matrix(NA, nrow=length(names_vec), ncol=1+length(timePointVec)))
#	colnames(timepointDF) <- c("subject", timePointVec)
#	timepointDF$subject <- names_vec
#	# iterate over all of the rows of dataframe check off if that data exists
#	for (i in 1:numrows) {
#		subj <- as.character(dataframe[i, subjColumn])
#		subj
#		time <- as.character(dataframe[i, timeColumn])
#		time
#		timepointDF[timepointDF$subject == subj, time] <- 1
#	}
#	timepointDF <- timepointDF[,c("subject", timePointVec)]
#	# create a vector to filter by
#	filter_vec <- c()
#	for (i in 1:nrow(timepointDF)) {
#		row <- as.vector(timepointDF[i, ])
#		names(row) <- NULL
#		row <- as.vector(row)
#		if (!(NA %in% row)) {
#			subj <- timepointDF[i, "subject"]
#			filter_vec <- c(filter_vec, subj)
#		}
#	}
#	dataframe <- dataframe[ , which(dataframe[ , subjColumn]         )] #ERB: Pick up here
#	return(dataframe)
#}


## Declare a function to "plot mean and standard error in Boxplot in R"
# input to fun.data in stat_summary
MinMeanSEMMax <- function(x) {
  v <- c(min(x), mean(x) - sd(x)/sqrt(length(x)), mean(x), mean(x) + sd(x)/sqrt(length(x)), max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}


## Declare a function to create a column of scanning sites
scanningSite_NASAAntartica <- function(dataframe, wintercol="winterover", subject_1col="subject_1", Timecol="Time") {
	dataframe$scanner <- NA
	num_rows <- nrow(dataframe)
	for (row in 1:num_rows) {
		if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_001") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_001") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_002") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_003") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_004") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_004") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"		
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_004") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_005") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_005") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_005") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_006") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_007") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_007") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_007") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_008") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_008") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_008") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_009") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_009") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_009") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_010") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_011") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_011") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_011") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_012") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_012") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_013") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_013") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "DLR_013") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_101") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_101") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_102") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_102") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_102") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_103") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_103") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_103") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_104") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_104") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_105") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_105") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_105") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_106") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_106") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_107") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_107") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_107") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_108") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_108") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_108") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_110") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_110") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_111") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_111") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_111") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_112") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_112") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_112") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_113") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_113") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "DLR_113") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_001") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_001") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "HOB"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_001") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_002") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_002") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "HOB"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_002") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_003") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_003") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "HOB"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_003") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_004") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_004") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CHR"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_004") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_005") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_005") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "HOB"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_005") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_006") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_006") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CHR"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_006") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_007") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CHR"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_007") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CHR"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_007") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_008") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_008") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "HOB"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_008") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_009") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_009") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "HOB"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_009") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_010") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "HOB"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_010") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CHR"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_010") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_011") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_011") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "HOB"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_011") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_012") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_012") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "HOB"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_012") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_013") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_013") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "HOB"
		} else if ((dataframe[row, wintercol] == "wo_2015") && (dataframe[row, subject_1col] == "concordia_013") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_101") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_101") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CHR"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_102") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_102") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "HOB"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_103") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_103") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CHR"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_104") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_104") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "HOB"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_104") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_105") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_105") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CHR"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_105") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_106") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_107") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "HOB"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_107") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CHR"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_108") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_108") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CHR"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_108") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_109") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "HOB"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_109") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_110") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_110") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "HOB"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_111") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_111") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "HOB"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_111") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_112") && (dataframe[row, Timecol] == "t0")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_112") && (dataframe[row, Timecol] == "t12")) {
			dataframe[row, "scanner"] <- "CHR"
		} else if ((dataframe[row, wintercol] == "wo_2016") && (dataframe[row, subject_1col] == "concordia_112") && (dataframe[row, Timecol] == "t18")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "phantoms") && (dataframe[row, subject_1col] == "BJ") && (dataframe[row, Timecol] == "t1")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "phantoms") && (dataframe[row, subject_1col] == "BJ") && (dataframe[row, Timecol] == "t4")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "phantoms") && (dataframe[row, subject_1col] == "BM") && (dataframe[row, Timecol] == "t1")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "phantoms") && (dataframe[row, subject_1col] == "BM") && (dataframe[row, Timecol] == "t2")) {
			dataframe[row, "scanner"] <- "CHR"
		} else if ((dataframe[row, wintercol] == "phantoms") && (dataframe[row, subject_1col] == "BM") && (dataframe[row, Timecol] == "t3")) {
			dataframe[row, "scanner"] <- "HOB"
		} else if ((dataframe[row, wintercol] == "phantoms") && (dataframe[row, subject_1col] == "BM") && (dataframe[row, Timecol] == "t4")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "phantoms") && (dataframe[row, subject_1col] == "EA") && (dataframe[row, Timecol] == "t1")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "phantoms") && (dataframe[row, subject_1col] == "EA") && (dataframe[row, Timecol] == "t2")) {
			dataframe[row, "scanner"] <- "CHR"
		} else if ((dataframe[row, wintercol] == "phantoms") && (dataframe[row, subject_1col] == "EA") && (dataframe[row, Timecol] == "t3")) {
			dataframe[row, "scanner"] <- "HOB"
		} else if ((dataframe[row, wintercol] == "phantoms") && (dataframe[row, subject_1col] == "EA") && (dataframe[row, Timecol] == "t4")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "phantoms") && (dataframe[row, subject_1col] == "GR") && (dataframe[row, Timecol] == "t1")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "phantoms") && (dataframe[row, subject_1col] == "GR") && (dataframe[row, Timecol] == "t2")) {
			dataframe[row, "scanner"] <- "CHR"
		} else if ((dataframe[row, wintercol] == "phantoms") && (dataframe[row, subject_1col] == "GR") && (dataframe[row, Timecol] == "t3")) {
			dataframe[row, "scanner"] <- "HOB"
		} else if ((dataframe[row, wintercol] == "phantoms") && (dataframe[row, subject_1col] == "GR") && (dataframe[row, Timecol] == "t4")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "phantoms") && (dataframe[row, subject_1col] == "PK") && (dataframe[row, Timecol] == "t1")) {
			dataframe[row, "scanner"] <- "CGN"
		} else if ((dataframe[row, wintercol] == "phantoms") && (dataframe[row, subject_1col] == "PK") && (dataframe[row, Timecol] == "t2")) {
			dataframe[row, "scanner"] <- "CHR"
		} else if ((dataframe[row, wintercol] == "phantoms") && (dataframe[row, subject_1col] == "PK") && (dataframe[row, Timecol] == "t4")) {
			dataframe[row, "scanner"] <- "CGN"
		} 
	}
	dataframe$scanner <- as.factor(dataframe$scanner)
	return(dataframe)
}







