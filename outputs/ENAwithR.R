## rm(list=ls(all=TRUE))

## loading libraries required to run the ENA scripts

## library("XLConnect")

##
##

##################################################################################################
###																							   ###
### Function that imports an Excel file and extracts three vectors and the matrix of exchanges ###
###																							   ###
##################################################################################################

##	fromEXCEL <- function(filename, num){
##		##	filename <- "amb_HW_ENA.xlsx"
##		##	num <- 1
##		theData <- readWorksheet(loadWorkbook(filename), sheet = num)
##		##
##		input <- theData$Z
##		export <- theData$E
##		respiration <- theData$R
##		exchanges <- theData[,-c(1:3)]
##		##
##		nomi <- colnames(exchanges)
##		rownames(exchanges) <- nomi
##		##
##		names(input) <- names(export) <- names(respiration) <- nomi
##		##
##		list_ret <- as.list(rep(NA,4))
##		list_ret[[1]] <- input
##		list_ret[[2]] <- exchanges
##		list_ret[[3]] <- export
##		list_ret[[4]] <- respiration
##		names(list_ret) <- c("Z", "T", "E", "R")
##		return(list_ret)
##	}

##
##

######################################################################################################
###												   												   ###
### Functions to convert dataframe input data (in a similar SCOR format) into vectors and matrices ###
###												   												   ###
######################################################################################################

## Conversion of SCOR vectors (available as input dataframes) into whole vectors. It is useful for (B), (Z), (E) and (R).

vec.ena <- function(df,x){
	V <- rep(0,x)
	for(i in 1:nrow(df))V[df[i,1]] <- df[i,2]
	return(V)
}


## Conversion of SCOR matrices (available as input dataframes) into whole matrix. 
## Mainly useful for the matrix of inter-compartmental[T] exchanges.

matr.ena <- function(df,x){
	tot_el <- x^2 
	M <- matrix(rep(0,tot_el), nrow=x)
	for(i in 1:nrow(df))M[df[i,1],df[i,2]] <- df[i,3]
	return(M)
}


## Conversion of the dataframe input (.csv) file into vectors and matrices for ENA.
## Filname should be "quoted". Empty cells are converted into zeros.

df.ena <- function(filename){
	ch <- nchar(filename)
	n.ch <- ch - 2
	nch.frag <- substring(filename,n.ch)
	#
	{
	if(nch.frag == "txt")AA <- read.table(filename,header=TRUE)
	#
	else AA <- read.csv2(filename,header=TRUE)
	}
	BB <- as.matrix(AA)
	#
	if(length(which(is.na(AA)==TRUE))!=0){
		aa <- which(is.na(AA)==TRUE)
		BB[aa] <- 0
		}
	#			
	# We create a list including information on: 1) total number of nodes; 2) number of living nodes;
	# 3) number of non-living nodes; 4) name of the nodes; 5) biomasses; 6) imports; 7) exports;
	# 8) respirations; 9) intercompartmental exchanges; 10) exchanges between living compartments
	#
	BB.L <- as.list(rep(NA,10))
	#
	BB.L[[1]] <- nrow(BB)
	tt <- BB.L[[1]]
	BB.L[[4]] <- t(t(colnames(BB)[6:(5+tt)]))
	{
	if(length(which(BB[,1]>0))!=0){
		BB.L[[2]] <- length(which(BB[,1]>0))
		}
	#
	else BB.L[[2]] <- 0
	}
	BB.L[[5]] <- BB[,2]
	BB.L[[6]] <- BB[,3]
	BB.L[[7]] <- BB[,4]
	BB.L[[8]] <- BB[,5]
	BB.L[[9]] <- BB[c(1:tt),c(6:(5+tt))]
	#
	{
	if(length(which(BB[,1]<0))!=0){
		BB.L[[3]] <- length(which(BB[,1]<0))
		neg <- which(BB[,1]<0)
		#
		BB.L[[10]] <- BB.L[[9]][-c(neg),-c(neg)]
		}
	#
	else {
		BB.L[[3]] <- 0
		BB.L[[10]] <- BB.L[[9]]
		}
	}
	#
	return(BB.L)
}

##
##

########################################################################################################################################
###																    																 ###
### Functions to depict all the flows (import, export, dissipations and inter-compartmentmental exchanges) as one big matrix [T_tot] ###
### In this case the input considered are vector and matrices (not dataframes)							     						 ###
### In the first function ("bigmatr.ena") the output will also show the name of each compartment, while this addional information    ###
### is not available in the case of the raw result of "bigmatrNOnames.ena"							    							 ###
###																    																 ###
########################################################################################################################################

bigmatr.ena <- function(cnames,inp,inter,outt,diss){
	tot_row <- length(inp)
	zeros_up <- rep(0, tot_row)
	core <- cbind(zeros_up, inter, outt, diss)
	inputs <- c(0, inp, 0, 0)
	zeros_down <- rep(0, tot_row + 3)
	final <- rbind(inputs, core, zeros_down, zeros_down)
	#
	names_tot <- c("Z", cnames, "E", "R")
	rownames(final) <- colnames(final) <- names_tot
	return(final)
}


bigmatrNOnames.ena <- function(inp,inter,outt,diss){
	tot_row <- length(inp)
	zeros_up <- rep(0, tot_row)
	core <- cbind(zeros_up, inter, outt, diss)
	inputs <- c(0, inp, 0, 0)
	zeros_down <- rep(0, tot_row+3)
	final <- rbind(inputs, core, zeros_down, zeros_down)
	colnames(final) <- NULL
	rownames(final) <- NULL
	#
	return(final)
}

##
##

#####################################################
###						 						  ###
### Balancing procedure with different algorithms ###
###						 						  ###
#####################################################

## "ETM" function represents an intermediate step for performing all the balancing procedure listed below

ETM <- function(inp,inter,outt,diss){
	tot_row <- length(inp)
	tot_el <- tot_row^2
	zeros_side <- rep(0, tot_row)
	core <- cbind(inter, outt, diss, zeros_side)
	#
	inputs <- c(inp, 0, 0, 0)
	zeros_down <- rep(0, tot_row+3)
	#
	final <- rbind(core, zeros_down, zeros_down, inputs)
	#
	return(final)
}


## Input Based Approach

IN.balance <- function(inp,inter,outt,diss){
			##
			## data input as: import vector (inp); intercompartmental exchange matrix [inter]; export (outt) 
			## and respiration vectors (diss); total number of nodes x
			## STEP 1
			##
			x <- nrow(inter)
			matrice.step1 <- ETM(inp,inter,outt,diss)
			INP <- apply(matrice.step1,2,sum)
			OUT <- apply(matrice.step1,1,sum)
			#
			DIFF <- rep(0,x)
			N_INP <- rep(0,x)
			N_OUT <- rep(0,x)
			#
			for(i in 1:x){
				N_INP[i] <- INP[i]
				N_OUT[i] <- OUT[i]
				DIFF[i] <- INP[i] - OUT[i]
				}
			#
			DIFF_IN <- (DIFF*100)/N_INP
			diff.df <- data.frame(cbind(N_INP,N_OUT,DIFF,DIFF_IN))
			##
			## STEP 2 - Partial Host Matrix [F]
			##
			ele <- x + 3
			F <- matrix(rep(0,ele^2),nrow=ele)
			#
			for(i in 1:ele)
				for(j in 1:ele){
				if(OUT[i]!=0) F[i,j] <- matrice.step1[i,j]/OUT[i]
				}
			##
			## STEP 3 - Matrix [R] = t[cF] - [I]
			## Transpose the N x N part of matrix F and substract the identity matrix to get matrix [R]
			##
			cF <- matrix(rep(0,x^2),nrow=x)
			#
			for(i in 1:x)
				for(j in 1:x){
				cF[i,j] <- F[i,j]
				}
			#
			I <- diag(1,x)
			R <- t(cF) - I
			##
			## STEP 4 - Inversion of matrix [R]
			##
			ddd <- det(R)
			if(ddd == 0)return(NA)
			RI <- solve(R)
			##
			## STEP 5 - Matrix [R] manipulation
			## Multiply every r_ij coefficient by the corresponding jth input in the [ETM] matrix, and change its sign
			## r_ij <- -r_ij * (t_N+1,j + t_N+2,j + t_N+3,j)
			##
			R.man <- matrix(rep(0,x^2),nrow=x)
			#
			for(i in 1:x)
				for(j in 1:x){
					R.man[i,j] <- (-RI[i,j]*inp[j])
					}
			#
			##
			## STEP 6 - Vector (U)
			## Sum i's row of the [R.man] matrix to build vector U(u_i)
			##
			U <- apply(R.man,1,sum)
			##
			## STEP 7 - [T.balIN]
			## Multiply each f_ij coefficient by the corresponding u_i, to obtain a balanced form [T.balIN] of the [inter] matrix
			##
			T.balIN <- matrix(rep(0,ele^2),nrow=ele)
			#
			for(i in 1:x)
				for(j in 1:ele){
					T.balIN[i,j] <- F[i,j] * U[i]
					}
			
			for(i in (x+1):ele)
				for(j in 1:ele){
						T.balIN[i,j] <- matrice.step1[i,j]
						}
			#
			return(T.balIN)
			}


## Output Based Approach

OUT.balance <- function(inp,inter,outt,diss){
				##
				## data input as: import vector (inp); intercompartmental exchange matrix [inter]; export (outt) 
				## and respiration vectors (diss); total number of nodes x
				## STEP 1
				##
				x <- nrow(inter)
				matrice.step1 <- ETM(inp,inter,outt,diss)
				matrice.step1 <- t(matrice.step1)
				#
				INP <- apply(matrice.step1,2,sum)
				OUT <- apply(matrice.step1,1,sum)
				##
				## STEP 2 - Partial Host Matrix [F]
				##
				ele <- x + 3
				F <- matrix(rep(0,ele^2),nrow=ele)
				#
				for(i in 1:ele)
					for(j in 1:ele){
							if(sum(matrice.step1[i,])!=0) F[i,j] <- matrice.step1[i,j]/sum(matrice.step1[i,])
							}
				##
				## STEP 3 - Matrix [R] = t[cF] - [I]
				## Transpose the N x N part of matrix F and substract the identity matrix to get matrix [R]
				##
				cF <- matrix(rep(0,x^2),nrow=x)
				#
				for(i in 1:x)
					for(j in 1:x){
						cF[i,j] <- F[i,j]
						}
				#
				I <- diag(1,x)
				R <- t(cF) - I
				##
				## STEP 4 - Inversion of matrix [R]
				##
				ddd <- det(R)
				if(ddd == 0)return(NA)
				RI <- solve(R)
				##
				## STEP 5 - Matrix [R] manipulation
				## Multiply every r_ij coefficient by the corresponding jth input in the [ETM] matrix, and change its sign
				## r_ij <- -r_ij * (t_N+1,j + t_N+2,j + t_N+3,j)
				##
				R.man <- matrix(rep(0,x^2),nrow=x)
				ER <- outt + diss
				#
				for(i in 1:x)
					for(j in 1:x){
						R.man[i,j] <- (-RI[i,j]*ER[j])
						}
				##
				## STEP 6 - Vector (U)
				## Sum i's row of the [R.man] matrix to build vector U(u_i)
				##
				U <- apply(R.man,1,sum)
				##
				## STEP 7 - [T.balOUT]
				## Multiply each f_ij coefficient by the corresponding u_i, to obtain a balanced form [T.balIN]
				## of the [inter] matrix and finally transpose it to calculate [T.balOUT]
				##
				T.balIN <- matrix(rep(0,ele^2),nrow=ele)
				#
				for(i in 1:x)
					for(j in 1:ele){
						T.balIN[i,j] <- F[i,j] * U[i]
						}
				#				
				for(i in (x+1):ele)
					for(j in 1:ele){
							T.balIN[i,j] <- matrice.step1[i,j]
							}
				#
				T.balOUT <- t(T.balIN)
				#
				return(T.balOUT)
				}


## Average Approach (AVG derived algorithm)

AVG.balance <- function(inp,inter,outt,diss){
				##
				## it calculates average coefficients using the corresponding values obtained 
				## with input-based and output-based approach
				##
				if(length(IN.balance(inp,inter,outt,diss)) == 1)return(NA)
				if(length(OUT.balance(inp,inter,outt,diss)) == 1)return(NA)
				#
				matrix.one <- IN.balance(inp,inter,outt,diss)
				matrix.two <- OUT.balance(inp,inter,outt,diss)
				#
				matrix.sum <- matrix.one + matrix.two
				#
				matrix.end <- 0.5 * matrix.sum
				#
				return(matrix.end)
				}


## Input - Output Based Approach

IO.balance <- function(inp,inter,outt,diss){
				## 
				## starting from unbalanced matrix [inter], the input based approach is applied 
				## to (1/2)[inter], the result is summed with (1/2)[inter]. The resulting matrix is
				## completely balanced using the output-based algorithm
				##
				x <- nrow(inter)
				if(length(IN.balance(inp,inter,outt,diss)) == 1)return(NA)
				#
				matrix.one <- 0.5 * IN.balance(inp,inter,outt,diss)
				matrix.two <- 0.5 * (ETM(inp,inter,outt,diss))
				#
				matrix.interm <- matrix.one + matrix.two
				#
				scambi<-matrix(rep(0,x^2),nrow=x)
				inpor <- rep(0,x)
				expor <- rep(0,x)
				res <- rep(0,x)
				#
				for(i in 1:x)
					for(j in 1:x)scambi[i,j] <- matrix.interm[i,j]
				#
				for(i in 1:x){ 
					inpor[i] <- matrix.interm[(x+3),i]
					expor[i] <- matrix.interm[i,(x+1)]
					res[i] <- matrix.interm[i,(x+2)]
					}
				#
				if(length(OUT.balance(inp,inter,outt,diss)) == 1)return(NA)
				matrix.final <- OUT.balance(inpor,scambi,expor,res)
				#
				return(matrix.final)
				}


## Output - Input Based Approach

OI.balance <- function(inp,inter,outt,diss){
				## 
				## is the same of IO but first output-based approach is implemented, and then input-based
				##
				x <- nrow(inter)
				if(length(OUT.balance(inp,inter,outt,diss)) == 1)return(NA)
				#
				matrix.one <- 0.5 * OUT.balance(inp,inter,outt,diss)
				matrix.two <- 0.5 * (ETM(inp,inter,outt,diss))
				#
				matrix.interm <- matrix.one + matrix.two
				#
				scambi<-matrix(rep(0,x^2),nrow=x)
				inpor <- rep(0,x)
				expor <- rep(0,x)
				res <- rep(0,x)
				#
				for(i in 1:x)
					for(j in 1:x)scambi[i,j] <- matrix.interm[i,j]
				#
				for(i in 1:x){ 
					inpor[i] <- matrix.interm[(x+3),i]
					expor[i] <- matrix.interm[i,(x+1)]
					res[i] <- matrix.interm[i,(x+2)]
					}
				#
				if(length(IN.balance(inp,inter,outt,diss)) == 1)return(NA)
				matrix.final <- IN.balance(inpor,scambi,expor,res)
				#
				return(matrix.final)
				}


## Average-2 Approach (AVG2 derived algorithm)
## This is the same of AVG, but using "Input - Output" and "Output - Input" based approaches

AVG2.balance <- function(inp,inter,outt,diss){
				##
				## the same of AVG but with [T]_IO*(bal) and [T]_OI*(bal) matrices
				##
				if(length(IO.balance(inp,inter,outt,diss)) == 1)return(NA)
				if(length(OI.balance(inp,inter,outt,diss)) == 1)return(NA)
				#
				matrix.one <- IO.balance(inp,inter,outt,diss)
				matrix.two <- OI.balance(inp,inter,outt,diss)
				#
				matrix.sum <- matrix.one + matrix.two
				#
				matrix.end <- 0.5 * matrix.sum
				#
				return(matrix.end)
}


## This function takes as an output the bigmatrix generated by one of the various balancing algorithms
## and creates an Excel file to be used for preparing the input to WAND

BALtoWAND <- function(BM_bal, nome1){
	nome <- paste(nome1, ".xlsx", sep = "")
	N <- nrow(BM_bal) - 3
	Z <- BM_bal[(N+3),c(1:N)]
	T <- BM_bal[c(1:N),c(1:N)]
	E <- BM_bal[c(1:N),(N+1)]
	R <- BM_bal[c(1:N),(N+2)]
	##
	vectors <- cbind(Z,E,R)
	colnames(vectors) <- c("Z", "E", "R")
	##
	wb <- loadWorkbook(nome, create = TRUE)
	createSheet(wb, name = "vectors")
	createSheet(wb, name = "exchanges")
	writeWorksheet(wb, vectors, sheet = "vectors", startRow = 1, startCol = 1)
	writeWorksheet(wb, T, sheet = "exchanges", startRow = 1, startCol = 1)
	saveWorkbook(wb)
}

##
##

#####################################################################################
###										 										  ###
### Functions to calculate Partial Feeding Matrix [G] and Partial Host Matrix [F] ###
###										 										  ###
#####################################################################################

## The function "partial.feeding" convert data from (Z) and [T] into partial feeding matrix [G].

partial.feeding <- function(inp,inter){
			tot_row <- length(inp)
			tot_el <- tot_row^2
			G <- matrix(rep(0,tot_el), nrow = tot_row)
			IN <- rep(0,tot_row)	
			#
			for (j in 1:tot_row)IN[j] <- (sum(inter[,j]) + sum(inp[j]))
           		#
           			for (i in 1:tot_row)
               				for (j in 1:tot_row)if(sum(IN[j])!=0)G[i,j] <- inter[i,j]/sum(IN[j])
      			return(G)
			}


## The function "partial.host" convert data from [T], (E) and (R) into partial host matrix [F].

partial.host <- function(inter,outt,diss){
			tot_row <- length(outt)
			tot_el <- tot_row^2
			F <- matrix(rep(NA,tot_el), nrow = tot_row)
			OUT <- rep(0,tot_row)
			#
			for (i in 1:tot_row)OUT[i]<-(sum(inter[i,]) + sum(outt[i]) + sum(diss[i]))
          		#
          			for (i in 1:tot_row)
              				for (j in 1:tot_row)if(sum(OUT[i])!=0)F[i,j] <- inter[i,j]/sum(OUT[i])
         		return(F)
              		}





#######################################################################
###								    ###
### Functions to normalize vectors (that could be (Z), (E) and (R)) ###
###								    ###
#######################################################################

## The function to normalize input vector is called "normalize.invector" and it is used for (Z).

normalize.invector <- function(inp,inter){
					tot_row <- length(inp)
					NZ <- rep(0,tot_row)
					IN <- rep(0,tot_row)
					#
					for (j in 1:tot_row)IN[j] <- (sum(inter[,j]) + sum(inp[j]))
					#
					for (j in 1:tot_row)if(sum(IN[j])!=0)NZ[j] <- inp[j]/sum(IN[j])
					return(NZ)
					}

					
## The function to normalize both export and dissipation vectors is called "normalize.outvector" and it is used for (E) and (R).

normalize.outvector <- function(inter,outt,diss){
					tot_row <- length(outt)
					NO <- matrix(rep(0,tot_row*2), nrow = tot_row)
					OUT <- rep(0,tot_row)
					#
					for (i in 1:tot_row)OUT[i]<-(sum(inter[i,]) + sum(outt[i]) + sum(diss[i]))
					#
					for (i in 1:tot_row)if(sum(OUT[i])!=0){NO[i,1] <- outt[i]/sum(OUT[i])
									NO[i,2] <- diss[i]/sum(OUT[i])
									}
					colnames(NO) <- c("NE", "NR")
					return(NO)
					}





########################################################################################################
###												     ###
### Functions to calculate Inverse Matrix of Leontief [LEO] and Inverse Matrix of Augustinovic [AUG] ###
###												     ###
########################################################################################################

## The function "leontief" estimates [LEO] when (Z) and [T] are available.

leontief <- function(inp,inter){
				tot_row <- length(inp)
				G <- partial.feeding(inp,inter)
				#
				ddd <- det(diag(1,tot_row) - G)
				if(ddd == 0)return(NA)
				#
				LEO <- solve(diag(1,tot_row) - G)
				return(LEO)
				}


## The function "augustinovic" estimates [AUG] when [T], (E) and (Z) are available.

augustinovic <- function(inter,outt,diss){
					tot_row <- length(outt)
					F <- partial.host(inter,outt,diss)
					#
					ddd <- det(diag(1,tot_row) - t(F))
					if(ddd == 0)return(NA)
					#
					AUG <- solve(diag(1,tot_row) - t(F))
					return(AUG)
					}





#####################################################################################
###										  ###
### Functions to build Total Dependency [IG] and Total Contribution [IF] matrices ###
###										  ###
#####################################################################################

## The function "total.dependency" convert data from (Z), [T], (E) and (R) into [IG].

total.dependency <- function(inp,inter,outt,diss){
						tot_row <- length(inp)
						tot_el <- tot_row^2
						#
						IN <- rep(0,tot_row)	
							for (j in 1:tot_row)IN[j] <- (sum(inter[,j]) + inp[j])
						OUT <- rep(0,tot_row)
							for (i in 1:tot_row)OUT[i]<-(sum(inter[i,]) + outt[i] + diss[i])
						#
						inv.leontief <- leontief(inp,inter)
						if(length(inv.leontief) == 1)return(NA)
						#
						leo <- (inv.leontief - diag(1,tot_row))
						interm1 <- rep(NA,tot_row)
							for(i in 1:tot_row)interm1[i] <- IN[i]/(inv.leontief[i,i] * OUT[i])
						#
						DEP <- matrix(rep(NA,tot_el), nrow = tot_row)
							for (i in 1:tot_row)
								for(j in 1:tot_row)DEP[i,j] <- leo[i,j] * interm1[i]
						return(DEP)
						}


## The function "total.contribution" convert data from (Z), [T], (E) and (R) into [IF].

total.contribution <- function(inp,inter,outt,diss){
						tot_row <- length(inp)
						tot_el <- tot_row^2
						#
						IN <- rep(0,tot_row)	
							for (j in 1:tot_row)IN[j] <- (sum(inter[,j]) + inp[j])
						OUT <- rep(0,tot_row)
							for (i in 1:tot_row)OUT[i]<-(sum(inter[i,]) + outt[i] + diss[i])
						#
						inv.augustinovic <- augustinovic(inter,outt,diss)
						if(length(inv.augustinovic) ==1)return(NA)
						#
						aug <- (inv.augustinovic - diag(1,tot_row))
						interm1 <- rep(NA,tot_row)
							for(i in 1:tot_row)interm1[i] <- IN[i]/(inv.augustinovic[i,i] * OUT[i])
						#
						CONTR <- matrix(rep(NA,tot_el), nrow = tot_row)
							for (i in 1:tot_row)
								for(j in 1:tot_row)CONTR[i,j] <- aug[i,j] * interm1[i]
						return(t(CONTR))
						}





########################################################
###						     ###
### Functions to carry out INPUT and OUTPUT analyses ###
###						     ###
########################################################

## The function "input.analysis".

input.analysis <- function(inp,inter,outt,diss){
						tot_row <- length(inp)
						tot_el <- tot_row^2
						F <- partial.host(inter,outt,diss)
						NE <- normalize.outvector(inter,outt,diss)[,1]
						NR <- normalize.outvector(inter,outt,diss)[,2]
						#
						ddd <- det(diag(1,tot_row) - t(F))
						if(ddd == 0)return(NA)
						IFT <- solve(diag(1,tot_row) - t(F))						
						#
						# non-zero elements in the input vector (Z)
						#
						count<-0
							for (i in 1:tot_row)if(inp[i]!=0)count <- count+1
						#
						# list for each input different from zero (Z)_i != 0
						#
						final_input <- vector("list", count)
							for(i in 1:count){final_input[[i]] <- vector("list", 4)
									final_input[[i]][[1]] <- rep(0,tot_row)
									final_input[[i]][[2]] <- matrix(rep(0,tot_el), nrow = tot_row)
									final_input[[i]][[3]] <- rep(0,tot_row)
									final_input[[i]][[4]] <- rep(0,tot_row)
									}
						#
						step <- 1
						for(i in 1:tot_row)if(inp[i]!=0){final_input[[step]][[1]][i] <- 1
										step <- step + 1
										}
						#
						intermediate_step <- vector("list", count)
							for (i in 1:count)intermediate_step[[i]] <- IFT %*% final_input[[i]][[1]]
						#
						# defining normalized inter-compartmental exchanges, export and dissipations
						#
						for(k in 1:count){
							for(i in 1:tot_row){
								final_input[[k]][[3]][i] <- intermediate_step[[k]][i] * NE[i]
								final_input[[k]][[4]][i] <- intermediate_step[[k]][i] * NR[i]
								for(j in 1:tot_row)final_input[[k]][[2]][i,j] <- intermediate_step[[k]][i] * F[i,j]
								}
							k <- k+1
							}
						#
						return(final_input)			
						}


## The function "output.analysis".

output.analysis <- function(inp,inter,outt){
					tot_row <- length(inp)
					tot_el <- tot_row^2
					G <- partial.feeding(inp,inter)
					NZ <- normalize.invector(inp,inter)
					#
					ddd <- det(diag(1,tot_row) - G)
					if(ddd == 0)return(NA)
					IGT <- solve(diag(1,tot_row) - G)
					#
					# non-zero elements in the export vector (E)
					#
					count<-0
						for (i in 1:tot_row)if(outt[i]!=0)count <- count+1
					#
					# list for each output different from zero (E)_i != 0
					#
					final_output <- vector("list", count)
						for(i in 1:count){final_output[[i]] <- vector("list", 4)
								final_output[[i]][[1]] <- rep(0,tot_row)
								final_output[[i]][[2]] <- matrix(rep(0,tot_el), nrow = tot_row)
								final_output[[i]][[3]] <- rep(0,tot_row)
								final_output[[i]][[4]] <- rep(0,tot_row)
								}
					#
					step <- 1
					for(i in 1:tot_row)if(outt[i]!=0){final_output[[step]][[3]][i] <- 1
									step <- step + 1
									}
					#
					intermediate_step <- vector("list", count)
						for (i in 1:count)intermediate_step[[i]] <- IGT %*% final_output[[i]][[3]]
					#
					# defining normalized import and inter-compartmental exchanges
					#
					for(k in 1:count){
						for(i in 1:tot_row){
							final_output[[k]][[1]][i] <- intermediate_step[[k]][i] * NZ[i]
							for(j in 1:tot_row)final_output[[k]][[2]][j,i] <- intermediate_step[[k]][i] * G[j,i]
							}
						k <- k+1
						}
					#
					return(final_output)
					}





####################################################################################################################################
###																 ###
### Trophic analysis both in the classical framework of Canonical Trophic Aggregation (CTA) and in the "extended" version of CTA ###
###																 ###
####################################################################################################################################

## "Extended" version of the CTA. In this case the effect of imports with TPs far from 0 (i.e. prey migration) can be considered.
## It is also allowed to modulate the TL of "non-living imports" (TL = 0 or TL = 1).

TP.extendedCTA <- function(inp,inter,outt,diss,nl,TPimports,TPnonliving){
	xx <- ncol(inter)
	#
	# defining "living" imports (Z_living) with internal flows (T_living) 
	# and computing input from non_living nodes to living (K)
	#
	living <- xx - nl
	Z_living <- inp[c(1:living)]
	TTT <- inter[c(1:living),c(1:living)]
	E_living <- outt[c(1:living)]
	R_living <- diss[c(1:living)]
	#
	# exctracting the matrix of acyclic flows between living compartments
	#
	TTT.s1 <- cycling.analysis(Z_living,TTT,E_living,R_living)
	#
	{
	if(TTT.s1[[1]] == 0) T_living <- TTT.s1[[2]]
	#
	else T_living <- TTT.s1[[8]]
	}
	#
	{
	if(nl!=0){
		LVSNL <- matrix(inter[c((living+1):xx),c(1:living)],nrow=nl,byrow=FALSE)
		}
	#
	else LVSNL <- matrix(rep(0,living),nrow=1)
	}
	#
	K <- apply(LVSNL,2,sum)
	#
	# vector for the non_living compartments TLs (TK). In this version all the TLs are set to 1
	#
	TK <- rep(0,living)
		for(i in 1:living){
			if(K[i]!=0)TK[i]<-1
			}
	#
	# checking non-living input with TLs equal to 0 (flows from non-living to primary producers);
	# the vector (TPnonliving) is composed of "living" elements
	#
	for(i in 1:living){
		if(K[i]!=0 & TPnonliving[i]==0){
			Z_living[i] <- Z_living[i] + K[i]
			K[i] <- 0
			}
	}
	#
	# building partial feeding matrix G_living [G] and normalized import vectors N_living (N) and NK_living (NK)
	#
	IN <- rep(NA,living)
	G <- matrix(rep(NA,living^2), nrow=living)
	N <- rep(NA,living)
	NK <- rep(NA,living)
	for (j in 1:living) IN[j] <- (sum(T_living[,j]) + Z_living[j] + K[j])
	#
	for (i in 1:living){
		for (j in 1:living){
			if(IN[j]==0) G[i,j] <- 0
	     		#
	     		else G[i,j] <- T_living[i,j]/IN[j]
	     	}
	}
	#
	for (j in 1:living){
		if(IN[j]==0) N[j] <- 0
		#
		else N[j] <- Z_living[j]/IN[j]
	}
	#
	for (j in 1:living){
		if(IN[j]==0) NK[j] <- 0
		#     	
	     	else NK[j] <- K[j]/IN[j]
	}
	#
	# calculation of [G] matrix powers (computation interrupted after "living" steps)
	#
	G.lista <- as.list(rep(NA,living))
	#
	G.lista[[1]] <- diag(rep(1,living))
	G.lista[[2]] <- G
	for(q in 3:living)G.lista[[q]] <- G%*%G.lista[[q-1]]
	#
	# checking external migratory imports (with TPs far from 0);
	# the vector (TPimports) is composed of "living" elements
	#
	levels <- 0
	for(i in 1:living){
		if((Z_living[i]!=0) & (TPimports[i]!=0)) levels <- levels+1
	}
	#
	# vectors containing TPs of the migratory imports only
	#
	EFFECTIVE <- rep(NA,levels)
	pat <- 1
	for(i in 1:living){
		if(TPimports[i]!=0){
			EFFECTIVE[pat] <- TPimports[i]
			pat <- pat+1
			}
	}
	#
	k <- 1
	initialize <- rep(0,living)
	NP <- rep(0,living)
	TO <- rep(0,living)
	NPC <- rep(0,living)
	TONE <- rep(0,living)
	ooo <- rep(NA,levels)
	eee <- rep(NA,levels)
	normLIST <- as.list(ooo)
	trophicLIST <- as.list(eee)
	#
	if(levels!=0){
		for(i in 1:levels){
				normLIST[[i]] <- initialize
			   	trophicLIST[[i]] <- initialize
			  	}
	for(i in 1:living){
        	if(Z_living[i]!=0){
        	        if(TPimports[i]==0)NP[i] <- N[i]
			#
			else{
				normLIST[[k]][i] <- N[i]
				trophicLIST[[k]][i] <- TPimports[i]
				k <- k + 1
				}
          		}
		}
	}
	#
	# list containing Transformation Matrices with a Normalized Import Vectors for each migratory import
	#
	TEN <- matrix(rep(0,living^2),nrow = living)
	tot_lev <- levels+2
	el <- rep(NA,tot_lev)
	WHOLE <- as.list(el)
	INPUTS <- as.list(el)		
	#		
	for (i in 1:tot_lev) WHOLE[[i]] <- TEN
	#
	for (i in 1:tot_lev){
			if(i==1)INPUTS[[i]] <- t(NP)
			#
			else{
				if(i==2)INPUTS[[i]] <- t(NK)
				#
				else INPUTS[[i]] <- t(normLIST[[i-2]])
				}
			}		
	#	
	# preparation of a Transformation Matrix for each migratory import
	# and its storage into the list called WHOLE
	#		
	for(i in 1:tot_lev){
		for(q in 1:living)WHOLE[[i]][q,] <- INPUTS[[i]]%*%G.lista[[q]]	
	}
	#
	# preparation of vectors used to calculate TPs of living nodes
	#
	TROPHICVEC <- as.list(el)
	for(i in 1:tot_lev){
			if(i==1)TROPHICVEC[[i]] <- c(1:living)
			#
			else{
				if(i==2)TROPHICVEC[[i]] <- c(2:(living+1))
				#
				else TROPHICVEC[[i]]<- EFFECTIVE[i-2] * rep(1,living) + c(1:living)
			}
	}		
	#
	# computation of partial TPs and their storage into a list called PARTIAL.list
	#
	PARTIAL.list <- as.list(el)	
	for(i in 1:tot_lev)PARTIAL.list[[i]] <- TROPHICVEC[[i]]%*%WHOLE[[i]]
	#
	# TPs of living nodes
	#
	tr.vec <- rep(0,living)
	for(i in 1:tot_lev)tr.vec <- tr.vec + PARTIAL.list[[i]]
	#
	final.effectiveTPs <- t(tr.vec)
	#
	spl <- pat + 1
	AAA <- as.list(rep(NA,pat))
	DDD <- as.list(rep(NA,pat))
	#
	for(i in 2:spl){
		infe <- trunc(TROPHICVEC[[i]][1]-1)
		supe <- trunc(TROPHICVEC[[i]][1])
		supe.el <- TROPHICVEC[[i]][1] - supe
		#
		{
		if(supe.el==0){
			ZEROs <- matrix(rep(0,(living*infe)),nrow = infe)
			cuts <- c((living + 1):(living+infe))
			AAA.s1 <- rbind(ZEROs,WHOLE[[i]])
			AAA[[i-1]] <- AAA.s1[-cuts,]
			DDD[[i-1]] <- matrix(rep(0,living^2),nrow = living)
			}
		else	{
			infe.el <- 1 - supe.el
			#
			M.infe <- WHOLE[[i]] * infe.el
			M.supe <- WHOLE[[i]] * supe.el
			#
			ZEROs.infe <- matrix(rep(0,(living * infe)),nrow = infe)
			cuts.infe <- c((living + 1):(living + infe))
			AAA.s1 <- rbind(ZEROs.infe,M.infe)
			AAA[[i-1]] <- AAA.s1[-cuts.infe,]
			#
			ZEROs.supe <- matrix(rep(0,(living * supe)),nrow = supe)
			cuts.supe <- c((living + 1):(living + supe))
			DDD.s1 <- rbind(ZEROs.supe,M.supe)
			DDD[[i-1]] <- DDD.s1[-cuts.supe,]
			}
		}
	}
	#
	ONE.LIND <- WHOLE[[1]]
	for(i in 1:pat)ONE.LIND <- ONE.LIND + AAA[[i]] + DDD[[i]]
	#
	EL <- as.list(rep(NA,12))
	EL[[1]] <- ONE.LIND									## Lindeman Transformation Matrix
	#
	if(nl!=0)final.effectiveTPs.w <- c(final.effectiveTPs,rep(1,nl))
	#
	EL[[2]] <- final.effectiveTPs.w								## Effective Trophic Levels of each species
	#											## in presence of migratory imports
	{
	if(nl!=0){
		EL[[3]] <- c(EL[[1]] %*% t(t(E_living)),outt[c(living+1):xx])			## Canonical Exports
		R.det <- sum(diss[c(living+1):xx])
		R.det.v <- rep(0,nl)
		R.det.v[nl] <- R.det
		EL[[4]] <- c(EL[[1]] %*% t(t(R_living)),R.det.v)				## Canonical Respirations
		#
		DetPool <- matrix(inter[c(1:living),c((living+1):xx)],nrow=living,byrow=FALSE)
		EL6.s1 <- apply(DetPool,1,sum)
		EL[[6]] <- EL[[1]] %*% t(t(EL6.s1))						## Returns to Detrital Pool
		#
		EL[[8]] <- sum(inp[c((living+1):xx)])						## Input to Detrital Pool
		EL[[9]] <- sum(inter[c((living+1):xx),c((living+1):xx)])			## Circulation within Detrital Pool
		}
	else	{
		 DetPool <- 0
		 EL[[3]] <- EL[[1]] %*% t(t(E_living))
		 EL[[4]] <- EL[[1]] %*% t(t(R_living))
		 EL[[6]] <- 0
		 EL[[8]] <- 0
		 EL[[9]] <- 0
		}
	}
	#
	detritivory.s1 <- apply(LVSNL,2,sum)
	EL[[7]] <- sum(detritivory.s1[which(TPnonliving==1)])					## Detritivory
	#
	EL5.s1 <- apply(T_living,1,sum)
	EL5.s2 <- EL[[1]] %*% t(t(EL5.s1))
	lll <- living + 1
	IN.liv <- apply(T_living,2,sum)
	prim <- which(Z_living != 0 & IN.liv == 0)
	Z.prim <- sum(Z_living[prim])
	EL5.s3 <- c(Z.prim,EL5.s2)
	detr.to.pp <- sum(detritivory.s1[which(TPnonliving==0)])
	EL5.s3 <- EL5.s3[-lll]
	EL5.s3[1] <- EL5.s3[1] + detr.to.pp
	EL[[5]] <- EL5.s3									## The Grazing Chain
	#
	Det.not.p <- DetPool
	if(nl!=0)Det.not.p[prim,] <- rep(0,nl)
	lind.s1 <- EL[[5]]
	pr.pro <- which(detritivory.s1 != 0 & TPnonliving == 0)
	{
	if(length(pr.pro)==0)lind.s1[1] <- lind.s1[1] + sum(Det.not.p) + EL[[8]]
	#
	else lind.s1[1] <- lind.s1[1] + sum(Det.not.p) + sum(detritivory.s1[pr.pro])
	}
	lind.s1[2] <- lind.s1[2] + EL[[7]]
	EL[[10]] <- lind.s1									## Lindeman Spine
	#
	lind.eff <- rep(0,(living-1))
	for(q in 2:living){
		if(EL[[10]][q-1] > 0.0001){
			lind.eff[q-1] <- EL[[10]][q]/EL[[10]][q-1]
			}
	}
	EL[[11]] <- lind.eff									## Trophic Efficiencies
	#
	mig <- which(Z_living != 0 & IN.liv != 0)
	if(length(mig)==0)mig <- 0
	EL[[12]] <- mig										## Detecting Migratory Imports
	return(EL)
}


## Canonical Trophic Aggregation (CTA) with all the imports from outside with Hypothetical TP = 0 and non living nodes set to TP = 1 
## (as primary producers). No migratory imports are considered and only a warning message is supplied in this case.

TP.CTA <- function(inp,inter,outt,diss,nl){
	xx <- ncol(inter)
	#
	# defining "living" imports (Z_living) with internal flows (T_living) 
	# and computing input from non_living nodes to living (K)
	#
	living <- xx - nl
	Z_living <- inp[c(1:living)]
	TTT <- inter[c(1:living),c(1:living)]
	E_living <- outt[c(1:living)]
	R_living <- diss[c(1:living)]
	#
	# exctracting the matrix of acyclic flows between living compartments
	#
	TTT.s1 <- cycling.analysis(Z_living,TTT,E_living,R_living)
	#
	{
	if(TTT.s1[[1]] == 0) T_living <- TTT.s1[[2]]
	#
	else T_living <- TTT.s1[[8]]
	}
	#
	{
	if(nl!=0){
		LVSNL <- matrix(inter[c((living+1):xx),c(1:living)],nrow=nl,byrow=FALSE)
		}
	#
	else LVSNL <- matrix(rep(0,living),nrow=1)
	}
	#
	K <- apply(LVSNL,2,sum)
	#
	# vector for the non_living compartments TLs (TK). In this version all the TLs are set to 1
	#
	TK <- rep(0,living)
	for(i in 1:living){
		if(K[i]!=0)TK[i]<-1
		}
	#
	# building partial feeding matrix G_living [G] and normalized import vectors N_living (N) and NK_living (NK)
	#
	IN <- rep(NA,living)
	G <- matrix(rep(NA,living^2), nrow=living)
	N <- rep(NA,living)
	NK <- rep(NA,living)
	for (j in 1:living) IN[j] <- (sum(T_living[,j]) + Z_living[j] + K[j])
	#
	for (i in 1:living){
		for (j in 1:living){
			if(IN[j]==0) G[i,j] <- 0
	     		#
	     		else G[i,j] <- T_living[i,j]/IN[j]
	     	}
	}
	#
	for (j in 1:living){
		if(IN[j]==0) N[j] <- 0
		#
		else N[j] <- Z_living[j]/IN[j]
	}
	#
	for (j in 1:living){
		if(IN[j]==0) NK[j] <- 0
		#     	
	     	else NK[j] <- K[j]/IN[j]
	}
	#
	G.lista <- as.list(rep(NA,living))
	#
	G.lista[[1]] <- diag(rep(1,living))
	G.lista[[2]] <- G
	for(q in 3:living)G.lista[[q]] <- G%*%G.lista[[q-1]]
	#
	# Trophic Transformation Matrices: [A] for input set to TL = 0 
	# and [B] for non-living inputs with TL = 1 (called AAA e BBB)
	#
	EL <- as.list(rep(NA,12))
	TROPHIC <- c(1:living)
	ONE <- rep(1,living)
	#
	# calculation of [A]
	#
	AAA <- matrix(rep(0,living^2),nrow=living)
	for(q in 1:living)AAA[q,] <- N%*%G.lista[[q]]
	#
	# calculation of [B]
	#
	BBB <- matrix(rep(0,living^2),nrow=living)
	for(q in 1:living)BBB[q,] <- NK%*%G.lista[[q]]
	#
	# multiplication of [A] and [B] matrices with the corresponding (TROPHIC) and (ONE) vectors of TLs
	#
	AAAPOS <- TROPHIC%*%AAA
	BBBPOS <- (ONE+TROPHIC)%*%BBB
	#
	# saving [A] and [B]
	#
	BBB <- rbind(rep(0,living),BBB)
	BBB <- BBB[-(living+1),]
	EL[[1]] <- AAA + BBB									## Lindeman Transformation Matrix
	#
	# final trophic positions vector for the living compartments (finalTPs)
	#
	intermTPs <- AAAPOS+BBBPOS
	finalTPs <- t(intermTPs)
	EL[[2]] <- c(finalTPs,rep(1,nl))							## Effective Trophic Levels of each species
	#
	{
	if(nl!=0){
		EL[[3]] <- c(EL[[1]] %*% t(t(E_living)),outt[c(living+1):xx])			## Canonical Exports
		R.det <- sum(diss[c(living+1):xx])
		R.det.v <- rep(0,nl)
		R.det.v[nl] <- R.det
		EL[[4]] <- c(EL[[1]] %*% t(t(R_living)),R.det.v)				## Canonical Respirations
		#
		DetPool <- matrix(inter[c(1:living),c((living+1):xx)],nrow=living,byrow=FALSE)
		EL6.s1 <- apply(DetPool,1,sum)
		EL[[6]] <- EL[[1]] %*% t(t(EL6.s1))					## Returns to Detrital Pool
		#
		EL[[8]] <- sum(inp[c((living+1):xx)])						## Input to Detrital Pool
		EL[[9]] <- sum(inter[c((living+1):xx),c((living+1):xx)])			## Circulation within Detrital Pool
		}
	else	{
		 DetPool <- 0
		 EL[[3]] <- EL[[1]] %*% t(t(E_living))
		 EL[[4]] <- EL[[1]] %*% t(t(R_living))
		 EL[[6]] <- 0
		 EL[[8]] <- 0
		 EL[[9]] <- 0
		}
	}
	#
	EL[[7]] <- sum(LVSNL)									## Detritivory
	#
	EL5.s1 <- apply(T_living,1,sum)
	EL5.s2 <- EL[[1]] %*% t(t(EL5.s1))
	lll <- living + 1
	IN.liv <- apply(T_living,2,sum)
	prim <- which(Z_living != 0 & IN.liv == 0)
	Z.prim <- sum(Z_living[prim])
	EL5.s3 <- c(Z.prim,EL5.s2)
	EL5.s3 <- EL5.s3[-lll]
	EL[[5]] <- EL5.s3									## The Grazing Chain
	#
	Det.not.p <- DetPool
	if(nl!=0)Det.not.p[prim,] <- rep(0,nl)
	lind.s1 <- EL[[5]]
	lind.s1[1] <- lind.s1[1] + sum(Det.not.p) + EL[[8]]
	lind.s1[2] <- lind.s1[2] + EL[[7]]
	EL[[10]] <- lind.s1									## Lindeman Spine
	#
	lind.eff <- rep(0,(living-1))
	for(q in 2:living){
		if(EL[[10]][q-1] > 0.0001){
			lind.eff[q-1] <- EL[[10]][q]/EL[[10]][q-1]
			}
	}
	EL[[11]] <- lind.eff									## Trophic Efficiencies
	#
	mig <- which(Z_living != 0 & IN.liv != 0)
	{
	if(length(mig)==0)mig <- 0
	#
	else	{
		TP.migr <- rep(0,living)
		TP.detr <- rep(0,living)
		#
		TP.migr[mig] <- trunc(EL[[2]])[mig]
		ONEs <- which(apply(LVSNL,2,sum)!=0)
		TP.detr[ONEs] <- 1
		EL <- TP.extendedCTA(inp,inter,outt,diss,nl,TP.migr,TP.detr)
		#
		EL[[13]] <- EL[[1]]%*%t(t(inp[1:living]))					## Canonical Distribution of inputs
		EL[[13]] <- c(EL[[13]],inp[(living+1):xx])					## (only in case of migrations)
		names(EL[[13]]) <- NULL
		}
	}
	EL[[12]] <- mig										## Detecting Migratory Imports
	#
	return(EL)
}


## The function "Higashi" returns the trophic positions (TPs) of all the compartments (both living and non-living).
## This algorithm is based on the original matrix-based unfolding method (Higashi et al. 1989, 1991, 1992).
## No matrix transformation is required for extracting the compartments' TPs (both in the case of cyclic and acyclic networks).
## TPs are computed without the need of an arbitrary trophic level assignment to non-living nodes.
## This method considers any passage of energy as incrementing the TP, also in the presence of non-trophic flows that
## represent a change of state without a transformation of energy ("trophic level inflation").

Higashi <- function(inp,inter){
	G <- partial.feeding(inp,inter)
	tutt <- nrow(G)
	I <- diag(1,tutt)
	#
	ddd <- det(I-G)
	if(ddd == 0)return(NA)
	IG <- solve(I-G)
	HIGA <- apply(IG,2,sum)
	#
	return(HIGA)
	}


## With the "Whipple" function is implemented the path-based network unfolding analysis (Whipple 1998).
## The whole set of flows in the matrix is considered (decyclization algorithms are not required) and this avoids that detritus
## is assigned to an arbitrary trophic level.
## The path-based method does not increase the trophic level in presence of non-trophic interactions which occur in the pathway such as,
## for example, egestive transfers from living to non-living compartments. In this way, the inconsistencies originated when applying
## the original matrix-based unfolding method ("trophic level inflation") are solved.

Whipple <- function(inp,inter,outt,diss,nl){
	t.el <- nrow(inter)
	l.el <- t.el - nl
	#
	ZZZ <- inp
	EEE <- outt + diss 
	TTT <- inter
	#
	##
	##
	##
	#
	if(nl!=0){
		#
		# Constructing a whole non-living node for which connections with living compartments are detailed as the
		# last row and column in a DET.TTT matrix and last element in DET.ZZZ & DET.EEE vectors
		#
		{if(nl>1){
			#
			DET.TTT <- matrix(rep(0,(l.el+1)^2),nrow=(l.el+1))
			DET.ZZZ <- rep(NA,(l.el+1))
			DET.EEE <- rep(NA,(l.el+1))
			#
			# Constructing non-living vectors and matrix
			#
			FROM.det <- matrix(rep(0,t.el*nl),nrow=nl)
			ESP.det <- rep(NA,nl)
			TO.det <- matrix(rep(0,t.el*nl),nrow=t.el)
			IMP.det <- rep(NA,nl)
			#
			# Sub-matrix of flows from non-living to all the system nodes
			#
			for(i in (l.el+1):t.el)
				for(j in 1:t.el)FROM.det[i-l.el,j] <- TTT[i,j]
			#
			# Vector of total outflows from non-living compartments (R) + (E)
				#
			for(i in (l.el+1):t.el)ESP.det[i-l.el] <- EEE[i]
			#
			# Sub-matrix of flows to non-living from all the system nodes
			#
			for(i in 1:t.el)
				for(j in (l.el+1):t.el)TO.det[i,j-l.el] <- TTT[i,j]
			#
			# Vector of import to non-living compartments
			#
			for(i in (l.el+1):t.el)IMP.det[i-l.el] <- ZZZ[i]
			#
			for(i in 1:l.el)
				for(j in 1:l.el)DET.TTT[i,j] <- TTT[i,j]
			#
			for(j in 1:l.el){
					DET.TTT[l.el+1,j] <- sum(FROM.det[,j])
					DET.TTT[j,l.el+1] <- sum(TO.det[j,])
					}
			for(i in (l.el+1):t.el)
				for(j in (l.el+1):t.el)DET.TTT[l.el+1,l.el+1] <- DET.TTT[l.el+1,l.el+1] + TTT[i,j]
			#
			for(i in 1:l.el){
					DET.ZZZ[i] <- ZZZ[i]
					DET.EEE[i] <- EEE[i]
					}
			DET.ZZZ[l.el+1] <- sum(IMP.det)
			DET.EEE[l.el+1] <- sum(ESP.det)
			#
			}
		else{DET.TTT <- TTT
		DET.ZZZ <- ZZZ
		DET.EEE <- EEE}
		}
		#
		###
		###
		###
		###
		###
		#
		tot.FROM.DET <- sum(DET.TTT[l.el+1,]) + DET.EEE[l.el+1]
		tot.TO.DET <- sum(DET.TTT[,l.el+1]) + DET.ZZZ[l.el+1]
		#
		# Assembling normalized contribution to and from the whole non-living node obtained from
		# the previous steps. Also external import and export are normalized. Also the possible
		# self-loop (from non-living to non-living) is normalized
		#
		NORM.FROM.det <- rep(NA,l.el)
		NORM.TO.det <- rep(NA,l.el)
		NORM.ZZZ.det <- DET.ZZZ[l.el+1]/tot.TO.DET
		NORM.EEE.det <- DET.EEE[l.el+1]/tot.FROM.DET
		NORM.DET.det <- DET.TTT[l.el+1,l.el+1]/tot.TO.DET
		ALL.WITH.DET <- DET.TTT[l.el+1,l.el+1]
		#
		for(i in 1:l.el){
				NORM.FROM.det[i] <- DET.TTT[l.el+1,i]/tot.FROM.DET
				NORM.TO.det[i] <- DET.TTT[i,l.el+1]/tot.TO.DET
				}
		#
		# Defining new matrix and vectors, only with living flows (e.g. effective energy transfer
		# with ingestion/assimilation)
		#
		TO <- DET.TTT[,nrow(DET.TTT)]
		FROM <- DET.TTT[nrow(DET.TTT),]
		TO <- TO[-length(TO)]
		FROM <- FROM[-length(FROM)]
		#
		## 
		##
		##
		#
		new.TTT <- DET.TTT[-nrow(DET.TTT),-nrow(DET.TTT)]
		new.ZZZ <- DET.ZZZ[-length(DET.ZZZ)]
		new.EEE <- DET.EEE[-length(DET.EEE)]
		#
		# Redistributing direct flows from non-living whole compartment to other living nodes
		#
		for(i in 1:l.el)
			for(j in 1:l.el)if(NORM.TO.det[i]!=0 & NORM.FROM.det[j]!=0)new.TTT[i,j] <- new.TTT[i,j] + TO[i]*NORM.FROM.det[j]/(sum(NORM.FROM.det)+NORM.EEE.det)
		#
		# Redistributing to export the flows to non-living whole compartment from the living
		# nodes showing direct connections to it
		#	
		for(i in 1:l.el){
			new.EEE[i] <- new.EEE[i] + TO[i]*NORM.EEE.det/(sum(NORM.FROM.det)+NORM.EEE.det)
			}		
		#
		# Allocating input to non-living node in terms of outflows from living compartments: exports
		# and intercompartmental exchanges
		#
		if(DET.ZZZ[length(DET.ZZZ)]!=0){
						for(i in 1:l.el)
							for(j in 1:l.el)if(NORM.TO.det[i]!=0 & NORM.FROM.det[j]!=0){
								new.TTT[i,j] <- new.TTT[i,j] + (DET.ZZZ[length(DET.ZZZ)]*(NORM.TO.det[i]/sum(NORM.TO.det))*NORM.FROM.det[j]/(sum(NORM.FROM.det)+NORM.EEE.det))
								}
		#
		#
		for(i in 1:l.el){
			new.EEE[i] <- new.EEE[i] + (DET.ZZZ[length(DET.ZZZ)]*(NORM.TO.det[i]/sum(NORM.TO.det))*NORM.EEE.det/(sum(NORM.FROM.det)+NORM.EEE.det))
			}
		#
		# Setting external input to living nodes with direct connections to non-living whole compartment
		#
		for(i in 1:l.el){
			new.ZZZ[i] <- new.ZZZ[i] + (DET.ZZZ[length(DET.ZZZ)]*(NORM.TO.det[i]/sum(NORM.TO.det)))
			}
		#
		# Checking for steady state after reassembling, when whole non-living node
		# is deleted from the system
		#
		IN.new <- apply(new.TTT,2,sum)+new.ZZZ
		OUT.new <- apply(new.TTT,1,sum)+new.EEE
		DIFF.corr <- rep(0,l.el)
		for(u in 1:l.el)if((IN.new[u]-OUT.new[u])<0.0001) DIFF.corr[u] <- 0
		#
		}
	}
	G <- matrix(rep(0,l.el^2),nrow=l.el)
	for(e in 1:l.el)
		for(f in 1:l.el){
				if(IN.new[f]!=0)G[e,f] <- new.TTT[e,f]/IN.new[f]
				}
	#
	I <- diag(1,l.el)
	ddd <- det(I-G)
	if(ddd == 0)return(NA)
	IG <- solve(I-G)
	WHIP <- apply(IG,2,sum)
	#
	return(WHIP)
	}





############################################################################################
###											 ###
### System level index for the assessment of the whole ecosystem functioning		 ###
### We also included the analysis of some sinle-species features that are extracted from ###
### global indices (e.g., AMI on input and output)					 ###
###											 ###
############################################################################################

## The function "TST" is used to calculate the Total System Throughput. The inputs are vectors of import (Z), 
## export (E) and dissipations (R) with the matrix of intercompartmental exchanges [T].

TST <- function(inp,inter,outt,diss){
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				B <- sum(A)
				return(B)
				}


## The function "AMI" calculates the Average Mutual Information of the ecosystem. Also in this case, as in the following,
## the only needed inputs are given by a matrix of internal exchanges [T] with vectors of import (Z), export (E) and respirations (R).

AMI <- function(inp,inter,outt,diss){
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				tot <- TST(inp,inter,outt,diss)
				xx <- ncol(A)
				ami <- 0
				for (i in 1:xx){
					for (j in 1:xx){
						if(A[i,j]==0)ami <- (ami+0)
							else if (sum(A[,j])==0)ami <- (ami+0)
								else if (sum(A[i,])==0)ami <- (ami+0)
									else ami <- ami+(A[i,j]/tot)*log2(A[i,j]*tot/(sum(A[i,])*sum(A[,j])))
							}
						}
				return(ami)
				}


## The function "DC" calculates the Development Capacity (the upper limit to the development of an ecosystem).

DC <- function(inp,inter,outt,diss){
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				tot <- TST(inp,inter,outt,diss)
				xx <- ncol(A)
				dc <- 0
				for (i in 1:xx){
					for (j in 1:xx){
						if(A[i,j]==0)dc <- (dc+0)
							else dc <- dc+A[i,j]*log2(A[i,j]/tot)
							}
						}
				dev.cap<-(-1*dc)
				return(dev.cap)
				}


## The function "AMI.max" estimates the AMI corresponding to the Development Capacity.

AMI.max <- function(inp,inter,outt,diss){
					dev <- DC(inp,inter,outt,diss)
					tot <- TST(inp,inter,outt,diss)
					a.m <- dev/tot
					return(a.m)
					}


## The function "ASC" to compute the Ascendency of the system (TST * AMI).

ASC <- function(inp,inter,outt,diss){
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				tot <- TST(inp,inter,outt,diss)
				develop <- DC(inp,inter,outt,diss)
				xx <- ncol(A)
				asc <- 0
				for (i in 1:xx){
					for (j in 1:xx){
						if(A[i,j]==0)asc <- (asc+0)
							else if (sum(A[,j])==0)asc <- (asc+0)
								else if (sum(A[i,])==0)asc <- (asc+0)
									else asc <- asc+A[i,j]*log2(A[i,j]*tot/(sum(A[i,])*sum(A[,j])))
							}
						}
				perc <- asc/develop
				end <- matrix(c(asc,perc),nrow=1,byrow=TRUE)
				rownames(end) <- c("Ascendency")
				colnames(end) <- c("Value", "Ratio")				
				return(end)
				}


## The function "OVERHEAD.imports" determines the Overhead on Imports.

OVERHEAD.imports <- function(inp,inter,outt,diss){
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				develop <- DC(inp,inter,outt,diss)
				kk <- ncol(A)-2
				io <- 0
				for (j in 1:kk){
					if(A[1,j]==0)io <- (io+0)
						else io <- io+A[1,j]*log2((A[1,j]^2)/(sum(A[1,])*sum(A[,j])))
						}
				#
				{
				if(sum(io!=0))import.o <- (-1*io)
				else import.o <- 0
				}
				#
				perc <- import.o/develop
				end <- matrix(c(import.o,perc),nrow=1,byrow=TRUE)
				rownames(end) <- c("Overhead on Imports")
				colnames(end) <- c("Value", "Ratio")
				return(end)
				}


## The function "OVERHEAD.exports" determines the Overhead on Exports.

OVERHEAD.exports <- function(inp,inter,outt,diss){
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				develop <- DC(inp,inter,outt,diss)
				kk <- ncol(A)-2
				uu <- ncol(A)-1
				eo <- 0
				for (i in 1:kk){
					if(A[i,uu]==0)eo <- (eo+0)
						else eo <- eo+A[i,uu]*log2((A[i,uu]^2)/(sum(A[i,])*sum(A[,uu])))
						}
				#
				{
				if(sum(eo!=0))export.o <- (-1*eo)
				else export.o <- 0
				}
				#
				perc <- export.o/develop
				end <- matrix(c(export.o,perc),nrow=1,byrow=TRUE)
				rownames(end) <- c("Overhead on Exports")
				colnames(end) <- c("Value", "Ratio")
				return(end)
				}


## The function "OVERHEAD.resp" determines the Dissipative Overhead.

OVERHEAD.resp <- function(inp,inter,outt,diss){
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				develop <- DC(inp,inter,outt,diss)
				xx <- ncol(A)
				kk <- ncol(A)-2
				do <- 0
				for (i in 1:kk){
					if(A[i,xx]==0)do <- (do+0)
						else do <- do+A[i,xx]*log2((A[i,xx]^2)/(sum(A[i,])*sum(A[,xx])))
						}
				#
				{
				if(sum(do!=0))dissipative.o <- (-1*do)
				else dissipative.o <- 0
				}
				#
				perc <- dissipative.o/develop
				end <- matrix(c(dissipative.o,perc),nrow=1,byrow=TRUE)
				rownames (end) <- c("Dissipative Overhead")
				colnames(end) <- c("Value", "Ratio")
				return(end)
				}


## The function "REDUNDANCY" determines the Redundancy (Overhead on internal flows).

REDUNDANCY <- function(inp,inter,outt,diss){
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				develop <- DC(inp,inter,outt,diss)
				kk <- ncol(A)-2
				r <- 0
				for (i in 2:kk){
					for (j in 2:kk){
						if(A[i,j]==0)r <- (r+0)
							else r <- r+A[i,j]*log2((A[i,j]^2)/(sum(A[i,])*sum(A[,j])))
							}
						}
				#
				{
				if(sum(r)!=0)redun <- (-1*r)
				else redun <- 0
				}
				#
				perc <- redun/develop
				end <- matrix(c(redun,perc),nrow=1,byrow=TRUE)
				rownames(end) <- c("Redundancy")
				colnames(end) <- c("Value", "Ratio")
				return(end)
				}


## The function "OVERHEAD" determines the total amount of System Overhead.

OVERHEAD <- function(inp,inter,outt,diss){
					o1 <- OVERHEAD.imports(inp,inter,outt,diss)[1,1]
					o2 <- OVERHEAD.exports(inp,inter,outt,diss)[1,1]
					o3 <- OVERHEAD.resp(inp,inter,outt,diss)[1,1]
					o4 <- REDUNDANCY(inp,inter,outt,diss)[1,1]
					develop <- DC(inp,inter,outt,diss)
					o <- o1 + o2 + o3 + o4
					perc <- o/develop
					end <- matrix(c(o,perc),nrow=1,byrow=TRUE)
					rownames(end) <- c("Total Overhead")
					colnames(end) <- c("Value", "Ratio")
					return(end)
					}


## The function "internal.DC" is used to calculate the Internal Development Capacity.

internal.DC <- function(inp,inter,outt,diss){
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				tot <- TST(inp,inter,outt,diss)
				kk <- ncol(A)-2
				ic <- 0
				for (i in 2:kk){
					for (j in 2:kk){
						if(A[i,j]==0)ic <- (ic+0)
							else ic <- ic+A[i,j]*log2(A[i,j]/tot)
							}
						}
				#
				{
				if(sum(ic)!=0)intc <- (-1*ic)
				else intc <- 0
				}
				#
				return(intc)
				}


## The function "internal.ASC" allows the computation of the Internal Ascendency.

internal.ASC <- function(inp,inter,outt,diss){
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				tot <- TST(inp,inter,outt,diss)
				red <- REDUNDANCY(inp,inter,outt,diss)[1,1]
				idc <- internal.DC(inp,inter,outt,diss)
				kk <- ncol(A)-2
				iasc <- 0
				for (i in 2:kk){
					for(j in 2:kk){
						if(A[i,j]!=0)iasc <- iasc+A[i,j]*log2(A[i,j]*tot/(sum(A[i,])*sum(A[,j])))
							}
						}
				perc <- iasc/idc
				perc1 <- red/idc
				end <- matrix(c(iasc,perc,red,perc1),nrow=2,byrow=TRUE)
				rownames(end) <- c("Internal Ascendency", "Redundancy")
				colnames(end) <- c("Value", "Ratio")
				return(end)
				}


## The function for computing the overall connectance also includes the effects of exogenous transfers.

O.C <- function(inp,inter,outt,diss){
				tot <- TST(inp,inter,outt,diss)
				el <- length(inp)
				#
				in.v <- apply(inter,2,sum) + inp
				out.v <- apply(inter,1,sum) + outt + diss
				#
				C.M <- matrix(rep(0,el^2),nrow=el)
				V.I <- rep(0,el)
				V.E <- rep(0,el)
				V.R <- rep(0,el)
				#
				for(i in 1:el)
					for(j in 1:el)if(inter[i,j]!=0)C.M[i,j] <- (inter[i,j]/tot) * log(inter[i,j]^2/(in.v[j]*out.v[i]))
				#
				for(i in 1:el){
					if(inp[i]!=0)V.I[i] <- (inp[i]/tot) * log(inp[i]/in.v[i])
					#
					if(outt[i]!=0)V.E[i] <- (outt[i]/tot) * log(outt[i]/out.v[i])
					#
					if(diss[i]!=0)V.R[i] <- (diss[i]/tot) * log(diss[i]/out.v[i])
				}
				#
				TT <- sum(-1 * sum(C.M))
				II <- sum(-2 * sum(V.I))
				EE <- sum(-2 * sum(V.E))
				RR <- sum(-2 * sum(V.R))
				#
				fi <- sum(TT + II + EE + RR)
				#
				m <- exp(fi/2)
				#
				return(m)
				}


## The intercompartmental connectance characterizes only the endogenous exchanges.

I.C <- function(inter){
		tot <- sum(inter)
		el <- nrow(inter)
		#
		in.v <- apply(inter,2,sum)
		out.v <- apply(inter,1,sum)
		#
		C.M <- matrix(rep(0,el^2),nrow=el)
		#
		for(i in 1:el)
			for(j in 1:el)if(inter[i,j]!=0)C.M[i,j] <- (inter[i,j]/tot) * log(inter[i,j]^2/(in.v[j]*out.v[i]))
		#
		fi <- sum(-1 * sum(C.M))
		#
		m <- exp(fi/2)
		#
		return(m)
		}


## The foodweb connectance pertains only to transfers among the living compartments.

FW.C <- function(inter,nl){
		el <- nrow(inter) - nl
		inter.l <- inter[c(1:el),c(1:el)]
		tot <- sum(inter.l)
		#
		in.v <- apply(inter.l,2,sum)
		out.v <- apply(inter.l,1,sum)
		#
		C.M <- matrix(rep(0,el^2),nrow=el)
		#
		for(i in 1:el)
			for(j in 1:el)if(inter.l[i,j]!=0)C.M[i,j] <- (inter.l[i,j]/tot) * log(inter.l[i,j]^2/(in.v[j]*out.v[i]))
		#
		fi <- sum(-1 * sum(C.M))
		#
		m <- exp(fi/2)
		#
		return(m)
		}


## Weighted connectance computed with the formula described by Zorach and Ulanowicz (2003).
## It is equivalent to the effective connectance (m) introduced by Ulanowicz and Wolff (1991).

effconn.m <- function(inp,inter,outt,diss){
				tot <- TST(inp,inter,outt,diss)
				in.v <- apply(inter,2,sum) + inp
				out.v <- apply(inter,1,sum) + outt + diss
				el <- nrow(inter)
				#
				C.M <- matrix(rep(1,el^2),nrow=el)
				#
				for(i in 1:el)
					for(j in 1:el)if(inter[i,j]!=0)C.M[i,j] <- (inter[i,j]^2/(in.v[j]*out.v[i]))^(-0.5*inter[i,j]/tot)
				#
				m <- prod(C.M)
				#
				return(m)
				}


## Link density (LD) estimated with the Bersier et al. (2002) formulation.
## It is slightly different from the Ulanowicz's effective connectance (m).
## With LD, the weighting applies directly to the taxa's equivalent numbers of prey and consumers.

LD.B <- function(inp,inter,outt,diss){
				tot <- TST(inp,inter,outt,diss)
				in.v <- apply(inter,2,sum) + inp
				out.v <- apply(inter,1,sum) + outt + diss
				el <- nrow(inter)
				#
				H.N.M <- matrix(rep(0,el^2),nrow=el)
				H.P.M <- matrix(rep(0,el^2),nrow=el)
				N.v <- rep(0,el)
				P.v <- rep(0,el)
				#
				for(i in 1:el)
					for(j in 1:el){
						if(inter[i,j]!=0){
							H.N.M[i,j] <- (-inter[i,j]/in.v[j]) * log2(inter[i,j]/in.v[j])
							#
							H.P.M[i,j] <- (-inter[i,j]/out.v[i]) * log2(inter[i,j]/out.v[i])
							}
						}
				#
				H.N <- apply(H.N.M,2,sum)
				H.P <- apply(H.P.M,1,sum)
				#
				for(i in 1:el){
					P.v[i] <- (out.v[i]/tot) * 2^H.P[i]
					#
					N.v[i] <- (in.v[i]/tot) * 2^H.N[i]
					}
				#
				LD <- 0.5 * (sum(P.v) + sum(N.v))
				#
				return(LD)
				}


## Biomass inclusive Ascendency (Ulanowicz and Abarca-Arenas, 1997).
## It is with intercompartmental exchanges only (equivalent to the Internal Ascendency).

ASC.B <- function(inp,inter,outt,diss,bi){
				tot <- TST(inp,inter,outt,diss)
				xx <- ncol(inter)
				B <- sum(bi)
				#
				asc <- 0
				for (i in 1:xx){
					for (j in 1:xx)if(inter[i,j]!=0)asc <- asc + (inter[i,j] * log2((inter[i,j] * B^2)/(tot * bi[i] * bi[j])))
						}
				#
				return(asc)
				}


## AMI on input for each compartment. This function estimates two versions of the index:
## 1) the first takes into account all the flows (intercompartmental exchanges both for living and non-living compartments
## as well as vectors of import, export and respiration); 2) the second includes living exchanges only. 
## These indices are described in Scotti et al. (2009).

AMI.in <- function(inp,inter,outt,diss,nl){
				tot <- TST(inp,inter,outt,diss)
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				A.in <- apply(A,2,sum)
				A.ou <- apply(A,1,sum)
				ww <- nrow(A)
				ll <- nrow(inter)
				#
				A.w <- matrix(rep(0,ww^2),nrow=ww)
				#
				for(i in 1:ww)
					for(j in 1:ww)if(A[i,j]!=0)A.w[i,j] <- A[i,j] * log2((A[i,j] * tot)/(A.in[j] * A.ou[i]))
				#
				A.w.c <- apply(A.w,2,sum)
				AMI.w.c <- A.w.c[2:(ll+1)]/A.in[2:(ll+1)]
				#
				#
				live <- ll - nl
				L <- inter[c(1:live),c(1:live)]
				tot.l <- sum(L)
				A.l.in <- apply(L,2,sum)
				A.l.ou <- apply(L,1,sum)
				#
				A.l <- matrix(rep(0,live^2),nrow=live)
				#
				for(i in 1:live)
					for(j in 1:live)if(L[i,j]!=0)A.l[i,j] <- L[i,j] * log2((L[i,j] * tot.l)/(A.l.in[j] * A.l.ou[i]))
				#
				A.l.c <- apply(A.l,2,sum)
				#
				AMI.l.c <- rep(0,live)
				#
				for(i in 1:live){
					if(A.l.c[i]!=0)AMI.l.c[i] <- A.l.c[i]/A.l.in[i]
					}
				#
				# this function returns two outcomes for the AMI on inflows:
				# 1) the first is estimated using the "classical" ENA approach (including
				# non-living nodes and exchanges with surroundings (i.e., import,
				# export and respiration vectors);
				# 2) the second is computed by taking into account the flows between
				# living compartments only
				#
				AMI.IN.FINAL <- as.list(rep(NA,2))
				names(AMI.IN.FINAL) <- c("ENA","living")
				AMI.IN.FINAL[[1]] <- AMI.w.c
				AMI.IN.FINAL[[2]] <- AMI.l.c
				#
				return(AMI.IN.FINAL)
				}


## AMI on export for each compartment. This index is measured:
## 1) for the whole flows, with the "classical" network analysis approach; 2) for living exchanges only.
## Further details are provided with the description of the function "AMI.in".

AMI.out <- function(inp,inter,outt,diss,nl){
				tot <- TST(inp,inter,outt,diss)
				A <- bigmatrNOnames.ena(inp,inter,outt,diss)
				A.in <- apply(A,2,sum)
				A.ou <- apply(A,1,sum)
				ww <- nrow(A)
				ll <- nrow(inter)
				#
				A.w <- matrix(rep(0,ww^2),nrow=ww)
				#
				for(i in 1:ww)
					for(j in 1:ww)if(A[i,j]!=0)A.w[i,j] <- A[i,j] * log2((A[i,j] * tot)/(A.in[j] * A.ou[i]))
				#
				A.w.r <- apply(A.w,1,sum)
				AMI.w.r <- A.w.r[2:(ll+1)]/A.ou[2:(ll+1)]
				#
				#
				live <- ll - nl
				L <- inter[c(1:live),c(1:live)]
				tot.l <- sum(L)
				A.l.in <- apply(L,2,sum)
				A.l.ou <- apply(L,1,sum)
				#
				A.l <- matrix(rep(0,live^2),nrow=live)
				#
				for(i in 1:live)
					for(j in 1:live)if(L[i,j]!=0)A.l[i,j] <- L[i,j] * log2((L[i,j] * tot.l)/(A.l.in[j] * A.l.ou[i]))
				#
				A.l.r <- apply(A.l,1,sum)
				#
				AMI.l.r <- rep(0,live)
				#
				for(i in 1:live){
					if(A.l.r[i]!=0) AMI.l.r[i] <- A.l.r[i]/A.l.ou[i]
					}
				#
				# this function returns two outcomes for the AMI on outflows:
				# 1) AMI on outflows using the "classical" ENA approach (including
				# non-living nodes and exchanges with surroundings (i.e., import,
				# export and respiration vectors);
				# 2) AMI on outflows computed by taking into account the flows between
				# living compartments only
				#
				AMI.OU.FINAL <- as.list(rep(NA,2))
				names(AMI.OU.FINAL) <- c("ENA","living")
				AMI.OU.FINAL[[1]] <- AMI.w.r
				AMI.OU.FINAL[[2]] <- AMI.l.r
				#
				return(AMI.OU.FINAL)
				}


## The Finn Cycling Index (FCI) is computed as the ratio between cycled throughflows and total system throughput.
## Using the "FCI" function, three different versions of the index are calculated. In the first case ("finn80"),
## diagonal elements of the Leontief matrix are used to assess the amount of material cycling within the ecosystem,
## - (Sii - 1)/Sii - and TST stands for the sum of internal transfers plus imports - TST* = sum(Z) + sum([T]).
## In the Ulanowicz's (1986) version ("ulan86"), TST* is replaced by the more widespread formulation (with total system
## throughput corresponding to the total amount of flows occurring in the system) - TST = sum(Z) + sum([T]) + sum(E) + sum(R).
## Finally, "szyr87" is measured using the diagonal elements of the total dependency matrix, with the TST computed as the sum
## of all the ecosystem flows (Szyrmer & Ulanowicz, 1987).

FCI <- function(inp,inter,outt,diss){
			tot_flows <- TST(inp,inter,outt,diss)
			tot_finn <- sum(inter) + sum(inp)
			xx <- ncol(inter)
			#
			IN <- apply(inter,2,sum) + inp
			#
			LEONTIEF <- leontief(inp,inter)
			if(length(LEONTIEF) == 1)return(NA)
			LEONT <- diag(LEONTIEF)
			#
			DEPENDENCY <- total.dependency(inp,inter,outt,diss)
			if(length(DEPENDENCY) == 1)return(NA)
			#
			DEPEN <- diag(DEPENDENCY)
			#
			finn80 <- rep(0,xx)
			szyr87 <- rep(0,xx)
			#
			for(i in 1:xx){
				if(LEONT[i]!=0)finn80[i] <- IN[i] * (LEONT[i]-1)/LEONT[i]
				#
				szyr87[i] <- IN[i] * DEPEN[i]
				}
			#
			ulan86 <- finn80	
			#
			finn.OK <- rep(0,3) 
			finn.OK[1] <- sum(finn80)/tot_finn		## Finn, 1980
			finn.OK[2] <- sum(ulan86)/tot_flows		## Ulanowicz, 1986
			finn.OK[3] <- sum(szyr87)/tot_flows		## Szyrmer & Ulanowicz, 1987
			#
			return(finn.OK)
			}


## Comprehensive Cycling Index (CCI), the updated and improved version of the FCI. It contains the contribution to cysling given by:
## (a) simple cycles; (b)compound paths and (c) and compound cycles. As CCI is cumbersome to compute, but correlated with FCI, we used
## a constant coefficient (1.142) for its rough estimation from FCI.

CCI <- function(inp,inter,outt,diss){
			k <- 1.142
			finn <- FCI(inp,inter,outt,diss)[2]
			if(is.na(finn) == TRUE)return(NA)
			allesina <- finn * k
			return(allesina)
			}





###########################################################################################################################################
###																	###
### The following functions are inspired to the pubblication of Puccia and Ulanowicz (1990) that analyzed indirect effects in ecosystem	###
### This is not a typical tool of network analysis and is also available at the Ulanowicz's website as a program called IMPACTS		###
### The function "IMA" was introduced by Vasas and Jordan (2006).									###
###																	###
###########################################################################################################################################

## The output of the function "Q.TI" is the Matrix of Direct Trophic Interactions [Q].

Q.TI <- function(inp,inter,outt,nl){
			pf <- partial.feeding(inp,inter)
			#
			tot.l <- length(inp)
			living <- tot.l - nl
			#
			out.l <- apply(inter,1,sum)
			OOO.l <- out.l + outt
			f.ast <- matrix(rep(0,tot.l^2),nrow=tot.l)
			#
			for(i in 1:tot.l)
				for(j in 1:living){
					if(OOO.l[i]!=0)f.ast[i,j] <- inter[i,j]/OOO.l[i]
					}
			#
			pht <- t(f.ast)
			Q <- (pf - pht)
			#
			return(Q)
			}	


## The output of the function "M.TI" is the Matrix of Total Trophic Impacts [M].

M.TI <- function(inp,inter,outt,nl){
				xx <- ncol(inter)
				living <- xx - nl
				#
				I <- diag(1,xx)
				Q <- Q.TI(inp,inter,outt,nl)
				#
				ddd <- det(I-Q)
				if(ddd == 0)return(NA)
				#
				IQinv <- solve(I-Q)
				M <- (IQinv - I)
				#
				return(M)
				}


## The "IMA" function is used for summing the absolute values of effects as measured by partial feeding and partial host matrices.
## In this way, we investigate the absolute effects of species and tropho-species (both positive and negative), estimating their
## global interacting power (key species). Using absolute values does not lead to the loss of information on sign structure.

IMA <- function(inp,inter,outt,nl){
				if(length(M.TI(inp,inter,outt,nl)) == 1)return(NA)
				TTI.m <- M.TI(inp,inter,outt,nl)
				P.TI <- abs(TTI.m)
				IMA.v <- apply(P.TI,1,sum)
				#
				return(IMA.v)
				}





##########################################################
###		     				       ###
### Cycling Analysis - The structure of network cycles ###
###		     				       ###
##########################################################

## These functions ("coppia.bis") and ("LINK.to.NODES") are used into the algorithm for decyclization.
## They do not return any ENA outcome, although implemented within the cycling analysis framework.

coppia.bis <- function(V,M){
		rr <- nrow(M)
		A <- ((V[2]-1)*rr) + V[1]
		for(i in 3:length(V)) A <- c(A,(((V[i]-1)*rr)+V[i-1]))
		#
		return(A)
		}


LINK.to.NODES <- function(VET){
	nnn <- length(VET)
	NOD <- c(VET,NA) %% nnn
	NOD[which(NOD==0)] <- nnn
	mmm <- min(which(is.na(NOD)==TRUE))
	NOD[mmm] <- NOD[1]
	#
	return(NOD)
	}


## Ulanowicz (1983, 1986) devised a procedure to extract cycles from the network and compare their activities with the remaining unidirectional flows.
## He adopted a backtracking algorithm to find all the simple cycles (Mateti and Deo, 1976), increasing its efficiency with a suitable pruning method
## (Knuth, 1973). The whole process can be summarized as follows: (a) a preliminary depth-first search to count the number of cycle arcs incident
## to each node; (b) ordering nodes by decreasing number of incoming cycle arcs (compartments with no cycle arcs are eliminated from further analysis);
## (c) scrutinizing each cycle in search for the smallest cycle arc; (d) identifying if this cycle link is common to other cycles and grouping cycles
## sharing it (nexus); (e) starting backtracking algorithm to delete the nexus and distribute its magnitude among member cycles (in proportion to circuit
## probabilities and magnitudes of links). Finally, the acyclic residual oncethrough flow web (tree, in graph theory) and the aggregated network of
## cycled medium result as two separated parts of the whole network.

cycling.analysis <- function(inp,inter,outt,diss){
	#
	nt.system <- nrow(inter)
	CRCC.FULL <- inter
	CRGG <- partial.feeding(inp,inter)
	CRFF <- partial.host(inter,outt,diss)
	CRCC.ONES <- matrix(rep(0,nt.system^2),nrow=nt.system)
	#
	for(i in 1:nt.system)
		for(j in 1:nt.system)if(inter[i,j]!=0)CRCC.ONES[i,j] <- 1
	#
	ROW.CRCC <- apply(CRCC.ONES,1,sum)
	COL.CRCC <- apply(CRCC.ONES,2,sum)
	#
	SUM.CRCC <- ROW.CRCC + COL.CRCC
	names(SUM.CRCC) <- c(1:nt.system)
	#
	FROM.ROW <- as.list(rep(NA,nt.system))
	for(i in 1:nt.system)FROM.ROW[[i]] <- which(CRCC.ONES[i,]!=0)
	#
	#
	## A) Defining nodes showing both inflows and outflows. The others are excluded
	## (nodes with only input or output are excluded) --> "BLACK" nodes
	## B) Ordering nodes for their total number of links (IN + OUT)
	#
	color <- rep("WHITE",nt.system)
	#
	for(i in 1:nt.system)if(ROW.CRCC[i]==0 | COL.CRCC[i]==0){
							color[i] <- "BLACK"
							SUM.CRCC[i] <- 0
							}
	#
	SRT <- sort(SUM.CRCC,decreasing=TRUE)
	EFF.SRT <- as.integer(names(which(SRT>0)))
	tr <- length(EFF.SRT)
	p <- c(1:nt.system)
	MASSI <- rep(0,nt.system)
	#
	for(i in 1:nt.system){
			if(length(FROM.ROW[[i]])!=0){
				gt <- length(FROM.ROW[[i]])
				#
				W <- which(EFF.SRT==FROM.ROW[[i]][1])
				#
				if(gt > 1){
					for(j in 2:gt)W <- c(W,which(EFF.SRT==FROM.ROW[[i]][j]))
					}
				if(length(W) > 0)MASSI[i] <- max(W)
				}
			}
	#
	#
	## Looking for self loops and deleting them, when found
	#
	self.cyc <- 0
	inds <- 1
	#
	for(i in 1:nt.system)if(CRCC.ONES[i,i]==1){
				self.cyc <- self.cyc + 1
				CRCC.ONES[i,i] <- 0
				inter[i,i] <- 0
				{
				if(inds==1){
					SC <- c(i,i,1,inter[i,i])
					inds <- 2
					}
				else SC <- rbind(SC,c(i,i,1,inter[i,i]))
				}
				#
				}
	#
	#
	## Selecting the more connected node as the starting point for the
	## DFS algorithm.
	#
	cycles.count <- rep(0,nt.system)
	cycles <- 0
	index <- 1
	smile <- 1
	STOPPA <- 0
	one <- 1
	WEAK <- 0
	#
	#
	## Backtracking algorithm including procedure to save and collect weakest arcs
	#
	{
	if(inds!=1)LOO <- 1
	else LOO <- 0
	}
	#
	while(length(color=="WHITE")>1 & sum(CRCC.ONES)!=LOO){
	#
	starting <- EFF.SRT[one]
	color[starting] <- "GRAY"									## Discover vertex starting <- u
	#
	u <- starting
	CIC <- u
	#
	entra <- 0
	#
	#
	counter <- 1
	STOPPA <- 0
	#
	#
	while(STOPPA != 1){
		for(v in counter:tr){
			#
			if(STOPPA == 1)break
			#
			if(entra==2){
				entra <- 0
				aaa <- length(CIC)
				bbb <- length(CIC) - 1
				nott <- CIC[aaa]
				meas <- CIC[bbb]
				#
				if(bbb==0){
					STOPPA <- 1
					break
					}
				#
				{
				if(which(EFF.SRT==nott)<MASSI[meas]){
						 counter <- which(EFF.SRT==nott) + 1
						 CIC <- CIC[-aaa]
						 }
				#
				else {
					CIC <- CIC[-aaa]
					aaa <- length(CIC)
					bbb <- length(CIC) - 1
					#
					if(bbb==0){
						STOPPA <- 1
						break
						}
					#
					nott <- CIC[aaa]
					meas <- CIC[bbb]
					#
					while(which(EFF.SRT==nott)==MASSI[meas]){
									CIC <- CIC[-aaa]
									aaa <- length(CIC)
									bbb <- length(CIC) - 1
									#
									if(bbb==0){
										STOPPA <- 1
										break
										}
									#
									nott <- CIC[aaa]
									meas <- CIC[bbb]
									}
					#
					counter <- which(EFF.SRT==nott) + 1
					CIC <- CIC[-aaa]
					}
				}
				#
				u <- CIC[length(CIC)]
				break
			}
			#
			if(CRCC.ONES[u,EFF.SRT[v]]!=0){							## Examine edge (u,EFF.SRT[v])
				{
				if(color[EFF.SRT[v]]=="WHITE" & length(CIC[CIC==EFF.SRT[v]])==0){	## (u,EFF.SRT[v]) is a tree edge
					#
					CIC <- c(CIC,EFF.SRT[v])
					u <- EFF.SRT[v]
					counter <- 1
					entra <- 0
					break
					}
				else {
					if(color[EFF.SRT[v]]=="GRAY"){					## (u,EFF.SRT[v]) is a back edge
						#
						entra <- 0
						CIC <- c(CIC,EFF.SRT[v])
						names(CIC) <- NULL
						cycles <- cycles + 1
						#
						{
						if(smile==1){				
							LI <- coppia.bis(CIC,inter)
							LL <- length(LI)
							#
							{
							if(index == 1)G.prod <- prod(CRGG[LI])
							#
							else G.prod[index] <- prod(CRGG[LI])
							}
							#
							cycles.count[LL+1] <- cycles.count[LL+1] + 1
							{
							if(LL < nt.system)WEAK <- c(LI,rep(0,nt.system - LL))
							#
							else WEAK <- LI
							}
							index <- index + 1
							}
						else{
							LI <- coppia.bis(CIC,inter)
							LL <- length(LI)
							#
							{
							if(index == 1)G.prod <- prod(CRGG[LI])
							#
							else G.prod[index] <- prod(CRGG[LI])
							}
							#
							cycles.count[LL+1] <- cycles.count[LL+1] + 1
							#
							{
							if(LL < nt.system){
								WEAK <- cbind(WEAK,c(LI,rep(0,nt.system - LL)))
								index <- index + 1
								}
							else	{
								WEAK <- cbind(WEAK,LI)
								index <- index + 1
								}
						     	}
						     }
						}
						smile <- smile + 1
						#
						aaa <- length(CIC)
						bbb <- length(CIC) - 1
						nott <- CIC[aaa]
						meas <- CIC[bbb]
						#
						{
						if(which(EFF.SRT==nott)<MASSI[meas]){
								counter <- which(EFF.SRT==nott) + 1
								CIC <- CIC[-aaa]
								}
						#
						else {	CIC <- CIC[-aaa]
							aaa <- length(CIC)
							bbb <- length(CIC) - 1
							nott <- CIC[aaa]
							meas <- CIC[bbb]
							#					
							while(which(EFF.SRT==nott)==MASSI[meas]){
											CIC <- CIC[-aaa]
											aaa <- length(CIC)
											bbb <- length(CIC) - 1
											#
											if(bbb==0){
												STOPPA <- 1
												break
												}
											#
											nott <- CIC[aaa]
											meas <- CIC[bbb]
											}
							#
							counter <- which(EFF.SRT==nott) + 1
							CIC <- CIC[-aaa]
							#						
							}
						}
						#
						u <- CIC[length(CIC)]
						break
						#
						}
					}
				}
			}
		}	
	#
	entra <- entra + 1
	#
	}
	#
	color[starting] <- "BLACK"
	one <- one + 1
	#
	for(h in 1:nt.system){
		CRCC.ONES[h,starting] <- 0
		CRCC.ONES[starting,h] <- 0
		}
	#
	if(one==length(EFF.SRT))break
	}
	#
	##
	##
	##
	#
	br <- 0
	if(length(which(diag(inter)!=0)) > 0){
				l.sc <- length(which(diag(inter)!=0))
				#
				{
				if(sum(WEAK)==0)
					#
					{
					if(l.sc == 1)
						{
						WEAK <- c(which(diag(inter)!=0)[1],rep(0,nt.system - 1))
						break
						}
					#
					else	
						{
						WEAK <- c(which(diag(inter)!=0)[1],rep(0,nt.system - 1))
						for(i in 2:l.sc)WEAK <- cbind(WEAK,c(which(diag(inter)!=0)[i],rep(0, nt.system - 1)))
						break
						}
					}
				#
				else {
					for(i in 1:l.sc)WEAK <- cbind(WEAK,c(which(diag(inter)!=0)[i],rep(0, nt.system - 1)))
					}
				}
	}
	#
	if(is.matrix(WEAK)==FALSE & sum(WEAK)!=0){
		WEAK <- matrix(WEAK,ncol=1,byrow=TRUE)
		}
	#
	rownames(WEAK) <- NULL
	FULL.CYCLES <- ncol(WEAK)
	#
	##
	##
	##
	#
	CRCC.NA <- matrix(rep(0,nrow(inter)^2),nrow=nrow(inter))
	for(i in 1:nrow(inter))
		for(j in 1:nrow(inter)){
			if(inter[i,j]>0)CRCC.NA[i,j] <- inter[i,j]
			}
	#
	CLASS.cyc <- as.list(NA)
	WEAK.copy <- WEAK
	indaco <- 1
	#
	#
	while(sum(WEAK)!=0){
			#
			CRCC.vec <- as.vector(CRCC.NA)
			CRCC.min <- CRCC.vec
			WEAK.vec <- as.vector(WEAK)
			#
			least <- which.min(CRCC.vec)
			nexus.elements <- which(WEAK.vec==least)
			#
			while(length(nexus.elements)==0){
				CRCC.min[least] <- NA
				least <- which.min(CRCC.min)
				nexus.elements <- which(WEAK.vec==least)
				}
			#
			nexus.vec <- trunc((nexus.elements-1)/nt.system) + 1
			lun <- length(nexus.vec)
			#
			CLASS.cyc[[indaco]] <- matrix(WEAK[,c(nexus.vec)], nrow=lun, byrow=TRUE)
			#
			#
			PRD <- G.prod[c(nexus.vec)]/sum(G.prod[c(nexus.vec)])
			WEAK[,nexus.vec] <- rep(0,nt.system)
			#
			for(i in 1:lun){
				abc <- CLASS.cyc[[indaco]][i,]
				abc.p <- abc[abc > 0]
				CRCC.NA[abc.p] <- CRCC.NA[abc.p] - (CRCC.vec[least] * PRD[i])
				}
			#
				{
				if(indaco==1) nxs.val <- CRCC.NA[least]
				#
				else nxs.val <- c(nxs.val,CRCC.NA[least])
				}
			#
			CRCC.NA[least] <- 0
			indaco <- indaco + 1
		#
	}
	#
	#
	#
	inter <- CRCC.FULL
	#
	if(self.cyc!=0){
		cycles.count[2] <- self.cyc
		FULL.CYCLES <- sum(cycles.count)
		ZER <- matrix(rep(0,nt.system^2),nrow=nt.system)
		diag(ZER) <- diag(inter)
		dia <- which(ZER!=0)
		dia.l <- length(dia)
		WEAK.sl <- matrix(rep(0,(dia.l * nt.system)), nrow = nt.system)
		#
		for(t in 1:dia.l){
			WEAK.sl[,t] <- c(dia[t],rep(0,(nt.system) - 1))
		}
		#
		{
		if(sum(WEAK.copy)!=0)WEAK.copy <- cbind(WEAK.copy,WEAK.sl)
		#
		else WEAK.copy <- WEAK.sl
		}
	}	
	#
	#
	#
	{
	if(sum(FULL.CYCLES)!=0){
		deci <- 1
		#
		WEAK.copy.NAs <- WEAK.copy
		#
		CRCC.LINE <- inter
		WEAK.copy.NAs[which(WEAK.copy.NAs==0)] <- NA
		jump <- 1
		CYCLE.distr <- rep(0,nt.system)
		#
		while(sum(na.omit(WEAK.copy.NAs))!=0){
			#
			if(jump == 1){
				CRCC.NAs <- CRCC.LINE	
				CRCC.NAs[which(CRCC.NAs <= 0)] <- NA
				}
			#
			un.elem <- unique(c(WEAK.copy.NAs))
			un.elem.n <- un.elem[which(is.na(un.elem)==FALSE)]
			#
			nex.min <- which.min(CRCC.NAs)
			nex.an <- which(WEAK.copy.NAs==nex.min)
			#
			{
			if(length(nex.an)!=0){
					rig <- trunc(nex.an/nt.system,1) + 1
					{
					if(deci == 1){
						WHOLE.cyc <- as.list(c(NA))
						WHOLE.FFF <- as.list(c(NA))
						WHOLE.FFF.v <- as.list(c(NA))
						WHOLE.nodes <- as.list(c(NA))
						#
						WHOLE.cyc[[deci]] <- t((WEAK.copy.NAs[,rig]))
						WHOLE.FFF[[deci]] <- matrix(CRFF[WHOLE.cyc[[deci]]], byrow = FALSE, nrow = length(nex.an))
						WHOLE.FFF[[deci]][which(is.na(WHOLE.FFF[[deci]])==TRUE)] <- 1
						WHOLE.FFF.v <- apply(WHOLE.FFF[[deci]],1,prod)
						WHOLE.FFF.s <- sum(WHOLE.FFF.v)
						WHOLE.FFF.r <- WHOLE.FFF.v/WHOLE.FFF.s
						#
						WHOLE.nodes[[deci]] <- matrix(rep(NA,((nt.system + 1)*length(nex.an))),nrow = length(nex.an))
						#
						for(q in 1:length(nex.an)){
							ve <- WHOLE.cyc[[deci]][q,which(is.na(WHOLE.cyc[[deci]][q,])==FALSE)]
							nexi <- CRCC.NAs[nex.min]
							rela <- WHOLE.FFF.r[q]
							dele <- nexi * rela
							CRCC.LINE[ve] <- CRCC.LINE[ve] - (dele)
							#
							WHOLE.nodes[[deci]][q,] <- LINK.to.NODES(WHOLE.cyc[[deci]][q,])
							#
							ll.cc <- length(which(is.na(WHOLE.cyc[[deci]][q,])==FALSE))
							CYCLE.distr[ll.cc] <- CYCLE.distr[ll.cc] + (dele * ll.cc)
							#
						}
						#
						WEAK.ARCS <- inter[nex.min]
						#
						aj <- nex.min%%nt.system
						if(aj == 0) aj <- nt.system
						WHOLE.nexus <- c(aj,trunc(nex.min/nt.system)+1)
						#
						deci <- deci + 1
					}
					#
					else	{
						WHOLE.cyc[[deci]] <- t((WEAK.copy.NAs[,rig]))
						WHOLE.FFF[[deci]] <- matrix(CRFF[WHOLE.cyc[[deci]]], byrow = FALSE, nrow = length(nex.an))
						WHOLE.FFF[[deci]][which(is.na(WHOLE.FFF[[deci]])==TRUE)] <- 1
						WHOLE.FFF.v <- apply(WHOLE.FFF[[deci]],1,prod)
						WHOLE.FFF.s <- sum(WHOLE.FFF.v)
						WHOLE.FFF.r <- WHOLE.FFF.v/WHOLE.FFF.s
						#
						WHOLE.nodes[[deci]] <- matrix(rep(NA,((nt.system + 1)*length(nex.an))),nrow = length(nex.an))
						#
						for(q in 1:length(nex.an)){
							ve <- WHOLE.cyc[[deci]][q,which(is.na(WHOLE.cyc[[deci]][q,])==FALSE)]
							nexi <- CRCC.NAs[nex.min]
							rela <- WHOLE.FFF.r[q]
							dele <- nexi * rela
							CRCC.LINE[ve] <- CRCC.LINE[ve] - (dele)
							#
							WHOLE.nodes[[deci]][q,] <- LINK.to.NODES(WHOLE.cyc[[deci]][q,])
							#
							ll.cc <- length(which(is.na(WHOLE.cyc[[deci]][q,])==FALSE))
							CYCLE.distr[ll.cc] <- CYCLE.distr[ll.cc] + (dele * ll.cc)
							#
						}
						#
						WEAK.ARCS <- c(WEAK.ARCS,inter[nex.min])
						#
						aj <- nex.min%%nt.system
						if(aj == 0) aj <- nt.system
						WHOLE.nexus <- rbind(WHOLE.nexus,c(aj,trunc(nex.min/nt.system)+1))
						#
						deci <- deci + 1
						#
						}
					}
					jump <- 1
					WEAK.copy.NAs <- as.matrix(WEAK.copy.NAs[,-rig])
					CRCC.NAs[nex.min] <- NA
				}
			#
			else 	{
				CRCC.NAs[nex.min] <- NA
				jump <- 2
				}
			}
		}
		#
		#
		#
		M.LINEAR <- CRCC.LINE
		diag(M.LINEAR) <- 0
		M.LINEAR[which(M.LINEAR < 1e-12)] <- 0
		M.CYCLES <- inter - CRCC.LINE
		names(WEAK.ARCS) <- NULL
		rownames(WHOLE.nexus) <- NULL
		N.CYCLE.distr <- CYCLE.distr/TST(inp,inter,outt,diss)
		cycles.count <- cycles.count[-1]
		#
		CY <- as.list(rep(NA,9))
		CY[[1]] <- FULL.CYCLES					## Total number of cycles removed
		CY[[2]] <- cycles.count					## Cycle analysis (number of cycles classified by cycle length)
		CY[[3]] <- WHOLE.nodes					## Full cycle analysis (classified by nexus)
		CY[[4]] <- WEAK.ARCS					## Link strength of the removed nexus
		CY[[5]] <- WHOLE.nexus					## Weak arcs representing nexus (ranked by removal order)
		CY[[6]] <- CYCLE.distr					## Cycle distribution
		CY[[7]] <- N.CYCLE.distr				## Normalized distribution
		CY[[8]] <- M.LINEAR					## Residual flows
		CY[[9]] <- M.CYCLES					## Aggregated cycles
		#
		return(CY)
		}
	#
	else{
		CY <- as.list(rep(NA,2))
		CY[[1]] <- 0
		CY[[2]] <- inter
		return(CY)
		}
	}
}





#####################################################################################################################################################
#####################################################################################################################################################


toTXT <- function(cnames,l,inp,inter,outt,diss,bi,filename_out){
		xx <- length(cnames)
		nl <- xx - l
		{
		if(nl!=0)LVSNL.v <- c(rep(1,l),rep(-1,nl))
		#
		else LVSNL.v <- rep(1,l)
		}
		#
		ZZZ <- vec.ena(inp,xx)
		TTT <- matr.ena(inter,xx)
		EEE <- vec.ena(outt,xx)
		RRR <- vec.ena(diss,xx)
		#
		M.TXT <- cbind(LVSNL.v,bi,ZZZ,EEE,RRR,TTT)
		colnames(M.TXT) <- c("Living","Biomasses","Imports","Exports","Respirations",cnames)
		#
		write.table(M.TXT, file = filename_out, row.names = FALSE, quote = TRUE)
		return(M.TXT)
		}



ENA.c <- function(filename_in, filename_out, routes, output_form){
	#
	#					## output_form <- c("text","files")
	#
	if(length(which(output_form == 0)) == 1)print("No stored results", quote = FALSE)
	#
	if(length(which(routes == "all")) == 1){
			routes <- c("cd","io","tp","ei","ca")
			}
	#
	#					# "cd" = Partial host and partial feeding matrices. Dependency and contribution coefficients
	#					# "io" = Input/Output analysis
	#					# "tp" = Trophic aggregation (after applyind de-cyclization algorithm to living exchanges)
	#					# "ei" = Whole system indices and FCI in case of cycling analyses performed
	#					# "ca" = Cycling analysis for the whole system
	#
	step.1 <- df.ena(filename_in)	
	#
	TN <- step.1[[1]]			## total number of compartments (living + non-living)
	LN <- step.1[[2]]			## living compartments
	NLN <- step.1[[3]]			## non-living compartments
	#
	ZZZ <- step.1[[6]]
	EEE <- step.1[[7]]
	RRR <- step.1[[8]]
	TTT <- step.1[[9]]
	#
	step.1[[4]] <- as.vector(step.1[[4]])
	rownames(step.1[[9]]) <- step.1[[4]]
	rownames(step.1[[10]]) <- step.1[[4]][1:LN]
	names(step.1[[5]]) <- step.1[[4]]
	names(step.1[[6]]) <- step.1[[4]]
	names(step.1[[7]]) <- step.1[[4]]
	names(step.1[[8]]) <- step.1[[4]]
	#
	#
	if(length(which(output_form == "text")) == 1){
	cat("NUMBER OF COMPARTMENTS IS",TN, file = filename_out, append = TRUE)
	cat("\nNUMBER OF LIVING COMPARTMENTS IS",LN, file = filename_out, append = TRUE)
	cat("\n\nLIVING COMPARTMENTS\n\n", file = filename_out, append = TRUE)
	#
	for(i in 1:LN){
		if(i < 10)cat("   ",i," ", file = filename_out, append = TRUE, sep = "")
		#
		if(i >= 10 & i < 100)cat("  ",i," ", file = filename_out, append = TRUE, sep = "")
		#
		if(i >= 100)cat(" ",i," ", file = filename_out, append = TRUE, sep = "")
		#
		write(step.1[[4]][i], file = filename_out, append = TRUE)
		}
	if(NLN!=0){
		cat("\nNON-LIVING COMPARTMENTS\n\n", file = filename_out, append = TRUE)
		for(i in 1:NLN){
			if(i < 10)cat("   ",LN+i," ", file = filename_out, append = TRUE, sep = "")
			#
			if(i >= 10 & i < 100)cat("  ",LN+i," ", file = filename_out, append = TRUE, sep = "")
			#
			if(i >= 100)cat(" ",LN+i," ", file = filename_out, append = TRUE, sep = "")
			#
			write(step.1[[4]][LN+i], file = filename_out, append = TRUE)
			}
	}
	cat("\n\nINPUT VECTOR\n\n","   ",
			formatC(c(1:TN),digits=1,width=10), "\n\n",
			"      ",formatC(ZZZ,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
	cat("\n\n\nEXPORT VECTOR\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(EEE,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
	cat("\n\n\nRESPIRATION VECTOR\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(RRR,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
	cat("\n\n\nEXCHANGE MATRIX\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n", file = filename_out, append = TRUE, sep="")
			for(i in 1:TN){
				for(j in 1:TN){
					if(j!=TN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(TTT[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
							#
							if(i >=10 & i < 100)cat("  ",i,"  ",formatC(TTT[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
							#
							if(i > 100)cat(" ",i,"  ",formatC(TTT[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
							}
						#
						else cat(formatC(TTT[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(TTT[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
		}
	#
	step.2 <- bigmatrNOnames.ena(ZZZ,TTT,EEE,RRR)
	INP <- apply(step.2,2,sum)[2:(TN+1)]
	OUT <- apply(step.2,1,sum)[2:(TN+1)]
	DIFF_I <- INP - OUT
	DIFF_O <- OUT - INP
	DIFF_IP <- 100*(INP - OUT)/INP
	DIFF_OP <- 100*(OUT - INP)/OUT
	BAL <- cbind(INP,OUT,DIFF_I,DIFF_O,DIFF_IP,DIFF_OP)
	#
	if(length(which(output_form == "text")) == 1){
	{
	if(length(which(abs(BAL[,5]) > 0.0001))==0){
		cat("\n\nCOMPARTMENTAL THROUGHPUTS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10), "\n\n",
			"      ",formatC(INP,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep ="")
		}
	else	{
		cat("\n\n\nWARNING! - NETWORK IS NON STEADY-STATE!\n\nCOMPARISON OF THROUGHPUTS\n\n", file = filename_out, append = TRUE)
		cat("        INPUTS     OUTPUTS    IMBALANCE\n\n", file = filename_out, append = TRUE)
			for(i in 1:TN){
				for(j in 1:3){
					#
					if(j==1){
						if(i < 10)cat("   ",i,"   ",formatC(BAL[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
						#
						if(i >= 10 & i < 100)cat("  ",i,"   ",formatC(BAL[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
						#
						if(i >= 100)cat(" ",i,"   ",formatC(BAL[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
						}
					#
					if(j!=1 & j!=3)cat(" ",formatC(BAL[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
					#
					if(j==3)cat(" ",formatC(BAL[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE, sep = "")
				}
			}
		}
	}
	#
	}
	#
	#
	#
	##
	##
	##
	#
	#
	#
	if(length(which(routes == "cd"))==1){
	#
	if(length(which(output_form == "text")) == 1)cat("\n\n\n\n\n*** STRUCTURE ANALYSES ***\n", file = filename_out, append = TRUE)
	#
	GGG <- partial.feeding(ZZZ,TTT)
	FFF <- partial.host(TTT,EEE,RRR)
	TCTC <- total.contribution(ZZZ,TTT,EEE,RRR)
	TDTD <- total.dependency(ZZZ,TTT,EEE,RRR)
	#
	if(length(which(output_form == "text")) == 1){
	cat("\n\nPARTIAL FEEDING MATRIX\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n", file = filename_out, append = TRUE, sep = "")
			for(i in 1:TN){
				for(j in 1:TN){
					if(j!=TN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(GGG[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 10 & i < 100)cat("  ",i,"  ",formatC(GGG[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 100)cat(" ",i,"  ",formatC(GGG[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							}
						#
						else cat(formatC(GGG[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(GGG[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
	#
	cat("\n\n\nPARTIAL HOST MATRIX\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n", file = filename_out, append = TRUE, sep = "")
			for(i in 1:TN){
				for(j in 1:TN){
					if(j!=TN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(FFF[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 10 & i < 100)cat("  ",i,"  ",formatC(FFF[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 100)cat(" ",i,"  ",formatC(FFF[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							}
						#
						else cat(formatC(FFF[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(FFF[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
	#
	if(length(TCTC) > 1){
	cat("\n\n\nTOTAL CONTRIBUTION COEFFICIENTS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n", file = filename_out, append = TRUE, sep = "")
			for(i in 1:TN){
				for(j in 1:TN){
					if(j!=TN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(TCTC[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 10 & i < 100)cat("  ",i,"  ",formatC(TCTC[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 100)cat(" ",i,"  ",formatC(TCTC[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							}
						#
						else cat(formatC(TCTC[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(TCTC[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
		}
	#
	if(length(TDTD) > 1){
	cat("\n\n\nTOTAL DEPENDENCY COEFFICIENTS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n", file = filename_out, append = TRUE, sep = "")
			for(i in 1:TN){
				for(j in 1:TN){
					if(j!=TN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(TDTD[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 10 & i < 100)cat("  ",i,"  ",formatC(TDTD[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 100)cat(" ",i,"  ",formatC(TDTD[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							}
						#
						else cat(formatC(TDTD[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(TDTD[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
		}
	if(length(which(output_form == "text")) == 1)cat("\n\n", file = filename_out, append = TRUE, sep = "")
	}
	#
	}
	#
	if(length(which(routes == "io"))==1){
	options(warn = -1)
	IAN <- input.analysis(ZZZ,TTT,EEE,RRR)
	OAN <- output.analysis(ZZZ,TTT,EEE)
	step.7 <- as.list(rep(NA,length(IAN)))
	step.8 <- as.list(rep(NA,length(OAN)))
	#
	if(is.na(IAN) != TRUE){
	for(k in 1:length(IAN)){
		#
		IAN.1 <- which(IAN[[k]][[1]]!=0)		
		IAN.2 <- rep(0,TN)
		IAN.2[IAN.1] <- ZZZ[IAN.1]
		IAN.3 <- IAN[[k]][[3]] * ZZZ[IAN.1]
		IAN.4 <- IAN[[k]][[4]] * ZZZ[IAN.1]
		IAN.5 <- IAN[[k]][[2]] * ZZZ[IAN.1]
		#
		if(length(which(output_form == "text")) == 1){
		cat("\n\n\nANALYSIS OF INPUT ",IAN.1, file = filename_out, append = TRUE, sep = "")
		cat("\n\n\nINPUT VECTOR\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(IAN.2,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
		cat("\n\n\nEXPORT VECTOR\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(IAN.3,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
		cat("\n\n\nRESPIRATIONS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(IAN.4,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
		cat("\n\n\nEXCHANGE MATRIX\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n", file = filename_out, append = TRUE, sep="")
			for(i in 1:TN){
				for(j in 1:TN){
					if(j!=TN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(IAN.5[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 10 & i < 100)cat("  ",i,"  ",formatC(IAN.5[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 100)cat(" ",i,"  ",formatC(IAN.5[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							}
						#
						else cat(formatC(IAN.5[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(IAN.5[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
		}
		step.7[[k]] <- bigmatr.ena(step.1[[4]],IAN.2,IAN.5,IAN.3,IAN.4)
		}
	}
	#
	if(length(which(output_form == "text")) == 1)cat("\n\n", file = filename_out, append = TRUE, sep = "")
	#
	if(is.na(OAN) != TRUE){
	for(k in 1:length(OAN)){
		#
		OAN.1 <- which(OAN[[k]][[3]]!=0)		
		OAN.2 <- OAN[[k]][[1]] * EEE[OAN.1]
		OAN.3 <- rep(0,TN)
		OAN.3[OAN.1] <- EEE[OAN.1]
		OAN.4 <- OAN[[k]][[4]] * EEE[OAN.1]
		OAN.5 <- OAN[[k]][[2]] * EEE[OAN.1]
		#
		if(length(which(output_form == "text")) == 1){
		cat("\n\n\nANALYSIS OF OUTPUT ",OAN.1, file = filename_out, append = TRUE, sep = "")
		cat("\n\n\nINPUT VECTOR\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(OAN.2,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
		cat("\n\n\nEXPORT VECTOR\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(OAN.3,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
		cat("\n\n\nRESPIRATIONS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(OAN.4,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
		cat("\n\n\nEXCHANGE MATRIX\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n", file = filename_out, append = TRUE, sep="")
			for(i in 1:TN){
				for(j in 1:TN){
					if(j!=TN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(OAN.5[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 10 & i < 100)cat("  ",i,"  ",formatC(OAN.5[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 100)cat(" ",i,"  ",formatC(OAN.5[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							}
						#
						else cat(formatC(OAN.5[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(OAN.5[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
		}
		step.8[[k]] <- bigmatr.ena(step.1[[4]],OAN.2,OAN.5,OAN.3,OAN.4)
		}
	}
	#
	if(length(which(output_form == "files")) == 1){
	#
	if(is.na(IAN) != TRUE){
	NAME.IA <- paste("InputAnalysis_",filename_out,sep="")
	for(k in 1:length(IAN)){		
		write.table(step.7[[k]], file = NAME.IA, row.names = FALSE, quote = FALSE, append = TRUE)
		}
		#
		}
	#
	if(is.na(OAN) != TRUE){
	NAME.OA <- paste("OutputAnalysis_",filename_out,sep="")
	for(k in 1:length(OAN)){
		write.table(step.8[[k]], file = NAME.OA, row.names = FALSE, quote = FALSE, append = TRUE)
		}
		#
		}
	}
	#
	}
	#
	#
	##
	##
	##
	#
	#
	#
	options(warn = 0)
	if(length(which(routes == "tp"))==1){
	#
	ZZZ.L <- ZZZ[1:LN]
	TTT.L <- TTT[c(1:LN),c(1:LN)]
	EEE.L <- EEE[1:LN]
	RRR.L <- RRR[1:LN]
	#
	step.9 <- cycling.analysis(ZZZ.L,TTT.L,EEE.L,RRR.L)		## Full cycling analysis (living)
	#
	if(length(which(output_form == "text")) == 1)cat("\n\n\n\nSEARCH FOR FEEDING CYCLES\n", file = filename_out, append = TRUE)
	{
	if(step.9[[1]]==0){
		if(length(which(output_form == "text")) == 1)cat("\nNO CYCLES DETECTED", file = filename_out, append = TRUE)
		RFL.c <- step.9[[2]]
		FCI.c <- 0
		}
	#
	else	{
		TOT.c <- step.9[[1]]		## Total number of living cycles
		COM.c <- step.9[[2]]		## Total living cycles classified by their length
		EAC.c <- step.9[[3]]		## Full cycle analysis (living)
		NEX.c <- step.9[[4]]		## Nexus (living)
		WAR.c <- step.9[[5]]		## Weak arc intensity (living)
		CYD.c <- step.9[[6]]		## Living cycle distributions
		NCD.c <- step.9[[7]]		## Normalized living cycle distributions
		RFL.c <- step.9[[8]]		## Residual living flows
		AFL.c <- step.9[[9]]		## Aggregated living cycles
		#
		if(length(WAR.c)==2)WAR.c <- matrix(c(WAR.c),nrow=1)
		#
		FCI.c <- sum(AFL.c)/TST(ZZZ,TTT,EEE,RRR)		## Finn cycling index for living compartments
		#
		if(length(which(output_form == "text")) == 1){
		#
		all.cyc <- 1
		#
		for(k in 1:nrow(WAR.c)){
			for(q in 1:nrow(EAC.c[[k]])){
				if(q == 51){
					cat("\n       . . .", file = filename_out, append = TRUE, sep="")
					all.cyc <- all.cyc + (nrow(EAC.c[[k]]) - 50)
					break
					}
					
				#
				if(q <= 50){
					t <- 1
					while(is.na(EAC.c[[k]][q,t])==FALSE){
								if(t==1){
									cat("\n     ", all.cyc, ".   ",EAC.c[[k]][q,t],"-", file = filename_out, append = TRUE, sep="")
									t <- t + 1
									all.cyc <- all.cyc + 1
									}
								else	{
									cat(EAC.c[[k]][q,t],"-", file = filename_out, append = TRUE, sep="")
									t <- t + 1
									}
								}
					}
				}
			cat("\nNEXUS ", k, " CONSISTS OF ", nrow(EAC.c[[k]]), " CYCLES. WEAK ARC (", WAR.c[k,1], ", ", WAR.c[k,2], ") = ", 
			formatC(NEX.c[k],digits=3,width=10,format="e"), "\n", file = filename_out, append = TRUE, sep="")
			}
		cat("\n\n\nA TOTAL OF", TOT.c, "CYCLES REMOVED", file = filename_out, append = TRUE)
		cat("\n\n\n\nFINN CYCLING INDEX IS ", formatC(FCI.c,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
		cat("\n\n\n\nCYCLE-FREE FEEDING TRANSFERS\n","   ",
		formatC(c(1:LN),digits=1,width=10),"\n", file = filename_out, append = TRUE, sep = "")
			for(i in 1:LN){
				for(j in 1:LN){
					if(j!=LN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(RFL.c[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 10 & i < 100)cat("  ",i,"  ",formatC(RFL.c[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 100)cat(" ",i,"  ",formatC(RFL.c[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							}
						#
						else cat(formatC(RFL.c[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(RFL.c[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
		#
		}
		#
		}
	}
	#
	#
	#
	##
	##
	##
	#
	#
	#
	step.11 <- TP.CTA(ZZZ,TTT,EEE,RRR,NLN)
	#
	LT.cta <- step.11[[1]]				## Lindeman Transformation Matrix
	TL.cta <- step.11[[2]]				## Effective Trophic Levels of each species
	CE.cta <- step.11[[3]]				## Canonical Exports
	CR.cta <- step.11[[4]]				## Canonical Respirations
	GC.cta <- step.11[[5]]				## The Grazing Chain
	DP.cta <- step.11[[6]]				## Returns to Detrital Pool
	DT.cta <- step.11[[7]]				## Detritivory
	ID.cta <- step.11[[8]]				## Input to Detrital Pool
	CD.cta <- step.11[[9]]				## Circulation within Detrital Pool
	LS.cta <- step.11[[10]]				## Lindeman Spine
	TE.cta <- step.11[[11]]				## Trophic Efficiencies
	TE.cta <- c(TE.cta,0)
	MI.cta <- step.11[[12]]				## Detecting Migratory Imports
	if(sum(MI.cta!=0)){
		CI.cta <- step.11[[13]]			## Canonical Distribution of inputs
		}					## (only in case of migrations)
	#
	#
	if(length(which(output_form == "text")) == 1){
	cat("\n\n\n\n\n*** LINDEMAN TROPHIC AGGREGATIONS ***\n\n", file = filename_out, append = TRUE)
	#
	if(sum(MI.cta!=0)){
		cat("\nTHE FOLLOWING INPUTS ARE SUSPECTED TO BE MIGRATORY INPUTS.", file = filename_out, append = TRUE)
		cat("\nTHEY WILL BE TREATED AS NON-PRIMARY INFLOWS.\n", file = filename_out, append = TRUE)
 		for(z in 1:length(MI.cta)){
 			cat("\nMIGRATORY INPUT TO ", MI.cta[z], " = ", 
 			formatC(ZZZ[MI.cta[z]],digits=4,width=10,format="e"), file = filename_out, append = TRUE, sep ="")
 			}
 		cat("\n\n", file = filename_out, append = TRUE)
 		}
	#
	cat("\nLINDEMAN TRANSFORMATION MATRIX\n\n","   ",
		formatC(c(1:LN),digits=1,width=10),"\n\n", file = filename_out, append = TRUE, sep = "")
			for(i in 1:LN){
				for(j in 1:LN){
					if(j!=LN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(LT.cta[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 10 & i < 100)cat("  ",i,"  ",formatC(LT.cta[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 100)cat(" ",i,"  ",formatC(LT.cta[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							}
						#
						else cat(formatC(LT.cta[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(LT.cta[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
	#
	cat("\n\nEFFECTIVE TROPHIC LEVELS OF EACH SPECIES\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(TL.cta,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
	#
	if(sum(MI.cta!=0)){
		cat("\n\n\nCANONICAL DISTRIBUTION OF INPUTS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(CI.cta,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
		}
	#
	cat("\n\n\nCANONICAL EXPORTS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(CE.cta,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
	#
	cat("\n\n\nCANONICAL RESPIRATIONS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(CR.cta,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
	#
	cat("\n\n\nTHE GRAZING CHAIN\n\n","   ",
			formatC(c(1:LN),digits=1,width=10),"\n\n",
			"      ",formatC(GC.cta,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
	#
	cat("\n\n\nRETURNS TO DETRITAL POOL\n\n","   ",
			formatC(c(1:LN),digits=1,width=10),"\n\n",
			"      ",formatC(DP.cta,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
	#
	cat("\n\n\nDETRITIVORY =", formatC(DT.cta,digits=4,width=10,format="e"), "\n", file = filename_out, append = TRUE)
	#
	cat("\nINPUT TO DETRITAL POOL =", formatC(ID.cta,digits=4,width=10,format="e"), "\n", file = filename_out, append = TRUE)
	#
	cat("\nCIRCULATION WITHIN DETRITAL POOL =", formatC(CD.cta,digits=4,width=10,format="e"), "\n\n", file = filename_out, append = TRUE)
	#
	cat("\nLINDEMAN SPINE\n\n","   ",
			formatC(c(1:LN),digits=1,width=10),"\n\n",
			"      ",formatC(LS.cta,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
	#
	cat("\n\n\nTROPHIC EFFICIENCIES\n\n","   ",
			formatC(c(1:LN),digits=1,width=10),"\n\n",
			"      ",formatC(TE.cta,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
	}
	#
	}
	#
	#
	#
	##
	##
	##
	#
	#
	#
	if(length(which(routes == "ei"))==1){
	#
	tst.i <- TST(ZZZ,TTT,EEE,RRR)
	ami.i <- AMI(ZZZ,TTT,EEE,RRR)
	dc.i <- DC(ZZZ,TTT,EEE,RRR)
	amidc.i <- AMI.max(ZZZ,TTT,EEE,RRR)
	asc.i <- ASC(ZZZ,TTT,EEE,RRR)
	oi.i <- OVERHEAD.imports(ZZZ,TTT,EEE,RRR)
	oe.i <- OVERHEAD.exports(ZZZ,TTT,EEE,RRR)
	or.i <- OVERHEAD.resp(ZZZ,TTT,EEE,RRR)
	r.i <- REDUNDANCY(ZZZ,TTT,EEE,RRR)
	o.i <- OVERHEAD(ZZZ,TTT,EEE,RRR)
	idc.i <- internal.DC(ZZZ,TTT,EEE,RRR)
	ia.i <- internal.ASC(ZZZ,TTT,EEE,RRR)
	oc.i <- O.C(ZZZ,TTT,EEE,RRR)
	ic.i <- I.C(TTT)
	fwc.i <- FW.C(TTT,NLN)
	#
	#
	if(length(which(output_form == "text")) == 1){
	cat("\n\n\n\n\n\n*** INFORMATION INDICES ***\n\n", file = filename_out, append = TRUE)
	#
	cat("\nTOTAL SYSTEM THROUGHPUT =", formatC(tst.i,digits=4,width=10,format="e"), "\n\n", file = filename_out, append = TRUE)
	#
	cat("\nDEVELOPMENT CAPACITY = ", formatC(dc.i,digits=4,width=10,format="e"), "\n", file = filename_out, append = TRUE, sep="")
	cat("\nASCENDENCY = ", formatC(asc.i[1],digits=4,width=10,format="e"), "              (",
		formatC(asc.i[2],digits=3,width=0,format="f"), ")\n", file = filename_out, append = TRUE, sep="")
	cat("\nOVERHEAD ON IMPORTS = ", formatC(oi.i[1],digits=4,width=10,format="e"), "     (",
		formatC(oi.i[2],digits=3,width=0,format="f"), ")\n", file = filename_out, append = TRUE, sep="")
	cat("\nOVERHEAD ON EXPORTS = ", formatC(oe.i[1],digits=4,width=10,format="e"), "     (",
		formatC(oe.i[2],digits=3,width=0,format="f"), ")\n", file = filename_out, append = TRUE, sep="")
	cat("\nDISSIPATIVE OVERHEAD = ", formatC(or.i[1],digits=4,width=10,format="e"), "    (",
		formatC(or.i[2],digits=3,width=0,format="f"), ")\n", file = filename_out, append = TRUE, sep="")
	cat("\nREDUNDANCY = ", formatC(r.i[1],digits=4,width=10,format="e"), "              (",
		formatC(r.i[2],digits=3,width=0,format="f"), ")\n", file = filename_out, append = TRUE, sep="")
	#
	cat("\n\nINTERNAL CAPACITY = ", formatC(idc.i,digits=4,width=10,format="e"), "\n", file = filename_out, append = TRUE, sep="")
	cat("\nINTERNAL ASCENDENCY = ", formatC(ia.i[1,1],digits=4,width=10,format="e"), "     (",
		formatC(ia.i[1,2],digits=3,width=0,format="f"), ")\n", file = filename_out, append = TRUE, sep="")
	cat("\nREDUNDANCY = ", formatC(r.i[1],digits=4,width=10,format="e"), "              (",
		formatC((1-ia.i[1,2]),digits=3,width=0,format="f"), ")\n", file = filename_out, append = TRUE, sep="")
	#
	#
	cat("\n\n\n*** CONNECTANCE INDICES ***\n\n", file = filename_out, append = TRUE)
	#
	cat("\nOVERALL CONNECTANCE =", formatC(oc.i,digits=3,width=10,format="f"), "\n", file = filename_out, append = TRUE)
	cat("\nINTERCOMPARTMENTAL CONNECTANCE =", formatC(ic.i,digits=3,width=10,format="f"), "\n", file = filename_out, append = TRUE)
	cat("\nFOODWEB CONNECTANCE =", formatC(fwc.i,digits=3,width=10,format="f"), "\n", file = filename_out, append = TRUE)
	#
	}
	#
	}
	#
	#
	#
	##
	##
	##
	#
	#
	#
	if(length(which(routes == "ca"))==1){
	#
	step.12 <- cycling.analysis(ZZZ,TTT,EEE,RRR)
	#
	if(length(which(output_form == "text")) == 1)cat("\n\n\n\n\n\n*** FULL CYCLE ANALYSIS ***\n\n", file = filename_out, append = TRUE)
	{
	if(step.12[[1]]==0){
		if(length(which(output_form == "text")) == 1)cat("\nNO CYCLES DETECTED", file = filename_out, append = TRUE)
		RFL.wc <- step.12[[2]]
		FCI.wc <- 0
		}
	#
	else	{
		TOT.wc <- step.12[[1]]		## Total number of cycles
		COM.wc <- step.12[[2]]		## Total cycles classified by their length
		EAC.wc <- step.12[[3]]		## Full cycle analysis
		NEX.wc <- step.12[[4]]		## Nexus
		WAR.wc <- step.12[[5]]		## Weak arc intensity
		CYD.wc <- step.12[[6]]		## Cycle distributions
		NCD.ww <- step.12[[7]]		## Normalized cycle distributions
		RFL.wc <- step.12[[8]]		## Residual flows
		AFL.wc <- step.12[[9]]		## Aggregated cycles
		#
		FCI.wc <- sum(AFL.wc)/TST(ZZZ,TTT,EEE,RRR)
		#
		if(length(WAR.wc)==2)WAR.wc <- matrix(c(WAR.wc),nrow=1)
		#
		if(length(which(output_form == "text")) == 1){
		#
		all.cyc.w <- 1
		#
		for(k in 1:nrow(WAR.wc)){
			for(q in 1:nrow(EAC.wc[[k]])){
				if(q == 51){
					cat("\n       . . .", file = filename_out, append = TRUE, sep="")
					all.cyc.w <- all.cyc.w + (nrow(EAC.wc[[k]]) - 50)
					break
					}
				#
				if(q <= 50){
					t <- 1
					while(is.na(EAC.wc[[k]][q,t])==FALSE){
								if(t==1){
									cat("\n     ", all.cyc.w, ".   ",EAC.wc[[k]][q,t],"-", file = filename_out, append = TRUE, sep="")
									t <- t + 1
									all.cyc.w <- all.cyc.w + 1
									}
								else	{
									cat(EAC.wc[[k]][q,t],"-", file = filename_out, append = TRUE, sep="")
									t <- t + 1
									}
								}
					}
				}
			cat("\nNEXUS ", k, " CONSISTS OF ", nrow(EAC.wc[[k]]), " CYCLES. WEAK ARC (", WAR.wc[k,1], ", ", WAR.wc[k,2], ") = ", 
			formatC(NEX.wc[k],digits=3,width=10,format="e"), "\n", file = filename_out, append = TRUE, sep="")
			}
		cat("\n\nA TOTAL OF", TOT.wc, "CYCLES REMOVED\n", file = filename_out, append = TRUE)
		#
		cat("\n\n\nCYCLE DISTRIBUTIONS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(CYD.wc,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
		#
		cat("\n\n\nNORMALIZED DISTRIBUTIONS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(NCD.ww,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
		#
		cat("\n\n\n\nFINN CYCLING INDEX IS ", formatC(FCI.wc,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
		cat("\n\n\n\nRESIDUAL FLOWS\n\n","   ",
		formatC(c(1:TN),digits=1,width=10),"\n", file = filename_out, append = TRUE, sep = "")
			for(i in 1:TN){
				for(j in 1:TN){
					if(j!=TN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(RFL.wc[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 10 & i < 100)cat("  ",i,"  ",formatC(RFL.wc[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 100)cat(" ",i,"  ",formatC(RFL.wc[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							}
						#
						else cat(formatC(RFL.wc[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(RFL.wc[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
		#
		cat("\n\n\nAGGREGATED CYCLES\n\n","   ",
		formatC(c(1:TN),digits=1,width=10),"\n", file = filename_out, append = TRUE, sep = "")
			for(i in 1:TN){
				for(j in 1:TN){
					if(j!=TN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(AFL.wc[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 10 & i < 100)cat("  ",i,"  ",formatC(AFL.wc[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 100)cat(" ",i,"  ",formatC(AFL.wc[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							}
						#
						else cat(formatC(AFL.wc[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(AFL.wc[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
		#
		}
	#
	}
	#
	}
	#
	}
	output_ENA.c <- as.list(rep(NA,13))
	#
	#
	##
	##
	#
	#
	output_ENA.c[[1]] <- step.1
	names(output_ENA.c[[1]]) <- c("Compartments","Living","NonLiving","Names",
	"Biomasses","Input","Export","Respiration","Exchange","LivingExchange")
	#
	output_ENA.c[[2]] <- bigmatr.ena(step.1[[4]],ZZZ,TTT,EEE,RRR)
	if(length(which(output_form == "files")) == 1){
		NAME.BM <- paste("BigMatrix_",filename_out,sep="")
		write.table(output_ENA.c[[2]], file = NAME.BM, row.names = FALSE, quote = FALSE)
		}
	#
	rownames(BAL) <- step.1[[4]]
	colnames(BAL) <- c("Inp","Out","Diff(I)","Diff(O)","Diff(I%)","Diff(O%)")
	output_ENA.c[[3]] <- BAL
	if(length(which(output_form == "files")) == 1){
		NAME.SC <- paste("SteadyStateCheck_",filename_out,sep="")
		write.table(BAL, file = NAME.SC, row.names = FALSE, quote = FALSE)
		}
	#
	if(length(which(routes == "cd"))==1){
	#
	rownames(GGG) <- step.1[[4]]
	colnames(GGG) <- step.1[[4]]
	output_ENA.c[[4]] <- GGG
	if(length(which(output_form == "files")) == 1){
		NAME.PF <- paste("PartialFeeding_",filename_out,sep="")
		write.table(output_ENA.c[[4]], file = NAME.PF, row.names = FALSE, quote = FALSE)
		}
	#
	rownames(FFF) <- step.1[[4]]
	colnames(FFF) <- step.1[[4]]
	output_ENA.c[[5]] <- FFF
	if(length(which(output_form == "files")) == 1){
		NAME.PH <- paste("PartialHost_",filename_out,sep="")
		write.table(output_ENA.c[[5]], file = NAME.PH, row.names = FALSE, quote = FALSE)
		}
	#
	rownames(TCTC) <- step.1[[4]]
	colnames(TCTC) <- step.1[[4]]
	output_ENA.c[[6]] <- TCTC
	if(length(which(output_form == "files")) == 1 & length(TCTC) > 1){
		NAME.TC <- paste("Contribution_",filename_out,sep="")
		write.table(output_ENA.c[[6]], file = NAME.TC, row.names = FALSE, quote = FALSE)
		}
	#
	rownames(TDTD) <- step.1[[4]]
	colnames(TDTD) <- step.1[[4]]
	output_ENA.c[[7]] <- TDTD
	if(length(which(output_form == "files")) == 1 & length(TDTD) > 1){
		NAME.TD <- paste("Dependency_",filename_out,sep="")
		write.table(output_ENA.c[[7]], file = NAME.TD, row.names = FALSE, quote = FALSE)
		}
	}
	#
	if(length(which(routes == "io"))==1){
	output_ENA.c[[8]] <- step.7
	#
	output_ENA.c[[9]] <- step.8
	}
	#
	#
	##
	##
	#
	#
	if(length(which(routes == "tp"))==1){
	output_ENA.c[[10]] <- step.9
	{
	if(length(output_ENA.c[[10]])==2){
		names(output_ENA.c[[10]]) <- c("FCIL","LinearExchL")
		rownames(output_ENA.c[[10]][[2]]) <- step.1[[4]][1:LN]
		colnames(output_ENA.c[[10]][[2]]) <- step.1[[4]][1:LN]
		if(length(which(output_form == "files")) == 1){
			NAME.LC <- paste("AcyclicFlowsL_",filename_out,sep="")
			write.table(output_ENA.c[[10]][[2]], file = NAME.LC, row.names = FALSE, quote = FALSE)
			}
		}
	#
	else	{
		output_ENA.c[[10]][[10]] <- FCI.c
		names(output_ENA.c[[10]]) <- c("TotalCyclesL","CyclesLengthL","FullCycleAnalysisL",
		"NexusL","WeakArcIntensityL","CycleDistrL","NormDistrL","ResidualFlowsL",
		"AggregatedCyclesL","FCIL")
		rownames(output_ENA.c[[10]][[8]]) <- step.1[[4]][1:LN]
		colnames(output_ENA.c[[10]][[8]]) <- step.1[[4]][1:LN]
		rownames(output_ENA.c[[10]][[9]]) <- step.1[[4]][1:LN]
		colnames(output_ENA.c[[10]][[9]]) <- step.1[[4]][1:LN]
		#
		if(length(which(output_form == "files")) == 1){
			NAME.CLL <- paste("CyclesLengthL_",filename_out,sep="")
			CYC.CLASS <- cbind(c(1:(LN-1)),output_ENA.c[[10]][[2]])
			colnames(CYC.CLASS) <- c("CycleLength","CyclesNumber")
			write.table(CYC.CLASS, file = NAME.CLL, row.names = FALSE, quote = FALSE)
			#
			NAME.RFL <- paste("ResidualFlowsL_",filename_out,sep="")
			write.table(output_ENA.c[[10]][[8]], file = NAME.RFL, row.names = FALSE, quote = FALSE)
			#
			NAME.ACL <- paste("AggregatedCyclesL_",filename_out,sep="")
			write.table(output_ENA.c[[10]][[9]], file = NAME.ACL, row.names = FALSE, quote = FALSE)
			}
		}
	}
	#
	output_ENA.c[[11]] <- step.11
	{
	if(sum(output_ENA.c[[11]][[12]])!=0){
		output_ENA.c[[11]][[6]] <- as.vector(output_ENA.c[[11]][[6]])
		output_ENA.c[[11]][[13]] <- as.vector(output_ENA.c[[11]][[13]])
		names(output_ENA.c[[11]][[13]]) <- step.1[[4]]
		names(output_ENA.c[[11]]) <- c("LindemanTransfM","EffTLs","CanonicalExports","CanonicalResp",
		"GrazingChain","ReturnsDetPool","Detritivory","InputDetPool","CircDetPool","LindemanSpine",
		"TrophicEff","MigratoryInputs","CanonicalInputs")
		#
		if(length(which(output_form == "files")) == 1){
			NAME.LTM <- paste("LindemanTransfM_",filename_out,sep="")
			write.table(output_ENA.c[[11]][[1]], file = NAME.LTM, row.names = FALSE, quote = FALSE)
			tp.s1 <- cbind(step.1[[4]],output_ENA.c[[11]][[2]],output_ENA.c[[11]][[13]],output_ENA.c[[11]][[3]],output_ENA.c[[11]][[4]])
			colnames(tp.s1) <- c("compartments","TPs","CanonicalINP","CanonicalEXP","CanonicalRESP")
			#
			NAME.ETL <- paste("EffTLs_",filename_out,sep="")
			write.table(tp.s1, file = NAME.ETL, row.names = FALSE, quote = FALSE)
			#
			lind.last <- cbind(output_ENA.c[[11]][[5]],output_ENA.c[[11]][[6]],output_ENA.c[[11]][[10]],c(output_ENA.c[[11]][[11]],0))
			NAME.LLV <- paste("LindemanSpine_",filename_out,sep="")
			colnames(lind.last) <- c("GrazingChain","RetDetPool","LindemanSpine","TrophicEff")
			write.table(lind.last, file = NAME.LLV, row.names = FALSE, quote = FALSE)
			}
		}
	#
	else	{
		output_ENA.c[[11]][[6]] <- as.vector(output_ENA.c[[11]][[6]])
		rownames(output_ENA.c[[11]][[1]]) <- step.1[[4]][1:LN]
		colnames(output_ENA.c[[11]][[1]]) <- step.1[[4]][1:LN]
		names(output_ENA.c[[11]][[2]]) <- step.1[[4]]
		names(output_ENA.c[[11]]) <- c("LindemanTransfM","EffTLs","CanonicalExports","CanonicalResp",
		"GrazingChain","ReturnsDetPool","Detritivory","InputDetPool","CircDetPool","LindemanSpine",
		"TrophicEff","MigratoryInputs")
		#
		if(length(which(output_form == "files")) == 1){
			NAME.LTM <- paste("LindemanTransfM_",filename_out,sep="")
			write.table(output_ENA.c[[11]][[1]], file = NAME.LTM, row.names = FALSE, quote = FALSE)
			tp.s1 <- cbind(step.1[[4]],output_ENA.c[[11]][[2]],ZZZ,output_ENA.c[[11]][[3]],output_ENA.c[[11]][[4]])
			colnames(tp.s1) <- c("compartments","TPs","INP","CanonicalEXP","CanonicalRESP")
			#
			NAME.ETL <- paste("EffTLs_",filename_out,sep="")
			write.table(tp.s1, file = NAME.ETL, row.names = FALSE, quote = FALSE)
			#
			lind.last <- cbind(output_ENA.c[[11]][[5]],output_ENA.c[[11]][[6]],output_ENA.c[[11]][[10]],c(output_ENA.c[[11]][[11]],0))
			NAME.LLV <- paste("LindemanSpine_",filename_out,sep="")
			colnames(lind.last) <- c("GrazingChain","RetDetPool","LindemanSpine","TrophicEff")
			write.table(lind.last, file = NAME.LLV, row.names = FALSE, quote = FALSE)
			}
		}
	}
	#
	}
	#
	#
	##
	##
	#
	#
	if(length(which(routes == "ei"))==1){
	output_ENA.c[[12]] <- matrix(rep(NA,32),nrow = 16)
	output_ENA.c[[12]][1,1] <- tst.i
	output_ENA.c[[12]][2,1]	<- dc.i
	output_ENA.c[[12]][3,] <- asc.i
	output_ENA.c[[12]][4,] <- oi.i
	output_ENA.c[[12]][5,] <- oe.i
	output_ENA.c[[12]][6,] <- or.i
	output_ENA.c[[12]][7,] <- r.i
	output_ENA.c[[12]][8,] <- o.i
	output_ENA.c[[12]][9,1] <- idc.i
	output_ENA.c[[12]][10,] <- ia.i[1,]
	output_ENA.c[[12]][11,] <- ia.i[2,]
	output_ENA.c[[12]][12,1] <- oc.i
	output_ENA.c[[12]][13,1] <- ic.i
	output_ENA.c[[12]][14,1] <- fwc.i
	#
	if(length(which(routes == "tp"))==1){
		output_ENA.c[[12]][15,1] <- FCI.c
		}
	#
	if(length(which(routes == "ca"))==1){
		output_ENA.c[[12]][16,1] <- FCI.wc
		}
	#
	colnames(output_ENA.c[[12]]) <- c("Value","Ratio")
	n.wsi <- c("TST","DC","ASC","OI","OE","DO","Redundancy","Overhead","IntDC",
	"IntAsc","Redundancy","OverallConn","IntercompConn","FWConn","LivingFCI","WholeFCI")
	rownames(output_ENA.c[[12]]) <- n.wsi
	#
	if(length(which(output_form == "files")) == 1){
		gl.att <- cbind(n.wsi,output_ENA.c[[12]])
		NAME.WSI <- paste("GlobalAttributes_",filename_out,sep="")
		colnames(gl.att) <- c("Index","Value","Ratio")
		write.table(gl.att, file = NAME.WSI, row.names = FALSE, quote = FALSE)
		}
	#
	}
	#
	#
	##
	##
	#
	#
	if(length(which(routes == "ca"))==1){
	output_ENA.c[[13]] <- step.12
	{
	if(length(output_ENA.c[[13]])==2){
		names(output_ENA.c[[13]]) <- c("FCI","LinearExch")
		rownames(output_ENA.c[[13]][[2]]) <- step.1[[4]]
		colnames(output_ENA.c[[13]][[2]]) <- step.1[[4]]
		#
		if(length(which(output_form == "files")) == 1){
			NAME.CY <- paste("AcyclicFlows_",filename_out,sep="")
			write.table(output_ENA.c[[13]][[2]], file = NAME.CY, row.names = FALSE, quote = FALSE)	
			}
		}
	#
	else	{
		output_ENA.c[[13]][[10]] <- FCI.wc
		names(output_ENA.c[[13]]) <- c("TotalCycles","CyclesLength","FullCycleAnalysis",
		"Nexus","WeakArcIntensity","CycleDistr","NormDistr","ResidualFlows",
		"AggregatedCycles","FCI")
		rownames(output_ENA.c[[13]][[8]]) <- step.1[[4]]
		colnames(output_ENA.c[[13]][[8]]) <- step.1[[4]]
		rownames(output_ENA.c[[13]][[9]]) <- step.1[[4]]
		colnames(output_ENA.c[[13]][[9]]) <- step.1[[4]]
		#
		if(length(which(output_form == "files")) == 1){
			NAME.CYL <- paste("CyclesLength_",filename_out,sep="")
			CYCW.CLASS <- cbind(c(1:(TN-1)),output_ENA.c[[13]][[2]])
			colnames(CYCW.CLASS) <- c("CycleLength","CyclesNumber")
			#
			write.table(CYCW.CLASS, file = NAME.CYL, row.names = FALSE, quote = FALSE)
			#
			NAME.RFY <- paste("ResidualFlows_",filename_out,sep="")
			write.table(output_ENA.c[[13]][[8]], file = NAME.RFY, row.names = FALSE, quote = FALSE)
			#
			NAME.ACY <- paste("AggregatedCycles_",filename_out,sep="")
			write.table(output_ENA.c[[13]][[9]], file = NAME.ACY, row.names = FALSE, quote = FALSE)		
			}
		}
	}
	#
	}
	#
return(output_ENA.c)
}





#####################################################################################################################################################
#####################################################################################################################################################





ENA.e <- function(filename_in, MIG_in, DET_in, filename_out, routes, output_form){
	#
	#					## output_form <- c("text","files")
	#
	if(length(which(output_form == 0)) == 1)print("No stored results", quote = FALSE)
	#
	if(length(which(routes == "all")) == 1){
			routes <- c("cd","io","tp","ei","ca","nf")
			}
	#
	#					# "cd" = Partial host and partial feeding matrices. Dependency and contribution coefficients
	#					# "io" = Input/Output analysis
	#					# "tp" = Trophic aggregation (after applyind de-cyclization algorithm to living exchanges)
	#					# "ei" = Whole system indices, including also: FCI (Ulanowicz, 1986), CCI, AMI, AMI of DC,
	#					#	 biomass inclusive ascendency, effective connectance (m) and Bersier's link density (LD)
	#					# "ca" = Cycling analysis for the whole system
	#					# "nf" = Indirect effects ---> [Q] [M]
	#
	step.1 <- df.ena(filename_in)	
	#
	TN <- step.1[[1]]			## total number of compartments (living + non-living)
	LN <- step.1[[2]]			## living compartments
	NLN <- step.1[[3]]			## non-living compartments
	BIOM <- step.1[[5]]
	#
	ZZZ <- step.1[[6]]
	EEE <- step.1[[7]]
	RRR <- step.1[[8]]
	TTT <- step.1[[9]]
	#
	ZZZ.init <- ZZZ
	EEE.init <- EEE
	RRR.init <- RRR
	TTT.init <- TTT
	#
	ZZZ.L <- ZZZ[1:LN]
	TTT.L <- TTT[c(1:LN),c(1:LN)]
	EEE.L <- EEE[1:LN]
	RRR.L <- RRR[1:LN]
	#
	step.1[[4]] <- as.vector(step.1[[4]])
	rownames(step.1[[9]]) <- step.1[[4]]
	rownames(step.1[[10]]) <- step.1[[4]][1:LN]
	names(step.1[[5]]) <- step.1[[4]]
	names(step.1[[6]]) <- step.1[[4]]
	names(step.1[[7]]) <- step.1[[4]]
	names(step.1[[8]]) <- step.1[[4]]
	#
	#
	if(length(which(output_form == "text")) == 1){
	cat("NUMBER OF COMPARTMENTS IS",TN, file = filename_out, append = TRUE)
	cat("\nNUMBER OF LIVING COMPARTMENTS IS",LN, file = filename_out, append = TRUE)
	cat("\n\nLIVING COMPARTMENTS\n\n", file = filename_out, append = TRUE)
	#
	for(i in 1:LN){
		if(i < 10)cat("   ",i," ", file = filename_out, append = TRUE, sep = "")
		#
		if(i >= 10 & i < 100)cat("  ",i," ", file = filename_out, append = TRUE, sep = "")
		#
		if(i >= 100)cat(" ",i," ", file = filename_out, append = TRUE, sep = "")
		#
		write(step.1[[4]][i], file = filename_out, append = TRUE)
		}
	if(NLN!=0){
		cat("\nNON-LIVING COMPARTMENTS\n\n", file = filename_out, append = TRUE)
		for(i in 1:NLN){
			if(i < 10)cat("   ",LN+i," ", file = filename_out, append = TRUE, sep = "")
			#
			if(i >= 10 & i < 100)cat("  ",LN+i," ", file = filename_out, append = TRUE, sep = "")
			#
			if(i >= 100)cat(" ",LN+i," ", file = filename_out, append = TRUE, sep = "")
			#
			write(step.1[[4]][LN+i], file = filename_out, append = TRUE)
			}
	}
	cat("\n\nINPUT VECTOR\n\n","   ",
			formatC(c(1:TN),digits=1,width=10), "\n\n",
			"      ",formatC(ZZZ,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
	cat("\n\n\nEXPORT VECTOR\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(EEE,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
	cat("\n\n\nRESPIRATION VECTOR\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(RRR,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
	cat("\n\n\nEXCHANGE MATRIX\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n", file = filename_out, append = TRUE, sep="")
			for(i in 1:TN){
				for(j in 1:TN){
					if(j!=TN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(TTT[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
							#
							if(i >=10 & i < 100)cat("  ",i,"  ",formatC(TTT[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
							#
							if(i > 100)cat(" ",i,"  ",formatC(TTT[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
							}
						#
						else cat(formatC(TTT[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(TTT[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
		}
	#
	step.2 <- bigmatrNOnames.ena(ZZZ,TTT,EEE,RRR)
	INP <- apply(step.2,2,sum)[2:(TN+1)]
	OUT <- apply(step.2,1,sum)[2:(TN+1)]
	DIFF_I <- INP - OUT
	DIFF_O <- OUT - INP
	DIFF_IP <- 100*(INP - OUT)/INP
	DIFF_OP <- 100*(OUT - INP)/OUT
	BAL <- cbind(INP,OUT,DIFF_I,DIFF_O,DIFF_IP,DIFF_OP)
	#				
	if(length(which(abs(BAL[,5]) > 0.0001))!=0){
			BMB <- AVG2.balance(ZZZ,TTT,EEE,RRR)
			if(length(BMB) > 1){
				ZZZ <- as.vector(BMB[nrow(BMB),c(1:TN)])
				EEE <- as.vector(BMB[c(1:TN),(TN+1)])
				RRR <- as.vector(BMB[c(1:TN),(TN+2)])
				#
				TTT <- BMB[c(1:TN),c(1:TN)]
				#
				if(NLN!=0){
					ZZZ.L <- ZZZ[c(1:LN)]
					EEE.L <- EEE[c(1:LN)]
					RRR.L <- RRR[c(1:LN)]
					TTT.L <- TTT[c(1:LN),c(1:LN)]
					}
				}
			}
	#
	if(length(which(output_form == "text")) == 1){
	{
	if(length(which(abs(BAL[,5]) > 0.0001))==0){			## Applying balancing algorithm in case of difference (IN - OUT) > 0.001%
		cat("\n\nCOMPARTMENTAL THROUGHPUTS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10), "\n\n",
			"      ",formatC(INP,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep ="")
		}
	else	{
		cat("\n\n\nWARNING! - NETWORK IS NON STEADY-STATE!\n\nCOMPARISON OF THROUGHPUTS\n\n", file = filename_out, append = TRUE)
		cat("        INPUTS     OUTPUTS    IMBALANCE\n\n", file = filename_out, append = TRUE)
			for(i in 1:TN){
				for(j in 1:3){
					#
					if(j==1){
						if(i < 10)cat("   ",i,"   ",formatC(BAL[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
						#
						if(i >= 10 & i < 100)cat("  ",i,"   ",formatC(BAL[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
						#
						if(i >= 100)cat(" ",i,"   ",formatC(BAL[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
						}
					#
					if(j!=1 & j!=3)cat(" ",formatC(BAL[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
					#
					if(j==3)cat(" ",formatC(BAL[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE, sep = "")
				}
			}
		#
		if(length(which(abs(BAL[,5]) > 0.0001))!=0 & length(BMB) > 1){
		cat("\n\n\n*** STEADY-STATE FLOWS AFTER BALANCING PROCEDURE ***\n", file = filename_out, append = TRUE)
		#
		cat("\n\nINPUT VECTOR\n\n","   ",
			formatC(c(1:TN),digits=1,width=10), "\n\n",
			"      ",formatC(ZZZ,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
		cat("\n\n\nEXPORT VECTOR\n\n","   ",
				formatC(c(1:TN),digits=1,width=10),"\n\n",
				"      ",formatC(EEE,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
		cat("\n\n\nRESPIRATION VECTOR\n\n","   ",
				formatC(c(1:TN),digits=1,width=10),"\n\n",
				"      ",formatC(RRR,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
		cat("\n\n\nEXCHANGE MATRIX\n\n","   ",
				formatC(c(1:TN),digits=1,width=10),"\n\n", file = filename_out, append = TRUE, sep="")
				for(i in 1:TN){
					for(j in 1:TN){
						if(j!=TN){
							if(j==1){
								if(i < 10)cat("   ",i,"  ",formatC(TTT[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
								#
								if(i >=10 & i < 100)cat("  ",i,"  ",formatC(TTT[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
								#
								if(i > 100)cat(" ",i,"  ",formatC(TTT[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
								}
							#
							else cat(formatC(TTT[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
							}
						#
						else cat(formatC(TTT[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
					}
				}
			}
		}
	}
	#
	}
	#
	#
	#
	##
	##
	##
	#
	#
	#
	if(length(which(routes == "cd"))==1){
	#
	if(length(which(output_form == "text")) == 1)cat("\n\n\n\n\n*** STRUCTURE ANALYSES ***\n", file = filename_out, append = TRUE)
	#
	GGG <- partial.feeding(ZZZ,TTT)
	FFF <- partial.host(TTT,EEE,RRR)
	TCTC <- total.contribution(ZZZ,TTT,EEE,RRR)
	TDTD <- total.dependency(ZZZ,TTT,EEE,RRR)
	#
	if(length(which(output_form == "text")) == 1){
	cat("\n\nPARTIAL FEEDING MATRIX\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n", file = filename_out, append = TRUE, sep = "")
			for(i in 1:TN){
				for(j in 1:TN){
					if(j!=TN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(GGG[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 10 & i < 100)cat("  ",i,"  ",formatC(GGG[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 100)cat(" ",i,"  ",formatC(GGG[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							}
						#
						else cat(formatC(GGG[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(GGG[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
	#
	cat("\n\n\nPARTIAL HOST MATRIX\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n", file = filename_out, append = TRUE, sep = "")
			for(i in 1:TN){
				for(j in 1:TN){
					if(j!=TN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(FFF[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 10 & i < 100)cat("  ",i,"  ",formatC(FFF[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 100)cat(" ",i,"  ",formatC(FFF[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							}
						#
						else cat(formatC(FFF[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(FFF[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
	#
	if(length(TCTC) > 1){
	cat("\n\n\nTOTAL CONTRIBUTION COEFFICIENTS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n", file = filename_out, append = TRUE, sep = "")
			for(i in 1:TN){
				for(j in 1:TN){
					if(j!=TN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(TCTC[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 10 & i < 100)cat("  ",i,"  ",formatC(TCTC[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 100)cat(" ",i,"  ",formatC(TCTC[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							}
						#
						else cat(formatC(TCTC[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(TCTC[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
		}
	#
	if(length(TDTD) > 1){
	cat("\n\n\nTOTAL DEPENDENCY COEFFICIENTS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n", file = filename_out, append = TRUE, sep = "")
			for(i in 1:TN){
				for(j in 1:TN){
					if(j!=TN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(TDTD[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 10 & i < 100)cat("  ",i,"  ",formatC(TDTD[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 100)cat(" ",i,"  ",formatC(TDTD[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							}
						#
						else cat(formatC(TDTD[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(TDTD[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
		}
	if(length(which(output_form == "text")) == 1)cat("\n\n", file = filename_out, append = TRUE, sep = "")
	}
	#
	}
	#
	if(length(which(routes == "io"))==1){
	options(warn = -1)
	IAN <- input.analysis(ZZZ,TTT,EEE,RRR)
	OAN <- output.analysis(ZZZ,TTT,EEE)
	step.7 <- as.list(rep(NA,length(IAN)))
	step.8 <- as.list(rep(NA,length(OAN)))
	if(is.na(IAN) != TRUE){
	for(k in 1:length(IAN)){
		
		IAN.1 <- which(IAN[[k]][[1]]!=0)		
		IAN.2 <- rep(0,TN)
		IAN.2[IAN.1] <- 1
		IAN.3 <- IAN[[k]][[3]]
		IAN.4 <- IAN[[k]][[4]]
		IAN.5 <- IAN[[k]][[2]]
		#
		if(length(which(output_form == "text")) == 1){
		cat("\n\n\nANALYSIS OF INPUT ",IAN.1, file = filename_out, append = TRUE, sep = "")
		cat("\n\n\nINPUT VECTOR\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(IAN.2,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
		cat("\n\n\nEXPORT VECTOR\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(IAN.3,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
		cat("\n\n\nRESPIRATIONS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(IAN.4,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
		cat("\n\n\nEXCHANGE MATRIX\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n", file = filename_out, append = TRUE, sep="")
			for(i in 1:TN){
				for(j in 1:TN){
					if(j!=TN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(IAN.5[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 10 & i < 100)cat("  ",i,"  ",formatC(IAN.5[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 100)cat(" ",i,"  ",formatC(IAN.5[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							}
						#
						else cat(formatC(IAN.5[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(IAN.5[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
		}
		step.7[[k]] <- bigmatr.ena(step.1[[4]],IAN.2,IAN.5,IAN.3,IAN.4)
		}
	}
	#
	if(length(which(output_form == "text")) == 1)cat("\n\n", file = filename_out, append = TRUE, sep = "")
	#
	if(is.na(OAN) != TRUE){
	for(k in 1:length(OAN)){
		options(warn = -1)
		OAN.1 <- which(OAN[[k]][[3]]!=0)		
		OAN.2 <- OAN[[k]][[1]]
		OAN.3 <- rep(0,TN)
		OAN.3[OAN.1] <- 1
		OAN.4 <- OAN[[k]][[4]]
		OAN.5 <- OAN[[k]][[2]]
		#
		if(length(which(output_form == "text")) == 1){
		cat("\n\n\nANALYSIS OF OUTPUT ",OAN.1, file = filename_out, append = TRUE, sep = "")
		cat("\n\n\nINPUT VECTOR\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(OAN.2,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
		cat("\n\n\nEXPORT VECTOR\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(OAN.3,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
		cat("\n\n\nRESPIRATIONS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(OAN.4,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
		cat("\n\n\nEXCHANGE MATRIX\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n", file = filename_out, append = TRUE, sep="")
			for(i in 1:TN){
				for(j in 1:TN){
					if(j!=TN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(OAN.5[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 10 & i < 100)cat("  ",i,"  ",formatC(OAN.5[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 100)cat(" ",i,"  ",formatC(OAN.5[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							}
						#
						else cat(formatC(OAN.5[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(OAN.5[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
		}
		step.8[[k]] <- bigmatr.ena(step.1[[4]],OAN.2,OAN.5,OAN.3,OAN.4)
		}
	}
	#
	if(length(which(output_form == "files")) == 1){
	#
	if(is.na(IAN) != TRUE){
	NAME.IA <- paste("InputAnalysis_",filename_out,sep="")
	for(k in 1:length(IAN)){		
		write.table(step.7[[k]], file = NAME.IA, row.names = FALSE, quote = FALSE, append = TRUE)
		}
		#
		}
	#
	if(is.na(OAN) != TRUE){
	NAME.OA <- paste("OutputAnalysis_",filename_out,sep="")
	for(k in 1:length(OAN)){
		write.table(step.8[[k]], file = NAME.OA, row.names = FALSE, quote = FALSE, append = TRUE)
		}
		#
		}
	}
	#
	}
	#
	#
	options(warn = 0)
	#
	if(length(which(routes == "nf"))==1){
		ind.eff.Q <- Q.TI(ZZZ,TTT,EEE,NLN)
		ind.eff.M <- M.TI(ZZZ,TTT,EEE,NLN)
		#
		step.9 <- ind.eff.Q
		step.10 <- ind.eff.M
		#
		BEN.PRED <- which(ind.eff.Q < 0 & ind.eff.M > 0)	# "Beneficial predators" (compartments with
		#							# negative direct impacts, but positive total trophic effect)
		#
		MAL.PREY <- which(ind.eff.Q > 0 & ind.eff.M < 0)	# "Malefic prey" (compartments with
		#							# positive direct impacts, but negative total trophic effect)
		#
		if(length(which(output_form == "text")) == 1){
			if(length(step.9) > 1 | length(step.10) > 1 ){
			cat("\n\n\n\n\n*** INDIRECT EFFECTS ***\n", file = filename_out, append = TRUE)
			}
			#
			if(length(step.9) > 1){
			cat("\n\nMATRIX OF DIRECT TROPHIC INTERACTIONS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n", file = filename_out, append = TRUE, sep="")
				for(i in 1:TN){
					for(j in 1:TN){
						if(j!=TN){
							if(j==1){
								if(i < 10)cat("   ",i,"  ",formatC(step.9[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
								#
								if(i >=10 & i < 100)cat("  ",i,"  ",formatC(step.9[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
								#
								if(i > 100)cat(" ",i,"  ",formatC(step.9[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
								}
							#
							else cat(formatC(step.9[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
							}
						#
						else cat(formatC(step.9[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
					}
				}
			}
			#
			#
			if(length(step.10) > 1){
			cat("\n\nMATRIX OF TOTAL TROPHIC IMPACTS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n", file = filename_out, append = TRUE, sep="")
				for(i in 1:TN){
					for(j in 1:TN){
						if(j!=TN){
							if(j==1){
								if(i < 10)cat("   ",i,"  ",formatC(step.10[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
								#
								if(i >=10 & i < 100)cat("  ",i,"  ",formatC(step.10[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
								#
								if(i > 100)cat(" ",i,"  ",formatC(step.10[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep="")
								}
							#
							else cat(formatC(step.10[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
							}
						#
						else cat(formatC(step.10[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
					}
				}
			}			
			#
			#
			if(length(step.9) > 1 & length(step.10) > 1 ){
			cat("\n\n\nSEARCH FOR BENEFICIAL PREDATORS\n\n", file = filename_out, append = TRUE)
			#
			{
			if(length(BEN.PRED)==0)cat("(NO BENEFICIAL PREDATORS FOUND)\n\n\n", file = filename_out, append = TRUE)
			else	{
				kk <- 1
				for(p in 1:length(BEN.PRED)){
					prey <- BEN.PRED[p] %% TN
					if(prey == 0) prey <- TN
					predat <- trunc(BEN.PRED[p]/TN)+1
					#
					if(kk < 10){
						cat("   ",kk,". PREDATOR ", predat, " IS BENEFICIAL TO PREY ", prey,"\n\n", 
						file = filename_out, append = TRUE, sep = "")
						kk <- kk + 1
						}
					#
					if(kk >= 10 & kk < 100){
						cat("  ",kk,". PREDATOR ", predat, " IS BENEFICIAL TO PREY ", prey,"\n\n",
						file = filename_out, append = TRUE, sep = "")
						kk <- kk + 1
						}
					if(kk >= 100){
						cat(" ",kk,". PREDATOR ", predat, " IS BENEFICIAL TO PREY ", prey,"\n\n",
						file = filename_out, append = TRUE, sep = "")
						kk <- kk + 1
						}
					}
				}
				
			}
			#
			#
			#
			cat("\nSEARCH FOR MALEFIC PREY\n\n", file = filename_out, append = TRUE)
			#
			{
			if(length(MAL.PREY)==0)cat("(NO MALEFIC PREY FOUND)\n\n\n", file = filename_out, append = TRUE)
			else	{
				kk <- 1
				for(p in 1:length(MAL.PREY)){
					prey <- MAL.PREY[p] %% TN
					if(prey == 0) prey <- TN
					predat <- trunc(MAL.PREY[p]/TN)+1
					#
					if(kk < 10){
						cat("   ",kk,". PREY ", prey, " IS MALEFIC TO PREDATOR ", predat,"\n\n", 
						file = filename_out, append = TRUE, sep = "")
						kk <- kk + 1
						}
					#
					if(kk >= 10 & kk < 100){
						cat("  ",kk,". PREY ", prey, " IS MALEFIC TO PREDATOR ", predat,"\n\n",
						file = filename_out, append = TRUE, sep = "")
						kk <- kk + 1
						}
					if(kk >= 100){
						cat(" ",kk,". PREY ", prey, " IS MALEFIC TO PREDATOR ", predat,"\n\n",
						file = filename_out, append = TRUE, sep = "")
						kk <- kk + 1
						}
					}
				}
			}
			#
			}
		#
		}
	}
	#
	#
	##
	##
	##
	#
	#
	#
	if(length(which(routes == "tp"))==1){
	#
	step.11 <- cycling.analysis(ZZZ.L,TTT.L,EEE.L,RRR.L)		## Full cycling analysis (living)
	#
	if(length(which(output_form == "text")) == 1)cat("\n\n\n\nSEARCH FOR FEEDING CYCLES\n", file = filename_out, append = TRUE)
	{
	if(step.11[[1]]==0){
		if(length(which(output_form == "text")) == 1)cat("\nNO CYCLES DETECTED", file = filename_out, append = TRUE)
		RFL.c <- step.11[[2]]
		FCI.c <- 0
		CCI.c <- 0
		}
	#
	else	{
		TOT.c <- step.11[[1]]		## Total number of living cycles
		COM.c <- step.11[[2]]		## Total living cycles classified by their length
		EAC.c <- step.11[[3]]		## Full cycle analysis (living)
		NEX.c <- step.11[[4]]		## Nexus (living)
		WAR.c <- step.11[[5]]		## Weak arc intensity (living)
		CYD.c <- step.11[[6]]		## Living cycle distributions
		NCD.c <- step.11[[7]]		## Normalized living cycle distributions
		RFL.c <- step.11[[8]]		## Residual living flows
		AFL.c <- step.11[[9]]		## Aggregated living cycles
		#
		if(length(WAR.c)==2)WAR.c <- matrix(c(WAR.c),nrow=1)
		#
		FCI.c <- FCI(ZZZ.L,TTT.L,EEE.L,RRR.L)[2]		## Finn cycling index for living compartments
		CCI.c <- CCI(ZZZ.L,TTT.L,EEE.L,RRR.L)
		#
		if(length(which(output_form == "text")) == 1){
		#
		all.cyc <- 1
		#
		for(k in 1:nrow(WAR.c)){
			for(q in 1:nrow(EAC.c[[k]])){
				if(q == 51){
					cat("\n       . . .", file = filename_out, append = TRUE, sep="")
					all.cyc <- all.cyc + (nrow(EAC.c[[k]]) - 50)
					break
					}
					
				#
				if(q <= 50){
					t <- 1
					while(is.na(EAC.c[[k]][q,t])==FALSE){
								if(t==1){
									cat("\n     ", all.cyc, ".   ",EAC.c[[k]][q,t],"-", file = filename_out, append = TRUE, sep="")
									t <- t + 1
									all.cyc <- all.cyc + 1
									}
								else	{
									cat(EAC.c[[k]][q,t],"-", file = filename_out, append = TRUE, sep="")
									t <- t + 1
									}
								}
					}
				}
			cat("\nNEXUS ", k, " CONSISTS OF ", nrow(EAC.c[[k]]), " CYCLES. WEAK ARC (", WAR.c[k,1], ", ", WAR.c[k,2], ") = ", 
			formatC(NEX.c[k],digits=3,width=10,format="e"), "\n", file = filename_out, append = TRUE, sep="")
			}
		cat("\n\n\nA TOTAL OF", TOT.c, "CYCLES REMOVED", file = filename_out, append = TRUE)
		if(is.na(FCI.c)!=TRUE){
			cat("\n\n\n\nFINN CYCLING INDEX IS ", formatC(FCI.c,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
			cat("\n\nCOMPREHENSIVE CYCLING INDEX IS ", formatC(CCI.c,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
			}
		cat("\n\n\n\nCYCLE-FREE FEEDING TRANSFERS\n","   ",
		formatC(c(1:LN),digits=1,width=10),"\n", file = filename_out, append = TRUE, sep = "")
			for(i in 1:LN){
				for(j in 1:LN){
					if(j!=LN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(RFL.c[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 10 & i < 100)cat("  ",i,"  ",formatC(RFL.c[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 100)cat(" ",i,"  ",formatC(RFL.c[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							}
						#
						else cat(formatC(RFL.c[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(RFL.c[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
		#
		}
		#
		}
	}
	#
	#
	#
	##
	##
	##
	#
	#
	#
	{
	if(length(MIG_in) == 1 & length(DET_in) == 1){
		step.12 <- TP.CTA(ZZZ,TTT,EEE,RRR,NLN)
		MI.cta <- step.12[[12]]				## Detecting Migratory Imports
		if(sum(MI.cta!=0)){
			CI.cta <- step.12[[13]]			## Canonical Distribution of inputs
			}					## (only in case of migrations)
		}
	#
	else	{
		if((length(MIG_in) > 1) & (length(DET_in) > 1)){
			step.12 <- TP.extendedCTA(ZZZ,TTT,EEE,RRR,NLN,MIG_in,DET_in)
			MI.cta <- which(MIG_in!=0)
			{
			if(NLN!=0){
				CI.step1 <- step.12[[1]]%*%t(t(ZZZ.L))
				CI.cta <- c(CI.step1,ZZZ[(LN+1):TN])
				}
			else CI.sta <- step.12[[1]]%*%t(t(ZZZ))
			}
			#						
			}
		#
		if((length(MIG_in) > 1) & (length(DET_in) == 1)){
			if(NLN == 1){
				DET_in <- rep(0,LN)
				DET_in[which(TTT[TN,(1:LN)] != 0)] <- 1
				}
			if(NLN > 1){
				DET_in <- rep(0,LN)
				NLTOL <- apply(TTT[((LN+1):TN),(1:LN)],2,sum)
				DET_in[which(NLTOL != 0)] <- 1
				}
			#
			step.12 <- TP.extendedCTA(ZZZ,TTT,EEE,RRR,NLN,MIG_in,DET_in)
			MI.cta <- which(MIG_in!=0)
			{
			if(NLN!=0){
				CI.step1 <- step.12[[1]]%*%t(t(ZZZ.L))
				CI.cta <- c(CI.step1,ZZZ[(LN+1):TN])
				}
			else CI.sta <- step.12[[1]]%*%t(t(ZZZ))
			}
			#
			}
		#
		if((length(MIG_in) == 1) & (length(DET_in) > 1)){
			MIG_in <- rep(0,LN)
			step.12 <- TP.extendedCTA(ZZZ,TTT,EEE,RRR,NLN,MIG_in,DET_in)
			}
		#
		}
	}
	#
	LT.cta <- step.12[[1]]				## Lindeman Transformation Matrix
	TL.cta <- step.12[[2]]				## Effective Trophic Levels of each species
	CE.cta <- step.12[[3]]				## Canonical Exports
	CR.cta <- step.12[[4]]				## Canonical Respirations
	GC.cta <- step.12[[5]]				## The Grazing Chain
	DP.cta <- step.12[[6]]				## Returns to Detrital Pool
	DT.cta <- step.12[[7]]				## Detritivory
	ID.cta <- step.12[[8]]				## Input to Detrital Pool
	CD.cta <- step.12[[9]]				## Circulation within Detrital Pool
	LS.cta <- step.12[[10]]				## Lindeman Spine
	TE.cta <- step.12[[11]]				## Trophic Efficiencies
	TE.cta <- c(TE.cta,0)
	#
	#
	if(length(which(output_form == "text")) == 1){
	cat("\n\n\n\n\n*** LINDEMAN TROPHIC AGGREGATIONS ***\n\n", file = filename_out, append = TRUE)
	#
	if(sum(MI.cta!=0)){
		cat("\nTHE FOLLOWING INPUTS ARE SUSPECTED TO BE MIGRATORY INPUTS.", file = filename_out, append = TRUE)
		cat("\nTHEY WILL BE TREATED AS NON-PRIMARY INFLOWS.\n", file = filename_out, append = TRUE)
 		for(z in 1:length(MI.cta)){
 			cat("\nMIGRATORY INPUT TO ", MI.cta[z], " = ", 
 			formatC(ZZZ[MI.cta[z]],digits=4,width=10,format="e"), file = filename_out, append = TRUE, sep ="")
 			}
 		cat("\n\n", file = filename_out, append = TRUE)
 		}
	#
	cat("\nLINDEMAN TRANSFORMATION MATRIX\n\n","   ",
		formatC(c(1:LN),digits=1,width=10),"\n\n", file = filename_out, append = TRUE, sep = "")
			for(i in 1:LN){
				for(j in 1:LN){
					if(j!=LN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(LT.cta[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 10 & i < 100)cat("  ",i,"  ",formatC(LT.cta[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 100)cat(" ",i,"  ",formatC(LT.cta[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							}
						#
						else cat(formatC(LT.cta[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(LT.cta[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
	#
	cat("\n\nEFFECTIVE TROPHIC LEVELS OF EACH SPECIES\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(TL.cta,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
	#
	if(sum(MI.cta!=0)){
		cat("\n\n\nCANONICAL DISTRIBUTION OF INPUTS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(CI.cta,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
		}
	#
	cat("\n\n\nCANONICAL EXPORTS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(CE.cta,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
	#
	cat("\n\n\nCANONICAL RESPIRATIONS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(CR.cta,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
	#
	cat("\n\n\nTHE GRAZING CHAIN\n\n","   ",
			formatC(c(1:LN),digits=1,width=10),"\n\n",
			"      ",formatC(GC.cta,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
	#
	cat("\n\n\nRETURNS TO DETRITAL POOL\n\n","   ",
			formatC(c(1:LN),digits=1,width=10),"\n\n",
			"      ",formatC(DP.cta,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
	#
	cat("\n\n\nDETRITIVORY =", formatC(DT.cta,digits=4,width=10,format="e"), "\n", file = filename_out, append = TRUE)
	#
	cat("\nINPUT TO DETRITAL POOL =", formatC(ID.cta,digits=4,width=10,format="e"), "\n", file = filename_out, append = TRUE)
	#
	cat("\nCIRCULATION WITHIN DETRITAL POOL =", formatC(CD.cta,digits=4,width=10,format="e"), "\n\n", file = filename_out, append = TRUE)
	#
	cat("\nLINDEMAN SPINE\n\n","   ",
			formatC(c(1:LN),digits=1,width=10),"\n\n",
			"      ",formatC(LS.cta,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
	#
	cat("\n\n\nTROPHIC EFFICIENCIES\n\n","   ",
			formatC(c(1:LN),digits=1,width=10),"\n\n",
			"      ",formatC(TE.cta,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
	}
	#
	}
	#
	#
	#
	##
	##
	##
	#
	#
	#
	if(length(which(routes == "ei"))==1){
	#
	tst.i <- TST(ZZZ,TTT,EEE,RRR)
	ami.i <- AMI(ZZZ,TTT,EEE,RRR)
	dc.i <- DC(ZZZ,TTT,EEE,RRR)
	amidc.i <- AMI.max(ZZZ,TTT,EEE,RRR)
	asc.i <- ASC(ZZZ,TTT,EEE,RRR)
	asc.bio <- ASC.B(ZZZ,TTT,EEE,RRR,BIOM)
	oi.i <- OVERHEAD.imports(ZZZ,TTT,EEE,RRR)
	oe.i <- OVERHEAD.exports(ZZZ,TTT,EEE,RRR)
	or.i <- OVERHEAD.resp(ZZZ,TTT,EEE,RRR)
	r.i <- REDUNDANCY(ZZZ,TTT,EEE,RRR)
	o.i <- OVERHEAD(ZZZ,TTT,EEE,RRR)
	idc.i <- internal.DC(ZZZ,TTT,EEE,RRR)
	ia.i <- internal.ASC(ZZZ,TTT,EEE,RRR)
	oc.i <- O.C(ZZZ,TTT,EEE,RRR)
	ic.i <- I.C(TTT)
	fwc.i <- FW.C(TTT,NLN)
	eff.c <- effconn.m(ZZZ,TTT,EEE,RRR)
	link.den <- LD.B(ZZZ,TTT,EEE,RRR)
	#
	AMI.in.w <- AMI.in(ZZZ,TTT,EEE,RRR,NLN)[[1]]
	AMI.in.l <- AMI.in(ZZZ,TTT,EEE,RRR,NLN)[[2]]
	if(NLN!=0) AMI.in.l <- c(AMI.in.l,rep(0,NLN))
	AMI.out.w <- AMI.out(ZZZ,TTT,EEE,RRR,NLN)[[1]]
	AMI.out.l <- AMI.out(ZZZ,TTT,EEE,RRR,NLN)[[2]]
	if(NLN!=0) AMI.out.l <- c(AMI.out.l,rep(0,NLN))
	step.13 <- cbind(AMI.in.w,AMI.in.l,AMI.out.w,AMI.out.l)
	#
	#
	if(length(which(output_form == "text")) == 1){
	cat("\n\n\n\n\n\n*** INFORMATION INDICES ***\n\n", file = filename_out, append = TRUE)
	#
	cat("\nTOTAL SYSTEM THROUGHPUT =", formatC(tst.i,digits=4,width=10,format="e"), "\n\n", file = filename_out, append = TRUE)
	#
	cat("\nDEVELOPMENT CAPACITY = ", formatC(dc.i,digits=4,width=10,format="e"), "\n", file = filename_out, append = TRUE, sep="")
	cat("\nASCENDENCY = ", formatC(asc.i[1],digits=4,width=10,format="e"), "                        (",
		formatC(asc.i[2],digits=3,width=0,format="f"), ")\n", file = filename_out, append = TRUE, sep="")
	cat("\nOVERHEAD ON IMPORTS = ", formatC(oi.i[1],digits=4,width=10,format="e"), "               (",
		formatC(oi.i[2],digits=3,width=0,format="f"), ")\n", file = filename_out, append = TRUE, sep="")
	cat("\nOVERHEAD ON EXPORTS = ", formatC(oe.i[1],digits=4,width=10,format="e"), "               (",
		formatC(oe.i[2],digits=3,width=0,format="f"), ")\n", file = filename_out, append = TRUE, sep="")
	cat("\nDISSIPATIVE OVERHEAD = ", formatC(or.i[1],digits=4,width=10,format="e"), "              (",
		formatC(or.i[2],digits=3,width=0,format="f"), ")\n", file = filename_out, append = TRUE, sep="")
	cat("\nREDUNDANCY = ", formatC(r.i[1],digits=4,width=10,format="e"), "                        (",
		formatC(r.i[2],digits=3,width=0,format="f"), ")\n", file = filename_out, append = TRUE, sep="")
	#
	cat("\n\nMAX AVERAGE MUTUAL INFORMATION = ", formatC(amidc.i,digits=4,width=10,format="e"), "\n", file = filename_out, append = TRUE, sep="")
	cat("\nAVERAGE MUTUAL INFORMATION = ", formatC(ami.i,digits=4,width=10,format="e"), "\n", file = filename_out, append = TRUE, sep="")
	cat("\nBIOMASS-INCLUSIVE ASCENDENCY = ", formatC(asc.bio,digits=4,width=10,format="e"), "\n", file = filename_out, append = TRUE, sep="")
	#
	cat("\n\nINTERNAL CAPACITY = ", formatC(idc.i,digits=4,width=10,format="e"), "\n", file = filename_out, append = TRUE, sep="")
	cat("\nINTERNAL ASCENDENCY = ", formatC(ia.i[1,1],digits=4,width=10,format="e"), "               (",
		formatC(ia.i[1,2],digits=3,width=0,format="f"), ")\n", file = filename_out, append = TRUE, sep="")
	cat("\nREDUNDANCY = ", formatC(r.i[1],digits=4,width=10,format="e"), "                        (",
		formatC((1-ia.i[1,2]),digits=3,width=0,format="f"), ")\n", file = filename_out, append = TRUE, sep="")
	#
	#
	cat("\n\n\n*** CONNECTANCE INDICES ***\n\n", file = filename_out, append = TRUE)
	#
	cat("\nOVERALL CONNECTANCE = ", formatC(oc.i,digits=3,width=10,format="f"), "\n", file = filename_out, append = TRUE, sep="")
	cat("\nINTERCOMPARTMENTAL CONNECTANCE = ", formatC(ic.i,digits=3,width=10,format="f"), "\n", file = filename_out, append = TRUE, sep="")
	cat("\nFOODWEB CONNECTANCE = ", formatC(fwc.i,digits=3,width=10,format="f"), "\n", file = filename_out, append = TRUE, sep="")
	cat("\nEFFECTIVE CONNECTANCE = ", formatC(eff.c,digits=3,width=10,format="f"), "\n", file = filename_out, append = TRUE, sep="")
	cat("\nLINK DENSITY = ", formatC(link.den,digits=3,width=10,format="f"), "\n", file = filename_out, append = TRUE, sep="")
	#
	cat("\n\n\n*** AVERAGE MUTUAL INFORMATION OF EACH COMPARTMENT ***\n\n\n", file = filename_out, append = TRUE)
		cat("        AMI IN W   AMI IN L   AMI OUT W  AMI OUT L\n\n", file = filename_out, append = TRUE)
			for(i in 1:TN){
				for(j in 1:4){
					#
					if(j==1){
						if(i < 10)cat("   ",i,"   ",formatC(step.13[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
						#
						if(i >= 10 & i < 100)cat("  ",i,"   ",formatC(step.13[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
						#
						if(i >= 100)cat(" ",i,"   ",formatC(step.13[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
						}
					#
					if(j!=1 & j!=4)cat(" ",formatC(step.13[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
					#
					if(j==4)cat(" ",formatC(step.13[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE, sep = "")
				}
			}
	}
	#
	}
	#
	#
	#
	##
	##
	##
	#
	#
	#
	if(length(which(routes == "ca"))==1){
	#
	step.14 <- cycling.analysis(ZZZ,TTT,EEE,RRR)
	#
	if(length(which(output_form == "text")) == 1)cat("\n\n\n\n\n*** FULL CYCLE ANALYSIS ***\n\n", file = filename_out, append = TRUE)
	{
	if(step.14[[1]]==0){
		if(length(which(output_form == "text")) == 1)cat("\nNO CYCLES DETECTED", file = filename_out, append = TRUE)
		RFL.wc <- step.14[[2]]
		FCI.wc <- 0
		CCI.wc <- 0
		}
	#
	else	{
		TOT.wc <- step.14[[1]]		## Total number of cycles
		COM.wc <- step.14[[2]]		## Total cycles classified by their length
		EAC.wc <- step.14[[3]]		## Full cycle analysis
		NEX.wc <- step.14[[4]]		## Nexus
		WAR.wc <- step.14[[5]]		## Weak arc intensity
		CYD.wc <- step.14[[6]]		## Cycle distributions
		NCD.ww <- step.14[[7]]		## Normalized cycle distributions
		RFL.wc <- step.14[[8]]		## Residual flows
		AFL.wc <- step.14[[9]]		## Aggregated cycles
		#
		FCI.wc <- FCI(ZZZ,TTT,EEE,RRR)[2]
		CCI.wc <- CCI(ZZZ,TTT,EEE,RRR)
		#
		if(length(WAR.wc)==2)WAR.wc <- matrix(c(WAR.wc),nrow=1)
		#
		if(length(which(output_form == "text")) == 1){
		#
		all.cyc.w <- 1
		#
		for(k in 1:nrow(WAR.wc)){
			for(q in 1:nrow(EAC.wc[[k]])){
				if(q == 51){
					cat("\n       . . .", file = filename_out, append = TRUE, sep="")
					all.cyc.w <- all.cyc.w + (nrow(EAC.wc[[k]]) - 50)
					break
					}
				#
				if(q <= 50){
					t <- 1
					while(is.na(EAC.wc[[k]][q,t])==FALSE){
								if(t==1){
									cat("\n     ", all.cyc.w, ".   ",EAC.wc[[k]][q,t],"-", file = filename_out, append = TRUE, sep="")
									t <- t + 1
									all.cyc.w <- all.cyc.w + 1
									}
								else	{
									cat(EAC.wc[[k]][q,t],"-", file = filename_out, append = TRUE, sep="")
									t <- t + 1
									}
								}
					}
				}
			cat("\nNEXUS ", k, " CONSISTS OF ", nrow(EAC.wc[[k]]), " CYCLES. WEAK ARC (", WAR.wc[k,1], ", ", WAR.wc[k,2], ") = ", 
			formatC(NEX.wc[k],digits=3,width=10,format="e"), "\n", file = filename_out, append = TRUE, sep="")
			}
		cat("\n\nA TOTAL OF", TOT.wc, "CYCLES REMOVED\n", file = filename_out, append = TRUE)
		#
		cat("\n\n\nCYCLE DISTRIBUTIONS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(CYD.wc,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
		#
		cat("\n\n\nNORMALIZED DISTRIBUTIONS\n\n","   ",
			formatC(c(1:TN),digits=1,width=10),"\n\n",
			"      ",formatC(NCD.ww,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
		#
		if(is.na(FCI.wc)!=TRUE){
			cat("\n\n\n\nFINN CYCLING INDEX IS ", formatC(FCI.wc,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
			cat("\n\nCOMPREHENSIVE CYCLING INDEX IS ", formatC(CCI.wc,digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
			}
		cat("\n\n\n\nRESIDUAL FLOWS\n\n","   ",
		formatC(c(1:TN),digits=1,width=10),"\n", file = filename_out, append = TRUE, sep = "")
			for(i in 1:TN){
				for(j in 1:TN){
					if(j!=TN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(RFL.wc[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 10 & i < 100)cat("  ",i,"  ",formatC(RFL.wc[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 100)cat(" ",i,"  ",formatC(RFL.wc[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							}
						#
						else cat(formatC(RFL.wc[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(RFL.wc[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
		#
		cat("\n\n\nAGGREGATED CYCLES\n\n","   ",
		formatC(c(1:TN),digits=1,width=10),"\n", file = filename_out, append = TRUE, sep = "")
			for(i in 1:TN){
				for(j in 1:TN){
					if(j!=TN){
						if(j==1){
							if(i < 10)cat("   ",i,"  ",formatC(AFL.wc[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 10 & i < 100)cat("  ",i,"  ",formatC(AFL.wc[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							#
							if(i >= 100)cat(" ",i,"  ",formatC(AFL.wc[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE, sep = "")
							}
						#
						else cat(formatC(AFL.wc[i,j],digits=3,width=10,format="e"), file = filename_out, append = TRUE)
						}
					#
					else cat(formatC(AFL.wc[i,j],digits=3,width=10,format="e"),"\n", file = filename_out, append = TRUE)
				}
			}
		#
		}
	#
	}
	#
	}
	#
	}
	output_ENA.e <- as.list(rep(NA,17))
	#
	#
	##
	##
	#
	#
	output_ENA.e[[1]] <- step.1
	names(output_ENA.e[[1]]) <- c("Compartments","Living","NonLiving","Names",
	"Biomasses","Input","Export","Respiration","Exchange","LivingExchange")
	#
	output_ENA.e[[2]] <- bigmatr.ena(step.1[[4]],ZZZ.init,TTT.init,EEE.init,RRR.init)
	if(length(which(output_form == "files")) == 1){
		NAME.BM <- paste("BigMatrix_",filename_out,sep="")
		write.table(output_ENA.e[[2]], file = NAME.BM, row.names = FALSE, quote = FALSE)
		}
	#
	rownames(BAL) <- step.1[[4]]
	colnames(BAL) <- c("Inp","Out","Diff(I)","Diff(O)","Diff(I%)","Diff(O%)")
	output_ENA.e[[3]] <- BAL
	if(length(which(output_form == "files")) == 1){
		NAME.SC <- paste("SteadyStateCheck_",filename_out,sep="")
		write.table(BAL, file = NAME.SC, row.names = FALSE, quote = FALSE)
		}
	#
	if(length(which(abs(BAL[,5]) > 0.0001))!=0){
		output_ENA.e[[4]] <- bigmatr.ena(step.1[[4]],ZZZ,TTT,EEE,RRR)
		if(length(which(output_form == "files")) == 1 & length(output_ENA.e[[4]]) > 1){
			NAME.BMBB <- paste("BigMatrixBal_",filename_out,sep="")
			write.table(output_ENA.e[[4]], file = NAME.BMBB, row.names = FALSE, quote = FALSE)
		}
	}
	#
	if(length(which(routes == "cd"))==1){
	#
	rownames(GGG) <- step.1[[4]]
	colnames(GGG) <- step.1[[4]]
	output_ENA.e[[5]] <- GGG
	if(length(which(output_form == "files")) == 1){
		NAME.PF <- paste("PartialFeeding_",filename_out,sep="")
		write.table(output_ENA.e[[5]], file = NAME.PF, row.names = FALSE, quote = FALSE)
		}
	#
	rownames(FFF) <- step.1[[4]]
	colnames(FFF) <- step.1[[4]]
	output_ENA.e[[6]] <- FFF
	if(length(which(output_form == "files")) == 1){
		NAME.PH <- paste("PartialHost_",filename_out,sep="")
		write.table(output_ENA.e[[6]], file = NAME.PH, row.names = FALSE, quote = FALSE)
		}
	#
	rownames(TCTC) <- step.1[[4]]
	colnames(TCTC) <- step.1[[4]]
	output_ENA.e[[7]] <- TCTC
	if(length(which(output_form == "files")) == 1 & length(TCTC) > 1){
		NAME.TC <- paste("Contribution_",filename_out,sep="")
		write.table(output_ENA.e[[7]], file = NAME.TC, row.names = FALSE, quote = FALSE)
		}
	#
	rownames(TDTD) <- step.1[[4]]
	colnames(TDTD) <- step.1[[4]]
	output_ENA.e[[8]] <- TDTD
	if(length(which(output_form == "files")) == 1&  length(TDTD) > 1){
		NAME.TD <- paste("Dependency_",filename_out,sep="")
		write.table(output_ENA.e[[8]], file = NAME.TD, row.names = FALSE, quote = FALSE)
		}
	}
	#
	if(length(which(routes == "io"))==1){
	output_ENA.e[[9]] <- step.7
	#
	output_ENA.e[[10]] <- step.8
	}
	#
	#
	##
	##
	#
	#
	if(length(which(routes == "nf"))==1){
		rownames(step.9) <- step.1[[4]]
		colnames(step.9) <- step.1[[4]]
		output_ENA.e[[11]] <- step.9
		if(length(which(output_form == "files")) == 1 & length(step.9) > 1){
			NAME.QQQ <- paste("DirectTrophicInteractions_",filename_out,sep="")
			write.table(output_ENA.e[[11]], file = NAME.QQQ, row.names = FALSE, quote = FALSE)
			}
		#
		rownames(step.10) <- step.1[[4]]
		colnames(step.10) <- step.1[[4]]
		output_ENA.e[[12]] <- step.10
		if(length(which(output_form == "files")) == 1 & length(step.10) > 1){
			NAME.MMM <- paste("MatrixTotalImpacts_",filename_out,sep="")
			write.table(output_ENA.e[[12]], file = NAME.MMM, row.names = FALSE, quote = FALSE)
			}
	}
	#
	#
	##
	##
	#
	#
	if(length(which(routes == "tp"))==1){
	output_ENA.e[[13]] <- step.11
	{
	if(length(output_ENA.e[[13]])==2){
		names(output_ENA.e[[13]]) <- c("FCIL","LinearExchL")
		rownames(output_ENA.e[[13]][[2]]) <- step.1[[4]][1:LN]
		colnames(output_ENA.e[[13]][[2]]) <- step.1[[4]][1:LN]
		if(length(which(output_form == "files")) == 1){
			NAME.LC <- paste("AcyclicFlowsL_",filename_out,sep="")
			write.table(output_ENA.e[[13]][[2]], file = NAME.LC, row.names = FALSE, quote = FALSE)
			}
		}
	#
	else	{
		output_ENA.e[[13]][[10]] <- FCI.c
		output_ENA.e[[13]][[11]] <- CCI.c
		names(output_ENA.e[[13]]) <- c("TotalCyclesL","CyclesLengthL","FullCycleAnalysisL",
		"NexusL","WeakArcIntensityL","CycleDistrL","NormDistrL","ResidualFlowsL",
		"AggregatedCyclesL","FCIL","CCIL")
		rownames(output_ENA.e[[13]][[8]]) <- step.1[[4]][1:LN]
		colnames(output_ENA.e[[13]][[8]]) <- step.1[[4]][1:LN]
		rownames(output_ENA.e[[13]][[9]]) <- step.1[[4]][1:LN]
		colnames(output_ENA.e[[13]][[9]]) <- step.1[[4]][1:LN]
		#
		if(length(which(output_form == "files")) == 1){
			NAME.CLL <- paste("CyclesLengthL_",filename_out,sep="")
			CYC.CLASS <- cbind(c(1:(LN-1)),output_ENA.e[[13]][[2]])
			colnames(CYC.CLASS) <- c("CycleLength","CyclesNumber")
			write.table(CYC.CLASS, file = NAME.CLL, row.names = FALSE, quote = FALSE)
			#
			NAME.RFL <- paste("ResidualFlowsL_",filename_out,sep="")
			write.table(output_ENA.e[[13]][[8]], file = NAME.RFL, row.names = FALSE, quote = FALSE)
			#
			NAME.ACL <- paste("AggregatedCyclesL_",filename_out,sep="")
			write.table(output_ENA.e[[13]][[9]], file = NAME.ACL, row.names = FALSE, quote = FALSE)
			}
		}
	}
	#
	output_ENA.e[[14]] <- step.12
	{
	if(sum(output_ENA.e[[14]][[12]])!=0){
		output_ENA.e[[14]][[6]] <- as.vector(output_ENA.e[[14]][[6]])
		output_ENA.e[[14]][[13]] <- output_ENA.e[[14]][[1]]%*%t(t(ZZZ.L))
		if(NLN!=0)output_ENA.e[[14]][[13]] <- c(output_ENA.e[[14]][[13]],ZZZ[(LN+1):TN])
		#
		names(output_ENA.e[[14]][[13]]) <- step.1[[4]]
		names(output_ENA.e[[14]]) <- c("LindemanTransfM","EffTLs","CanonicalExports","CanonicalResp",
		"GrazingChain","ReturnsDetPool","Detritivory","InputDetPool","CircDetPool","LindemanSpine",
		"TrophicEff","MigratoryInputs","CanonicalInputs")
		#
		if(length(which(output_form == "files")) == 1){
			NAME.LTM <- paste("LindemanTransfM_",filename_out,sep="")
			write.table(output_ENA.e[[14]][[1]], file = NAME.LTM, row.names = FALSE, quote = FALSE)
			tp.s1 <- cbind(step.1[[4]],output_ENA.e[[14]][[2]],output_ENA.e[[14]][[13]],output_ENA.e[[14]][[3]],output_ENA.e[[14]][[4]])
			colnames(tp.s1) <- c("compartments","TPs","CanonicalINP","CanonicalEXP","CanonicalRESP")
			#
			NAME.ETL <- paste("EffTLs_",filename_out,sep="")
			write.table(tp.s1, file = NAME.ETL, row.names = FALSE, quote = FALSE)
			#
			lind.last <- cbind(output_ENA.e[[14]][[5]],output_ENA.e[[14]][[6]],output_ENA.e[[14]][[10]],c(output_ENA.e[[14]][[11]],0))
			NAME.LLV <- paste("LindemanSpine_",filename_out,sep="")
			colnames(lind.last) <- c("GrazingChain","RetDetPool","LindemanSpine","TrophicEff")
			write.table(lind.last, file = NAME.LLV, row.names = FALSE, quote = FALSE)
			}
		}
	#
	else	{
		output_ENA.e[[14]][[6]] <- as.vector(output_ENA.e[[14]][[6]])
		rownames(output_ENA.e[[14]][[1]]) <- step.1[[4]][1:LN]
		colnames(output_ENA.e[[14]][[1]]) <- step.1[[4]][1:LN]
		names(output_ENA.e[[14]][[2]]) <- step.1[[4]]
		names(output_ENA.e[[14]]) <- c("LindemanTransfM","EffTLs","CanonicalExports","CanonicalResp",
		"GrazingChain","ReturnsDetPool","Detritivory","InputDetPool","CircDetPool","LindemanSpine",
		"TrophicEff","MigratoryInputs")
		#
		if(length(which(output_form == "files")) == 1){
			NAME.LTM <- paste("LindemanTransfM_",filename_out,sep="")
			write.table(output_ENA.e[[14]][[1]], file = NAME.LTM, row.names = FALSE, quote = FALSE)
			tp.s1 <- cbind(step.1[[4]],output_ENA.e[[14]][[2]],ZZZ,output_ENA.e[[14]][[3]],output_ENA.e[[14]][[4]])
			colnames(tp.s1) <- c("compartments","TPs","INP","CanonicalEXP","CanonicalRESP")
			#
			NAME.ETL <- paste("EffTLs_",filename_out,sep="")
			write.table(tp.s1, file = NAME.ETL, row.names = FALSE, quote = FALSE)
			#
			lind.last <- cbind(output_ENA.e[[14]][[5]],output_ENA.e[[14]][[6]],output_ENA.e[[14]][[10]],c(output_ENA.e[[14]][[11]],0))
			NAME.LLV <- paste("LindemanSpine_",filename_out,sep="")
			colnames(lind.last) <- c("GrazingChain","RetDetPool","LindemanSpine","TrophicEff")
			write.table(lind.last, file = NAME.LLV, row.names = FALSE, quote = FALSE)
			}
		}
	}
	#
	}
	#
	#
	##
	##
	#
	#
	if(length(which(routes == "ei"))==1){
	#
	if(length(which(routes == "tp"))!=1){
		FCI.c <- FCI(ZZZ.L,TTT.L,EEE.L,RRR.L)[2]
		CCI.c <- CCI(ZZZ.L,TTT.L,EEE.L,RRR.L)
	}
	#
	if(length(which(routes == "ca"))!=1){
		FCI.wc <- FCI(ZZZ,TTT,EEE,RRR)[2]
		CCI.wc <- CCI(ZZZ,TTT,EEE,RRR)
	}
	#
	output_ENA.e[[15]] <- matrix(rep(NA,46),nrow = 23)
	output_ENA.e[[15]][1,1] <- tst.i
	output_ENA.e[[15]][2,1]	<- dc.i
	output_ENA.e[[15]][3,] <- asc.i
	output_ENA.e[[15]][4,] <- oi.i
	output_ENA.e[[15]][5,] <- oe.i
	output_ENA.e[[15]][6,] <- or.i
	output_ENA.e[[15]][7,] <- r.i
	output_ENA.e[[15]][8,] <- o.i
	output_ENA.e[[15]][9,1] <- amidc.i
	output_ENA.e[[15]][10,1] <- ami.i
	output_ENA.e[[15]][11,1] <- asc.bio
	output_ENA.e[[15]][12,1] <- idc.i
	output_ENA.e[[15]][13,] <- ia.i[1,]
	output_ENA.e[[15]][14,] <- ia.i[2,]
	output_ENA.e[[15]][15,1] <- oc.i
	output_ENA.e[[15]][16,1] <- ic.i
	output_ENA.e[[15]][17,1] <- fwc.i
	output_ENA.e[[15]][18,1] <- eff.c
	output_ENA.e[[15]][19,1] <- link.den
	output_ENA.e[[15]][20,1] <- FCI.c
	output_ENA.e[[15]][21,1] <- CCI.c
	output_ENA.e[[15]][22,1] <- FCI.wc
	output_ENA.e[[15]][23,1] <- CCI.wc
	#
	colnames(output_ENA.e[[15]]) <- c("Value","Ratio")
	n.wsi <- c("TST","DC","ASC","OI","OE","DO","Redundancy","Overhead",
	"AMImax","AMI","BiomASC",
	"IntDC","IntAsc","Redundancy",
	"OverallConn","IntercompConn","FWConn","m","LD",
	"LivingFCI","LivingCCI","WholeFCI","WholeCCI")
	rownames(output_ENA.e[[15]]) <- n.wsi
	#
	if(length(which(output_form == "files")) == 1){
		gl.att <- cbind(n.wsi,output_ENA.e[[15]])
		NAME.WSI <- paste("GlobalAttributes_",filename_out,sep="")
		colnames(gl.att) <- c("Index","Value","Ratio")
		write.table(gl.att, file = NAME.WSI, row.names = FALSE, quote = FALSE)
		}
	#
	rownames(step.13) <- step.1[[4]]
	colnames(step.13) <- c("AMI_IN_W","AMI_IN_L","AMI_OUT_W","AMI_OUT_L")
	output_ENA.e[[16]] <- step.13
	if(length(which(output_form == "files")) == 1){
		NAME.AMISC <- paste("SingleCompAMI_",filename_out,sep="")
		write.table(output_ENA.e[[16]], file = NAME.AMISC, row.names = FALSE, quote = FALSE)
		}
	#
	}
	#
	#
	##
	##
	#
	#
	if(length(which(routes == "ca"))==1){
	output_ENA.e[[17]] <- step.14
	{
	if(length(output_ENA.e[[17]])==2){
		names(output_ENA.e[[17]]) <- c("FCI","LinearExch")
		rownames(output_ENA.e[[17]][[2]]) <- step.1[[4]]
		colnames(output_ENA.e[[17]][[2]]) <- step.1[[4]]
		#
		if(length(which(output_form == "files")) == 1){
			NAME.CY <- paste("AcyclicFlows_",filename_out,sep="")
			write.table(output_ENA.e[[17]][[2]], file = NAME.CY, row.names = FALSE, quote = FALSE)	
			}
		}
	#
	else	{
		output_ENA.e[[17]][[10]] <- FCI.wc
		output_ENA.e[[17]][[11]] <- CCI.wc
		names(output_ENA.e[[17]]) <- c("TotalCycles","CyclesLength","FullCycleAnalysis",
		"Nexus","WeakArcIntensity","CycleDistr","NormDistr","ResidualFlows",
		"AggregatedCycles","FCI","CCI")
		rownames(output_ENA.e[[17]][[8]]) <- step.1[[4]]
		colnames(output_ENA.e[[17]][[8]]) <- step.1[[4]]
		rownames(output_ENA.e[[17]][[9]]) <- step.1[[4]]
		colnames(output_ENA.e[[17]][[9]]) <- step.1[[4]]
		#
		if(length(which(output_form == "files")) == 1){
			NAME.CYL <- paste("CyclesLength_",filename_out,sep="")
			CYCW.CLASS <- cbind(c(1:(TN-1)),output_ENA.e[[17]][[2]])
			colnames(CYCW.CLASS) <- c("CycleLength","CyclesNumber")
			#
			write.table(CYCW.CLASS, file = NAME.CYL, row.names = FALSE, quote = FALSE)
			#
			NAME.RFY <- paste("ResidualFlows_",filename_out,sep="")
			write.table(output_ENA.e[[17]][[8]], file = NAME.RFY, row.names = FALSE, quote = FALSE)
			#
			NAME.ACY <- paste("AggregatedCycles_",filename_out,sep="")
			write.table(output_ENA.e[[17]][[9]], file = NAME.ACY, row.names = FALSE, quote = FALSE)		
			}
		}
	}
	#
	}
	#
return(output_ENA.e)
}
