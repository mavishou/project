
#get gcc between 2 genes
#from R package rsgcc v1.0.4, edit by Hou Mei
gcc <- function( x, y ) {

  if (!is.vector(x))
    stop("x must be a vector.")
  len = length(x)

  if (!is.vector(y) || len != length(y))
    stop("x is a vector; y should be a vector of the same length")

  #sort x
  tmp <- sort.int(x, na.last = NA, decreasing = FALSE, index.return = TRUE)
  valuex_x <- tmp$x
  valuey_x <- y[tmp$ix] #y sorted by the rank of x

  #sort y
  tmp <- sort.int(y, na.last = NA, decreasing = FALSE, index.return = TRUE)
  valuey_y <- tmp$x
  valuex_y <- x[tmp$ix]  #x sorted by the rank of y

  weight <- t(2*seq(1,len,1) - len - 1)

  gccxy <- sum(weight*valuex_y)/sum(weight*valuex_x)
  gccyx <- sum(weight*valuey_x)/sum(weight*valuey_y)

  #edit by Hou Mei
  #return( data.frame(gccxy, gccyx) )
  gccs <- c(gccxy,gccyx)
  #return the abs max of gccs
  return(max(gccs[abs(gccs)== max(abs(gccs))])) 

}


#get the GCC between a given gene and the whole genes
#oneGene is the name of a given gene
#exprsMatrix is the exression matrix of selected samples, each row represent a gene, each column represent a sample
getGCCor <- function(oneGene, exprsMatrix){
	oneGene_exprs <- exprsMatrix[oneGene,]
	gcc_result <- rep(0,nrow(exprsMatrix))
	names(gcc_result) <- rownames(exprsMatrix)
	for(i in 1:nrow(exprsMatrix)){
		gcc_result[i] <- gcc(oneGene_exprs, exprsMatrix[i,])
	}
	return(gcc_result)
}

