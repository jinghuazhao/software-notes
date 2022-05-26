########################################################
# plink2metasoft_subroutine.R                                  
#   Subroutine of plink2metasoft.py
#   This routine replaces mathematical expressions with their values.
#   Free license -- you are free to use it in any ways 
#   Buhm Han (2012)                                    
########################################################

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outfile <- args[2]
data <- read.table(infile)
rsids <- as.matrix(data[,1])
x <- as.matrix(data[,2:ncol(data)])
out <- mat.or.vec(nrow(x), ncol(x))
for (i in 1:nrow(x)) {
  for (j in 1:ncol(x)) {
    if (is.na(x[i,j])) {
      out[i,j] <- NA
    } else {
      out[i,j] <- eval(parse(text=x[i,j]))
    }
  }
}
outframe <- cbind(data.frame(rsids), data.frame(out))
write.table(format(outframe, digits=5),quote=FALSE,
            file=outfile,sep="\t",row.names=FALSE,col.names=FALSE)
    

