
#http://ofmind.net/doc/r-tips

read.cb <- function(header=TRUE, ...)
    utils::read.table(file="clipboard", header=header, ...)


write.cb <- function(data, sep="\t", header=TRUE, row.names=FALSE,
        col.names=ifelse(header && row.names, NA, header), qmethod="double", ...)
    utils::write.table(data, file="clipboard", sep=sep, row.names=row.names,
        col.names=col.names, qmethod=qmethod, ...)
        
        
