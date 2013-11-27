#test.R
#
# Par Guillaume Lahaie
#
# Comment charge-t-on un fichier dans R?

microRNA <- read.csv("dataMirna.csv")
result <-list()
for(i in 2:11) {
    string = paste("Lib",i-1)
    temp <- list()
    temp[['name']] <- string
    temp[['mean']] <- mean(microRNA[[i]])
    temp[['min']] <- min(microRNA[[i]])
    temp[['max']] <- max(microRNA[[i]])
    temp[['median']] <- median(microRNA[[i]])
    temp[['sd']] <- sd(microRNA[[i]])
    result[[i-1]] <- temp
}

write.csv(result)
