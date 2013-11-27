#Rscript.exe --no-save
data = read.table("dataMirna.csv", header= TRUE, sep = ',')

stats = array(0, dim=c(5,10))
colnames(stats) = colnames(data[,2:11])
rownames(stats) = c('minimum', 'maximum', 'moyenne', 'ecart-type', 'mediane')

for(i in 1:10) {
    stats[1,i] = min(data[,i+1])
    stats[2,i] = max(data[,i+1])
    stats[3,i] = mean(data[,i+1])
    stats[4,i] = sd(data[,i+1])
    stats[5,i] = median(data[,i+1])
}

write.table(stats, 'resultat.txt', row.names=TRUE, col.names=TRUE, quote=FALSE)


#On essaye pour l'abondance
sums = array(0, dim=c(1,10))
colnames(sums) = colnames(data[,2:11])
rownames(sums) = c('somme')

for(i in 1:10) {
    sums[1,i] = sum(data[,i+1])
}

#On a les sommes, comment faire l'histogramme
   
postscript('histo.ps')
par(cex.axis=.9)
couleurs = c('blue','green','red','yellow','black','purple','orange','red','blue','yellow','green')
barplot(sums, main='Distribution des abondances totales des 10 librairies', xlab='librairies', ylab='abondance totale', col=couleurs)
dev.off()

postscript('camembert.ps')
pie(sums)
title('Distribution des abondances totales des 10 librairies')
dev.off()
