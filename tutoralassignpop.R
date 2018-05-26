##este script es un tutorial del programa assignpop
library(assignPOP)

##leer el archivo de texto de poblaciones simuladas
YourGenepop <- read.Genepop( "simGenepop.txt", pop.names=c("pop_A","pop_B","pop_C"), haploid = FALSE)

##remueve los loci de baja variacion (este paso es opcional, p = 0.95 indicates the removal of loci having the frequency of an allele
## greater than 0.95.)
YourGenepopRd <- reduce.allele( YourGenepop, p = 0.95)

##realiza Monte-Carlo cross-validation. El siguiente comando muestrea el 50, 70 y 90% de individuos de cada poblacion por los valores 10, 25 y 50% de Fst mayores
assign.MC( YourGenepopRd, dir="Result-folder/", train.inds=c(0.5,0.7,0.9), train.loci=c(0.1,0.25,0.5,1), loci.sample="fst", iterations=30, model="svm")

##se requieren dos paquetes para correr la cross-validation klaR y MASS

##Calcular la presicion de la asignacion de los resultados cross-validation
accuMC <- accuracy.MC( dir = "Result-folder/" )

##crear boxplot de asignamiento
library(ggplot2)
accuracy.plot( accuMC, pop=c("all", "pop_A", "pop_B", "pop_C")) +
  ylim(0, 1) + #Set y limit between 0 and 1
  annotate("segment",x=0.4,xend=3.6,y=0.33,yend=0.33,colour="red",size=1) +
  #Add a red horizontal line at y = 0.33 (null assignment rate for 3 populations)
ggtitle("Monte-Carlo cross-validation using genetic loci")+
  #Add a plot title
  theme(plot.title = element_text(size=20, face="bold"))
#Edit plot title text size

##Realiza k-fold cross-validation. El siguiente script divide individuos de cada poblacion en 3, 4 o 5 grupos
assign.kfold( YourGenepopRd, k.fold=c(3,4,5), train.loci=c(0.1,0.25,0.5,1),
              loci.sample="random", dir="Result-folder2/", model="lda" )

##crear grafica de probabilidad de membresia
membership.plot( dir = "Result-folder2/" )

##concatena datos geneticos y morfologicos
YourIntegrateData <- compile.data( YourGenepopRd, "morphData.csv" )

##repite paso III
assign.MC( YourIntegrateData, dir="Resultinteg-folder/", train.inds=c(0.5,0.7,0.9),
           train.loci=c(0.1,0.25,0.5,1), loci.sample="fst", iterations=30, model="svm")

##repite paso 4
accuMCint <- accuracy.MC( dir = "Resultinteg-folder/" )

##crear boxplot de precision de asignamiento con datos mixtos
library(ggplot2)
accuracy.plot( accuMCint, pop=c("all", "pop_A", "pop_B", "pop_C")) +
  ylim(0, 1) + #Set y limit between 0 and 1
  annotate("segment",x=0.4,xend=3.6,y=0.33,yend=0.33,colour="red",size=1) +
  #Add a red horizontal line at y = 0.33 (null assignment rate for 3 populations)
ggtitle("Monte-Carlo cross-validation using genetic loci")+
  #Add a plot title
  theme(plot.title = element_text(size=20, face="bold"))
#Edit plot title text size

##realiza k-fold cross validation
assign.kfold( YourIntegrateData, k.fold=c(3,4,5), train.loci=c(0.1,0.25,0.5,1),
              loci.sample="random", dir="Resultinteg-folder3/", model="lda" )

##crear grafica de probabilidad de membresia
membership.plot( dir = "Resultinteg-folder3/" )

##los resultados de probabilidad de membresia son diferentes en comparación con datos geneticos unicamente.
