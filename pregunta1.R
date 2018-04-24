##script pregunta 1
getwd()

##1)cargar la dataframe (por cierto no carga en R windows)
fullmat<-read.delim(file="../meta/maizteocintle_SNP50k_meta_extended.txt")
class(fullmat)

##2) este archivo es una dataframe

##3)como se ven las primeras 6 lineas

head(fullmat, n=6)

##4)cuantas muestras hay?? (si las muestras son Nsiembra=9107)
summary(fullmat)
fullmat$NSiembra

##5)de cuantos estados hay muestras (en 19 estados)
fullmat$Estado.1

##6)cuantas muestras fueron colectadas antes de 1980 (8 muestras)
fullmat$A.o._de_colecta
colecta<-select(fullmat, A.o._de_colecta)
colecta
sort(colecta$A.o._de_colecta, decreasing = TRUE)

##7)cuantas muestras hay de cada raza (47 muestras)
class(fullmat$Raza)
levels(fullmat$Raza)

##8)en promedio a que altitud fueron colectadas las muestras (1519.242)
altura<-select(fullmat, Altitud) %>% unlist
altura
mean(altura)

##9)a que altitud maxima y minima fueron colectados (max=2769, min=5)
max(altura)
min(altura)

##10)crea una nueva df de datos solo con muestras de la raza olotillo
Olotillo<-subset(fullmat, Raza=="Olotillo")
View(Olotillo)
class(Olotillo)

##11)crea una df solo con los datos Reventador, Jala y ancho
RJA<-subset(fullmat, Raza ==c("Reventador", "Jala", "Ancho"))
levels(RJA$Raza)

##12)Escribe la matriz anterior a un archivo llamado "submat.cvs" en /meta

write.csv(RJA,"../meta/submat.csv")
dir()


