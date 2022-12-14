### Skrypt umozliwia przeporowadzenie dwoch roznych filtracji. 
### Pierwsza polega na usunieciu 50% odczytow o najnizszej sredniej intensywnosci sygnalu. Druga to filtracja adaptacyjna oparta o dekompozycje mieszaniny Gaussa


setwd("C:/Users/marty/Desktop/Praca magisterska/analiza_w_R/dane")

## Wgranie danych z kroku 1
BEpozostal<-as.matrix(readRDS("ekspresja_wszystkich.rds"))    #dane(informacje o ekspresji genow) po preprocessingu BEZ usuniecia efektu paczki 
BEusuniety<-readRDS("z_korekta_batch_effect.rds")  #dane po preprocessingu plikow Z usunietym efektem paczki


## FILTRACJA 1- usuniecie 50% odczytow o najnizszej sredniej intensywnosci sygnalu; na danych zlogarytmizowanych 

#Srednia po cechach
BEpozostal_srednia <- rowMeans(BEpozostal)  
BEusuniety_srednia <- rowMeans(BEusuniety)

#Wyznaczenie progu filtracji
BEpozostal_med<- summary(BEpozostal_srednia)[3]
BEusuniety_med<-summary(BEusuniety_srednia)[3]

#Dane po filtracji  (BEpozostal: 19,764-> 9 882; BEusuniety: 19,764->9,882) 
BEpozostal_wynikiFiltracji50<- BEpozostal[which(BEpozostal_srednia>= BEpozostal_med),] 
BEusuniety_wynikiFilttracji50 <- BEusuniety[which(BEusuniety_srednia>= BEusuniety_med),]

#Zapis wynikow
saveRDS(BEpozostal_wynikiFiltracji50,file="BEpozostal_filtracja50.rds")
saveRDS(BEusuniety_wynikiFilttracji50,file="BEusuniety_filtracja50.rds")


## FILTRACJA 2 - filtracja adaptacyjna oparta o dekompozycje mieszaniny Gaussa (progi wyznaczane za pomoca programu GaMRed, ktory dziala w srodowisku Matlab)

#Zapis srednich do txt -> aby moc wgrac te pliki do programu GaMRed
write.table(BEpozostal_srednia, file="srednia_BEpozostal.txt",append = FALSE, sep="\t", dec = ".",
            row.names = TRUE, col.names = TRUE)

write.table(BEusuniety_srednia, file="srednia_BEusuniety.txt",append = FALSE, sep="\t", dec = ".",
            row.names = TRUE, col.names = TRUE)

#Pierwszy punkt przeciecia 
BEpozostal_pp<- 4.17
BEusuniety_pp<- 4.17

#Usuniecie danych na podstawie wyznaczonego punktu przeciecia (pozostaja dane po prawej od punktu)
etap1BEpozostal<- BEpozostal[which(BEpozostal_srednia>= BEpozostal_pp),] #18288 
etap1BEusuniety<- BEusuniety[which(BEusuniety_srednia>= BEusuniety_pp),] #18291

#Wariancja wyliczona na danych, ktore pozostaly po pierwszym etapie filtracji
BEpozostal_war<- rowVars(etap1BEpozostal)
BEusuniety_war<- rowVars(etap1BEusuniety)

#Zapis wariancji do txt -> aby moc wgrac te pliki do programu GaMRed
write.table(BEpozostal_war, file="wariancja_BEpozostal.txt",append = FALSE, sep="\t", dec = ".",
            row.names = TRUE, col.names = TRUE)

write.table(BEusuniety_war, file="wariancja_BEusuniety.txt",append = FALSE, sep="\t", dec = ".",
            row.names = TRUE, col.names = TRUE)

#Drugi punkt przeciecia 
BEpozostal_pp2<-0.0807
BEusuniety_pp2<-0.1799

#Usuniecie danych na podstawie wyznaczonego punktu przeciwcia 
etap2BEpozostal<- etap1BEpozostal[which(BEpozostal_war>= BEpozostal_pp2),] #15683
etap2BEusuniety<- etap1BEusuniety[which(BEusuniety_war>= BEusuniety_pp2),] #9531

#Zapis wynikow
saveRDS(etap2BEpozostal,'BEpozostal_GaMRed.rds')
saveRDS(etap2BEusuniety,'BEusuniety_GaMRed.rds')


