#Metoda opiera sie na przeprowadzeniu testu Kruskala-Wallisa a następnie testu Dunna (testu post-hoc). Tak wyselekcjonowane cechy przekazywane są do SVM. 

# Tym kodem analizowany jest zbior:
# 3: dane z filtracja typu 50 % i bez korekty batch effect   (dane <- readRDS("BEpozostal_filtracja50.rds"))
# 4: dane z filtracja typu 50 % i z korekta batch effect     (dane <- readRDS("BEusuniety_filtracja50.rds"))
# 5: dane z filtracja adaptacyjna i bez korekty batch effect (dane <- readRDS("BEpozostal_GaMRed.rds"))  
# 6: dane z filtracja adaptacyjna i z korekta batch effect   (dane <- readRDS("BEusuniety_GaMRed.rds"))

#Pakiety 
library(stats)
library(e1071)
library(caret)
library(dunn.test)

#Wybranie konkretengo zbioru danych (z zakresu od 3 do 6)
setwd("C:/Users/marty/Desktop/Praca magisterska/analiza_w_R/dane") 
dane <- as.data.frame(readRDS("BEpozostal_GaMRed.rds"))
dane <-as.data.frame(dane)

#Podzial danych wzgledem 18 grup
plik<-readTargets("MILEdataInformation.txt")
grupy <- paste(plik$Characteristics.leukemia.class,sep='')
nazwy <- paste(plik$GSM,sep='')
unikalne<-unique(grupy) #nazwy wszyskich grup
pogrupowane <- lapply(1:length(unikalne), function(i) dane[,which(grupy==unikalne[i])])

#Detekcja liczby corów potrzebna do operacji zrównoleglania 
no_cores <- detectCores() - 1 

#Ustawienie ziarna w celu wygenerowania powtarzalnego losowego pobierania próbek 
set.seed(130)

#Zmienne do zapisywania wynikow
przewidziane_klasy<-list()
model_a_prawidlowa_klasa<-list()
poprawnosc_modelu<-list()
cechy_koncowe_do_svm<-list()

#Podzial na podgrupy
gg<-t(dane)
fff<-as.data.frame(grupy)
nowe<-cbind(gg,fff)
folds1 <- createFolds(nowe$grupy,30)

#Zapis podgrup
saveRDS(folds1,'FOLDS130.rds')

#Petla w ktorej wykonywana jest pelna analiza statystyczna oraz klasyfikacja danych
for (i in 1:30){
  #Podział na zbiór uczący i testowy
  zbior_uczacy<-dane[(-folds1[[i]])] #uczący
  zbior_testowy<-dane[folds1[[i]]] #testowy
  
  #Na zbiorze uczaczym wykonuje testy statystyczne: KW, DUNNA i na ich podstawie uzyskuje cechy do SVM
  
  ### ----- Test1: KRUSKAL-WALLIS
  # Testowane hipotezy: nie sprawdzamy, która grupa się różni a jedynie czy wogóle się jakaś wyróżnia!
  #HO:  występuje równość dystrybuant rozkładów w porównywanych populacjach (grupy nie różnią się od siebie w sposób istotny statystycznie)
  #HA:  nie ma równości dystrybuant rozkładów w porównywanych populacjach (co najmniej jedna z grup różni się istotnie od innych )
  
  #1. Usuwam z 'pogrupowane' pliki ktore są w zbiorze testowym 
  grupy_nowe<-grupy[-folds1[[i]]] 
  unikalne_nowe<-unique(grupy_nowe) 
  pogrupowane_nowe <- lapply(1:length(unikalne_nowe), function(k) zbior_uczacy[,which(grupy_nowe==unikalne_nowe[k])] ) 
  
  grupy_usuniete<-grupy[folds1[[i]]] 
  
  #2.Porównuje cechę 1 pomiędzy wszystkimi grupami, 2,3....
  testKW<- function(m) {
    cecha <- lapply(1:length(unikalne), function(k) pogrupowane_nowe[[k]][m,] ) #wszystkie grupy, cecha 1... RUSZAM TYLKO CECHY BO GRUPY SPRAWDZAM WSZYSTKIE JEDNOCZESNIE !!!
    testKW_wyniki<- kruskal.test(cecha)
    p_wartoscKW<- testKW_wyniki$p.value 
  }

  #Proces zrównoleglania 
  start.time <- proc.time()
  registerDoParallel(cores=no_cores)  
  cl <- makePSOCKcluster(no_cores)  
  
  clusterExport(cl,c("testKW","unikalne","pogrupowane_nowe"))
  clusterEvalQ(cl,library(stats))
  
  p_wartosciKW<-parLapply(cl,1:dim(zbior_uczacy)[1], testKW)
  
  stop.time <- proc.time()
  run.time <- stop.time - start.time
  print(run.time)
  
  stopCluster(cl)  
  
  #3. Korekta na wielokrotne testowanie -> Benjamini Hochberg
  p_wartosciKW_PoKorekcji<- p.adjust(p_wartosciKW, method='hb')
  
  #4. Sprawdzenie które cechy mają pwartość większa niż alfa -> one nas nie interesują bo są nieróżnicujące
  ktoreH0<-which(p_wartosciKW_PoKorekcji>=0.05) #nie różnicują; brak podstaw do odrzucenia H0; p-wartośc> 0.05
  
  #5. Usunięcie ich
  skorygowane_KW<-zbior_uczacy[-ktoreH0,] #usuwam wiersze
  
  #6. Na nowo grupuje wszystkie dane ale nie będziemy już w nich mieli cech, które nie różnicowały (nie różniły się między sobą w grupie)
  pogrupowane_KW <- lapply(1:length(unikalne_nowe), function(i) skorygowane_KW[,which(grupy_nowe==unikalne_nowe[i])] ) #TERAZA PRACUJE NA TEJ ZMIENNEJ !!!
  
  
  ### ------- Test2: DUNN -> porównanie każdej pary z każdą 
  #H0: prawdopodobieństwo zaobserwowania losowej wartości w pierwszej grupie, która jest większa niż losowa wartość w drugiej grupie, wynosi połowę

  testDunn<- function(j) {
    cecha <- lapply(1:length(pogrupowane_KW), function(k) pogrupowane_KW[[k]][j,] )   #z każdej frupy biorę ceche równą il k to numer grupy
    wynikiDunn <- dunn.test(cecha,method="bh",list = "FALSE",table="FALSE") #tutaj wkładam np. cecha1 dla wszystkich grup 
  }
  
  #1. Proces zrównoleglania 
  start.time <- proc.time()
  registerDoParallel(cores=no_cores)  
  cl <- makePSOCKcluster(no_cores)  
  clusterExport(cl,c("testDunn","pogrupowane_KW"))
  clusterEvalQ(cl,library(dunn.test))
  
  wynikiDunna<- parLapply(cl,1:dim(skorygowane_KW)[1],testDunn) 
  
  stop.time <- proc.time()
  run.time <- stop.time - start.time
  print(run.time)
  stopCluster(cl)  
  
  #2. Proces intepretacji wyników (dokłanie opisany w pracy)
  kombinacje<-combn(1:18,2)
  
  roznicujacy_gen_Dunn<-list()
  
  for (d in 1:length(wynikiDunna)){  #do to iteracja po cechach 
    dana_cecha<-wynikiDunna[[d]]$P.adjusted
    cecha_decyzje<-ifelse(dana_cecha>0.05/2,0,1)  #H0=0, HA=1;  dostaje 153: 0 albo 1
    
    #teraz schodzimy niżej - do każdej kombinacji ale DLA TEJ SAMEJ CECHY 
    jest<- 0
    for (c in 1:length(unikalne)){ # po każdej grupie czyli 1:18 
      if (jest==0){ #to pomaga nie robić kolejnyh sprawdzen jesli gen jest już w liście różnicujacych 
        
        #szukam kombinacji danej grupy w obu wierszach
        numerV1<-as.numeric(which(kombinacje[1,]==c))
        numerV2<-as.numeric(which(kombinacje[2,]==c))
        numery<-c(numerV1,numerV2)
        
        #teraz przechodze do p-wartosci pod 'numery' i sprawdzam czy mamy Ha dla kazdej pary 
        wybrane_pwartosci<-cecha_decyzje[numery] #wycigam tylko kombinacje z intersujcą mnie grupą 
        ileHA<-sum(wybrane_pwartosci)
        
        #jesli roznucuje to mam sume taką jak ilość kombinacji z tą grupą 
        if (ileHA==length(numery)){
          roznicujacy_gen_Dunn[[d]]<-rownames(skorygowane_KW)[d] 
          jest<-1
        }
        
      }
      
    }
  }
  
  #3. Wyciagamy pelne dane dla wybranych genow
  pelneD<-skorygowane_KW[unlist(roznicujacy_gen_Dunn),]
  cechy_koncowe_do_svm[[i]]<-unlist(roznicujacy_gen_Dunn)
  
  
  ###---- WYBRANE CECHY IDĄ DO SVM JAKO ZBIÓR UCZĄCY 
  
  #1. Budujemy model SVM
  g<-as.data.frame(cbind(t(pelneD),grupy_nowe))
  
  f<-factor(grupy_nowe)
  fit_svm<-svm(f~.,data=t(pelneD),kernel="radial")
  
  #2. Pedykcja klas
  test_pred<-predict(fit_svm, t(zbior_testowy))
  przewidziane_klasy[[i]]<-test_pred
  
  #3. Poprawne klasy (liczone zeby zobaczyć ile klas predykowane jest poprawnie)
  model_a_prawidlowa_klasa[[i]]<-cbind(grupy_usuniete,as.data.frame(przewidziane_klasy[[i]]))
  colnames(model_a_prawidlowa_klasa[[i]])<-c("Klasy z SVM","Klasy poprawne")

  poprawne=0
  niepoprawne=0
  for (m in 1:dim(model_a_prawidlowa_klasa[[i]])[1]){
    prawidlowa<-as.character(model_a_prawidlowa_klasa[[i]][m,1])
    przewidziana<-as.character(model_a_prawidlowa_klasa[[i]][m,2])
    if (prawidlowa==przewidziana){
      poprawne<-poprawne+1
    } else {
      niepoprawne<-niepoprawne+1
    }
  }
  poprawnosc_modelu[[i]]<-cbind(poprawne,niepoprawne)
  
  closeAllConnections() 
}


#Zapis wyników 
saveRDS(cechy_koncowe_do_svm,'cechy_koncowe_do_svm.rds')
saveRDS(model_a_prawidlowa_klasa,'model_a_prawidlowa_klasa.rds')
saveRDS(poprawnosc_modelu,'poprawnosc_modelu.rds')
saveRDS(przewidziane_klasy,'przewidziane_klasy.rds')



