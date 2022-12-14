### Metoda opiera sie na przeprowadzeniu testu t-Studenta a nastepnie wybraniu 100 najczestszych genow. Tak wyselekcjonowane cechy przekazywane sa do SVM. 

# Tym kodem analizowany jest zbior:
# 1: dane z filtracją typu 50 % i bez korekty batch effect   
# 2: dane z filtracją typu 50 % i z korekta batch effect     

## Pakiety 
library(e1071)
library(caret)
library(dunn.test)
library(limma)
library(sva)
library(foreach)
library(doParallel)
library(stats)
library(GA)

## Wczytanie danych
dane <- readRDS("ekspresja_wszystkich.rds") 
dane <-as.data.frame(dane)

## Podzial na grupy
plik <- readTargets("MILEdataInformation.txt")
grupy <- paste(plik$Characteristics.leukemia.class,sep='')
nazwy <- paste(plik$GSM,sep='')

unikalne<-unique(grupy) #nazwy wszyskich grup
pogrupowane <- lapply(1:length(unikalne), function(i) dane[,which(grupy==unikalne[i])] )

## Liczba watkow uzyta do zrownoleglania (kod realizowano na klastrze obliczeniowym Ziemowit)
no_cores <- 25

## Ustawienie ziarna w celu wygenerowania powtarzalnego losowego pobierania probek 
set.seed(123)

## Zmienne do zapisywania wynikow
cechy_testT_all<-list()
cechy_testT_uniq<-list()

wynik_GA<-list()
cechy_GA<-list()
cechy_koncowe_do_svm<-list()

przewidziane_klasy<-list()
model_a_prawidlowa_klasa<-list()
poprawnosc_modelu<-list()

## Podzial danych na podgrupy
gg<-t(dane)
fff<-as.data.frame(grupy)
nowe<-cbind(gg,fff)
folds <- createFolds(nowe$grupy,30)

saveRDS(folds,'fold_newT123.rds')

## Funkcja dla testu t
kombinacje<- t(combn(1:length(unikalne), 2))
#kombinacje<- t(combn(1:4, 2))
testT<- function(m) {
  wynikiTestT<- lapply(1:dim(pogrupowane_nowe[[1]])[1], function(h) t.test(pogrupowane_nowe[[kombinacje[m,1]]][h,],pogrupowane_nowe[[kombinacje[m,2]]][h,],var.equal = T))#testy dla każdej sondy porownane między porównywanymi grupami
  statystykaT<- lapply(1:dim(pogrupowane_nowe[[1]])[1], function(j) wynikiTestT[[j]]$statistic) #t-statistic
  statystykaT_abs<-abs(unlist(statystykaT))#statistic w jednej liscie 
  znaczace100<-sort(statystykaT_abs,decreasing=T)[1:100]
  names<-unlist(lapply(1:100, function(f) rownames(zbior_uczacy[which(unlist(statystykaT_abs)==znaczace100[f]),]))) #szukam który to numer sondy
}

## Petla dla wszystkich procesow
## Na zbiorze uczacym wykonywany jest test T, redukowana ilosc cech za pomoca wyboru 100 najczestszych genow ktore nastepnie przekazywane sa do SVM 

for (i in 1:30){  
  
  #Podzial na zbior uczacy i testowy 
  zbior_uczacy<-dane[(-folds[[i]])] #uczący
  zbior_testowy<-dane[folds[[i]]] #testowy
  
  #Usuwanie z 'pogrupowane' plikow ktore sa w zbiorze testowym 
  grupy_nowe<-grupy[-folds[[i]]] 
  grupy_usuniete<-grupy[folds[[i]]] 
  unikalne_nowe<-unique(grupy_nowe) 
  pogrupowane_nowe <- lapply(1:length(unikalne_nowe), function(k) zbior_uczacy[,which(grupy_nowe==unikalne_nowe[k])] ) 
  
  ## Przeprowadzenie testu t 
  
  registerDoParallel(cores=no_cores)  
  cl <- makePSOCKcluster(no_cores)  
  clusterExport(cl,c("testT","pogrupowane_nowe","kombinacje","zbior_uczacy"))
  clusterEvalQ(cl,library(stats))
  clusterEvalQ(cl,library(doParallel))
  
  #wynikiTestT<- mclapply(1:dim(kombinacje)[1],testT,mc.cores=25)
  wynikiTestT<- mclapply(1:dim(kombinacje)[1],testT,mc.cores=1) 
  
  stopCluster(cl) 
  closeAllConnections()
  
  # Wybieranie 100 najczęesciej wystepujacych cech 
  tabela<-table(unlist(wynikiTestT)) #tabela zlicza ilosci genow 
  tabela_sort<-sort(tabela,decreasing=T)[1:100,]
  unikalne100<-names(tabela_sort)
  
  cechy_testT_all[[i]]<- unlist(wynikiTestT) 
  cechy_testT_uniq[[i]]<- unikalne100 
  
  ## SVM
  
  # Zbudowanie modelu SVM
  calosc<-zbior_uczacy[unikalne100,]
  
  f<-factor(grupy_nowe)
  fit_svm<-svm(f~.,data=t(calosc),kernel="radial")
  
  # Predykcja
  test_pred<-predict(fit_svm, t(zbior_testowy))
  przewidziane_klasy[[i]]<-test_pred
  
  #Poprawne klasy 
  model_a_prawidlowa_klasa[[i]]<-cbind(grupy_usuniete,as.data.frame(przewidziane_klasy[[i]]))
  colnames(model_a_prawidlowa_klasa[[i]])<-c("Klasy poprawne","Klasy z SVM")
  
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
  
 print(i) #aby wiedziec w ktorym momencie realziacji jest aktualnie kod 
 
  ## Zapis wynikow do .rds
  saveRDS(cechy_testT_all,'cechy_testT_all123.rds')
  saveRDS(cechy_testT_uniq,'cechy_testT_uniq123.rds')
  
  saveRDS(przewidziane_klasy,'przewidziane_klasy123.rds')
  saveRDS(model_a_prawidlowa_klasa,"model_a_prawidlowa_klasa123.rds")
  saveRDS(poprawnosc_modelu,"poprawnosc_modelu123.rds")
}


