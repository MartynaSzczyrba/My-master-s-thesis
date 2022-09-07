# Tym kodem analizowany jest zbior:
# 1: dane z filtracja typu 50 % i bez korekty batch effect   
# 2: dane z filtracja typu 50 % i z korektą batch effect     

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

#Wczytanie danych
dane <- readRDS("ekspresja_wszystkich.rds") 
dane <-as.data.frame(dane)

## Podzial na grupy
plik<-readTargets("MILEdataInformation.txt")
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

testT<- function(m) {
  wynikiTestT<- lapply(1:dim(pogrupowane_nowe[[1]])[1], function(h) t.test(pogrupowane_nowe[[kombinacje[m,1]]][h,],pogrupowane_nowe[[kombinacje[m,2]]][h,],var.equal = T))#testy dla każdej sondy porownane między porównywanymi grupami
  statystykaT<- lapply(1:dim(pogrupowane_nowe[[1]])[1], function(j) wynikiTestT[[j]]$statistic) #t-statistic
  statystykaT_abs<-abs(unlist(statystykaT))#statistic w jednej liscie 
  znaczace100<-sort(statystykaT_abs,decreasing=T)[1:100]
  names<-unlist(lapply(1:100, function(f) rownames(zbior_uczacy[which(unlist(statystykaT_abs)==znaczace100[f]),]))) #szukam który to numer sondy
}

## Funkcje dla GA

#Funkcja inicjujaca populacje poczatkowa   
initial_population <- function(object) {
  mtx<-matrix(ncol=object@nBits,nrow=object@popSize,0)
  for (i in 1: object@popSize){
    numery<-sample(1:object@nBits,200,rep=FALSE) 
    mtx[i,numery]<-1
  }
  return(mtx)
}

#Funkcja dopasowania 
FitFun<- function(population){
  reg<-glm(data$status~. , data=data[,c(1,1+which(population==1))])
  family =binomial(logit)
  BIC.f<-BIC(reg)
  return(-BIC.f)
}

## Petla dla wszystkich procesow
## Na zbiorze uczaczym wykonywany jest test T, redukukowana ilosc cech za pomoca GA ktore nastepnie przekazywane sa do SVM 

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
  
  wynikiTestT<- mclapply(1:2,testT,mc.cores=1) 
  
  stopCluster(cl) 
  closeAllConnections()
  
  tylko_unikalne <-unique(unlist(wynikiTestT)) 
  
  cechy_testT_all[[i]]<- unlist(wynikiTestT) #czyli bede miala wszytkie nazwy w jednym wieszu wpisane dla każdej z 30 iteracji
  cechy_testT_uniq[[i]]<- tylko_unikalne    #tu mam uniklalne wybrane 
  
  ## Przeprowadzadzenie algorytmu genetycznego (GA)
  wyciagniete<- zbior_uczacy[tylko_unikalne,] #wybrane cechy z ich ekspresja 
  wyciagniete<-t(wyciagniete) #bo geny musze miec w kolumnach
  klasy<-make.names(grupy_nowe)
  status<-as.factor(klasy)
  data<-as.data.frame(cbind(wyciagniete,status))
  
  #Parametry do GA
  param_nBits=ncol(data)-1
  col_names=colnames(data[1:length(data)-1])
  
  #Przeprowadzenie GA 
  wynik<-ga(type="binary",fitness = FitFun,population = initial_population, nBits = param_nBits, maxiter = 1000, popSize = 100, names=col_names, pcrossover = 0.9, pmutation = 0.02,parallel = T)
  
  #Wyciągniecie cech i zapis wyników 
  wybrane<- col_names[wynik@solution[1,]==1]
  wybrane_do_SVM<- data[,wybrane] 
  
  wynik_GA[[i]]<-wynik
  cechy_koncowe_do_svm[[i]] <- wybrane_do_SVM
  
  ## SVM
  # Zbudowanie modelu SVM
  calosc<-zbior_uczacy[wybrane,]
  
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
  print(i) #aby wiedziec w ktorym momencie realizacji jest aktualnie kod
  
  ## Zapis wynikow do .rds 
  saveRDS(cechy_testT_all,'cechy_testT_all123.rds')
  saveRDS(cechy_testT_uniq,'cechy_testT_uniq123.rds')
  saveRDS(wynik_GA,'wynik_GA123.rds')
  saveRDS(cechy_GA,'cechy_GA123.rds')
  
  saveRDS(przewidziane_klasy,'przewidziane_klasy123.rds')
  saveRDS(model_a_prawidlowa_klasa,"model_a_prawidlowa_klasa123.rds")
  saveRDS(poprawnosc_modelu,"poprawnosc_modelu123.rds")
}


