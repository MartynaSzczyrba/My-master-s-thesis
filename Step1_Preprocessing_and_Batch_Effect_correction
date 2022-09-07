setwd("C:/Users/marty/Desktop/Praca magisterska/analiza_w_R/dane") 

## Instalacja pakietow
#install.packages("BiocManager")
#BiocManager::install("limma")
library("limma")

#BiocManager::install("GEOquery")
library("GEOquery")

#BiocManager::install("affy")
library(affy)

#BiocManager::install("hgu133plus2frmavecs") 
library(hgu133plus2frmavecs)

#BiocManager::install("ff")
library(ff)

#BiocManager::install("frma") 
library(frma)

#BiocManager::install("sva") 
library(sva)

#BiocManager::install("foreach") 
library(foreach)

#BiocManager::install("rlist") 
library(rlist)

#BiocManager::install("doParallel") 
library(doParallel)

#BiocManager::install("stats") 
library(stats)

## Wczytanie pliku zawierajacego miedzy innymi nazwy plikow CEL na ktorych bedzie opierala sie cala analiza 
plik<-readTargets("MILEdataInformation.txt")
pelneNazwyPlikow<-paste(plik$GSM,'.CEL',sep='')
nazwyPlikow<-paste(plik$GSM,'',sep='')

## Pobieranie i rozpakowywanie danych
for (i in nazwyPlikow){
  gse=getGEOSuppFiles(i, makeDirectory = FALSE)
  gunzip(paste(i,".CEL.gz",sep = ""),paste(i, ".CEL", sep=""), remove = TRUE)
  unlink(paste(i, ".CHP.gz", sep=""))
}

## Wgranie Custom CDF z BrainArrary http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp#v19
#install.packages("hgu133plus2hsentrezgprobe_19.0.0.tar.gz",repos=NULL,getwd=sciezka_k)
library("hgu133plus2hsentrezgprobe")

#install.packages("hgu133plus2hsentrezgcdf_19.0.0.tar.gz",repos=NULL,getwd=sciezka_k)
library("hgu133plus2hsentrezgcdf")

#install.packages("hgu133plus2hsentrezg.db_19.0.0.tar.gz",repos=NULL,getwd=sciezka_k)
library("hgu133plus2hsentrezg.db")

## Wgranie danych i przeprowadzenie preprocesingu za pomocą frma, ktore obejmuje korekte tla, normalizacje oraz sumaryzacje
zakres_danych <- 1:length(pelneNazwyPlikow)

preprocesing <- function(i) {
  dane<-ReadAffy(filenames = pelneNazwyPlikow[i])
  dane@annotation <-"hgu133plus2hsentrezg" 
  dane@cdfName <- "hgu133plus2hsentrezg" 
  frma(dane,background="rma", normalize="quantile", summarize="median_polish", target="probeset", input.vecs=NULL, output.param=NULL, verbose=FALSE)
}
dane_po_prepocesingu <- foreach(i=zakres_danych) %dopar% {preprocesing(i)} #umożliwia wykonanie petli bez pętli 

## Usuwanie efektu paczki
instytucje <- paste(plik$Institution,sep='') # instytucje potrzebne do usunięcia efektu paczki 
edane<-exprs(po_preprocesingu[[1]]) 
for (i  in 2:length(zakres_danych)) {
  nowa<- exprs(po_preprocesingu[[i]])
  edane<-cbind(edane, nowa)
}

batch_effect <- ComBat(dat=edane, batch=instytucje, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

## Zapisywanie wstepnnie przetworzonych danych 
saveRDS(edane, file = "ekspresja_wszystkich.rds")
saveRDS(batch_effect, file = "z_korekta_batch_effect.rds")

