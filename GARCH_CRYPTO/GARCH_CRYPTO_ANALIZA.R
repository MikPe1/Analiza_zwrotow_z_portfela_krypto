
library(xts)
library(fBasics) 
library(tseries)
library(FinTS)
library(flextable)
library(tidyverse)
library(car) 
library(rugarch) 
library(fGarch) 
source("functions07.R")

#dane
data <- read.csv("./Dane.csv")
data$date <- as.Date(data$date)
row.names(data) <- NULL
data <- subset (data, select = -X)
in_sample<-data[1:599,]
out_sample<-data[600:630,]
nrow(out_sample)
nrow(in_sample)

max(na.omit(in_sample$total_r))
min(na.omit(in_sample$total_r))
plot(in_sample$date, in_sample$total_r, type = "l", col = "red",
     main = "Dzienne zwroty z portfela",
     xlab = "", ylab = "zwroty")
abline(h = 0, col = "black")

#Wykres ACF dla kwadratow zwrotow z portfela
acf(in_sample$total_r^2, lag.max = 36, na.action = na.pass,
    ylim = c(-0.05, 1), # dopasowanie wartosci osi y dla wykresu
    col = "darkblue", lwd = 7,
    main = "ACF kwadratów zwrotów z portfela")
#histogram
hist(in_sample$total_r, prob = T, breaks = 100, main = "Rozkład dziennych zwrotów z portfela")
# porownanie z rozkladem normalnym
curve(dnorm(x,
            mean = mean(in_sample$total_r, na.rm = T),
            sd = sd(in_sample$total_r, na.rm = T)),
      col = "darkblue", lwd = 2, add = TRUE)
#badany rozklad ma cechy rozkladu leptokurtycznego

# wykres quantile-quantile
qqnorm(in_sample$total_r)
qqline(in_sample$total_r, col = 2)
#statystyki opisowe
cbind(Metric = rownames(basicStats(in_sample$total_r)), X = basicStats(in_sample$total_r)) %>%  regulartable() %>% autofit()

 
# Wykres dziennych zwrotow z wybranego portfela cechuje grupowanie wariancji. Szczegolnie wysoka wariancja obserwowana jest na przelomie w kwietniu 2021 roku, gdy kursy wielu kryptowalut zaczely sie oslabiac.
#  
# Na wykresie dziennych zwrotów z portfolio można zauważyć grupowanie sie wariancji i okresy podwyższonego niepokoju w okolicach pandemii na przełomie 1/2 kwartału 2020 oraz na końcu okresu w kwietniu 2021, który był czasem załamań na rynku kryptowalut.
# Na podstawie analizy wizulanej stwierdzono, że zwroty charakteryzowały się zróżnicowaniem wartości co potwierdza również wykres ACF, na którym obserwujemy opóźnienia przekraczające poziom istotności. Również ocena wizualna rozkładu zwrotów.
# Wykres quantile-quantile oraz wartości statystyk wskazują na brak normalnosci oraz leptokurtyczność rozkładu.
# Ponizej dla potwierdzenia wnioskow przeprowadzona zostanie dalsza diagnostyka.

##################################
# statystyka Jarque-Bera
jarque.bera.test(na.omit(in_sample$total_r))
# H0 o normalności badanego rozkladu jest silnie odrzucana.

# statystyka Durbina-Watsona
# funkcja durbinWatsonTest() wymaga modelu liniowego
# jako argumentu
# wykorzystamy tu model skladający się tylko ze stalej
durbinWatsonTest(lm(formula = in_sample$total_r ~ 1),
                 # sprawdzm autokorelacje dla 5 pierwszych opoznien
                 max.lag = 5)
# wystepuje autokorelacja 1 i 5 rzedu.
# występowanie efektów ARCH wśród zwrotów
ArchTest(in_sample$total_r, lags = 10)
# H0 o braku efektów ARCH jest silnie odrzucana.

# statystyka DW dla kwadratów zwrotów
durbinWatsonTest(lm(formula = in_sample$total_r ^ 2 ~ 1),
                 # sprawdźmy autokorelację dla 5 pierwszych opóźnień kwadratów
                 max.lag = 10)
# opóźnienia 1,2,5 są istotne na poziomie istotnosci 5%. 
# Na poziomie istotnosci 10% istotne jest jeszcze opoznienie 3 i 4.

# Modele GARCH

# Autorzy rozwazyli estymację 4 modeli z rodziny GARCH:
# GARCH(1,1), GARCH(1,2) i GARCH(2,1), GARCH(2,2). 

#GARCH(1,1)
k.garch11 <- garchFit(formula = ~ garch(1, 1),
                      data = na.omit(in_sample$total_r),
                      include.mean = F,
                      cond.dist = "norm",
                      trace = F)

#GARCH(1,2)
k.garch12 <- garchFit(formula = ~ garch(2, 1),
                      data = na.omit(in_sample$total_r),
                      include.mean = F,
                      cond.dist = "norm",
                      trace = F)

#GARCH(2,1)
k.garch21 <- garchFit(formula = ~ garch(1, 2),
                      data = na.omit(in_sample$total_r),
                      include.mean = F,
                      cond.dist = "norm",
                      trace = F)

#GARCH(2,2)
k.garch22 <- garchFit(formula = ~ garch(2, 2),
                      data = na.omit(in_sample$total_r),
                      include.mean = F,
                      cond.dist = "norm",
                      trace = F)
#GARCH(1,1) mean=T
k.garch11m <- garchFit(formula = ~ garch(1, 1),
                       data = na.omit(in_sample$total_r),
                       include.mean = T,
                       cond.dist = "norm",
                       trace = F)
#GARCH(1,1) mean=T
k.garch11mstd <- garchFit(formula = ~ garch(1, 1),
                       data = na.omit(in_sample$total_r),
                       include.mean = T,
                       cond.dist = "std",
                       trace = F)
#GARCH(1,2) mean=T
k.garch12m <- garchFit(formula = ~ garch(2, 1),
                       data = na.omit(in_sample$total_r),
                       include.mean = T,
                       cond.dist = "norm",
                       trace = F)

#GARCH(2,1) mean=T
k.garch21m <- garchFit(formula = ~ garch(1, 2),
                       data = na.omit(in_sample$total_r),
                       include.mean = T,
                       cond.dist = "norm",
                       trace = F)

#GARCH(2,2) mean=T
k.garch22m <- garchFit(formula = ~ garch(2, 2),
                       data = na.omit(in_sample$total_r),
                       include.mean = T,
                       cond.dist = "norm",
                       trace = F)
# Wyniki dla wszystkich modeli
k.garch11@fit$matcoef
k.garch12@fit$matcoef
k.garch21@fit$matcoef
k.garch22@fit$matcoef
k.garch11m@fit$matcoef
k.garch12m@fit$matcoef
k.garch21m@fit$matcoef
k.garch22m@fit$matcoef
# GARCH(2,2) odrzucam bo jest zlym modelem dla analizowanych danych na co wskazuje p-value dla parametrow oraz kryteria informacyjne.
plot(k.garch12, which = 11)
plot(k.garch21, which = 11)
plot(k.garch11, which = 11)
plot(k.garch12m, which = 11)
plot(k.garch21m, which = 11)
plot(k.garch11m, which = 11)
compare.ICs(c("k.garch12", "k.garch21","k.garch11","k.garch12m", "k.garch21m","k.garch11m"))


# GARCH(1,1) z rownaniem sredniej jest nalepszym modelem, jest to jedyny z wybranych modeli GARCH 
# w ktorym wszystkie parametry sa istotne a kryteria informacyjne najnizsze. 
# Nastepnie nalezy przeprowadzic diagnostyke wybranego modelu w celu 
# sprawdzenia czy efekt ARCH zostal wymodelowany oraz czy rozklad 
# wystandaryzowanych reszt jest rozkladem normalnym.

######################################################

# Badanie normalnosci rozkladu wystandaryzowanych reszt
stdres <- k.garch11@residuals / sqrt(k.garch11@h.t)
hist(stdres, breaks = 20, prob = T,
     main = "Histogram wystandaryzowanych reszt, model GARCH(1,1) ")
# ggestosc rozkladu normalnego
curve(dnorm(x, mean = mean(stdres, na.rm = T),
            sd = sd(stdres,na.rm = T)),
      col = "darkblue", lwd = 2, add = TRUE)
# test Jarque-Bera
jarque.bera.test(stdres)
# normalnosc rozkladu jest silnie odrzucana w tescie Jarque Bera, 
#co jest czeste w badaniach przeprowadzonych z uzyciem modeli GARCH.


# test Durbina-Watsona
durbinWatsonTest(lm(formula = stdres^2 ~ 1),
                 max.lag = 10) # dla 10 pierwszych opóźnień
# test na występowanie efektów ARCH wśród wystandaryzowanych reszt
ArchTest(stdres, lags = 10)
plot(k.garch11, which=11)

# Badanie normalnosci rozkladu wystandaryzowanych reszt
stdres <- k.garch11m@residuals / sqrt(k.garch11m@h.t)
hist(stdres, breaks = 20, prob = T,
     main = "Histogram wystandaryzowanych reszt, model GARCH(1,1) z rownaniem sredniej ")
# ggestosc rozkladu normalnego
curve(dnorm(x, mean = mean(stdres, na.rm = T),
            sd = sd(stdres,na.rm = T)),
      col = "darkblue", lwd = 2, add = TRUE)
# test Jarque-Bera
jarque.bera.test(stdres)
# normalnosc rozkladu jest silnie odrzucana w tescie Jarque Bera, 
#co jest czeste w badaniach przeprowadzonych z uzyciem modeli GARCH.


# test Durbina-Watsona
durbinWatsonTest(lm(formula = stdres ~ 1),
                 max.lag = 10) # dla 10 pierwszych opóźnień
# test na występowanie efektów ARCH wśród wystandaryzowanych reszt
ArchTest(stdres, lags = 10)
plot(k.garch11, which=11)


# Wykres ACF dla kwadratow wystandaryzowanych reszt oraz
# Test na wystepowanie efektow ARCH (cechuje sie wysokim P-value), 
# wskazuja, ze zjawisko zostalo wymodelowane - brak efektow ARCH w kwadratach
# wystandaryzowanych reszt.


# EGARCH

#############
# EGARCH(1,1)
spec <- ugarchspec(
  # równanie warunkowej wariancji
  variance.model = list(model = "eGARCH",
                        garchOrder = c(1, 1)),
  # sGARCH oznacza standardowy model GARCH
  # równanie warunkowej średniej
  mean.model = list(armaOrder = c(1,0), include.mean = T),
  # zakładany rozkład warunkowy reszt
  distribution.model = "norm")

k.egarch11 <- ugarchfit(spec = spec, data = na.omit(in_sample$total_r))
# wyniki
k.egarch11@fit$matcoef
(info_egarch =infocriteria(k.egarch11))
plot(k.egarch11, which=11)
plot(k.egarch11, which=12)
# Wyniki dla kryteriow informacyjnych w modelu EGARCH sa w wiekszosci lepsze niz 
#w modelu GARCH. W oszacowanym modelu alpha jest nieistotna na kazdym poziomie istotnosci
# wiec warunkowa wariancja reaguje symetrycznie
# na wiadomości napływające na rynek. Wykres News-Impact
# curve wskazuje, ze wniosek jest wlasciwy.
# Wykres ACF wskazuje, ze zjawisko zostalo wymodelowane.
# Dla pewnosci przeprowadzone zostaly dalsze testy DW i ArchTest

stdres <- k.egarch11@fit$residuals / sqrt(k.egarch11@fit$var)
# test Durbina-Watsona
durbinWatsonTest(lm(formula = stdres^2 ~ 1),
                 max.lag = 10) # dla 5 pierwszych opóźnień


# test na występowanie efektów ARCH wśród wystandaryzowanych reszt
ArchTest(stdres, lags = 10)

# Przeprowadzone testy pozwalaja wnioskowac, ze efekt ARCH nie wystepuje.

####################
# AR(1)GARCH-in-Mean(1,1)

spec <- ugarchspec(
  # równanie warunkowej wariancji
  variance.model = list(model = "sGARCH",
                        garchOrder = c(1, 1)),
  # równanie warunkowej średniej - dołączamy do równania wyraz wolny
  mean.model = list(armaOrder = c(1, 0),
                    include.mean = T,
                    # dodajemy jeden element do równania warunkowej średniej,
                    # może to być odch. standardowe (pow=1) lub wariancja (pow=2)
                    archm = T, archpow = 1),
  # zakładany rozkład warunkowy reszt
  distribution.model = "norm")

# funkcja nie akceptuje braków danych, więc musimy skorzystać
# z opcji na.omit()
k.garchm11 <- ugarchfit(spec = spec, data = na.omit(in_sample$total_r))

# wyniki
k.garchm11@fit$matcoef
infocriteria(k.garchm11)

# Kryteria informacyjne w modelu GARCH-in-Mean sa nieznacznie lepsze od zwyklego modelu GARCH(1,1).

stdres <- k.garchm11@fit$residuals / sqrt(k.garchm11@fit$var)
# test Durbina-Watsona
durbinWatsonTest(lm(formula = stdres^2 ~ 1),
                 max.lag =10) # dla 5 pierwszych opóźnień
hist(stdres, breaks = 20, prob = T,
     main = "Histogram wystandaryzowanych reszt \n model AR(1)GARCH-in-Mean-t(1,1) ")
# ggestosc rozkladu normalnego
curve(dnorm(x, mean = mean(stdres, na.rm = T),
            sd = sd(stdres,na.rm = T)),
      col = "darkblue", lwd = 2, add = TRUE)
# test Jarque-Bera
jarque.bera.test(stdres)
# test na występowanie efektów ARCH wśród wystandaryzowanych reszt
ArchTest(stdres, lags = 10)
plot(k.garchm11, which=11)

# Wykres ACF i wynik testu na wystepowanie efektu ARCH pokazuja,
# ze zostal on wymodelowany.

###################
#GARCH-in-Mean(1,1) rozklad reszt - t-Studenta

spec <- ugarchspec(
  # równanie warunkowej wariancji
  variance.model = list(model = "sGARCH",
                        garchOrder = c(1, 1)),
  # równanie warunkowej średniej - dołączamy do równania wyraz wolny
  mean.model = list(armaOrder = c(1, 0),
                    include.mean = T,
                    # dodajemy jeden element do równania warunkowej średniej,
                    # może to być odch. standardowe (pow=1) lub wariancja (pow=2)
                    archm = T, archpow = 1),
  # zakładany rozkład warunkowy reszt
  distribution.model = "std")

# funkcja nie akceptuje braków danych, więc musimy skorzystać
# z opcji na.omit()
k.garchmst11 <- ugarchfit(spec = spec, data = na.omit(in_sample$total_r))

# wyniki
k.garchmst11@fit$matcoef
infocriteria(k.garchmst11)

# Kryteria informacyjne w modelu AR(1)GARCH-in-Mean-t(1,1) sa znacznie lepsze.
# Wszystkie parametry w modelu sa istotne na poziomie istotnosci 10%. 

stdres <- k.garchmst11@fit$residuals / sqrt(k.garchmst11@fit$var)
# test Durbina-Watsona
durbinWatsonTest(lm(formula = stdres^2 ~ 1),
                 max.lag = 10) # dla 5 pierwszych opóźnień

# test na występowanie efektów ARCH wśród wystandaryzowanych reszt
ArchTest(stdres, lags = 10)
plot(k.garchmst11, which=11)

# Wykres ACF i wynik testu na wystepowanie efektu ARCH pokazuja,
# ze zostal on wymodelowany.

#####
#AR(1)GARCH(1,1)
spec <- ugarchspec(
  # równanie warunkowej wariancji
  variance.model = list(model = "sGARCH",
                        garchOrder = c(1, 1)),
  # równanie warunkowej wartości oczekiwanej
  mean.model = list(armaOrder = c(0, 0),
                    include.mean = T),
  # zakładany rozkład warunkowy reszt
  distribution.model = "std") # rozkład t-Studenta

# funkcja nie akceptuje braków danych, więc musimy skorzystać
# z opcji na.omit()
k.argarch11 <- ugarchfit(spec = spec, data = na.omit(in_sample$total_r))
k.argarch11@fit$matcoef
infocriteria(k.argarch11)

# GARCH-t(1,1) ma najlepsze kryteria informacyjne z badanych modeli.
# wszystkie parametry sa istotne na poziomie istotnosci 5%
# parametr shape wskazuje na 4 stopnie swobody w rozkladzie t-studenta.

stdres <- k.argarch11@fit$residuals / sqrt(k.argarch11@fit$var)
# test Durbina-Watsona
durbinWatsonTest(lm(formula = stdres^2 ~ 1),
                 max.lag = 10) # dla 5 pierwszych opóźnień

# test na występowanie efektów ARCH wśród wystandaryzowanych reszt
ArchTest(stdres, lags = 10)
plot(k.argarch11, which=11)



#############
#GARCH-t(1,1)

spec <- ugarchspec(
  # równanie warunkowej wariancji
  variance.model = list(model = "sGARCH",
                        garchOrder = c(1, 1)),
  # równanie warunkowej wartości oczekiwanej
  mean.model = list(armaOrder = c(1, 0),
                    include.mean = T),
  # zakładany rozkład warunkowy reszt
  distribution.model = "norm") # rozkład t-Studenta

# funkcja nie akceptuje braków danych, więc musimy skorzystać
# z opcji na.omit()
k.garcht11 <- ugarchfit(spec = spec, data = na.omit(in_sample$total_r))
k.garcht11@fit$matcoef
infocriteria(k.garcht11)

hist(stdres, breaks = 20, prob = T,
     main = "Histogram wystandaryzowanych reszt, model GARCH-t(1,1)")
# ggestosc rozkladu normalnego
curve(dnorm(x, mean = mean(stdres, na.rm = T),
            sd = sd(stdres,na.rm = T)),
      col = "darkblue", lwd = 2, add = TRUE)
# test Jarque-Bera
jarque.bera.test(stdres)
# GARCH-t(1,1) ma najlepsze kryteria informacyjne z badanych modeli.
# wszystkie parametry sa istotne na poziomie istotnosci 5%
# parametr shape wskazuje na 4 stopnie swobody w rozkladzie t-studenta.

stdres <- k.garcht11@fit$residuals / sqrt(k.garcht11@fit$var)
# test Durbina-Watsona
durbinWatsonTest(lm(formula = stdres^2 ~ 1),
                 max.lag = 10) # dla 5 pierwszych opóźnień

hist(stdres, breaks = 20, prob = T,
     main = "Histogram wystandaryzowanych reszt model AR(1)GARCH-t(1,1)")
# dodajmy jeszcze gęstość rozkładu normalnego
curve(dnorm(x, mean = mean(stdres, na.rm = T),
            sd = sd(stdres,na.rm = T)),
      col = "darkblue", lwd = 2, add = TRUE)

#statystyki opisowe
basicStats(stdres)

# test Jarque-Bery
jarque.bera.test(stdres)

# test na występowanie efektów ARCH wśród wystandaryzowanych reszt
ArchTest(stdres, lags = 10)
plot(k.garcht11, which=11)

# Wykres ACF i wynik testu na wystepowanie efektu ARCH pokazuja,
# ze zostal on wymodelowany.


############
#TGARCH(1,1) (reszty z rozkadu t-studenta)

spec <- ugarchspec(
  # równanie warunkowej wariancji
  variance.model = list(model = "fGARCH",
                        garchOrder = c(1, 1),
                        submodel = "TGARCH"),
  # model="fGARCH" (family GARCH) razem z submodel="TGARCH"
  # równanie warunkowej średniej
  mean.model = list(armaOrder = c(1, 0),
                    include.mean = T),
  # zakładany rozkład warunkowy reszt
  distribution.model = "norm")

# funkcja nie akceptuje braków danych, więc musimy skorzystać
# z opcji na.omit()
k.tgarcht11 <- ugarchfit(spec = spec, data = na.omit(in_sample$total_r))
# wyniki
k.tgarcht11@fit$matcoef
infocriteria(k.tgarcht11)

plot(k.tgarcht11, which=11)

stdres <- k.tgarcht11@fit$residuals / sqrt(k.tgarcht11@fit$var)
# test Durbina-Watsona
durbinWatsonTest(lm(formula = stdres^2 ~ 1),
                 max.lag = 5) # dla 5 pierwszych opóźnień

# test na występowanie efektów ARCH wśród wystandaryzowanych reszt
ArchTest(stdres, lags = 10)

# W modelu TGARCH(1,1) parametr eta jest nieistotny, wykres ACF dla kwadratow standaryzowanych reszt
# nie jest czysty, test na wystepowanie efektu ARCH takze wskazuje, ze zjawisko nie zostalo
# wymodelowane, poniewaz hipoteza o braku efektu ARCH jest silnie odrzucana.


in_sample = na.omit(in_sample)
plot(in_sample$date, in_sample$total_r, type = "l", col = "red",
     main = "Dzienne zwroty z portfela",
     xlab = "", ylab = "zwroty")
abline(h = 0, col = "black")
plot(in_sample$date, cumsum(in_sample$total_r), type = "l", col = "blue",
     main = "Skumulowane zwroty z portfela",
     xlab = "", ylab = "zwroty")



  plot(in_sample$date, k.garch11@h.t, type = "l", col = "darkgreen",
       ylab = "cvar", xlab="data", main = "Warunkowa wariancja w modelu \n GARCH(1,1)")
  plot(in_sample$date, k.garch11m@h.t, type = "l", col = "darkgreen",
       ylab = "cvar", xlab="data", main = "Warunkowa wariancja w modelu \n GARCH(1,1) mean=T")
  plot(in_sample$date, k.garch11mstd@h.t, type = "l", col = "darkgreen",
       ylab = "cvar", xlab="data", main = "Warunkowa wariancja w modelu \n GARCH(1,1) mean=T dist=std")
  plot(in_sample$date, k.garch11std@h.t, type = "l", col = "darkgreen",
       ylab = "cvar", xlab="data", main = "Warunkowa wariancja w modelu \n GARCH(1,1) dist=std")
  plot(in_sample$date, k.argarch11@h.t, type = "l", col = "darkgreen",
       ylab = "cvar", xlab="data", main = "Warunkowa wariancja w modelu \n GARCH(1,1) dist=std")
  
  
  
  
  plot(in_sample$date, k.egarch11@fit$var, type = "l", col = "darkgreen",
       ylab = "cvar", xlab="data", main = "Warunkowa wariancja w modelu \n EGARCH(1,1)")
  
  plot(in_sample$date, k.garchm11@fit$var, type = "l", col = "darkgreen",
       ylab = "cvar", xlab="data", main = "Warunkowa wariancja w modelu \n AR(1) GARCH-in-Mean(1,1)")
  
  plot(in_sample$date, k.garcht11@fit$var, type = "l", col = "darkgreen",
       ylab = "cvar", xlab="data", main = "Warunkowa wariancja w modelu \n AR(1) GARCH(1,1)")
    

#ANALIZA VAR
  
  data = na.omit(data)
  data$date = as.Date(data$date)
  data$nobs <- 1:length(data$total_r)
  data$total_rstd = (data$total_r - mean(data$total_r, na.rm=T))/sd(data$total_r, na.rm=T)
  in_samp= data[1:599,]
  in_samp$date = as.Date(in_samp$date)
  poczatek <- data$nobs[data$date == as.Date("2022-04-27")]
  koniec  <- data$nobs[data$date == as.Date("2022-05-27")]
  data_forecast <- data[poczatek:koniec, ]
  data_forecast$date = as.Date(data_forecast$date)
  basicStats(data$total_rstd)
  
  q01 <- quantile(in_samp$total_rstd, 0.01, na.rm = T)
  q01
  
  qnorm(0.01, 0, 1)
  
  spec <- ugarchspec(variance.model = list(model = "sGARCH",
                                           garchOrder = c(1, 1)),
                     mean.model = list(armaOrder = c(1, 0),
                                       include.mean = T),
                     distribution.model = "norm")
  
  in_samp.garch11 <- ugarchfit(spec = spec,
                               data = na.omit(in_samp$total_r))
  plot(in_samp.garch11, which = 10)
  
  in_samp$VaR <- q01 * in_samp.garch11@fit$sigma
  tail(in_samp)
  
  plot(in_samp$date, in_samp$total_r)
  plot(in_samp$date, in_samp$total_r, col = "red", lwd = 1, type = 'l',
       ylim = c(-0.4, 0.4))
  abline(h = 0, lty = 2)
  lines(in_samp$date, in_samp$VaR, type = 'l', col = "green")
  
  (sum(in_samp$total_r < in_samp$VaR) / length(in_samp$VaR))*100
############
# GARCH (1,1)


VaR <- rep(NA, times = koniec - poczatek + 1)
# 
for (k in poczatek:koniec) {
    tmp.data <- data[data$nobs <= (k - 1), ]
    tmp.data <- tmp.data[as.Date("2020-09-03") <= tmp.data$date, ]
    tmp.data$rstd <- (tmp.data$total_r - mean(tmp.data$total_r, na.rm = T)) /
        sd(tmp.data$total_r, na.rm = T)
    q01 <- quantile(tmp.data$rstd, 0.01, na.rm = T)
    spec <- ugarchspec(variance.model = list(model = "sGARCH",
                                             garchOrder = c(1, 1)),
                       mean.model = list(armaOrder = c(0, 0),
                                         include.mean = F),
                       distribution.model = "norm")
    tmp.garch11 <- ugarchfit(spec = spec, data = na.omit(tmp.data$total_r))
    sigma.forecast  <- ugarchforecast(tmp.garch11, n.ahead = 1)
    sigma.forecast2 <- sigma.forecast@forecast$sigmaFor[1, 1]
    VaR[k - poczatek + 1] <- q01 * sigma.forecast2
}

# przypisanie oszacowania modelu VaR do zbioru danych
data_forecast$VaR <- VaR

data_forecast = na.omit(data_forecast)
# wykres zwrotóW w porownaniu z wynikami modelu VAR w probie out of sample
plot(data_forecast$date, data_forecast$total_r, col = "red", lwd = 1, type = 'l', 
     xlab="data", 
     ylab="VaR i zwroty z portfela",
     main = 'VaR Out-Of-Sample GARCH(1,1)',
     ylim = c(-0.4, 0.4))
abline(h = 0, lty = 2)
lines(data_forecast$date, data_forecast$VaR, type = 'l', col = "green")
# w ilu przypadkach straty przekroczyły zakładany poziom VaR?
paste0("Wartość VaR przekroczona ",sum(data_forecast$total_r < data_forecast$VaR)," razy (",
       round(sum(data_forecast$total_r < data_forecast$VaR) / length(data_forecast$VaR),2),"%).")

#GARCH(1,1) mean=T
poczatek <- data$nobs[data$date == as.Date("2022-04-26")]
koniec  <- data$nobs[data$date == as.Date("2022-05-26")]
data_forecast <- data[poczatek:koniec, ]
VaR <- rep(NA, times = koniec - poczatek + 1)
# 
for (k in poczatek:koniec) {
  tmp.data <- data[data$nobs <= (k - 1), ]
  tmp.data <- tmp.data[as.Date("2020-09-03") <= tmp.data$date, ]
  tmp.data$rstd <- (tmp.data$total_r - mean(tmp.data$total_r, na.rm = T)) /
    sd(tmp.data$total_r, na.rm = T)
  q01 <- quantile(tmp.data$rstd, 0.01, na.rm = T)
  spec <- ugarchspec(variance.model = list(model = "sGARCH",
                                           garchOrder = c(1, 1)),
                     mean.model = list(armaOrder = c(0, 0),
                                       include.mean = T),
                     distribution.model = "norm")
  tmp.garch11 <- ugarchfit(spec = spec, data = na.omit(tmp.data$total_r))
  sigma.forecast  <- ugarchforecast(tmp.garch11, n.ahead = 1)
  sigma.forecast2 <- sigma.forecast@forecast$sigmaFor[1, 1]
  VaR[k - poczatek + 1] <- q01 * sigma.forecast2
}

# przypisanie oszacowania modelu VaR do zbioru danych
data_forecast$VaR <- VaR

data_forecast = na.omit(data_forecast)
# wykres zwrotóW w porownaniu z wynikami modelu VAR w probie out of sample
plot(data_forecast$date, data_forecast$total_r, col = "red", lwd = 1, type = 'l', 
     xlab="data", 
     ylab="VaR i zwroty z portfela",
     main = 'VaR Out-Of-Sample GARCH(1,1) mean = T',
     ylim = c(-0.4, 0.4))
abline(h = 0, lty = 2)
lines(data_forecast$date, data_forecast$VaR, type = 'l', col = "green")
# w ilu przypadkach straty przekroczyły zakładany poziom VaR?
paste0("Wartość VaR przekroczona ",sum(data_forecast$total_r < data_forecast$VaR)," razy (",
       round(sum(data_forecast$total_r < data_forecast$VaR) / length(data_forecast$VaR),2),"%).")

#GARCH(1,1) mean = T dist = std

poczatek <- data$nobs[data$date == as.Date("2022-04-26")]
koniec  <- data$nobs[data$date == as.Date("2022-05-26")]
data_forecast <- data[poczatek:koniec, ]
VaR <- rep(NA, times = koniec - poczatek + 1)
# 
for (k in poczatek:koniec) {
  tmp.data <- data[data$nobs <= (k - 1), ]
  tmp.data <- tmp.data[as.Date("2020-09-03") <= tmp.data$date, ]
  tmp.data$rstd <- (tmp.data$total_r - mean(tmp.data$total_r, na.rm = T)) /
    sd(tmp.data$total_r, na.rm = T)
  q01 <- quantile(tmp.data$rstd, 0.01, na.rm = T)
  spec <- ugarchspec(variance.model = list(model = "sGARCH",
                                           garchOrder = c(1, 1)),
                     mean.model = list(armaOrder = c(0, 0),
                                       include.mean = T),
                     distribution.model = "std")
  tmp.garch11 <- ugarchfit(spec = spec, data = na.omit(tmp.data$total_r))
  sigma.forecast  <- ugarchforecast(tmp.garch11, n.ahead = 1)
  sigma.forecast2 <- sigma.forecast@forecast$sigmaFor[1, 1]
  VaR[k - poczatek + 1] <- q01 * sigma.forecast2
}

# przypisanie oszacowania modelu VaR do zbioru danych
data_forecast$VaR <- VaR

data_forecast = na.omit(data_forecast)
# wykres zwrotóW w porownaniu z wynikami modelu VAR w probie out of sample
plot(data_forecast$date, data_forecast$total_r, col = "red", lwd = 1, type = 'l', 
     xlab="data", 
     ylab="VaR i zwroty z portfela",
     main = 'VaR Out-Of-Sample GARCH(1,1) mean = T, dist = std',
     ylim = c(-0.4, 0.4))
abline(h = 0, lty = 2)
lines(data_forecast$date, data_forecast$VaR, type = 'l', col = "green")
# w ilu przypadkach straty przekroczyły zakładany poziom VaR?
paste0("Wartość VaR przekroczona ",sum(data_forecast$total_r < data_forecast$VaR)," razy (",
       round(sum(data_forecast$total_r < data_forecast$VaR) / length(data_forecast$VaR),2),"%).")

#GARCH(1,1) mean = T AR(1) dist = std

poczatek <- data$nobs[data$date == as.Date("2022-04-26")]
koniec  <- data$nobs[data$date == as.Date("2022-05-26")]
data_forecast <- data[poczatek:koniec, ]
VaR <- rep(NA, times = koniec - poczatek + 1)
# 
for (k in poczatek:koniec) {
  tmp.data <- data[data$nobs <= (k - 1), ]
  tmp.data <- tmp.data[as.Date("2020-09-03") <= tmp.data$date, ]
  tmp.data$rstd <- (tmp.data$total_r - mean(tmp.data$total_r, na.rm = T)) /
    sd(tmp.data$total_r, na.rm = T)
  q01 <- quantile(tmp.data$rstd, 0.01, na.rm = T)
  spec <- ugarchspec(variance.model = list(model = "sGARCH",
                                           garchOrder = c(1, 1)),
                     mean.model = list(armaOrder = c(1, 0),
                                       include.mean = T),
                     distribution.model = "std")
  tmp.garch11 <- ugarchfit(spec = spec, data = na.omit(tmp.data$total_r))
  sigma.forecast  <- ugarchforecast(tmp.garch11, n.ahead = 1)
  sigma.forecast2 <- sigma.forecast@forecast$sigmaFor[1, 1]
  VaR[k - poczatek + 1] <- q01 * sigma.forecast2
}

# przypisanie oszacowania modelu VaR do zbioru danych
data_forecast$VaR <- VaR

data_forecast = na.omit(data_forecast)
# wykres zwrotóW w porownaniu z wynikami modelu VAR w probie out of sample
plot(data_forecast$date, data_forecast$total_r, col = "red", lwd = 1, type = 'l', 
     xlab="data", 
     ylab="VaR i zwroty z portfela",
     main = 'VaR Out-Of-Sample GARCH(1,1) mean = AR(1), dist = std',
     ylim = c(-0.4, 0.4))
abline(h = 0, lty = 2)
lines(data_forecast$date, data_forecast$VaR, type = 'l', col = "green")
# w ilu przypadkach straty przekroczyły zakładany poziom VaR?
paste0("Wartość VaR przekroczona ",sum(data_forecast$total_r < data_forecast$VaR)," razy (",
       round(sum(data_forecast$total_r < data_forecast$VaR) / length(data_forecast$VaR),2),"%).")

#GARCH(1,1) mean = T AR(1)

poczatek <- data$nobs[data$date == as.Date("2022-04-26")]
koniec  <- data$nobs[data$date == as.Date("2022-05-26")]
data_forecast <- data[poczatek:koniec, ]
VaR <- rep(NA, times = koniec - poczatek + 1)
# 
for (k in poczatek:koniec) {
  tmp.data <- data[data$nobs <= (k - 1), ]
  tmp.data <- tmp.data[as.Date("2020-09-03") <= tmp.data$date, ]
  tmp.data$rstd <- (tmp.data$total_r - mean(tmp.data$total_r, na.rm = T)) /
    sd(tmp.data$total_r, na.rm = T)
  q01 <- quantile(tmp.data$rstd, 0.01, na.rm = T)
  spec <- ugarchspec(variance.model = list(model = "sGARCH",
                                           garchOrder = c(1, 1)),
                     mean.model = list(armaOrder = c(1, 0),
                                       include.mean = T),
                     distribution.model = "norm")
  tmp.garch11 <- ugarchfit(spec = spec, data = na.omit(tmp.data$total_r))
  sigma.forecast  <- ugarchforecast(tmp.garch11, n.ahead = 1)
  sigma.forecast2 <- sigma.forecast@forecast$sigmaFor[1, 1]
  VaR[k - poczatek + 1] <- q01 * sigma.forecast2
}

# przypisanie oszacowania modelu VaR do zbioru danych
data_forecast$VaR <- VaR

data_forecast = na.omit(data_forecast)
# wykres zwrotóW w porownaniu z wynikami modelu VAR w probie out of sample
plot(data_forecast$date, data_forecast$total_r, col = "red", lwd = 1, type = 'l', 
     xlab="data", 
     ylab="VaR i zwroty z portfela",
     main = 'VaR Out-Of-Sample GARCH(1,1) mean = AR(1)',
     ylim = c(-0.4, 0.4))
abline(h = 0, lty = 2)
lines(data_forecast$date, data_forecast$VaR, type = 'l', col = "green")
# w ilu przypadkach straty przekroczyły zakładany poziom VaR?
paste0("Wartość VaR przekroczona ",sum(data_forecast$total_r < data_forecast$VaR)," razy (",
       round(sum(data_forecast$total_r < data_forecast$VaR) / length(data_forecast$VaR),2),"%).")

################
#EGARCH (1,1)

data$nobs <- 1:nrow(data)
poczatek <- data$nobs[data$date == as.Date("2022-04-26")]
koniec  <- data$nobs[data$date == as.Date("2022-05-26")]
data_forecast <- data[poczatek:koniec, ]
VaR <- rep(NA, times = koniec - poczatek + 1)


for (k in poczatek:koniec) {
    tmp.data <- data[data$nobs <= (k - 1), ]
    tmp.data <- tmp.data[as.Date("2020-09-03") <= tmp.data$date, ]
    tmp.data$rstd <- (tmp.data$total_r - mean(tmp.data$total_r, na.rm = T)) /
        sd(tmp.data$total_r, na.rm = T)
    q01 <- quantile(tmp.data$rstd, 0.01, na.rm = T)
    spec <- ugarchspec(
  # równanie warunkowej wariancji
  variance.model = list(model = "eGARCH",
                        garchOrder = c(1, 1)),
  # sGARCH oznacza standardowy model GARCH
  # równanie warunkowej średniej
  mean.model = list(armaOrder = c(0, 0), include.mean = F),
  # zakładany rozkład warunkowy reszt
  distribution.model = "norm")
    tmp.garch11 <- ugarchfit(spec = spec, data = na.omit(tmp.data$total_r))
    sigma.forecast  <- ugarchforecast(tmp.garch11, n.ahead = 1)
    sigma.forecast2 <- sigma.forecast@forecast$sigmaFor[1, 1]
    VaR[k - poczatek + 1] <- q01 * sigma.forecast2
}

# przypisanie oszacowania modelu VaR do zbioru danych
data_forecast$VaR <- VaR
# wykres zwrotóW w porownaniu z wynikami modelu VAR w probie out of sample
plot(data_forecast$date, data_forecast$total_r, col = "black", lwd = 1, type = 'l', 
     xlab="data", 
     ylab="VaR i zwroty z portfela",
     main = 'VaR Out-Of-Sample EGARCH(1,1)',
     ylim = c(-0.4, 0.4))
abline(h = 0, lty = 2)
lines(data_forecast$date, data_forecast$VaR, type = 'l', col = "blue")
# w ilu przypadkach straty przekroczyły zakładany poziom VaR?
paste0("Wartość VaR przekroczona ",sum(data_forecast$total_r < data_forecast$VaR)," razy (",
       round(sum(data_forecast$total_r < data_forecast$VaR) / length(data_forecast$VaR),2),"%).")
####################
#AR(1)GARCH-in-Mean-t (1,1)


for (k in poczatek:koniec) {
    tmp.data <- data[data$nobs <= (k - 1), ]
    tmp.data <- tmp.data[as.Date("2020-09-03") <= tmp.data$date, ]
    tmp.data$rstd <- (tmp.data$total_r - mean(tmp.data$total_r, na.rm = T)) /
        sd(tmp.data$total_r, na.rm = T)
    q01 <- quantile(tmp.data$rstd, 0.01, na.rm = T)
    spec <- ugarchspec(
      # równanie warunkowej wariancji
      variance.model = list(model = "sGARCH",
                            garchOrder = c(1, 1)),
      # równanie warunkowej średniej - dołączamy do równania wyraz wolny
      mean.model = list(armaOrder = c(1, 0),
                        include.mean = T,
                        # dodajemy jeden element do równania warunkowej średniej,
                        # może to być odch. standardowe (pow=1) lub wariancja (pow=2)
                        archm = T, archpow = 1),
      # zakładany rozkład warunkowy reszt
      distribution.model = "std")

    tmp.garch11 <- ugarchfit(spec = spec, data = na.omit(tmp.data$total_r))
    sigma.forecast  <- ugarchforecast(tmp.garch11, n.ahead = 1)
    sigma.forecast2 <- sigma.forecast@forecast$sigmaFor[1, 1]
    VaR[k - poczatek + 1] <- q01 * sigma.forecast2
}

# przypisanie oszacowania modelu VaR do zbioru danych
data_forecast$VaR <- VaR
# wykres zwrotóW w porownaniu z wynikami modelu VAR w probie out of sample
plot(data_forecast$date, data_forecast$total_r, col = "black", lwd = 1, type = 'l', 
     xlab="data", 
     ylab="VaR i zwroty z portfela",
     main = 'VaR Out-Of-Sample GARCH-in-Mean(1,1)',
     ylim = c(-0.4, 0.4))
abline(h = 0, lty = 2)
lines(data_forecast$date, data_forecast$VaR, type = 'l', col = "blue")
# w ilu przypadkach straty przekroczyły zakładany poziom VaR?
paste0("Wartość VaR przekroczona ",sum(data_forecast$total_r < data_forecast$VaR)," razy (",
       round(sum(data_forecast$total_r < data_forecast$VaR) / length(data_forecast$VaR),2),"%).")

################
# GARCH-t(1,1)

for (k in poczatek:koniec) {
  tmp.data <- data[data$nobs <= (k - 1), ]
  tmp.data <- tmp.data[as.Date("2020-09-03") <= tmp.data$date, ]
  tmp.data$rstd <- (tmp.data$total_r - mean(tmp.data$total_r, na.rm = T)) /
    sd(tmp.data$total_r, na.rm = T)
  q01 <- quantile(tmp.data$rstd, 0.01, na.rm = T)
  spec <- ugarchspec(
    # równanie warunkowej wariancji
    variance.model = list(model = "sGARCH",
                          garchOrder = c(1, 1)),
    # równanie warunkowej wartości oczekiwanej
    mean.model = list(armaOrder = c(0, 0),
                      include.mean = F),
    # zakładany rozkład warunkowy reszt
    distribution.model = "std") # rozkład t-Studenta
  
  tmp.garch11 <- ugarchfit(spec = spec, data = na.omit(tmp.data$total_r))
  sigma.forecast  <- ugarchforecast(tmp.garch11, n.ahead = 1)
  sigma.forecast2 <- sigma.forecast@forecast$sigmaFor[1, 1]
  VaR[k - poczatek + 1] <- q01 * sigma.forecast2
}

data_forecast$VaR <- VaR
# wykres zwrotóW w porownaniu z wynikami modelu VAR w probie out of sample
plot(data_forecast$date, data_forecast$total_r, col = "black", lwd = 1, type = 'l', 
     xlab="data", 
     ylab="VaR i zwroty z portfela",
     main = 'VaR Out-Of-Sample GARCH-t(1,1)',
     ylim = c(-0.4, 0.4))
abline(h = 0, lty = 2)
lines(data_forecast$date, data_forecast$VaR, type = 'l', col = "blue")
# w ilu przypadkach straty przekroczyły zakładany poziom VaR?
paste0("Wartość VaR przekroczona ",sum(data_forecast$total_r < data_forecast$VaR)," razy (",
       round(sum(data_forecast$total_r < data_forecast$VaR) / length(data_forecast$VaR),2),"%).")


###############
#TGARCH (1,1)
for (k in poczatek:koniec) {
    tmp.data <- data[data$obs <= (k - 1), ]
    tmp.data <- tmp.data[as.Date("2020-09-03") <= tmp.data$date, ]
    tmp.data$rstd <- (tmp.data$total_r - mean(tmp.data$total_r, na.rm = T)) /
        sd(tmp.data$total_r, na.rm = T)
    q01 <- quantile(tmp.data$rstd, 0.01, na.rm = T)
    spec <- ugarchspec(
  # równanie warunkowej wariancji
  variance.model = list(model = "fGARCH",
                        garchOrder = c(1, 1),
                        submodel = "TGARCH"),
  # model="fGARCH" (family GARCH) razem z submodel="TGARCH"
  # równanie warunkowej średniej
  mean.model = list(armaOrder = c(0, 0),
                    include.mean = F),
  # zakładany rozkład warunkowy reszt
  distribution.model = "std")

    tmp.garch11 <- ugarchfit(spec = spec, data = na.omit(tmp.data$total_r))
    sigma.forecast  <- ugarchforecast(tmp.garch11, n.ahead = 1)
    sigma.forecast2 <- sigma.forecast@forecast$sigmaFor[1, 1]
    VaR[k - poczatek + 1] <- q01 * sigma.forecast2
}

# przypisujemy oszacowania VaRów do zbioru
data_forecast$VaR <- VaR
# wyres zwrotóW vs. VaR w okresie OUT-OF-SAMPLE
plot(data_forecast$date, data_forecast$total_r, col = "black", lwd = 1, type = 'l', 
     xlab="data", 
     ylab="VaR i zwroty z portfela",
     main = 'VaR Out-Of-Sample TGARCH(1,1)',
     ylim = c(-0.4, 0.4))
abline(h = 0, lty = 2)
lines(data_forecast$date, data_forecast$VaR, type = 'l', col = "blue")
# w ilu przypadkach straty przekroczyły zakładany poziom VaR?
paste0("Wartość VaR przekroczona ",sum(data_forecast$total_r < data_forecast$VaR)," razy (",
       round(sum(data_forecast$total_r < data_forecast$VaR) / length(data_forecast$VaR),2),"%).")

#WNIOSKI
#W badaniu przeanalizowane zostały modele GARCH do parametrów p=2 i q = 2 zarówno z uwzględnieniem 
#równania średniej jak i bez tego równania. Wyestymowane zostały także modele GARCH-t, GARCH-in-Mean, 
#EGARCH oraz TGARCH a następnie wybrane te z nich, które cechowały się wymodelowaniem efektu 
#ARCH z badanego portfela oraz najlepszymi kryteriami informacyjnymi. Dla tak wyselekcjonowanych modeli 
#zostały przedstawione wykresy funkcji warunkowej wariancji dla dwóch okresów in sample – jednym 
#dłuższym i drugim krótszym. Wyniki uzyskane za pomocą modeli różnią się nieznacznie w przypadku 
#dłuższego okresu in sample, w krótszym różnice są zauważalne już na samym wykresie funkcji. 
#Dla najlepszych modeli – GARCH(1,1) z równaniem średniej, AR(1)GARCH-in-Mean-t(1,1), 
#AR(1)GARCH(1,1) oraz AR(1)GARCH-t(1,1) został oszacowany model Value at Risk. Zarówno dla dłuższego 
#jak i krótszego okresu in sample wyniki pokazują, że model niedoszacował straty w dwóch przypadkach,
#co pokazuje, że model jest dobry a przekroczenie wartości VaR było związane z wyjątkową losowością.
