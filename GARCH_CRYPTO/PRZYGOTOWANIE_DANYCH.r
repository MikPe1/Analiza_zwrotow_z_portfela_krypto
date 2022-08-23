library(tidyverse)
library(rvest)
library(xts)
library(urca)
library(lmtest)
library(fBasics)
library(tseries)
library(car)
library(FinTS) 
library(fGarch) 
library(rugarch) 
source("functions07.R")



#zaladowanie danych o poszczegolnych kryptowalutach

#MONERO
MONERO <- read.csv("MONERO.CSV",
                  header = TRUE,
                  sep = ";",
                  dec = ",",
                  stringsAsFactors = F)
MONERO$Date <- as.Date(MONERO$Date)
MONERO <- MONERO[, c("Date", "Close", "Market.Cap")]
MONERO <- MONERO[MONERO$Date >= as.Date("2020-09-05"), ]

#POLYGON
POLYGON <- read.csv("POLYGON.CSV",
                   header = TRUE,
                   sep = ";",
                   dec = ",",
                   stringsAsFactors = F)
POLYGON$Date <- as.Date(POLYGON$Date)
POLYGON <- POLYGON[, c("Date", "Close", "Market.Cap")]
POLYGON <- POLYGON[POLYGON$Date >= as.Date("2020-09-05"), ]

#MANA (DECENTRALAND)
MANA <- read.csv("MANA.CSV",
                   header = TRUE,
                   sep = ";",
                   dec = ",",
                   stringsAsFactors = F)
MANA$Date <- as.Date(MANA$Date)
MANA <- MANA[, c("Date", "Close", "Market.Cap")]
MANA <- MANA[MANA$Date >= as.Date("2020-09-05"), ]

#KUSAMA
KUSAMA <- read.csv("KUSAMA.CSV",
                   header = TRUE,
                   sep = ";",
                   dec = ",",
                   stringsAsFactors = F)
KUSAMA$Date <- as.Date(KUSAMA$Date)
KUSAMA <- KUSAMA[, c("Date", "Close", "Market.Cap")]
KUSAMA <- KUSAMA[KUSAMA$Date >= as.Date("2020-09-05"), ]

head(KUSAMA, n = 10)
tail(KUSAMA, n = 10)


colnames(MANA) <- paste0("MANA_",colnames(MANA))
colnames(POLYGON) <- paste0("POLYGON_",colnames(POLYGON))
colnames(MONERO) <- paste0("MONERO_",colnames(MONERO))
colnames(KUSAMA) <- paste0("KUSAMA_",colnames(KUSAMA))

#dopasowanie do dlugosci danych MONERO
range <- 1:nrow(MONERO)
#Numery kolumn Close i Market.Cap
kol <- c(2,3)
#Laczenie danych
data <- data.frame(id = range, date = MANA$MANA_Date[range])
data <- cbind(data,
      MANA[range,kol],
      POLYGON[range,kol],
      MONERO[range,kol],
      KUSAMA[range,kol])

data$total_cap = data$MANA_Market.Cap + data$POLYGON_Market.Cap + data$MONERO_Market.Cap + data$KUSAMA_Market.Cap

data$MANA_p = (data$MANA_Market.Cap/data$total_cap)
data$POLYGON_p = (data$POLYGON_Market.Cap/data$total_cap)
data$MONERO_p = (data$MONERO_Market.Cap/data$total_cap)
data$KUSAMA_p = (data$KUSAMA_Market.Cap/data$total_cap)

#Obliczanie logarytmicznych stóp zwrotu dla poszczegolnych kryptowalut

data$MANA_r = diff.xts(log(data$MANA_Close))
data$POLYGON_r = diff.xts(log(data$POLYGON_Close))
data$MONERO_r = diff.xts(log(data$MONERO_Close))
data$KUSAMA_r = diff.xts(log(data$KUSAMA_Close))

#Obliczanie stop zwrotu z calego portfela
data$total_r = (data$MANA_p*data$MANA_r
                + data$POLYGON_p*data$POLYGON_r
                + data$MONERO_p*data$MONERO_r
                + data$KUSAMA_p*data$KUSAMA_r)

#zapisuje
write.csv(data,"./Dane.csv")

#Rysowanie wykresow i stop zwrotu dla poszczegolnych zmiennych
par(mfrow=c(2,2)) 
plot(data$date,data$MANA_Close, type = "l")
plot(data$date,data$MANA_r, type = "l",ylim = c(-0.4, 0.4))

plot(data$date,data$POLYGON_Close, type = "l")
plot(data$date,data$POLYGON_r, type = "l",ylim = c(-0.4, 0.4))

plot(data$date,data$MONERO_Close, type = "l")
plot(data$date,data$MONERO_r, type = "l",ylim = c(-0.4, 0.4))

plot(data$date,data$KUSAMA_Close, type = "l")
plot(data$date,data$KUSAMA_r, type = "l",ylim = c(-0.4, 0.4))
par(mfrow=c(1,1)) 


#statystyki opisowe
basicStats(data$total_r)

#histogram
hist(data$total_r, prob = T, breaks = 100)
curve(dnorm(x,
            mean = mean(data$total_r, na.rm = T),
            sd = sd(data$MANA_r, na.rm = T)),
      col = "darkblue", lwd = 2, add = TRUE)

#badanie normalności
jarque.bera.test(na.omit(data$total_r))
#h0 o normalności silnie odrzucona


#wykres ACF zwrotów
acf(data$total_r, lag.max = 30, na.action = na.pass,
    ylim = c(-0.3, 0.3), 
    col = "darkblue", lwd = 7,
    main = "Wykres ACF zwrotów portfela")
#Istotne okazały się być opóźnienia: 1, 5, 23, 24, 28

# wykres ACF kwadratów zwrotów
acf(data$total_r^2, lag.max = 30, na.action = na.pass,
    ylim = c(-0.3, 0.3), 
    col = "darkblue", lwd = 7,
    main = "Wykres ACF kwadrat?w zwrot?w FTSE100")
#Wszystkie opóźnienia istotne
#Na podstawie wykresów można podejrzewać, że występują efekty ARCH

#statystyka Durbina-Watsona
durbinWatsonTest(lm(formula = data$total_r ~ 1),
                 max.lag = 5)
#sprawdzenie autokorelacji dla pierwszych 5 opóźnień

#występowanie efektów ARCH wśród zwrotów
ArchTest(data$total_r, lags = 5)
# H0 o braku efektów silnie odrzucana

#statystyka DW dla kwadratów zwrotów
durbinWatsonTest(lm(formula = data$total_r ^ 2 ~ 1),
                 max.lag = 5) # sprawdźmy autokorelację dla 5 pierwszych opóźnień kwadratów

# dla 95-procentowego przedziału ufności opóźnienia 1,2 i 5 są istotne

#test występowania efektów ARCH
ArchTest(data$total_r, lag = 5, demean=TRUE)$p.value
#na podstawie wyników testu można odrzucić h0 o braku efektów ARCH

