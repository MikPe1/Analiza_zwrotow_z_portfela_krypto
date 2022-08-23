##########################################################################################

import_data<-function(dataset)
	{
	 data_<- read.csv(paste("data/",
	                        dataset,
	                        ".txt",
	                        sep = ""),
	                  sep = ",",
	                  dec = ".",
	                  header = T)
	 data_$Date <- as.Date(as.character(data_$Date),
	                       "%Y%m%d")
	 data_ <- data_[, c(2,6)]
	 names(data_)[2] <- dataset
	 return(data_)
	}

##########################################################################################

compare.ICs <- function(models_list)
{
	n_ <- length(models_list)

	for(i in 1:n_)
		{
		ICs_ <- data.frame(t(get(models_list[i])@fit$ics))
		ICs_$model <- models_list[i]
		if(i == 1) ICs <- ICs_ else ICs <- rbind(ICs, ICs_)
		}
	return(ICs)
}