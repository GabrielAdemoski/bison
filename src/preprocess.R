library(readxl, warn.conflicts=FALSE)
library(dplyr, warn.conflicts=FALSE)
library(tidyr, warn.conflicts=FALSE)
library(comprehenr, warn.conflicts=FALSE)
library(docstring, warn.conflicts=FALSE)


split_catagory = function(x, col){
  #' Split by catagory
  #'
  #' @description Split a dataframe containing labled taxa counts into individual dataframes by a given catagory.
  #' @param x. A dataframe containing all taxa count and catagory data.
  #' @param col. A string of a column name for the column you want to split by.
  #' @return a list of dataframes corresponding to samples split by subcatagory.
  #' The list will be in the order i.e. winter, spring, summer, fall or 2018, 2019, 2020, 2021
  #' Each sample (row) contains taxa names, their corresponding OTU amounts, index metadata, and catagory metadata.

  x = split(x, x[,col])
  for (i in 1:(length(x))) {
    x[[i]] = x[[i]][, !names(x[[i]]) %in% col]
  }
  return(to_list(for (i in 1:(length(x))) x[[i]]))
}

load_data = function(path, cols, zero_frac=0.95) {
  #' Load data
  #'
  #' @description Load and parse data from an excel file.
  #' @details The excel file found at 'path' must contain an index column which has unique sample identifiers,
  #' individual columns for each taxa where the taxa names are in QIIME format
  #'
  #' @param path. The path to the excel file to be loaded. The path variable must be a string.
  #' @param cols. A vector of column names which you are looking for. Each column name must be a string.
  #' @param zero_frac. A fraction, [0,1] which represents the fraction of a taxa which may be zeros. For example, if a
  #' taxa has 100 samples and a zero_frac = 0.95, then the taxa that have 95 zeros or more will be filtered out.
  #' @return a list of dataframes corresponding to samples split by catagories.
  #' Each sample contains taxa names, their corresponding abundance amounts, index metadata, and catagory metadata.
  #' Sample taxa that do not start with k_ (QIIME code for kingdom) are filtered out.

  if (! file.exists(path)){stop(paste(path, 'does not exist.')}
  data <- read_excel(path)
  taxa = select(data, contains("k_"))
  filt = taxa[colSums(taxa != 0.0) >= dim(taxa)[1]*(1-zero_frac)]
  ret = list()
  for(col in cols){
    ret[[length(ret)+1]] = split_catagory(select(data, colnames(filt), matches(col), index), col)
  }
  return(ret)
}



