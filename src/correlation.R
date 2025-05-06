library(dplyr, warn.conflicts=FALSE)
library(tidyr, warn.conflicts=FALSE)
library(jsonlite, warn.conflicts=FALSE)
library(stringr, warn.conflicts=FALSE)
library(docstring, warn.conflicts=FALSE)

source("src/network.R")
source("src/preprocess.R")
source("src/phylo.R")


individual.corr <- function(categories, levels, xlsx.prefix){
  #' Network analysis of categories and levels
  #'
  #' @description Perform network analysis for categorical time series taxanomic level data.
  #' @details This function itterates over all possible combnations of category and taxanomic level. The requisite data
  #' is in data/. This data was generated using QIMME and represents all taxa which showed significant OTU abundance.
  #' A network is created for every combination of category and taxanomic level resulting in
  #' (number of levels * number of categories) network runs and resulting correlation matricies and network analyses.
  #' All network analyses are concatenated into one variable called 'info'. This variable is saved as a .csv file and
  #' saved in the same folder (corr.dir) as the correlation matracies at the end of the run.
  #'
  #' @param There are no function arguments.
  #' @return Nothing is returned.

  main.dir = paste0(getwd(), "/")
  data.dir = "data/"
  corr.dir = "data/corr/"

  # There is currently only one directory to be created, however, if more were needed then add them to the vector in
  # for loops argument
  for (dir in c(corr.dir)){
    dir.create(file.path(main.dir, dir), showWarnings = FALSE, recursive = TRUE)
  }

    info = data.frame(Level = character(), Category = character(), Variable = character(), Stability = double(),
                      OptimalLambda = double(), NumVert = double(), NumEdges = double(), AvgDegree = double(),
                      AvgPathLen = double(), ClusteringCoeff = double(), PositiveEdges = double(),
                      NegativeEdges = double(), NumLoneNodes = integer(), NumCuster = integer(),
                      BiggestCluster = integer(), SmallestCluster = integer())

  # zeallot must be called after other libraries and before load_data or you'll get an error
  library(zeallot, warn.conflicts=FALSE)

  for (ca in categories){
        for (lvl in levels){
          c(season, year) %<-% load_data(paste0(data.dir, xlsx.prefix, lvl, ".xlsx"), cols=catagories)

          if (ca == "Seasons"){
            splt = season
            hsh = c('Winter', 'Spring', "Summer", "Fall")
         ## copy this else if for more variables like STEC presence
          }else if(ca == "Years"){
            splt = year
            hsh = c(2018, 2019, 2020, 2021)
          }
          for (group in hsh){
            group.data = splt[[match(group, hsh)]]
            cat(paste0("\n\n", "l", lvl, "_", group, "\n\n"))
            csv.name = paste0(main.dir, corr.dir, "l", lvl, "_", group, ".csv")
            if (!file.exists(csv.name))
            {
              # Make phyloseq object from the data
              c(otu.mat, otu.hash) %<-% format_for_phylo(group.data)
              if (lvl == "2"){
                ps = lvl2(otu.mat, otu.hash)
              }else if (lvl == "6"){
                ps = lvl6(otu.mat, otu.hash)
              }
              else{stop("Only level 2 (Phylum) and level 6 (Genus) data are accepted")}

              c(corr, stab, lambda) %<-% network(ps, verbose = TRUE)
              # Save correlation matrix for later use
              write.table(corr, file=csv.name)
              stats = network_stats(corr)
              stats = modifyList(list(Level = lvl, Category = ca, Variable = group,
                                      Stability = stab, OptimalLambda = lambda), stats)
              info = rbind(info, stats)
            }
          }
        }
      }

  stat.csv = paste0(corr.dir, "individual_network_stats.csv")
  print(stat.csv)
  if (! file.exists(stat.csv)){write.csv(info, stat.csv, row.names = FALSE)}
  print("Done running individual.corr---------------------------------------------------")
}


delta.corr <- function(catagories, levels, xlsx.prefix){
  #' Change in categories over time
  #'
  #' @description Perform network analysis for categorical time series taxanomic level data.
  #' @details This function itterates over all possible combnations of category and taxanomic level. The requisite data
  #' is in data/{category}_json/. This data was generated using QIMME and represents all taxa which showed significant
  #' log fold change between time 1 to time 2 for a given category. This information is then used to filter out all taxa
  #' which were not considered to be significant. A network is created for both time 1 taxa OTU amounts and time 2 taxa
  #' OTU amounts per category resulting in (2 * number of levels * number of categories) network runs and resulting
  #' correlation matricies and network analyses. All network analyses are concatenated into one variable called 'info'.
  #' This variable is saved as a .csv file and saved in the same folder (corr.dir) as the correlation matracies at the
  #' end of the run.
  #'
  #' @param There are no function arguments.
  #' @return Nothing is returned.

  main.dir = paste0(getwd(), "/")
  data.dir = "data/"
  corr.dir = "data/delta_corr/"


  # There is currently only one directory to be created, however, if more were needed then add them to the vector in
  # for loops argument
  for (dir in c(corr.dir)){
    dir.create(file.path(main.dir, dir), showWarnings = FALSE, recursive = TRUE)
  }
  info = data.frame(Level = character(), Category = character(), Variable = character(), Stability = double(),
                      OptimalLambda = double(), NumVert = double(), NumEdges = double(), AvgDegree = double(),
                      AvgPathLen = double(), ClusteringCoeff = double(), PositiveEdges = double(),
                      NegativeEdges = double(), NumLoneNodes = integer(), NumCuster = integer(),
                      BiggestCluster = integer(), SmallestCluster = integer())

  # zeallot must be called after other libraries and before load_data or you'll get an error
  library(zeallot, warn.conflicts=FALSE)

  for (ca in categories){
    ca.folder = paste0(data.dir, ca, "_json")
    if (!(dir.exists(ca.folder))){stop(ca.folder, " does not exist!")}
    # Load all JSON file names
    jsons = list.files(ca.folder)
    for (json in jsons){
      # Extract the taxa names that changed significantly
      lfc = fromJSON(txt=paste0(ca.folder,"/", json))$datasets[[1]]$id
      lfc = lfc[grepl("k_", lfc)]
      # remove .json from the end of the file name
      json = str_sub(json, start=0, end=-6)
      # Split the file names based on _ (\\_) or (|) - (\\-), the () are grep patterns
      ss = strsplit(json, split="\\_|\\-")

      for (s in ss){
        lvl = s[3]
        from = s[4]
        to = s[5]

        # Make sure lvl is a valid input
        if(!(lvl %in% levels)){
          #ERROR
          stop("The level must be one of the following: ", list(levels))
        }else{
          c(season, year) %<-% load_data(paste0(data.dir, xlsx.prefix, str_sub(lvl, start=2), ".xlsx"),
                                         cols=catagories)
        }

        # Make sure from and to are valid inputs
        if (ca == "Seasons"){
          splt = season
          hsh = c('Winter', 'Spring', "Summer", "Fall")
          if (!(from %in% hsh)){
            #ERROR
            stop("The 'from' season must be in the following: ", list(hsh))
          } else if(!(to %in% hsh)){
            #ERROR
            stop("The 'to' season must be in the following: ", list(hsh))
          }
        }else if(ca == "Years"){
            splt = year
            hsh = c(2018, 2019, 2020, 2021)
          if (!(from %in% hsh)){
            #ERROR
            stop("The 'from' year must be in the following: ", list(hsh))
          } else if(!(to %in% hsh)){
            #ERROR
            stop("The 'to' year must be in the following: ", list(hsh))
          }
        }

        for (group in c(from, to)){
          group.data = splt[[match(group, hsh)]]
          group.data = select(group.data, all_of(lfc), index)
          var.name = paste0(lvl, "_", from,"-",to, "_", group)
          cat("\n\n", var.name, "\n\n")
          csv.name = paste0(main.dir, corr.dir, var.name, ".csv")
          if (!file.exists(csv.name))
          {
            # Make phyloseq object from the data
            c(otu.mat, otu.hash) %<-% format_for_phylo(group.data)
            if (lvl == "l2"){
              ps = lvl2(otu.mat, otu.hash)
            }else if (lvl == "l6"){
              ps = lvl6(otu.mat, otu.hash)
            }
            else{stop("Only level 2 (Phylum) and level 6 (Genus) data are accepted")}
            # Check to make sure there is more than one taxa
            if (dim(tax_table(ps))[1] > 1){
              # Save correlation matrix, network stability, and optimal StARS lambda for later use
              c(corr, stab, lambda) %<-% network(ps, verbose = FALSE)
            }
            # If there isn't, return a null matrix (all zeros)
            else{
              corr = matrix(0, 1, 1)
              tax = data.frame(tax_table(ps))
              last2 = tax[tail(colnames(tax), n=2)]
              # combine the first letter of the second to last column with the final column and remove from the final column
              tax.names = paste0(str_sub(last2[[1]], start=4, end=4), ". ", str_sub(last2[[2]], start=4))
              # last.col = tax_table(ps)[,ncol(tax_table(ps))][1]
              colnames(corr) = tax.names
              rownames(corr) = tax.names
              stab = -1
              lambda = -1
            }
            stats = network_stats(corr)
            stats = modifyList(list(Level = lvl, Category = ca, Variable = var.name,
                                    Stability = stab, OptimalLambda = lambda), stats)
            info = rbind(info, stats)
            write.table(corr, file=csv.name)
          }
        }
      }
    }
  }
  stat.csv = paste0(corr.dir, "delta_network_stats.csv")
  print(stat.csv)
  if (! file.exists(stat.csv)){write.csv(info, stat.csv, row.names = FALSE)}
  print("Done running delta.corr----------------------------------------------------")
}