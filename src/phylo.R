library(phyloseq, warn.conflicts=FALSE)
library(dplyr, warn.conflicts=FALSE)        # filter and reformat data frames
library(tibble, warn.conflicts=FALSE)       # Needed for converting column to row names
library(comprehenr, warn.conflicts=FALSE)
library(hash, warn.conflicts=FALSE)
library(docstring, warn.conflicts=FALSE)


format_for_phylo <- function(OTU.t){
  #' Format sample OTU data for phyloseq
  #'
  #' @description Produce the data needed to create a phyloseq object from OTU/taxa sample data provided. Two new data
  #' objects are created, otu.mat and otu.hash.
  #' @details
  #' otu.mat is a matrix where row names are index values, column names are OTU keys, and the values of the matrix are
  #' the OTU amounts.
  #' otu.hash is a hash of key-value pairs similar to a python dictionary. In this case the keys are arbitrary OTU
  #' integers and the values are the taxa names.
  #'
  #' @param OTU.t. A dataframe containing taxa names, their corresponding OTU amounts, index metadata, and catagory
  #' metadata.
  #' @return a list containing otu.mat and otu.hash

  OTU.t <- OTU.t %>% tibble::column_to_rownames("index")
  taxa = colnames(OTU.t)
  # create arbitrary OTU keys because phyloseq needs an OTU-to-taxa hash
  otus = to_list(for (i in 1:(length(taxa))) paste("OTU", i))
  otu.hash = hash(keys = otus, values = taxa)
  colnames(OTU.t) = otus
  # index is rownames and OTU keys is colnames
  OTU = t(OTU.t)
  otu.mat <- as.matrix(OTU)

  return (list(otu.mat, otu.hash))
}

lvl2 <- function(otu.mat, otu.hash) {
  #' Level 2 data to phyloseq
  #'
  #' @description Create a phyloseq object from level 2 (Phylum) OTU amount data and an OTU-to-taxa hash.
  #'
  #' @param otu.mat a matrix where row names are index values, column names are arbitrary OTU keys, and the values
  #' of the matrix are the OTU amounts.
  #' @param otu.hash OTU-to-taxa hash. The OTU keys are the same as in otu.mat.
  #' @return phy. A phyloseq object containing taxa names, OTU amounts, and arbitrary OTU keys.


  tax.table = data.frame(OTU = character(), Kingdom = character(), Phylum = character())
  library(zeallot, warn.conflicts=FALSE)
  for (key in keys(otu.hash)) {
    c(k, p) %<-% unlist(strsplit(otu.hash[[key]], ";"))
    tax.table = rbind(tax.table, list(OTU = key, Kingdom = k, Phylum = p))
  }
  tax.table <- tax.table %>% tibble::column_to_rownames("OTU")
  tax_mat = as.matrix(tax.table)
  phy = phyloseq(otu_table(otu.mat, taxa_are_rows = TRUE), tax_table(tax_mat))
  phy = subset_taxa(phy, Phylum != "p__")
  phy = subset_taxa(phy, Phylum != "__")

  return (phy)
}

lvl6 <- function(otu.mat, otu.hash) {
  #' Level 6 data to phyloseq
  #'
  #' @description Create a phyloseq object from level 6 (Genus) OTU amount data and an OTU-to-taxa hash.
  #'
  #' @param otu.mat a matrix where row names are index values, column names are arbitrary OTU keys, and the values
  #' of the matrix are the OTU amounts.
  #' @param otu.hash OTU-to-taxa hash. The OTU keys are the same as in otu.mat.
  #' @return phy. A phyloseq object containing taxa names, OTU amounts, and arbitrary OTU keys.

  tax.table = data.frame(OTU = character(), Kingdom = character(), Phylum = character(), Class = character(),
                       Order = character(), Family = character(), Genus = character())
  library(zeallot, warn.conflicts=FALSE)
  for (key in keys(otu.hash)) {
    c(k, p, c, o, f, g) %<-% unlist(strsplit(otu.hash[[key]], ";"))
    tax.table = rbind(tax.table, list(OTU = key, Kingdom = k, Phylum = p, Class = c, Order = o, Family = f, Genus = g))
  }
  tax.table <- tax.table %>% tibble::column_to_rownames("OTU")
  tax_mat = as.matrix(tax.table)
  phy = phyloseq(otu_table(otu.mat, taxa_are_rows = TRUE), tax_table(tax_mat))
  phy = subset_taxa(phy, Genus != "g__")
  phy = subset_taxa(phy, Genus != "__")

  return (phy)
}