getCorData <- function(m, sub = NULL){
  if(!is.null(sub)){
    return(dplyr::select(as.data.frame(t(m)), any_of(sub)))
  } else{
    return(as.data.frame(t(m)))
  }
}
buildESDF <- function(rds.file){
  enrichment.scores <- readRDS(rds.file)
  subset <- unlist(strsplit(readLines("/home/jan1/bioinf-tools/pipelines/ssGSVA/sub.txt"), ",\\s*"))
  
  # df <- getCorData(enrichment.scores[[fname]], sub = subset)
  
  es.df <- getCorData(enrichment.scores, sub = subset)
  
  return(es.df)
}
getSubset <- function(p, n){
  return(as.vector(n[grep(pattern = p, n)]))
}
main <- function(fpath, subset){
  rds.files <- list.files(path = fpath,pattern = "\\.RDS$", full.names = TRUE)
  
  if(length(rds.files) > 0){
    enrichment.scores <- list()
    es.df             <- list()
    
    for(i in seq_along(rds.files)){
      fname          <- mesocore::getFileName(rds.files[i])
      es.df[[fname]] <- buildESDF(rds.files[i])
      
    }
    res <- Reduce(cbind, es.df)
    
  } else{
    stop("Missing and enrichment score table. Please provide at least one enrichment scroe table to subset!")
  }
  
  outfolder <- here::here(fpath)
  write.csv(res, here::here(outfolder, "enrichment_experiment_df.csv"), row.names = F)
  
}

#---- Check Packages ----
pkgs <- c("here", "docopt")
mesocore::handleRequirements(pkgs)

#---- Header ----
"Mesothelioma AI Pipeline - Building the Correlations DF

Usage: survival_analysis.R [options]

Options:
  -h --help                     Show this screen.
  -e --enrichment_scores=<PATH> Path to the enrichment scores data frame folder.
  -s --subset=<PATH>            A .txt file that specifies what subset of ES to select.
"-> doc

#---- Flag Arguments ----
arguments <- docopt(doc, quoted_args = TRUE, help = TRUE)
print(arguments)

fpath <- arguments$enrichment_scores

if(!is.null(arguments$subset)){
  subset            <- unlist(strsplit(readLines(arguments$subset), ",\\s*"))
} else{
  subset <- NULL
}

#---- Using Msig.db.all database ----
# fpath  <- "/home/jan1/Documents/Cancer_Studies_PhD/Studies/Study_CONFIRM/results/ssGSEA/CONFIR_GROUP_NR"
# subset <- "/home/jan1/bioinf-tools/pipelines/ssGSVA/sub.txt"

main(fpath, subset)
