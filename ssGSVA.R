handleRequirements  <- function(pkgs){ # This function checks package requirements and install them if they are missing
  suppressMessages(if (!require("BiocManager", character.only = TRUE)) { # First check via R BioConductior
    install.packages("BiocManager")
    BiocManager::install()
  } else {
    ipkgs <- sapply(pkgs, function(...) require(..., character.only = TRUE))
    if (any(!ipkgs)) {
      BiocManager::install(pkgs[!ipkgs])
      install.packages(pkgs[!ipkgs])
    } else {
      message("\n\nCool! your machine has everything is needed.\n\n")
    }
  })
  
  print("Loading required packages...")
  library(pacman)
  pacman::p_load(pkgs, install = TRUE, character.only = TRUE) # Check via RCran and other repositories
  return(pacman::p_loaded()) # Return loaded packages
}
prepareGSVAMatrix <- function(m){

  fe <- tools::file_ext(m)

  if(fe == "csv"){
    x <- read.csv(m)
    
  } else if (fe %in% c("xls", "xlsx")) {
    # Install and load the readxl package if not already installed
    if (!require(readxl)) {
      install.packages("readxl")
      library(readxl)
    }
    x <- read_excel(fe)
    
  } else {
    stop("Unsupported file type.")
  }
  
  gene.symbols <- x[,1]
  x   <- as.matrix(x[,-1])
  sid <- colnames(x)
  
  rownames(x) <- gene.symbols
  
  return(
    list(x = x,
         sid = sid,
         genes  = gene.symbols)
  )
}

#---- Package Installation ----
pkgs <- c("GSVA", "here", "docopt", "GSEABase")
handleRequirements(pkgs)

#---- Header ----
"Mesothelioma AI Pipeline - ssGSVA Experiments Wrapper

Usage: ssGSEA_quickstart.R [options]

Options:
  -h --help                  Show this screen.
  --matrix=<string>          Specify the path to your gct file.
  --gene_signature=<string>  Path to the gene signature(s).
  --verbose=<value>          If set to T prints all messages [default: F].
  --version                  
"-> doc

#---- Arguments ----
arguments <- docopt(doc, quoted_args = TRUE, help = TRUE)
print(arguments)

#------------------ Load dataset and parameters into R environment -------------

# Command Line Arguments
gsva.obj <- prepareGSVAMatrix(here(arguments$matrix))
gs       <- getGmt(arguments$gene_signature)
verbose  <- as.integer(arguments$verbose)

# Local Arguments
# m   <- "E:/R-development/ssGSEA/medusa150/MEDUSA150_FPKM.csv"
# gmt <- "E:/R-development/ssGSEA/quickstart/db/h.all.v7.0.symbols.gmt"
# gsva.obj <- prepareGSVAMatrix(m)
# msigdb <- "E:/GitHubRepos/ssGSEA2/db/msigdb"

start_time <- Sys.time()
ES <- gsva(gsva.obj$x, gs, method = "ssgsea", verbose = T)


exp.name <- tools::file_path_sans_ext(basename((arguments$gene_signature)))
outdir <- here(dirname(here(arguments$matrix)), 'results')
dir.create(outdir, showWarnings = T, recursive = T)
saveRDS(ES, here(outdir, paste0(exp.name, ".RDS")))

# Timer
end_time <- Sys.time()
execution_time <- end_time - start_time
print(execution_time)
writeLines(paste0('Analysis ran for ', round(execution_time[[1]], 2), ' min'),
           con = here(outdir,'time.log'))