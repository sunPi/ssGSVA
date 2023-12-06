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
library("docopt")

#---- Header ----
"Mesothelioma AI Pipeline - ssGSVA Experiments Wrapper

Usage: ssGSEA_quickstart.R [options]

Options:
  -h --help                  Show this screen.
  --matrix=<string>          Specify the path to your gct file.
  --gene_signature=<string>  Path to the gene signature(s).
  --verbose=<value>          If set to T prints all messages [default: F].
  --version                  Show version.
"-> doc

#---- Arguments ----
arguments <- docopt(doc, quoted_args = TRUE, help = TRUE, version = 'ssGSVA Wrapper 1.0\n')
print(arguments)

#---- Package Installation ----
pkgs <- c("GSVA", "here", "GSEABase")
suppressMessages(mesocore::handleRequirements(pkgs))

#------------------ Load dataset and parameters into R environment -------------
# Command Line Arguments
gsva.obj <- prepareGSVAMatrix(here(arguments$matrix))
gs       <- getGmt(arguments$gene_signature)
verbose  <- as.integer(arguments$verbose)

# Local Arguments
start_time <- Sys.time()
ES <- gsva(gsva.obj$x, gs, method = "ssgsea", verbose = T)

exp.name <- tools::file_path_sans_ext(basename((arguments$gene_signature)))
outdir <- here(dirname(arguments$matrix), mesocore::getFileName(arguments$matrix))
dir.create(outdir, showWarnings = F, recursive = T)
saveRDS(ES, here(outdir, paste0(exp.name, ".RDS")))

# Timer
end_time <- Sys.time()
execution_time <- end_time - start_time
print(execution_time)
writeLines(paste0('Analysis ran for ', round(execution_time[[1]], 2), ' min'),
           con = here(outdir,'time.log'))