#---- Functions ----
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
    x <- read_excel(m)
    
  } else {
    stop("Unsupported file type.")
  }
  
  gene.symbols <- x[,1]
  x   <- as.matrix((x[,-1]))
  sid <- colnames(x)
  
  rownames(x) <- gene.symbols
  
  return(
    list(x = x,
         sid = sid,
         genes  = gene.symbols)
  )
}
transform.ES <- function(ES){
  res <- as.data.frame(ES)
  res <- tibble::rownames_to_column(res, "Pathway_Name")
  res <- as.data.frame(t(res))
  colnames(res) <- res[1,]
  res <- res[-1,]
  res <- as.data.frame(tibble::rownames_to_column(res, "SampleID"))
  return(res)
}

#---- Package Installation ----
pkgs <- c("GSVA", "here", "GSEABase", "docopt")
suppressMessages(mesocore::handleRequirements(pkgs))

#---- Header ----
"Mesothelioma AI Pipeline - ssGSVA Experiments Wrapper

Usage: ssGSEA_quickstart.R [options]

Options:
  -h --help                  Show this screen.
  --matrix=<string>          Specify the path to your gct file.
  --gene_signature=<string>  Path to the gene signature(s).
  --outfolder=<string>       Path to the results directory.
  --verbose=<value>          If set to T prints all messages [default: F].
  --version                  Show version.
"-> doc

#---- Arguments ----
arguments <- docopt(doc, quoted_args = TRUE, help = TRUE, version = 'ssGSVA Wrapper 1.0\n')
print(arguments)

#------------------ Load dataset and parameters into R environment -------------
# Command Line Arguments
gsva.obj  <- prepareGSVAMatrix(normalizePath(arguments$matrix, winslash = "/", mustWork = FALSE))
gs        <- getGmt(arguments$gene_signature)
# verbose   <- as.integer(arguments$verbose)
# gsva.obj <- prepareGSVAMatrix("/home/jr453/Documents/Cancer_Studies_PhD/Studies/Study_PARP/datasets/MEDUSA/RNAseq/M126_TPMpc.csv")
# gs        <- getGmt("~/bioinf-tools/db/molsigs/msigdb/h.all.v2024.1.Hs.symbols.gmt")

# Main
start_time <- Sys.time()

gsea.set <- ssgseaParam(gsva.obj$x, gs)
ES <- gsva(gsea.set, verbose = T)

exp.name <- tools::file_path_sans_ext(basename((arguments$gene_signature)))
if(arguments$outfolder == "NULL"){
  outdir <- here(dirname(arguments$matrix), mesocore::getFileName(arguments$matrix))
} else{
  outdir <- normalizePath(arguments$outfolder, winslash = "/", mustWork = FALSE)
}

dir.create(outdir, showWarnings = F, recursive = T)
saveRDS(ES, here(outdir, paste0(exp.name, ".RDS")))

ES.t <- transform.ES(ES)

write.csv(ES.t, here(outdir, paste0(exp.name, ".csv")), row.names = F)

# Timer
end_time <- Sys.time()
execution_time <- end_time - start_time
print(execution_time)
writeLines(paste0('Analysis ran for ', round(execution_time[[1]], 2), ' min'),
           con = here(outdir,'time.log'))