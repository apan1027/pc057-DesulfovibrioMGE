# Execute the same QMD R expressions with knitr, then write GitHub Markdown.
# This driver requires no network access and never invokes optional refit source.
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(getwd())
stopifnot(file.exists(file.path(root,"DESCRIPTION")))
if(dir.exists(file.path(root,"renv/library"))) {
  if(!requireNamespace("renv",quietly=TRUE))stop("Install renv before loading the restored project library.")
  renv::load(project=root,quiet=TRUE)
}
# Child Rscript checks use the same library selection as the notebook process.
Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
steps <- jsonlite::fromJSON(file.path(root,"metadata/steps.json"))$step
if(length(args)) {
  if(!all(args %in% steps))stop("Unknown requested notebook.")
  steps <- args
}
pandoc <- Sys.getenv("PC057_PANDOC",unset=Sys.which("pandoc"))
if(!nzchar(pandoc) || !file.exists(pandoc))stop("Pandoc is required; set PC057_PANDOC if it is not on PATH.")
stage <- file.path(root,"checks/render-sources");dir.create(stage,recursive=TRUE,showWarnings=FALSE)
convert <- function(lines) {
  out <- character();i <- 1L
  while(i<=length(lines)) {
    if(grepl("^```\\{r\\}$",lines[i])) {
      i<-i+1L;opts<-character()
      while(i<=length(lines) && grepl("^#\\| ",lines[i])) {
        pair<-strsplit(sub("^#\\| ","",lines[i]),": ",fixed=TRUE)[[1]]
        key<-pair[1];value<-paste(pair[-1],collapse=": ")
        if(key=="eval")value<-"FALSE"
        else if(value %in% c("true","false"))value<-toupper(value)
        else if(key=="label")value<-paste0('"',value,'"')
        opts<-c(opts,paste0(key,"=",value));i<-i+1L
      }
      out<-c(out,paste0("```{r ",paste(opts,collapse=", "),"}"))
    }else{out<-c(out,lines[i]);i<-i+1L}
  }
  out
}
for(name in steps) {
  message("Executing frozen workflow notebook: ",name)
  input<-file.path(root,"analyses",paste0(name,".qmd"));rmd<-file.path(stage,paste0(name,".Rmd"));md<-file.path(stage,paste0(name,".md"))
  writeLines(convert(readLines(input,warn=FALSE)),rmd)
  e<-new.env(parent=globalenv());e$params<-list(name=name)
  knitr::opts_knit$set(root.dir=file.path(root,"analyses"))
  knitr::knit(rmd,output=md,envir=e,quiet=TRUE)
  target<-file.path(root,"analyses",paste0(name,".md"))
  status<-system2(pandoc,c(shQuote(md),"--from=markdown","--to=gfm","--standalone","--output",shQuote(target)))
  if(status!=0)stop("Pandoc conversion failed for ",name)
  stopifnot(file.exists(target),file.exists(file.path(root,"analyses/data",name,"receipt.json")))
}
message("PASS: ",length(steps)," notebooks executed; no model fitting, bootstrap or annotation invoked.")
