# Workflow plumbing only. Scientific calculations are not invoked by these helpers.
pc057_root <- function() {
  p <- normalizePath(getwd())
  for (i in seq_len(8L)) {
    desc <- file.path(p,"DESCRIPTION")
    if(file.exists(desc) && any(grepl("^Package: pc057$",readLines(desc)))) return(p)
    parent <- dirname(p); if(identical(parent,p)) break; p <- parent
  }
  stop("Run from the PC057 repository or its analyses directory.")
}
pc057_python <- function() {
  p <- Sys.getenv("PC057_PYTHON",unset="python3")
  if(!nzchar(Sys.which(p)) && !file.exists(p)) stop("Python 3.11 or newer is required.")
  p
}
pc057_command <- function(command,args) {
  output <- system2(command,vapply(args,shQuote,character(1)),stdout=TRUE,stderr=TRUE)
  status <- attr(output,"status"); if(length(output))cat(paste(output,collapse="\n"),"\n")
  if(!is.null(status) && status!=0L)stop("Command failed: ",basename(command))
  invisible(output)
}
pc057_step_begin <- function(name,path_target,path_resource) {
  root <- pc057_root(); allowed <- jsonlite::fromJSON(file.path(root,"metadata/steps.json"))$step
  if(!name %in% allowed)stop("Unknown step: ",name)
  target <- normalizePath(path_target(),mustWork=FALSE)
  expected <- normalizePath(file.path(root,"analyses/data",name),mustWork=FALSE)
  if(target!=expected)stop("qproj output path does not match the selected step.")
  bundle <- normalizePath(file.path(path_resource,"pc057"),mustWork=FALSE)
  if(!file.exists(file.path(bundle,"inputs/models/revised_results.rds")))stop("Prepare inputs first: python3 dev/prepare_inputs.py")
  dir.create(target,recursive=TRUE,showWarnings=FALSE)
  stage <- tempfile(pattern=".staging-",tmpdir=target);dir.create(stage)
  list(name=name,root=root,target=target,stage=stage,bundle=bundle)
}
pc057_verify_inputs <- function(step) {
  lines <- pc057_command(pc057_python(),c(file.path(step$root,"dev/prepare_inputs.py"),"--verify",step$bundle))
  writeLines(lines,file.path(step$stage,"input-verification.json"))
}
pc057_export_tables <- function(step) {
  pc057_command(pc057_python(),c(file.path(step$root,"dev/export_supplementary_tables.py"),"--bundle",step$bundle,"--output",file.path(step$stage,"tables")))
}
pc057_verify_saved_results <- function(step) {
  lines <- pc057_command(file.path(R.home("bin"),"Rscript"),c("--vanilla",file.path(step$root,"dev/verify_saved_results.R"),"--bundle",step$bundle))
  writeLines(lines,file.path(step$stage,"saved-results-verification.txt"))
}
pc057_export_spatial_margins <- function(step) {
  pc057_command(pc057_python(),c(file.path(step$root,"dev/export_spatial_margins.py"),"--bundle",step$bundle,"--output",file.path(step$stage,"spatial-margins.csv")))
}
pc057_draw_figure <- function(step,figure) {
  fn <- get(paste0("pc057_plot_",figure),mode="function")
  fn(step$bundle,step$stage)
  stopifnot(file.exists(file.path(step$stage,paste0(figure,".png"))))
}
pc057_step_finish <- function(step) {
  files <- list.files(step$stage,recursive=TRUE,full.names=TRUE)
  if(!length(files))stop("The step produced no outputs.")
  relative <- substring(files,nchar(step$stage)+2L)
  for(i in seq_along(files)) {
    dest <- file.path(step$target,relative[i]);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
    if(!file.copy(files[i],dest,overwrite=TRUE))stop("Unable to publish a staged output.")
  }
  receipt <- list(step=step$name,status="pass",operation="frozen-result verification, export or redraw only",outputs=relative)
  jsonlite::write_json(receipt,file.path(step$target,"receipt.json"),pretty=TRUE,auto_unbox=TRUE)
  unlink(step$stage,recursive=TRUE)
  data.frame(step=step$name,status="pass",output=relative)
}
