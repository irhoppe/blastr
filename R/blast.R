
#' Nucleotide BLAST
#' 
#' Find similar nucleotide sequences using NCBI's BLAST+ API.
#' 
#' Function for running nucleotide BLAST queries, including traditional BLASTN and MegaBLAST tasks as well searches
#' optimized for short and discontiguous sequences. Accepts control parameters (requesting defaults when appropriate), 
#' generates a BLAST+ command, and passes command to system2() to interface with BLAST+ API. Accepts, parses, and returns
#' result of the query to the user.
#' 
#' @param db character; name of the BLAST database against which to compare query sequence(s)
#' @param query character vector giving the nuecleotide sequence(s) to submit to the BLAST+ API
#' @param options a `blastControl()`-derived list of arguments
blastn <- function( db, query, task=c("megablast", "blastn", "short", "dc"), options=blastControl() ){
  
  query <- parse_blast_query(query)
  task <- switch(match.arg(task), "megablast"="megablast", "blastn"="blastn", "short"="blastn-short", "dc"="dc-megablast")
  args <- c( dbflag="-db", db=db, taskflag="-task", task=task, options )
  blast_res <- system2( "blastn", args=args, input=query, stdout=TRUE )
  parsed_blast_res <- parse_blast_res( blast_res, args["outfmt"], args["fmtstr"] )
  
  return(parsed_blast_res)
  
}

# blastn.character <- function( db, query, task=c("megablast", "blastn", "short", "dc"), options=blastControl() ){
#   
#   
#   
# }
# 
# blastn.DNAStringSet <- function( db, query, task=c("megablast", "blastn", "short", "dc"), options=blastControl() ){
#   
#   query <- as.character(query)
#   query_names <- paste0(">",names(query))
#   query <- c(query_names, query)[order(c(seq_along(query_names),seq_along(query)))]
#   
#   
#   
#   
#   
#   task <- switch(match.arg(task), "megablast"="megablast", "blastn"="blastn", "short"="blastn-short", "dc"="dc-megablast")
#   args <- c( dbflag="-db", db=db, taskflag="-task", task=task, options )
#   blast_res <- system2( "blastn", args=args, input=query, stdout=TRUE )
#   parsed_blast_res <- parse_blast_res( blast_res, args["outfmt"], args["fmtstr"] )
#   
#   return(parsed_blast_res)
#   
# }

# blastp <- function( db, query, task=c("blastp", "short"))


blastnControl <- function(evalue=10, entrez_query=NULL, outfmt=6, fmtstr=NULL, remote=TRUE, arg_string=NULL ){
  
  if(!is.numeric(evalue))     stop("Control parameter 'evalue' must be a real number!")
  if(!outfmt%in%0:18)     stop("Not a valid output format.")
  if((!(outfmt%in%c(6,7,10))) & !is.null(fmtstr) ){
    message("Format string only meaningful for tabular formats (did you mean outfmt 6, 7, or 10?). Ignoring parameter 'fmtstr'.")
    fmtstr <- NULL
  }
  
  if( !is.null(fmtstr) ){                                                          # check validity of fmtstr fields
    fmtfieldlist <- strsplit(fmtstr, "\\s")[[1]]
    fmtfieldlist <- split(fmtfieldlist, fmtfieldlist %in% fmtfields())
    if( !is.null(fmtfieldlist$`FALSE`) ){
      message(sprintf("The following field codes are not valid format field strings: %s", paste(fmtfieldlist$`FALSE`, collapse=", ")))
      message("Ignoring those fields for now. Check `fmtspec` for a listing of valid field codes.")
    }
    fmtstr <- if( is.null(fmtfieldlist$`TRUE`) ) NULL else paste(fmtfieldlist$`TRUE`, collapse=" ")
  }
  
  dotargs <- as.character(as.list(match.call(expand.dots=FALSE)[-1])$`...`)
  
  controls <- c(
    evalueflag="-evalue", evalue=evalue, 
    outfmtflag="-outfmt", outfmt=outfmt, fmtstr=fmtstr, 
    arg_string, 
    remote=if(remote) "-remote" else NULL
  )
  
  return(controls)
  
}

blastControl <- function(...){
  
  app <- switch(as.character(sys.call(-1)[1]),        # use the name of calling function to identify which BLAST program is being invoked (influences which arguments are available)
                blastn="blastnControl", blastp="blastpControl", blastx="blastxControl", 
                tblastx="tblastxControl", tblastn="tblastnControl", psiblast="psiblastControl", 
                rpsblast="rpsblastControl", rpstblastn="rpsblastnControl")
  
  args <- as.list(match.call(expand.dots=TRUE)[-1])   # pull the arguments passed to the current function for passing on to the appropriate control function
  
  do.call(app, args)                                  # invoke the appropriate control function, passing on the present arguments
  
}

parse_blast_res <- function( res, outfmt, fmtstr ){
  
  if( outfmt %in% c("6","7","10") ){
    tsep <- switch( outfmt, "6"="\t", "7"="\t", "10"=",")
    if(is.na(fmtstr)) fmtstr <- "qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore"
    res_df <- read.table(textConnection(res), sep=tsep, quote="" )
    names(res_df) <- strsplit(fmtstr, split="\\s")[[1]]
  } else {
    message( "Parsing not yet implemented for non-tabular output formats. Returning un-parsed results. Sorry!" )
    res_df <- res
  }
  
  return(tibble::tibble(res_df))
  
}

parse_blast_input <- function(query){
  
  if(is.list(query)){
    query <- sapply(query,as.character)
  }
  if(class(query)=="DNAStringSet"){
    query <- as.character(query)
    if(!is.null(query_names <- names(query))){
      query_names <- paste0(">",names(query))
      query <- c(query_names, query)[order(c(seq_along(query_names), seq_along(query)))]
    }
  } else if(class(query)=="DNAString"){
    query <- as.character(query)
  } else if(is.character(query) & !is.null(query_names <- names(query))){
    query_names <- paste0(">",names(query))
    query <- c(query_names, query)[order(c(seq_along(query_names), seq_along(query)))]
  }
  
  return(query)
  
}


