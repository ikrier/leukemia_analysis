#This is the function from Gviz that I'm just using explicitly in the script
import.bam <- function (file, selection) 
{
  if (!file.exists(paste(file, "bai", sep = ".")) && !file.exists(paste(paste(head(strsplit("xxx.bam", 
                                                                                            ".", fixed = TRUE)[[1]], -1), collapse = "."), "bai", 
                                                                        sep = "."))) 
    stop("Unable to find index for BAM file '", file, "'. You can build an index using the following command:\n\t", 
         "library(Rsamtools)\n\tindexBam(\"", file, "\")")
  sinfo <- scanBamHeader(file)[[1]]
  if (parent.env(environment())[["._trackType"]] == "DataTrack") {
    res <- if (!as.character(seqnames(selection)[1]) %in% 
               names(sinfo$targets)) {
      mcols(selection) <- DataFrame(score = 0)
      selection
    }
    else {
      param <- ScanBamParam(what = c("pos", "qwidth"), 
                            which = selection, flag = scanBamFlag(isUnmappedQuery = FALSE))
      x <- scanBam(file, param = param)[[1]]
      cov <- coverage(IRanges(x[["pos"]], width = x[["qwidth"]]))
      if (length(cov) == 0) {
        mcols(selection) <- DataFrame(score = 0)
        selection
      }
      else {
        GRanges(seqnames = seqnames(selection), ranges = IRanges(start = start(cov), 
                                                                 end = end(cov)), strand = "*", score = runValue(cov))
      }
    }
  }
  else {
    res <- if (!as.character(seqnames(selection)[1]) %in% 
               names(sinfo$targets)) {
      mcols(selection) <- DataFrame(id = "NA", group = "NA")
      selection[0]
    }
    else {
      param <- ScanBamParam(what = c("pos", "qwidth", "strand", 
                                     "qname"), which = selection, flag = scanBamFlag(isUnmappedQuery = FALSE))
      x <- scanBam(file, param = param)[[1]]
      GRanges(seqnames = seqnames(selection), ranges = IRanges(start = x[["pos"]], 
                                                               width = x[["qwidth"]]), strand = x[["strand"]], 
              id = make.unique(x[["qname"]]), group = x[["qname"]])
    }
  }
  return(res)
}