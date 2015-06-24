# =======
#   License
# =======
#   This code is released under the GNU General Public License 3.0. A copy
# of this license is in the LICENSE.txt file.
# copyright Irina Krier 2015
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

strandedBamImport <- function (file, selection)
{
  if (!file.exists(paste(file, "bai", sep = ".")))
    stop("Unable to find index for BAM file '", file, "'. You can
         build an index using the following command:\n\t",
         "library(Rsamtools)\n\tindexBam(\"", file, "\")")
  sinfo <- scanBamHeader(file)[[1]]
  res <- if (!as.character(seqnames(selection)[1]) %in%
             names(sinfo$targets))
  {
    mcols(selection) <- DataFrame(score = 0)
    selection
  }else
  {
    param <- ScanBamParam(what = c("pos", "qwidth", "strand"),
                          which = selection, flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate=FALSE))
    x <- scanBam(file, param = param)[[1]]
    if(length(x$pos)==0)
    {
      mcols(selection) <- DataFrame(plus=0, minus=0)
      selection
    }else
    {
      gr <- GRanges(strand=x[["strand"]],ranges=IRanges(x[["pos"]],width = x[["qwidth"]]),seqnames=seqnames(selection)[1])
      grs <- split(gr, strand(gr))
      cov <- lapply(grs[c("+", "-")],function(y){coverage(ranges(y),width=end(selection))})
      pos <- sort(unique(unlist(lapply(cov, function(y) c(start(y),end(y))))))
      
      GRanges(seqnames = seqnames(selection)[1],
              ranges=IRanges(start=head(pos, -1), end=tail(pos, -1)),
              plus=as.numeric(cov[["+"]][head(pos, -1)]),
              minus=-as.numeric(cov[["-"]][head(pos, -1)]),
              both=as.numeric(cov[["+"]][head(pos, -1)])+as.numeric(cov[["-"]][head(pos, -1)]))
    }
  }
  return(res)
}