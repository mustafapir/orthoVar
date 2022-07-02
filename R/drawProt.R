#' Draw lollipop plot for matching variants
#'
#' This function takes matching variant table, protein ids to be plotted and species names of these protein ids and exports a
#' lolliplot.
#'
#' @param table data frame consisting of matching variants, the output of `orthoFind` function.
#' @param protein_id IDs of proteins to be drawn.
#' @param organisms Name of organisms that protein ids belong to.
#' @importFrom trackViewer lolliplot
#' @importFrom seqinr read.fasta
#'
#' @return plot
#' @export

drawProt<-function(table, protein_id, organisms){

  glist<-list()
  flist<-list()

  for(i in 1:length(protein_id)){

    info<-getinfo(protid = protein_id[i], organism = organisms[i])
    colors<-c("#FF1F5B","#009ADE")
    table<-table[table[[paste0(organisms[i], "_ID")]] == protein_id[i],]
    SNP <- table[[paste0(organisms[i], "_aapos")]]
    gene1 <- GRanges("chr1", IRanges(SNP, width=1, names=paste0(table[[paste0(organisms[i], "_from")]],
                                                                table[[paste0(organisms[i], "_aapos")]],
                                                                table[[paste0(organisms[i], "_to")]])))
    gene1$border <- sample(c("gray30"), length(SNP), replace=TRUE)

    gene1.rot <- gene1
    gene1.rot$label.parameter.rot <- 45
    gene1.rot$SNPsideID<-"top"

    features1 <- GRanges("chr1", IRanges(start=c(1, info$tbl$pfam_start, info$len),
                                         width=c(0, info$tbl$pfam_end-info$tbl$pfam_start, 0),
                                         names=c("", info$dom, "")))
    features1$fill <- c("#FF1F5B", colors[1:length(info$tbl$pfam_start)], "#009ADE")
    features1$height <- c(0, rep(0.2, length(info$tbl$pfam_start)), 0)

    glist[i]<-gene1.rot
    flist[i]<-features1

  }

  trackViewer::lolliplot(glist,
                         flist,
                         yaxis=FALSE, ylab = protein_id, type = "pie")
}
