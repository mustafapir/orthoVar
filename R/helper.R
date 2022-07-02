#'
#' Helper functions
#'
#'
#' @export

listSpecies<-function(...){
  names(martDatax)
}


#' Helper functions
#'
#' @param species get species from parent function
#' @noRd


shortnames<-function(species1, species){
  allsp<-c(species1, species)
  gsub('(\\b\\pL).*? ', '\\L\\1', allsp, perl = TRUE)
}

ort<-function(species1, species){

  spnames<-shortnames(species1, species)

  ensembl = biomaRt::useMart("ensembl", dataset=paste0(spnames[1], "_gene_ensembl"))
  lst<-lapply(spnames[2:length(spnames)], function(x) biomaRt::getBM(attributes = c("ensembl_gene_id",
                                                                                    paste0(x, "_homolog_ensembl_gene"),
                                                                                    paste0(x, "_homolog_orthology_confidence")),
                                                                     mart = ensembl))


  lst1<-Map(function(x, y) {
    a<-x[x[,2] != "" & x[,3] == 1, 1:2]
    colnames(a)<-c("Gene1Symbol", "Gene2Symbol")
    a$Gene1SpeciesName<-species1
    a$Gene2SpeciesName<-y
    a
  }, lst, species)

  return(lst1)

}



prot<-function(species1, species, martData){

  spnames<-shortnames(species1, species)
  allsp<-c(species1, species)

  martList<-lapply(spnames, function(x) {

    mlist<-biomaRt::useMart("ensembl", paste0(x, "_gene_ensembl"))
    mlist
  })

  martRefseq<-Map(function(x, y) {

    mrefseq<-biomaRt::getBM(c("ensembl_gene_id", martData), mart = x) %>%
      dplyr::filter(get(martData) != "" & !is.na(get(martData))) %>%
      dplyr::filter(ensembl_gene_id != "")

  }, martList, allsp)
  nms<-lapply(allsp, function(x) paste0(c("Gene_name_", "Refseq_"), x))
  martRefseq2<-Map(setNames, martRefseq, nms)

  return(martRefseq2)
}



getlinks<-function(species1, species, annot){

  allsp<-c(species1, species)
  allsp2<-tolower(substr(allsp, 1, 1))
  allsp21<-tolower(substr(allsp, 1, 100))
  allsp3<-stringr::str_split_fixed(allsp21, " ", 2)
  spnames<-gsub(" ", "+", allsp21)
  spnames2<-gsub(" ", "_", allsp21)

  if(annot == "ncbi"){

    all_links<-c()
    for (i in 1:length(allsp)){
      page<-rvest::read_html(paste0("https://www.ncbi.nlm.nih.gov/genome/?term=", spnames[i]))
      link<-page %>% rvest::html_nodes(xpath = '//*[@id="maincontent"]/div/div[1]/div[2]/span[1]/a[3]') %>% rvest::html_attr('href')
      all_links<-c(all_links, link)
    }

  } else {
    all_links<-c()
    for (i in 1:length(allsp)){
      page<-rvest::read_html(paste0("https://ftp.ensembl.org/pub/current_fasta/", spnames2[i], "/pep/"))
      link<-page %>% rvest::html_nodes(xpath = '/html/body/pre/a') %>% rvest::html_attr('href')
      link<-link[grepl("all",link)]
      link<-paste0("https://ftp.ensembl.org/pub/current_fasta/", spnames2[i], "/pep/", link)
      all_links<-c(all_links, link)
    }
  }


  return(all_links)

}

ortho_convert<-function(species1, species, ...){

  # filter data based on user input
  orr<-orthology %>%
    dplyr::filter(Gene1SpeciesName == species1, Gene2SpeciesName %in% species)

  # prepare organism names for use in biomart
  spnames<-shortnames(species1, species)

  # get annotations from ensembl
  ensembl<-lapply(spnames, function(x) biomaRt::useMart("ensembl", dataset=paste0(x, "_gene_ensembl")))
  lst<-lapply(ensembl, function(x) {
    biomaRt::getBM(attributes = c("external_gene_name", "ensembl_gene_id"), mart = x) %>%
      dplyr::filter(external_gene_name != "") %>%
      dplyr::rename(Gene2Symbol = external_gene_name) %>%
      dplyr::distinct()
  })

  colnames(lst[[1]])[1]<-"Gene1Symbol"

  # Join annotations to orthology table
  orr1<-dplyr::left_join(orr, lst[[1]], by = "Gene1Symbol")
  orr2<-Map(function(x,y) {

    y1<-dplyr::filter(orr1, Gene2SpeciesName == y)
    dplyr::left_join(y1, x, by = "Gene2Symbol") %>%
      dplyr::select(ensembl_gene_id.x, ensembl_gene_id.y, Gene1SpeciesName, Gene2SpeciesName) %>%
      dplyr::rename(Gene1Symbol = ensembl_gene_id.x, Gene2Symbol = ensembl_gene_id.y)
  },
  lst[-1], species)


  return(orr2)

}



