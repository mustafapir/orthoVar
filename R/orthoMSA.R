#' Create multiple sequence alignment table for use in orthoFind function
#'
#' This function takes species names and other optional parameters as argument and exports a data frame which is suitable to use as input in
#' orthoFind function to find orthologous variants between species.
#' By default, gene orthology data is taken from AllianceGenome website and custom orthology data can be used by specifying `customOrt`
#' parameter. For details, please check Readme file.
#'
#' @param species1 Name of the main species name which will be used in msa. Default is `Homo sapiens`.
#' @param species A character string or character vector specifying other species names used in msa
#' @param humanSeqFile Path of fasta file consisting of human protein sequences. Default is `NA`, which downloads file from NCBI.
#' @param seqFiles A character string or character vector specifying path of fasta files consisting of protein sequences of other
#' species. Default is `NA`, which downloads files from NCBI.
#' @param customOrt data frame consisting of gene orthology data for given species. Default is `NA`, which takes data from AllianceGenome.
#' @param annot source of annotation. Either "ncbi" or "ensembl" can be used. Default is ncbi.
#' @importFrom R.utils gunzip
#' @importFrom seqinr read.fasta
#'
#' @return dataframe object
#' @export


orthoMSA<-function(species1 = "Homo sapiens", species, humanSeqFile = NA, seqFiles = NA, customOrt = NA, annot = "ncbi"){

  if(!is.na(customOrt)){
    orthology<-customOrt
  }

  if(annot == "ncbi"){
    for(i in 1:7){
      martData[[i]][2]<-"refseq_peptide"
    }
    annot1<-"refseq_peptide"
  }
  else if(annot == "ensembl"){
    martData<-martData
    annot1<-"ensembl_peptide_id"
  }
  else {stop("`annot` must be either `ncbi` or `ensembl`")}

  if(is.na(humanSeqFile)){
    if(!dir.exists(file.path(getwd(), "human_sequence_file"))) {dir.create(file.path(getwd(), "human_sequence_file"), showWarnings = FALSE)}
    cat("\n Downloading Homo sapiens protein sequence fasta file.. \n")
    download.file(downloadLinks[["Homo sapiens"]], destfile = paste0(file.path(getwd(), "human_sequence_file"), "/Homo_sapiens_protein.faa.gz"),
                  method = "auto", quiet = TRUE)
  }

  if(is.na(seqFiles)){
    cat("\n Downloading other protein sequence fasta files.. \n")
    if(!dir.exists(file.path(getwd(), "other_sequence_files"))) {dir.create(file.path(getwd(), "other_sequence_files"), showWarnings = FALSE)}
    for(i in 1:length(species)){
      download.file(downloadLinks[[species[i]]], destfile = paste0(file.path(getwd(), "other_sequence_files"), "/", i, "_", species[i], "_protein.faa.gz"),
                    method = "auto", quiet = TRUE)
    }
  }

  cat("\n Reading fasta files.. \n")

  hfpath<-file.path(getwd(), "human_sequence_file")
  fpath<-file.path(getwd(), "other_sequence_files")
  seqHFiles1<-list.files(path = hfpath, full.names = TRUE)
  seqFiles1<-list.files(path = fpath, full.names = TRUE)

  if(tools::file_ext(seqHFiles1) == "gz") {R.utils::gunzip(seqHFiles1)}
  if(tools::file_ext(seqFiles1) == "gz") {Map(R.utils::gunzip, seqFiles1)}

  humanSeqFile<-list.files(path = "human_sequence_file", full.names = TRUE)
  seqFiles<-list.files(path = "other_sequence_files", full.names = TRUE)

  f<-function(aa, bb){
    eval(substitute(a <- b, list(a = as.name(aa), b = bb)))
  }
  seqList<-Map(f, paste0("file_", 1:length(seqFiles)), Map(function(a, x, y){seqinr::read.fasta(a, seqtype = x, as.string = y)},
                                                           a = seqFiles, x = "AA", y = TRUE))

  orthologyy<-orthology %>%
    dplyr::select(Gene1Symbol)
  martList<-list()
  martRefseq<-list()
  cat("\n Preparing tables.. \n")
  for(i in 1:length(species)){
    orthologyx<-orthology %>% dplyr::filter(Gene2SpeciesName == species[i]) %>% dplyr::distinct()
    orthologyy<-merge(orthologyy, orthologyx[,c(1,3)], by = "Gene1Symbol", all = TRUE, allow.cartesian = TRUE)
    colnames(orthologyy)[i+1]<-paste0("Gene_name_",i)
    martList<-c(martList, biomaRt::useMart("ENSEMBL_MART_ENSEMBL", martData[[species[i]]][1]))
    martRefseq<-c(martRefseq, list(biomaRt::getBM(c("external_gene_name", paste(martData[[species[i]]][2])), mart = martList[[i]]) %>%
                                     dplyr::filter(get(martData[[species[i]]][2]) != "" & !is.na(get(martData[[species[i]]][2]))) %>%
                                     dplyr::filter(external_gene_name != "")))
    colnames(martRefseq[[i]])[2]<-paste0("refseq_", i)
  }
  martList<-c(martList, biomaRt::useMart("ENSEMBL_MART_ENSEMBL", martData[["Homo sapiens"]][1]))
  martRefseq<-c(martRefseq, list(biomaRt::getBM(c("external_gene_name", paste(martData[["Homo sapiens"]][2])), mart = martList[[length(martList)]]) %>%
                                   dplyr::filter(get(martData[["Homo sapiens"]][2]) != "" & !is.na(get(martData[["Homo sapiens"]][2]))) %>%
                                   dplyr::filter(external_gene_name != "")))

  orthologyy<-unique(orthologyy) %>%
    dplyr::filter_at(dplyr::vars(-Gene1Symbol), dplyr::any_vars(!is.na(.)))

  df<-merge(orthologyy, martRefseq[[length(martRefseq)]], by.x = "Gene1Symbol", by.y = "external_gene_name")
  for(i in 1:length(species)){
    df<-merge(df, martRefseq[[i]], by.x = paste0("Gene_name_",i), by.y = "external_gene_name", allow.cartesian = TRUE, all = TRUE)
  }
  df<-df %>% dplyr::filter(!is.na(get(annot1)))
  df<-unique(df)

  humanSeq<-seqinr::read.fasta(humanSeqFile, seqtype = "AA", as.string = TRUE)
  humanSeq<-data.frame(Human_seq = unlist(humanSeq))
  humanSeq[[annot1]]<-rownames(humanSeq)
  suppressWarnings(humanSeq<-tidyr::separate(humanSeq, get(annot1), annot1, sep = "\\."))
  df<-merge(df, humanSeq, by = annot1)


  for(i in 1:length(species)){
    seq<-data.frame(seq = unlist(seqList[[i]]))
    seq[[paste0("refseq_", i)]]<-rownames(seq)
    if (species[i] != "Caenorhabditis elegans"){
      suppressWarnings(seq<-tidyr::separate(seq, 2, paste0("refseq_", i), sep = "\\."))
    }
    colnames(seq)[1]<-paste0("sequence_", i)
    df<-merge(df, seq, by = paste0("refseq_", i), all = TRUE)
  }

  df<-df %>% dplyr::filter_at(dplyr::vars(-Gene1Symbol), dplyr::any_vars(!is.na(.)))
  df<-df %>% dplyr::filter(!is.na(Human_seq)) %>%
    dplyr::filter_at(dplyr::vars(-!dplyr::matches("sequence_[0-9]")), dplyr::any_vars(!is.na(.)))
  cat("\n Aligning sequences.. \n\n")
  Sys.sleep(1)
  pb <- pbapply::timerProgressBar(min = 1, max = length(df[[1]]), style = 2)
  seqdf<-data.frame(matrix(NA, ncol = 2*(length(species)+1)))
  for(i in 1:length(df[[1]])){
    seqlength<-length(species) + 1
    seqchar<-as.character(df[i,(2*seqlength+1):length(df)])

    names(seqchar)<-df[i,seqlength:1]
    k<-which(!is.na(seqchar))
    invisible(capture.output(alignment<-msa(seqchar[k], type = "protein")))
    seqdf[i,2*k-1]<-BiocGenerics::rownames(alignment)
    for(j in 1:length(k)){seqdf[i,2*j]<-toString(unmasked(alignment)[j])}
    #seqdf[i,2*k]<-toString(unmasked(alignment))
    pbapply::setTimerProgressBar(pb, i)
  }
  speciesx<-c(species1, species)
  for(i in 1:length(speciesx)){
    colnames(seqdf)[2*i-1]<-paste0(speciesx[i],"_ID")
    colnames(seqdf)[2*i]<-paste0(speciesx[i], "_seq")
  }
  seqdf$index<-1:length(seqdf[[1]])
  return(seqdf)
}
