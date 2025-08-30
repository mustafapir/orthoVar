#' Create multiple sequence alignment table for use in orthoFind function
#'
#' This function takes species names and other optional parameters as argument and exports a data frame which is suitable to use as input in
#' orthoFind function to find orthologous variants between species.
#' By default, gene orthology data is taken from AllianceGenome website and custom orthology data can be used by specifying `customOrt`
#' parameter. For details, please check Readme file.
#'
#' @param species1 Name of the main species name which will be used in msa. Default is `Homo sapiens`.
#' @param species A character string or character vector specifying other species names used in msa
#' @param seqFile1 Path of fasta file consisting of protein sequence of first organism. Default is `NA`, which downloads the file from NCBI.
#' @param seqFiles A character string or character vector specifying path of fasta files consisting of protein sequences of
#' species indicated in `species` parameter. Default is `NA`, which downloads files from NCBI.
#' @param customOrt data frame consisting of gene orthology data for given species. Default is `NA`, which takes data from AllianceGenome.
#' If your species is not among the ones indicated in output of `listSpecies()` function, then use either custom orthology list or set this
#' to `customOrt = "ensembl"` to get orthology data from biomart.
#' @param annot source of annotation. Either "ncbi" or "ensembl" can be used. Default is ncbi.
#' @importFrom R.utils gunzip
#' @importFrom seqinr read.fasta
#'
#' @return dataframe object
#' @export


orthoMSA <- function(species1 = "Homo sapiens", species, seqFile1 = NA, seqFiles = NA, customOrt = NaN, annot = "ncbi", ...) {
  options(warn = -1)

  `%dopar%` <- foreach::`%dopar%`

  # Set orthology data ----
  Sys.sleep(1)
  cat(paste0("\r", "Downloading protein fasta files"))
  if (length(customOrt) > 0 && all(customOrt != "ensembl")) {
    orthology <- customOrt
    orthology <- ortho_convert(species1, species)
  } else if (customOrt == "ensembl") {
    orthology <- ort(species1, species)
  } else {
    orthology <- ortho_convert(species1, species)
  }


  # Check whether orthology data is a list ----
  checkort <- inherits(orthology, "list")


  # Set annotations ----
  martData <- ifelse(annot == "ncbi", "refseq_peptide", "ensembl_peptide_id")

  cat(paste0("\r", "Downloading protein fasta files.  "))

  # Get urls ----
  urls <- getlinks(species1, species, annot)

  # Download fasta files ----
  if (!dir.exists(file.path(getwd(), "sequence_files"))) {
    dir.create(file.path(getwd(), "sequence_files"), showWarnings = FALSE)
  }
  urlnumb <- ifelse(!is.na(seqFile1), 2, 1) # in the case where user enter file manually,
  urlnumb2 <- ifelse(!is.na(seqFiles), 1, (length(species) + 1)) # arrange the naming of files accordingly.
  if ((urlnumb2 - urlnumb >= 0)) {
    Sys.sleep(1)
    cat(paste0("\r", "Downloading protein fasta files..  "))

    dest <- paste0(file.path(getwd(), "sequence_files"), "/seq_", seq(urlnumb, urlnumb2, 1), ".faa.gz")
    urls_2 <- urls[urlnumb:urlnumb2]

    Map(function(u, d) download.file(u, d, method = "auto", quiet = TRUE), urls_2, dest)

    Sys.sleep(1)
    cat(paste0("\r", "Downloading protein fasta files... "))
  }

  # Read fasta files ----
  Sys.sleep(1)
  flush.console()
  Sys.sleep(1)
  cat(paste0("\r", "Reading fasta files               "))

  if (!is.na(seqFiles)) {
    file.copy(
      from = seqFiles,
      to = paste0("sequence_files/", "seq_", seq(2, length(species) + 1, 1), ".faa.gz")
    ) # copy manually entered files
  }

  if (!is.na(seqFile1)) {
    file.copy(
      from = seqFile1,
      to = paste0("sequence_files/seq_1.faa.gz")
    )
  }

  Sys.sleep(1)
  cat(paste0("\r", "Reading fasta files.              "))

  fpath <- file.path(getwd(), "sequence_files")
  seqFiles1 <- list.files(path = fpath, full.names = TRUE)
  purrr::map_if(seqFiles1, function(x) tools::file_ext(x) == "gz", function(x) Map(R.utils::gunzip, x)) # unzip all files in the directory
  sfs <- list.files(path = "sequence_files", full.names = TRUE)

  Sys.sleep(1)
  cat(paste0("\r", "Reading fasta files..             "))

  f <- function(aa, bb) {
    eval(substitute(a <- b, list(a = as.name(aa), b = bb)))
  }
  seqList <- Map(f, paste0("file_", 1:length(sfs)), Map(function(a, x, y) {
    seqinr::read.fasta(a, seqtype = x, as.string = y)
  },
  a = sfs, x = "AA", y = TRUE
  ))

  Sys.sleep(1)
  cat(paste0("\r", "Reading fasta files...            "))

  # Map annotations ----
  ortx <- data.table::rbindlist(orthology, use.names = FALSE)
  orty <- ortx %>%
    dplyr::filter(Gene1SpeciesName == species1) %>%
    dplyr::select(Gene1Symbol) %>%
    dplyr::distinct()

  Sys.sleep(1)
  flush.console()
  Sys.sleep(1)
  cat(paste0("\r", "Preparing tables                  "))

  allsp <- c(species1, species)
  nms2 <- paste0("Gene_name_", allsp)
  sseq <- seq(1, 3 * (length(species) - 1) + 1, 3)
  final_ort <- orthology %>%
    purrr::reduce(dplyr::full_join, by = "Gene1Symbol") %>% # merge orthology data of all species
    dplyr::select(1, all_of(sseq + 1)) %>%
    dplyr::distinct() %>%
    dplyr::filter_at(dplyr::vars(-Gene1Symbol), dplyr::any_vars(!is.na(.)))

  colnames(final_ort) <- c(paste0("Gene_name_", allsp))

  Sys.sleep(1)
  cat(paste0("\r", "Preparing tables.                 "))

  martref <- prot(species1, species, martData)

  Sys.sleep(1)
  cat(paste0("\r", "Preparing tables..                "))

  for (i in seq_along(allsp)) { # merge orthology and annotations
    final_ort <- merge(final_ort, martref[[i]],
      by = colnames(final_ort)[i],
      all = TRUE, allow.cartesian = TRUE
    )
  }


  for (i in seq_along(allsp)) {
    seq_names <- paste0("Refseq_", allsp[i])

    Seq1 <- data.frame(unlist(seqList[[i]]))
    names(Seq1) <- paste0("seq_", allsp[i])
    Seq1[[seq_names]] <- rownames(Seq1)
    suppressWarnings(Seq1 <- tidyr::separate(Seq1, get(seq_names), seq_names, sep = "\\."))
    final_ort <- merge(final_ort, Seq1, by = seq_names, all = TRUE)
  }

  rowAny <- function(x) rowSums(x) > 1
  final_ort <- final_ort %>% # remove rows having only one non-NA value
    dplyr::filter(rowAny(dplyr::across(
      .cols = dplyr::starts_with("seq"),
      .fns = ~ !is.na(.x)
    )))

  cat(paste0("\r", "Preparing tables...               "))


  # Sequence alignment ----
  Sys.sleep(2)
  flush.console()
  Sys.sleep(1)
  cat(paste0("\r", "Aligning sequences                "))
  Sys.sleep(2)

  seqlength <- length(species) + 1

  cl <- snow::makeCluster(parallel::detectCores()-1)
  doSNOW::registerDoSNOW(cl)
  pb <- pbapply::timerProgressBar(min = 1, max = length(final_ort[[1]]), style = 2)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)


  seqdf <- foreach::foreach(i = 1:length(final_ort[[1]]), .packages = c('msa', 'BiocGenerics', 'seqinr'), .combine = "rbind", .options.snow = opts, .errorhandling = "pass") %dopar% {
    seqchar <- as.character(final_ort[i, (2 * seqlength + 1):length(final_ort)])

    names(seqchar) <- final_ort[i, seqlength:1]
    k <- which(!is.na(seqchar))

    #seqinr::write.fasta(as.list(strsplit(seqchar[k], "\t")), names(seqchar), paste0("seqfile_", i, ".fasta"), open = "w", nbchar = 60, as.string = FALSE)
    fastafile<-Biostrings::AAStringSet(seqchar[k])
    Biostrings::writeXStringSet(fastafile, paste0("seqfile_", i, ".fasta"))

    invisible(capture.output(alignment <- msa(paste0("seqfile_", i, ".fasta"), type = "protein", order = "input")))
    unlink(paste0("seqfile_", i, ".fasta"))
    unlink(paste0("seqfile_", i, ".aln"))
    unlink(paste0("seqfile_", i, ".dnd"))

    lenx<- 2 * (length(species) + 1)
    seqdf_mid<-data.frame(matrix(NA, ncol = lenx))
    #seqdf_mid<-character(lenxL)
    #seqdf_mid[1, 2 * k - 1] <- BiocGenerics::rownames(alignment)
    seqdf_mid[1, 2 * k - 1] <- names(seqchar[k])
    for (j in k) {
      seqdf_mid[1, 2 * j] <- toString(unmasked(alignment)[which(k == j)])
    }

    return(seqdf_mid)
  }
  pbapply::closepb(pb)
  parallel::stopCluster(cl)

  for (i in 1:length(final_ort[[1]])) {
    unlink(paste0("seqfile_", i, ".fasta"))
    unlink(paste0("seqfile_", i, ".aln"))
    unlink(paste0("seqfile_", i, ".dnd"))
  }

  speciesx <- c(species1, species)
  for (i in 1:length(speciesx)) {
    colnames(seqdf)[2 * i - 1] <- paste0(speciesx[i], "_ID")
    colnames(seqdf)[2 * i] <- paste0(speciesx[i], "_seq")
  }
  seqdf$index <- 1:length(seqdf[[1]])
  return(seqdf)
}
