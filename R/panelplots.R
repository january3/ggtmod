
## a helper function to make names of gene sets from MSigDB easier to read
.cleanup_names <- function(ids) {
  ids <- gsub("_", " ", ids)

  #return(ids)

  ids <- strsplit(ids, " ")
  min.l <- min(sapply(ids, length))

  max_prefix <- 0
  for(i in 1:min.l) {

    if(length(unique(sapply(ids, function(x) x[i]))) == 1L) {
      max_prefix <- i
    }
  }

  if(max_prefix > 0) {
    ids <- lapply(ids, function(x) x[ -1:-max_prefix ])
  }

  ids <- unlist(lapply(ids, paste, collapse=" "))

  if(!any(grepl("[a-z]", ids))) {
    ids <- tolower(ids)
    substr(ids, 1, 1) <- toupper(substr(ids, 1, 1))
  }

  ids
}




#' Plot gene set enrichments as a panelplot
#'
#' Generate a ggplot panelplot from a data frame containing gene set enrichments.
#'
#' And here the gory details.
#' @param df data frame containing the gene set enrichments
#' @param es_col the column name with the effect sizes
#' @param pval_col column name with the (adjusted) p-values
#' @param cleanup_names if TRUE, the module names will be edited to make them more
#'                      legible
#' @return A ggplot object.
#' @importFrom ggplot2 ggplot xlab ylab geom_bar coord_flip aes_string
#' @examples
#' # no example yet. soon!
#' @export
plot_gsea <- function(df, es_col="E", pval_col="adj.P.Val",
                          cleanup_names=TRUE) {

  df <- df[ order(-df[[pval_col]]), ]

  df$log10_pval <- -log10(df[[pval_col]])

  if(cleanup_names) {
    df$Title <- .cleanup_names(df$Title)
  }
  df$Title <- factor(df$Title, levels=df$Title)

  ggplot(df,
         aes_string(x="Title", y="log10_pval", fill=es_col)) + geom_bar(stat="identity") +
    xlab("") +
    ylab("-log10(FDR)") +
    coord_flip()


}
