suppressPackageStartupMessages(library(PerfMeas))
suppressPackageStartupMessages(library(dplyr))


calculate_ranked_AUPRC <- function(rank, actual){
    df <- dplyr::data_frame(
        "rank" = rank,
        "actual" = actual)
    df <- calculate_precision_and_recall(df)
    lst <- list("lst" = list("precision" = df$precision, "recall" = df$recall))
    return(PerfMeas::AUPRC(lst))
}

calculate_precision_and_recall <- function(df){
    df %>% 
        dplyr::arrange(rank) %>% 
        dplyr::mutate(predicted = ifelse(is.na(rank), 0, 1)) %>% 
        dplyr::mutate(predicted = ifelse(is.nan(rank), 0, 1)) %>% 
        dplyr::mutate(true_positive = 
                          ifelse((actual == 1 & predicted == 1), 1, 0)) %>% 
        dplyr::mutate(false_positive = 
                          ifelse((actual == 0 & predicted == 1), 1, 0)) %>% 
        calculate_recall %>% 
        calculate_precision
}

calculate_recall <- function(df){
    dplyr::mutate(
        df,
        recall = cumsum(true_positive) / sum(actual))
}

calculate_precision <- function(df){
    dplyr::mutate(
        df, 
        precision = 
            cumsum(true_positive) / 
            (cumsum(true_positive) + cumsum(false_positive)))
}