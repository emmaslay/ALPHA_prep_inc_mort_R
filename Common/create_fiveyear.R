dat <- dat %>% 
    dplyr::mutate(fiveyear = cut(years_one, c(1989,seq(1995,2025,5)), include.lowest=TRUE,
                                 right = FALSE, labels = c(1:7)))
