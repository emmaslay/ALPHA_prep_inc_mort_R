dat <- dat %>%
    dplyr::mutate(fouryear = cut(years_one, c(1989,2000,seq(2005,2025,4)),
                                 include.lowest=TRUE, right=FALSE, labels=0:6))
