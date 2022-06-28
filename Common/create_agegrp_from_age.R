
dat <- dat %>%
    dplyr::mutate(agegrp = cut(age, c(seq(0,90,5),150), include.lowest=TRUE, right=FALSE,
                               labels = 0:18))
