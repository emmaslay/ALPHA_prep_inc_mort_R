dat$yob <- format(dat$dob,"%Y")

dat <- dat %>% 
    dplyr::mutate(birth_cohort = cut(as.numeric(yob), seq(1900,2030,5), include.lowest=TRUE,
                                     right=FALSE, labels=1:26)) %>%
    dplyr::select(-yob)
