library(tidyverse); library(haven)

dat <- haven::read_dta(paste0(alphapath,
                                      "/ALPHA/clean_data/",sitename,
                                      "/residency_",sitename,".dta"))
dat <- dat %>% dplyr::select(idno, study_name, entry_date, exit_date, residence)

## For Kisumu keep only Gem
lowsite <- tolower(sitename)
if(lowsite=="kisumu") {
  dat <- dat %>% dplyr::filter(residence==2)
}

dat <- dat %>% dplyr::select(-residence)

dat <- dat %>% dplyr::group_by(study_name,idno) %>%
  dplyr::arrange(idno, entry_date) %>%
  dplyr::mutate(order = row_number()) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(names_from=order,names_sep="",values_from=c(entry_date, exit_date)) %>%
  dplyr::arrange(study_name,idno)

# save the dataset
if(!dir.exists(paste0(alphapath,"/ALPHA/prepared_data/",sitename))) {
  dir.create(paste0(alphapath,"/ALPHA/prepared_data/",sitename),recursive=TRUE)
} 
saveRDS(dat, paste0(alphapath,"/ALPHA/prepared_data/",sitename,
               "/residency_summary_dates_",sitename,".rds"))
