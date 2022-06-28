######################## PREPARE THE IMPUTATIONS

library(tidyverse);library(rlang)

nimp <- 70

dat <- readRDS(paste0(alphapath,"/ALPHA/incidence_ready_data/", sitename,
                      "/incidence_ready_midpoint_",sitename,".rds"))

if(!dir.exists(paste0(alphapath,"/ALPHA/incidence_ready_data/",
                      sitename,"/mi_data"))) {
  dir.create(paste0(alphapath,"/ALPHA/incidence_ready_data/",
                    sitename,"/mi_data"),recursive=TRUE)  
}


### GENERATE LIST OF IDNO WITH WHICH TO GENERATE A DATASET OF RANDOM NUMBERS FOR
### EACH PERSON
### # this is important because we need a number for each person not each record
### # in the dataset

dat_sub <- select(dat,study_name,idno) %>%
  distinct()

for(i in 1:nimp) {
  set.seed(i)
  dat_sub <- dat_sub %>%
    dplyr::mutate("random{i}" := runif(1,0,1))
}

saveRDS(dat_sub, file=paste0(alphapath,"/ALPHA/incidence_ready_data/",
                             sitename,"/mi_data/random_numbers_",sitename,".rds"))


#### generate the mi files with the seroconversion date randomly allocated within
#### the interval
#### # EACH DATASET CONTAINS ONE IMPUTATION OF THE SEROCONVERSION DATE, SO WE END UP
#### # WITH AS MANY DATASETS AS IMPUTATIONS
#### # THE FILE NAMING IS IMPORTANT AND MUST FOLLOW A PATTERN WITH ONLY THE IMPUTATION
#### # NUMBER CHANGING FOR EACH FILE

for(x in 1:nimp) {
  set.seed(x)
  print(paste0("Imputation ",x))
  
  dat <- readRDS(paste0(alphapath,"/ALPHA/incidence_ready_data/", sitename,
                        "/incidence_ready_midpoint_",sitename,".rds"))
  
  dat <- dplyr::left_join(dat,dat_sub)
  
  dat <- dat %>% dplyr::select(-sero_conv_date)
  rnum <- paste0("random",x)
  dat <- dat %>%
    dplyr::mutate(sero_conv_date = case_when(!is.na(first_pos_date) & !is.na(last_neg_date) &
                                               first_pos_date>last_neg_date ~ 
                                               ((first_pos_date - last_neg_date)*!!rlang::parse_expr(rnum)) +
                                               last_neg_date,
                                             tested_first_opportunity==1 ~ 
                                               (first_pos_date - fifteen)*!!rlang::parse_expr(rnum) + fifteen,
                                             TRUE ~ as.Date(NA,origin="1970-01-01")))
  
  dat <- dat %>%
    dplyr::mutate(last_neg_date = as.Date(last_neg_date,origin="1970-01-01"),
                  first_pos_date = as.Date(first_pos_date,origin="1970-01-01"),
                  sero_conv_date = as.Date(sero_conv_date,origin="1970-01-01"),
                  start_ep_date = as.Date(start_ep_date,origin="1970-01-01"),
                  end_ep_date = as.Date(end_ep_date,origin="1970-01-01"))
  
  
  # put as not a failure if not resident at time
  dat <- dat %>%
    dplyr::mutate(serocon_fail = ifelse(sero_conv_date>start_ep_date & 
                                          sero_conv_date<=end_ep_date &
                                          !is.na(sero_conv_date), 1, 0))
  dat <- dat %>%
    dplyr::mutate(end_ep_date = ifelse(serocon_fail==1,sero_conv_date,end_ep_date))
  
  dat <- dat %>%
    dplyr::mutate(last_neg_date = as.Date(last_neg_date,origin="1970-01-01"),
                  first_pos_date = as.Date(first_pos_date,origin="1970-01-01"),
                  sero_conv_date = as.Date(sero_conv_date,origin="1970-01-01"),
                  start_ep_date = as.Date(start_ep_date,origin="1970-01-01"),
                  end_ep_date = as.Date(end_ep_date,origin="1970-01-01"))
  
  # drop observations that begin on or before DOB
  # Take out records whose DOB is on or later than episode end date
  dat <- dat %>% dplyr::filter(dob<end_ep_date)
  # drop observations that begin on or after failure
  dat <- dat %>% group_by(idno) %>%
    mutate(first_sero_conv_date = as.Date(ifelse(is.finite(min(sero_conv_date,na.rm=TRUE)),
                                             min(sero_conv_date, na.rm=TRUE),NA),origin="1970-01-01"),
           on_after_sero_conv = ifelse(start_ep_date>=first_sero_conv_date,1,0)) %>%
    filter(is.na(on_after_sero_conv) | on_after_sero_conv!=1) %>%
    select(-c(first_sero_conv_date,on_after_sero_conv))
  # drop observations that have a start_ep_date after sero_conv_date
  dat <- dat %>%
    filter(!(start_ep_date>sero_conv_date & !is.na(sero_conv_date)))
  
  dat <- dat %>%
    dplyr::arrange(study_name,idno,start_ep_date) %>%
    dplyr::group_by(study_name,idno) %>%
    dplyr::mutate(ep_num = row_number(),
           entry = start_ep_date,
           exit = end_ep_date,
           t0 = as.numeric((start_ep_date-dob)/365.25),
           t = as.numeric((end_ep_date - dob)/365.25),
           d = serocon_fail) %>%
    dplyr::select(-starts_with("random"))
  
  saveRDS(dat,paste0(alphapath,"/ALPHA/incidence_ready_data/", sitename,
                     "/mi_data/incidence_ready_mi_",sitename,"_",x,".rds"))
  
}

