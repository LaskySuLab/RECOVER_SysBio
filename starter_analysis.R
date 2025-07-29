repo_top_loc <- gsub("recover-nbr.+", "recover-nbr", getwd())
setwd(repo_top_loc)

get_folder_path <- function(loc = "", fld_str) {
  if(nchar(loc) > 20) error("Unable to find desired folder")
  all_fls <- list.files(loc)
  if(fld_str %in% all_fls) {
    return(file.path(loc, fld_str))
  } else {
    get_folder_path(paste0(loc, "../"), fld_str)
  }
}

pf_loc <- get_folder_path(fld_str="project-files")
of_loc <- get_folder_path(fld_str="output-files")

source("helper_script.R")

load_libs(c("multcomp", "broom", "flextable", "officer"))

bargs_in <- getArgs(defaults = list(add_log=0))

today <- format(Sys.Date(), "%Y%m%d")

output_id = gsub("user.name=", "", grep("^user.name=", system('git config --list', intern=T), value=T))

if(length(output_id) == 0) {
  output_id = "no_git_username"
} else if(output_id %in% c(NA, "")) {
  output_id = "no_git_username"
}
 
output_loc <- glue("{of_loc}/{output_id}/{today}")
dir.create(output_loc, recursive = T)

if(bargs_in$add_log %in% 1){
  log_loc <- file.path(output_loc, "an_init.log")
  cat(glue("-------- for run info check log file at {log_loc} ----- \n\n"))
  sf <- file(log_loc, open = "wt")
  sink(sf, split = T)
  sink(sf, type = "message")
  
  cat(glue("------------------- Haack init file - {format(Sys.time(), '%H:%M')} ------------------- \n\n"))
  
}


adult_env_list <- get_env_list("adult")

ds_dd <- adult_env_list$ds_dd()
core <- adult_env_list$core()

ps_combined_pheno <- adult_env_list$ps_combined_pheno()
ps_pasc_ds <- adult_env_list$ps_pasc_ds()
formds_list <- adult_env_list$formds_list("pasc_symptoms", "biospecimens", "new_covid_infection", "additional_tests_calculations")

g1 <- core %>% 
  summarise(n_day = n(),
            .by=enroll_dt) %>% 
  arrange(enroll_dt) %>% 
  mutate(cm = cumsum(n_day)) %>% 
  ggplot(aes(x=enroll_dt, y=cm)) + 
  geom_line() + 
  scale_x_date(date_breaks = "6 months",
               labels = scales::label_date_short(),
               minor_breaks = NULL)

# This saves in output-files
# when you him stop on the top right it'll be copied into project-files
ggsave(glue("{output_loc}/example_output_{today}.jpg"), g1, width=5, height=5, dpi=200)



