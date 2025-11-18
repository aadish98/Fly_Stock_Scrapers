# Section 1 -----------------------------------------------------------------
# Data Wrangling ------------------------------------------------------------

filter_df_list_on_freq <- function(df_list){
  for(df_name in names(df_list)){
    df <- df_list[[df_name]] %>%
      mutate(frequency = as.numeric(frequency)) %>%
      filter(frequency > 1) %>%
      rename(experiment_frequency = frequency) %>%
      arrange(desc(experiment_frequency))
    df_list[[df_name]] <- df
  }
  return(df_list)
}

get_todays_version <- function(){
  return(paste0("v", format(Sys.Date(), "%m%d%y")))
}

read_all_csvs_recursive <- function(folder_path) {
  # List all CSV files in the folder and its subfolders
  csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  # Initialize a list to store the data frames
  data_list <- list()
  
  # Loop through each CSV file and read it into a data frame
  for (csv_file in csv_files) {
    df <- read.csv(csv_file, check.names = FALSE, encoding= 'Latin-1')
    file_name <- tools::file_path_sans_ext(basename(csv_file))
    data_list[[file_name]] <- df
  }
  
  return(data_list)
}

read_all_csvs_as_dt_recursive <- function(folder_path) {
  # List all CSV files in the folder and its subfolders
  csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  # Initialize a list to store the data frames
  data_list <- list()
  
  # Loop through each CSV file and read it into a data frame
  for (csv_file in csv_files) {
    df <- read.csv(csv_file, check.names = FALSE, encoding= 'Latin-1')
    df <- as.data.table(df)
    file_name <- tools::file_path_sans_ext(basename(csv_file))
    data_list[[file_name]] <- df
  }
  
  return(data_list)
}

get_sorted_df <- function(df, sorting_invariant, descending= FALSE){
  if(descending){
    df <- df %>%
      arrange(desc(sorting_invariant))
    return(df)
  } 
  df <- df %>%
    arrange(sorting_invariant)
  return(df)
}

get_latest_version_folderpath <- function(directory_path){
  directory_path <- file.path(directory_path)
  dirs <- list.dirs(directory_path, recursive = FALSE, full.names = TRUE)
  v_dirs <- dirs[grepl("v[0-1][0-9][0-3][0-9][0-9][0-9]$", basename(dirs))]
  
  # Extract the dates and convert to Date class for comparison
  dates <- as.Date(sub("v","", basename(v_dirs)), format = "%m%d%y")
  
  # Get the directory with the latest date
  latest_dir <- v_dirs[which.max(dates)]
  return(latest_dir)
}

get_all_sleep_wake_gene_sets <- function(folder_path, file_keywords, path_keywords, blacklist_keywords) {
  # List all CSV files in the folder and its subfolders
  # Create a regex pattern to match all keywords
  file_keyword_pattern <- paste(file_keywords, collapse = "|")
  path_keyword_pattern <- paste(path_keywords, collapse = "|")
  blacklist_keyword_pattern <- paste(blacklist_keywords, collapse = "|")
  # Filter files that have all the keywords in their filename
  csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  csv_files <- csv_files[sapply(csv_files, function(x) any(grepl(path_keyword_pattern, x, ignore.case = FALSE)))]
  csv_files <- csv_files[sapply(csv_files, function(x) any(grepl(file_keyword_pattern, tools::file_path_sans_ext(x), ignore.case = FALSE)))]
  csv_files <- csv_files[sapply(csv_files, function(x) any(!grepl(blacklist_keyword_pattern, tools::file_path_sans_ext(x), ignore.case = FALSE)))]
  print(csv_files)
  data_list <- list()
  filepath_list <- list()
  # Loop through each CSV file and read it into a data frame
  for (csv_file in csv_files) {
    df <- read.csv(csv_file)
    file_name <- tools::file_path_sans_ext(basename(csv_file))
    data_list[[file_name]] <- df
    filepath_list[[file_name]] <- csv_file
  }
  return(data_list)
}

flatten_list_columns <- function(df) {
  for (col in names(df)) {
    if (is.list(df[[col]])) {
      df[[col]] <- sapply(df[[col]], paste, collapse = "; ")
    }
  }
  return(df)
}

get_version_name <- function(folder_path){
  
  # Extract all components of the path
  path_components <- strsplit(folder_path, "/|\\\\")[[1]]
  
  # Get the last non-empty component
  last_folder <- tail(path_components[nchar(path_components) > 0], 1)
  
  return(last_folder)
}

write_all_results_csv <- function(df_list, dir_path, row.names = TRUE){
  # Check if the directory exists; if not, create it
  dir_path <- file.path(dir_path)
  if (!dir.exists(dir_path)) {
    print(paste0("Creating directory ", dir_path, " ..."))
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Iterate over the list of dataframes
  for (df_name in names(df_list)) {
    # Construct the file path
    file_path <- file.path(dir_path, paste0(df_name, ".csv"))
    
    # Write the dataframe to the file
    print(df_list[[df_name]])
    write.csv(df_list[[df_name]], file_path, row.names = TRUE)
  }
  print(paste0("Results written in ", dir_path, "!"))
}


# Section 2 -----------------------------------------------------------------
# Selenium Client Server ------------------------------------------------------
setup_selenium_firefox_client <- function(){
  
  # Close the remote driver object to stop the Selenium server
  res <- list()
  res["remDr"] <- remDr
  res["rD"] <- rD
  return(res)
}

close_selenium_server <- function(remDr, rD){
  remDr$close()
  rD$server$stop()
}

prepare_BDSC_stock_search <- function(remDr){
  url = "https://bdsc.indiana.edu/Home/Search"
  # For navigating to the url 
  remDr$navigate(url)
  Sys.sleep(2)
}  

# Section 3 -----------------------------------------------------------------
# Querying on BDSC Website --------------------------------------------------

# Gets list of symbol id's for relevant gene driver lines
get_BDSC_driver_symbol_list <- function(gene_name, driver_line="GAL4"){
  # Select 'all' for gene constraints (open list, then select all)
  remDr$findElement(using = 'xpath', value = "/html/body/main/section[1]/div/div/div/div/table[1]/tbody/tr[1]/td/span/span[1]/span")$clickElement()
  remDr$findElement(using = 'xpath', value = "/html/body/div[2]/div/div/div/div[2]/ul/li[2]")$clickElement()
  
  # Select 'exactly matching' for gene name
  remDr$findElement(using = 'xpath', value = '/html/body/main/section[1]/div/div/div/div/table[1]/tbody/tr[2]/td[2]/span[1]/span[1]/span')$clickElement()
  remDr$findElement(using = 'xpath', value = '/html/body/div[4]/div/div/div/div[2]/ul/li[3]/span')$clickElement()
  
  # Add gene
  
  text_field1 <- remDr$findElement(using = "css", value = "#gene1")
  text_field1$clearElement()
  text_field1$sendKeysToElement(list(gene_name))
  
  # Add 2nd gene of interest
  if (driver_line != 'Mi{MIC}'){
    value_text_field <- '#gene2'
  } else {
    value_text_field <- '#componentSymbol'
  }
  
  text_field2 <- remDr$findElement(using = "css", value = value_text_field)
  text_field2$clearElement()
  text_field2$sendKeysToElement(list(driver_line))
  Sys.sleep(0.2)
  
  # Click search
  remDr$findElement(using = 'xpath', value = '//*[@id="AdvancedSearch"]')$clickElement()
  Sys.sleep(4)
  
  # Process symbol table
  table_element <- remDr$findElement(using = "class name", value = 'k-table')
  Sys.sleep(0.5)
  
  # Extract HTML content of the table
  table_html <- table_element$getElementAttribute("outerHTML")[[1]]
  
  # Convert HTML content to a dataframe (example using xml)
  table_data <- data.frame(readHTMLTable(table_html))
  col_name <- names(table_data)[grepl("Component.Symbol", names(table_data))][1]
  
  text_field1$clearElement()
  text_field2$clearElement()
  return(table_data[[col_name]])
}

get_BDSC_responder_RNAi_line_symbol_list <- function(gene_name, driver_line="GAL4"){
  # Select 'all' for gene constraints (open list, then select all)
  remDr$findElement(using = 'xpath', value = "/html/body/main/section[1]/div/div/div/div/table[1]/tbody/tr[1]/td/span/span[1]/span")$clickElement()
  remDr$findElement(using = 'xpath', value = "/html/body/div[2]/div/div/div/div[2]/ul/li[2]")$clickElement()
  
  # Select 'RNAi' for the categories field
  
  full_xpath <- '/html/body/main/section[1]/div/div/div/div/table[1]/tbody/tr[2]/td[1]/span'
  remDr$findElement(using = 'xpath', value = full_xpath)$clickElement()
  
  # Select 'exactly matching' for gene name

  full_xpath <- '/html/body/div[3]/div/div/div/div[2]/ul/li[16]'
  remDr$findElement(using = 'xpath', value = full_xpath)$clickElement()
  # Sys.sleep(1)
  
  remDr$findElement(using = 'xpath', value = '/html/body/main/section[1]/div/div/div/div/table[1]/tbody/tr[2]/td[2]/span[1]/span[1]/span')$clickElement()
  remDr$findElement(using = 'xpath', value = '/html/body/div[4]/div/div/div/div[2]/ul/li[3]/span')$clickElement()
  
  # Add gene
  
  text_field1 <- remDr$findElement(using = "css", value = "#gene1")
  text_field1$clearElement()
  text_field1$sendKeysToElement(list(gene_name))
  
  # Click search
  remDr$findElement(using = 'xpath', value = '//*[@id="AdvancedSearch"]')$clickElement()
  Sys.sleep(6)
  
  # Process symbol table
  table_element <- remDr$findElement(using = "class name", value = 'k-table')
  Sys.sleep(1)
  
  # Extract HTML content of the table
  table_html <- table_element$getElementAttribute("outerHTML")[[1]]
  
  # Convert HTML content to a dataframe (example using xml)
  table_data <- data.frame(readHTMLTable(table_html))
  col_name <- names(table_data)[grepl("Component.Symbol", names(table_data))][1]
  
  text_field1$clearElement()
  return(table_data[[col_name]])
}

get_stock_list <- function(symbol_id_list){
  
  stock_list <- list()
  comment_list <- list()
  genotype_list <- list()
  symbol_id_list <- gsub(" ", "", symbol_id_list, fixed = TRUE)
  symbol_id_format <- strsplit(symbol_id_list, ";")
  if(is.na(symbol_id_format)){
    return(NA)
  }
  if(length(symbol_id_format) == 0){
    return(NA)
  }
  for (symbol_id in symbol_id_format){
    if(all(is.na(symbol_id))){
      next
    }
    if(all(symbol_id == "")){
      next
    }
    for(real_sym_id in symbol_id){
      #Add symbol_id from list
      if((is.na(real_sym_id))){
        next
      }
      if((real_sym_id == "")){
        next
      }
      text_field1 <- remDr$findElement(using = "xpath", value = '//*[@id="presearch"]')
      text_field1$clearElement()
      text_field1$sendKeysToElement(list(real_sym_id))
      # Click search
      remDr$findElement(using = 'xpath', value = '//*[@id="SimpleSearch"]')$clickElement()
      Sys.sleep(6)
      # Process stock table
      table_element1 <- remDr$findElement(using = "css", value = '.k-master-row > td:nth-child(2)')
      stockid <- table_element1$getElementAttribute("innerHTML")[[1]]
      Sys.sleep(1)
      xpath <- "//div[contains(text(), 'Chr ')]"
      table_element2 <- remDr$findElement(using = "xpath", value = xpath)
      comments <- table_element2$getElementText()[[1]]
      Sys.sleep(1)
      table_element3 <- remDr$findElement(using = "css", value = 'td.k-table-td:nth-child(3)')
      # Extract HTML content of the table
      genotype <- table_element3$getElementAttribute("innerHTML")[[1]] # getElementText()[[1]]
      Sys.sleep(1)
      stock_list <- c(stock_list, stockid)
      comment_list <- c(comment_list, comments)
      genotype_list <- c(genotype_list, genotype)
      #stock_list <- c(stock_list, paste0("stock ID - ", stockid, " ",genotype, " - ", comments))
    }
    Sys.sleep(1)
  }
  stock_string <- paste(stock_list, collapse = ", ")
  comments_string <- paste(comment_list, collapse = '\n')
  genotype_string <- paste(genotype_list, collapse = '\n')
  res <- list(
    stock_id = stock_string,
    map = comments_string,
    genotypes = genotype_string
  )
  return(res)
} 

# Section 4 -----------------------------------------------------------------
# Querying on BDSC DB's -----------------------------------------------------

get_filtered_db <- function(dt, sub_str_list, col_name) {
  for (sub_str in sub_str_list) {
    pattern <- paste0("\\b", sub_str, "\\b")
    dt <- dt[str_detect(get(col_name), regex(pattern, ignore_case = FALSE))]
  }
  return(dt)
}

get_filtered_db_multiple_cols <- function(dt, sub_str_list, col_names) {
    dt_filtered <- dt
    for (sub_str in sub_str_list) {
      pattern <- paste0("\\b", sub_str, "\\b")
      dt_filtered <- dt_filtered[Reduce(`|`, lapply(col_names, function(col_name) {
        str_detect(get(col_name), regex(pattern, ignore_case = FALSE))
      }))]
    }
    return(dt_filtered)
}





