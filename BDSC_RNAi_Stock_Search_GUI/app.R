options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Install BiocManager if not already installed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# library(BiocManager)
# BiocManager::install("biomaRt")
# Install biomaRt using BiocManager
# Set a writable library path
#.libPaths(c("/usr/lib/R", .libPaths()))
# Install and load the checkpoint package for reproducible environments

# Load necessary libraries
library(shiny)
library(bslib)
library(shinythemes)
library(shinyjs)
library(DT)
library(readr)
library(devtools)
library(ggplot2)
library(conflicted)
library(data.table)
library(stringr)
library(dplyr)
library(biomaRt)
library(tidyr)
library(readxl)
library(DBI)
library(RPostgres)  
library(readr)

# global helper functions -----------------------------------------------------------------------



read_all_csvs_as_dt_recursive <- function(folder_path) {
  
  # List all CSV files in the folder and its subfolders
  csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  # Initialize a list to store the data frames
  data_list <- list()
  
  # Loop through each CSV file and read it into a data frame
  for (csv_file in csv_files) {
    df <- read.csv(csv_file, check.names = FALSE)#, encoding= 'Latin-1')
    df <- as.data.table(df)
    file_name <- tools::file_path_sans_ext(basename(csv_file))
    data_list[[file_name]] <- df
  }
  
  return(data_list)
}

read_all_tsvs_as_dt_recursive <- function(folder_path) {
  
  # List all CSV files in the folder and its subfolders
  tsv_files <- list.files(folder_path, pattern = "\\.tsv$", full.names = TRUE, recursive = TRUE)
  # Initialize a list to store the data frames
  data_list <- list()
  
  # Loop through each CSV file and read it into a data frame
  for (tsv_file in tsv_files) {
    df <- readr::read_tsv(tsv_file)
    df <- as.data.table(df)
    file_name <- tools::file_path_sans_ext(basename(tsv_file))
    data_list[[file_name]] <- df
  }
  
  return(data_list)
}

get_filtered_db_multiple_cols <- function(dt_filtered, sub_str_list, col_names) {
  
  # Check for each substring in the list
  for (sub_str in sub_str_list) {
    pattern <- sub_str
    
    # Use Reduce to apply a logical OR across results from multiple columns
    filter_conditions <- lapply(col_names, function(col_name) {
      if (col_name %in% names(dt_filtered)) {
        # Return a logical vector indicating rows where the substring is found
        str_detect(dt_filtered[[col_name]], regex(pattern, ignore_case = FALSE))
      } else {
        # If the column does not exist, return a logical vector of FALSE with the same length as the dataframe
        # This means that missing columns won't affect the filtering
        rep(FALSE, nrow(dt_filtered))
      }
    })
    
    # Apply the filtering conditions
    dt_filtered <- dt_filtered[Reduce(`|`, filter_conditions), ]
  }
  
  return(dt_filtered)
}
# Function to filter values that start with "FBgn"
filter_FBgn <- function(df, gene_col) {
  # Use grepl to create a logical vector
  fbgn_rows <- grepl("^FB", df[[gene_col]])
  # Subset the data frame based on the logical vector
  return(df[fbgn_rows, ])
}
# Define global functions for application:

# Converts official gene symbols to FBgnIDs
# *Switches mart based on availability
conflict_prefer("select", "dplyr")
conflicts_prefer(dplyr::filter)



# Reads the BDSC Stock Info DB
get_bdsc_db <- function(){
  bdsc_db_list <- list()
  #data_url <- "https://www.dropbox.com/scl/fi/a153dk5kn4tomlnwouiuq/stockcomps_map_comments.csv?rlkey=2bh1ygp2espra4egnbla1hzny&dl=0"
  #download.file(data_url, "Data/stockcomps_map_comments.csv")
  #readLines(url(data_url))
  bdsc_db_list["stockcomps_map_comments"] <- read.csv("data/stockcomps_map_comments.csv")
  return(bdsc_db_list)
}

get_RNAi_stock_df <- function(df_list, gene_list, df_name){
  print("HERE AT LEAST")
  new_rows <- data.frame()
  new_row <- list()
  for(gene in gene_list){
    next
    new_row = list()
    search_str_list <- c(gene, "RNAi")
    colList <- c("comment1", "comment2", "comment3")
    rnai_stockcomps <- get_filtered_db_multiple_cols(bdsc_db_list$stockcomps_map_comments, search_str_list, colList)
    num_stocks_for_gene <- nrow(rnai_stockcomps)
    if(num_stocks_for_gene > 0){
      for(i in 1:num_stocks_for_gene){
        new_row <- list()
        for(df_colname in names(df_list[[df_name]])){
          new_row[[df_colname]] <- df_list[[df_name]] %>%
            filter(flybase_gene_id == gene) %>%
            pull(df_colname)
          if(is.na(new_row[[df_colname]])){
            print("MAYBE?")
            new_row[[df_colname]] <- "NA"
          }
        }

        new_row[["BDSC_RNAi_stock_ids"]] <- rnai_stockcomps$`Stk #`[i]
        new_row[["genotype_list"]] <- rnai_stockcomps$Genotype[i]
        new_row[["symbol_list"]] <- rnai_stockcomps$component_symbol[i]
        new_row[["map_list"]] <- rnai_stockcomps$mapstatement[i]
        new_row[["comment1"]] <- rnai_stockcomps$comment1[i]
        new_row[["comment2"]] <- rnai_stockcomps$comment2[i]
        new_row[["comment3"]] <- rnai_stockcomps$comment3[i]
        print(names(new_rows))
        print(names(new_row))
        new_rows <- rbind(new_rows, new_row)
      }
    } else {
      new_row <- list()
      print(gene)
      for(df_colname in names(df_list[[df_name]])){
        new_row[[df_colname]] <- df_list[[df_name]] %>%
          filter(flybase_gene_id == gene) %>%
          pull(df_colname)
        print("THIS COND????")
        print(df_colname)
        print(new_row[[df_colname]])
        if(is.na(new_row[[df_colname]])){
          new_row[[df_colname]] <- "NA"
        }
      }
      print("HERE8")
      new_row[["BDSC_RNAi_stock_ids"]] <- "NA"
      new_row[["genotype_list"]] <-"NA"
      new_row[["symbol_list"]] <- "NA"
      new_row[["map_list"]] <- "NA"
      new_row[["comment1"]] <- "NA"
      new_row[["comment2"]] <- "NA"
      new_row[["comment3"]] <- "NA"
      print(length(names(new_rows)))
      print(length(names(new_row)))
      new_rows <- rbind(new_rows, new_row)
    }
    
  }
  if(nrow(new_rows) > 0) {
    print("HERE9")
    df <- data.frame()
    df[is.na(df) | df == ""] <- "NA"
    df <- df[order(df$BDSC_RNAi_stock_ids == "NA"), ]
    print("HERE10")
    
    return(df)
  } else {
    return(df_list[[df_name]])
  }
}

# TODO Query BDSC files / implement DB querying ####
get_Transgenic_construct_info <- function(df){
  original_names <- names(df)
  # joined_df <- df %>%
  #   left_join(db_list$insertion_mapping_fb_2024_04, 
  #             by = c("symbol_list" = "insertion_symbol"), keep = FALSE) %>% 
  #   mutate(transgenic_product_id = `FBti#`)
  return(df)
}

get_References_df <- function(df, entity_id_col){
  
  ref_df <- df %>%
    left_join(db_list$entity_publication_fb_2024_04 %>% 
                select("entity_id", "FlyBase_publication_id"), 
              by = setNames("entity_id", paste0(entity_id_col)), 
              keep = FALSE, relationship = "many-to-many") %>%
    left_join(db_list$fbrf_pmid_pmcid_doi_fb_2024_04, 
              by = c("FlyBase_publication_id" = "FBrf"),
              keep = FALSE) %>%
    filter(pub_type == "paper")
  return(ref_df)
}

get_Allele_stock_df <- function(df_list, gene_list, df_name, gene_col){
  df <- df_list[[df_name]] %>%
    left_join(db_list$Alleles_2_Gene_fbal_to_fbgn_fb_2024_04 %>%
                select("GeneID", "AlleleID", "AlleleSymbol"),
              by = setNames("GeneID", "flybase_gene_id"),
              keep = FALSE)
  if(nrow(df) > 0) {
    print(nrow(df))
    df[is.na(df) | df == ""] <- "NA"
    df <- df[df$AlleleID != "NA",] %>% 
      distinct(AlleleID, .keep_all = TRUE)
    return(df)
  } else {
    return(df_list[[df_name]])
  }
}


conflict_prefer("select", "dplyr")
# BDSC Stock Info
bdsc_db_list <- read_all_csvs_as_dt_recursive("data")
# FlyBase Info
db_list <- read_all_tsvs_as_dt_recursive("data")
# The main 'search' function:
# 1. Converts the gene name to FBgnID
# 2. Searches for stocks with comment containing '<FBgnID>' and 'RNAi'
# 3. Aggregates the map, comments, geneotype and insertion info
# 4. Appends the stock info to the uploaded file(s)
get_search_objects <- function(df_list, columnsList, flybase_gene_ids, progress_updater){
  
  total_rows <- (length(df_list) *100)
  current_progress <- 0
  
  print("BDSC Table: ")
  print(nrow(bdsc_db_list$stockcomps_map_comments))
  print(nrow(bdsc_db_list[["stockcomps_map_comments"]]))
  
  df_names <- names(df_list)
  search_summary <- list()
  total_stocks <- 0
  progress_updater(current_progress/total_rows, paste0("Processing files"))
  for(df_name in df_names){
    gene_col <- columnsList[[df_name]]
    clean_col <- df_list[[df_name]][[gene_col]]
    clean_col[clean_col == ""] <- NA
    clean_col <- na.omit(clean_col)
    num_fb_gn_id <- sum(grepl("^FB*", clean_col))
    total_genes <- length(unique(clean_col))
    FBgnIDcol <- (num_fb_gn_id/total_genes) > 0.5
    print("Selected col is FBgn?")
    print(FBgnIDcol)
    # Search for FBgnID
    if(any(grepl(",", df_list[[df_name]][[gene_col]]))){
      print("GO Analysis end point detected")
      df_list[[df_name]] <- df_list[[df_name]] %>%
        separate_longer_delim(gene_col, delim = ',') %>%
        mutate(!!gene_col := trimws(.data[[gene_col]])) %>%
        distinct(.data[[gene_col]], .keep_all = TRUE)
    }
    if(!FBgnIDcol){
      print("NOT A FBGnID col!!!")
      stop("NOT FBGnID COL!")
      print("Getting values$flybase ID")
      if ("flybase_gene_id" %in% names(df_list[[df_name]])) {
        df_list[[df_name]] <- df_list[[df_name]] %>%
          select(-flybase_gene_id)
      }
      # Following logic gets fbgnid's for readily available entries in current_symbol, else searches in the synonyms col
      
      df_list[[df_name]] <- df_list[[df_name]] %>%
        mutate(temp_gene_col = tolower(.data[[gene_col]])) %>%
        left_join(flybase_gene_ids %>%
                    mutate(temp_external_gene_name = tolower(current_symbol)),  # Create a lowercase column for flybase_gene_ids
                  by = c("temp_gene_col" = "temp_external_gene_name")) %>%
        mutate(flybase_gene_id = primary_FBid) %>%
        distinct(.data[[gene_col]], .keep_all = TRUE) %>%
        arrange(is.na(primary_FBid)) %>%
        select(-primary_FBid, -temp_gene_col, -current_symbol)  # Remove the temporary column
      
      unmatched_genes <- df_list[[df_name]] %>%
        filter(is.na(flybase_gene_id))
      
      if(nrow(unmatched_genes) > 0){
        print(unmatched_genes[[gene_col]])
        print(length(unmatched_genes[[gene_col]]))
        unmatched_genes_list <- unmatched_genes[[gene_col]]
        synonym_df <- flybase_gene_ids %>%
          select(primary_FBid, symbol_synonyms = `symbol_synonym(s)`) %>%
          separate_longer_delim(c("symbol_synonyms"), delim = "|") %>%         # Split synonyms into separate rows
          mutate(symbol_synonyms = str_trim(symbol_synonyms))
        
        unmatched_genes_updated <- unmatched_genes %>%
          mutate(temp_gene_col = tolower(.data[[gene_col]])) %>%
          left_join(
            synonym_df %>%
              mutate(temp_symbol_synonyms = tolower(symbol_synonyms)),
            by = setNames("temp_symbol_synonyms", "temp_gene_col")
          ) %>%
          mutate(flybase_gene_id = ifelse(
            !is.na(primary_FBid),
            primary_FBid,
            paste0("FBgnID for ", .data[[gene_col]], " not found!")
          )) %>%
          filter(!!sym(gene_col) %in% unmatched_genes_list) %>%
          select(-primary_FBid, -temp_gene_col, -symbol_synonyms)
        
        print(unmatched_genes_updated)

        df_list[[df_name]] <- df_list[[df_name]] %>%
          filter(!is.na(flybase_gene_id)) %>%  # Matched genes from the first join
          bind_rows(unmatched_genes_updated) %>%
          distinct(flybase_gene_id, .keep_all = TRUE)   # Add updated unmatched genes
      }
      df_list[[df_name]] <- df_list[[df_name]] %>%
        mutate(
          flybase_gene_id = ifelse(
            is.na(flybase_gene_id),
            paste0("FBgnID for ", .data[[gene_col]], " not found!"),
            flybase_gene_id
          )
        )%>%
        arrange(str_detect(flybase_gene_id, "not found!")) %>%
        select(-`symbol_synonym(s)`)
      
      gene_list <- df_list[[df_name]]$flybase_gene_id
      if(any(is.na(df_list[[df_name]]$flybase_gene_id))){
        stop("FEEL SOMETHING YOU SLIMY SNAIL")
      }
      print(names(df_list[[df_name]]))
    }else{ # Selected column consists of majority FBgnID's
      gene_list <- na.omit(df_list[[df_name]][[gene_col]])
    }
    gene_list <- unique(gene_list)
    gene_list_size <- length(gene_list)
    print(paste0("Length of gene list: ", gene_list_size))
    print("HERE20")
    print("UNIQUE")
    progress_updater(current_progress/total_rows, paste0("Getting RNAi stocks for genes in ", df_name))
    print(df_name)
    search_summary[[paste0(df_name, "_RNAi_Stock")]] <- get_RNAi_stock_df(df_list, gene_list, df_name)
    current_progress = current_progress + 1
    progress_updater(current_progress/total_rows, paste0("Getting References for RNAi stocks from ", df_name))
    search_summary[[paste0(df_name, "_RNAi_Stock_Refs")]] <- get_Transgenic_construct_info(get_RNAi_stock_df(df_list, gene_list, df_name))# get_References_df(search_summary[[paste0(df_name, "_RNAi_Stock")]], "transgenic_product_id")
    current_progress = current_progress + 1
    progress_updater(current_progress/total_rows, paste0("Getting Allele stocks for genes in ", df_name))
    
    
    
    print("HERE3993")
    search_summary[[paste0(df_name, "_Allele_Stock")]] <- get_Allele_stock_df(df_list, gene_list, df_name, gene_col)
    current_progress = current_progress + 1
    progress_updater(current_progress/total_rows, paste0("Getting References for Allele stocks from ", df_name))
    search_summary[[paste0(df_name, "_Allele_Stock_Refs")]] <- get_References_df(search_summary[[paste0(df_name, "_Allele_Stock")]], "AlleleID")
    current_progress = current_progress + 1
    progress_updater(current_progress/total_rows, paste0("Getting References for Allele stocks from ", df_name))
  }
  names(search_summary) <- gsub(".csv", "", names(search_summary))
  names(search_summary) <- gsub(".xlsx", "", names(search_summary))
  names(search_summary) <- gsub(".xls", "", names(search_summary))
  return(search_summary)
}

# ui ---------------------------------------------------------------------------

# Define UI for application:

ui <- fluidPage( useShinyjs(),
                 theme = shinytheme("simplex"),
                 
                 navbarPage(
                   "BDSC Stock Search", id = "tabs",
                   tabPanel("Upload",
                            h3("Select and upload gene set files in CSV/ Excel file format"),
                            fileInput("files", "Choose Gene Set Files", multiple = TRUE, accept = c(".csv", ".xlsx", ".xls")),
                            actionButton("reset", "Reset"),
                            actionButton("confirmUpload", "Confirm files"),
                            tags$hr(),
                            tableOutput("fileTable")
                   ),
                   tabPanel("Select Columns",
                            h3("Choose official gene symbol/ FBgnID column for uploaded files:"),
                            uiOutput("columnSelectors"),
                            actionButton("confirmColumns", "Confirm columns")
                   ),
                   tabPanel("Processing",
                            h3("Search for RNAi stock IDs completed!"),
                            uiOutput("progressUI"),
                            uiOutput("nextButtonUI")
                   ),
                   tabPanel("Results",
                            h3("Search Complete!"),
                            tableOutput('resultsTable'),  # Displays the data frame info
                            downloadButton("download_all", "Download all stock id files"),
                            actionButton("restart", "Restart")
                   )
                 )
)

# server -----------------------------------------------------------------------

options(shiny.maxRequestSize = 100*1024^2)
server <- function(input, output, session) {
  values <- reactiveValues(files = NULL, results_list = NULL, df_list = NULL, cols =NULL, flybase_gene_ids = NULL)
  deleteTempFiles <- function() {
    temp_directory <- tempdir()
    files_in_temp <- list.files(temp_directory, pattern = "csv$",full.names = TRUE)
    print(files_in_temp)
    # Remove specific files safely
    unlink(files_in_temp, recursive = TRUE)
  }
  observe({
    print("Session started for user ")
    print("Test deployment")
    deleteTempFiles()
    values$flybase_gene_ids <- as.data.frame(read_tsv("data/fb_synonym_fb_2024_04.tsv"))[c("current_symbol", "primary_FBid", 'symbol_synonym(s)')]
  })
  
  # Trigger only when input$files changes
  observeEvent(input$files, {
    output$fileTable <- renderTable({
      req(input$files)
      # Create a dataframe with file details
      data <- data.frame()
      if(any(!is.null(input$files))){
        data <- data.frame(
          Index  = 1:length(input$files$name),
          Name = input$files$name,
          Size = format_size(input$files$size),
          Rows = sapply(input$files$datapath, count_rows),
          stringsAsFactors = FALSE
        )
      } else {
        data <- data.frame(
          data.frame(Name = character(), Size = numeric(), Rows = numeric())
        )
      }
      data
    }, ignoreNULL = FALSE)  # Ensure that NULL values also trigger the observer
  })
  
  observeEvent(input$reset, {
    # Reset action, clear files
    # Clear the file input on the UI
    # Update the table to show no data
    values$files <- NULL
    values$results_list <- NULL
    values$df_list = NULL
    values$cols = NULL
    
    output$fileTable <- renderTable({
      data.frame(Name = character(), Size = numeric(), Rows = numeric())
    })
    shinyjs::reset("files")
  })
  
  # Upload event
  observeEvent(input$confirmUpload, {
    req(input$files)
    values$files <- input$files
    values$results_list <- NULL
    values$cols <- NULL
    values$df_list <- NULL
    
    updateNavbarPage(session, inputId = "tabs",selected = "Select Columns")
  })
  
  # Helper function to format file sizes into human-readable form
  format_size <- function(size_in_bytes) {
    sapply(size_in_bytes, function(size) {
      if (size >= 1048576) {
        paste0(format(round(size / 1048576, 2), nsmall = 2), " MB")
      } else {
        paste0(format(round(size / 1024, 2), nsmall = 2), " KB")
      }
    })
  }
  
  # Helper function to count rows in each file
  count_rows <- function(file_path) {
    ext <- tools::file_ext(file_path)
    if (ext == "csv") {
      return(nrow(read.csv(file_path)))
    } else if (ext == "xlsx" | ext == "xls") {
      return(nrow(read_excel(file_path)))
    }
  }
  
  # Dynamic UI for selecting columns from each uploaded file
  output$columnSelectors <- renderUI({
    req(values$files)
    lapply(1:nrow(values$files), function(i) {
      file <- values$files[i, ]
      if (grepl("xlsx", file$datapath) | grepl("xls", file$datapath)){
        colnames <- names(read_excel(file$datapath))
      } else {
        colnames <- names(read.csv(file$datapath, nrows = 1))
      }
      #colnames <- colnames[!colnames %in% c("flybase_gene_id") ]
      tags$div(class = "card", style = "margin: 0px; padding: 0px; width: 50%;", 
               tags$div(class = "card-body",
                        tags$h5(class = "card-title", HTML(paste0("Select official gene symbol/ FBgnID column in file: <br><strong>", 
                                                                  file$name, "</strong>"))),
                        selectInput(paste0("col_", file$name), "", colnames)
               ),
               tags$hr()
      )
    })
  })
  
  # Confirm col selection and read files
  observeEvent(input$confirmColumns, {
    req(values$files)  # Ensure files are available
    
    # First, store column selections
    values$cols <- lapply(1:nrow(values$files), function(i) {
      file <- values$files[i, ]
      input[[paste0("col_", file$name)]]
    })
    
    # Update the navigation to move to the "Processing" tab
    
    # Read files and apply the selected columns
    
    
    # Process data and update progress here
    withProgress(message = 'Searching for BDSC stocks...\n', value = 0, {
      for (i in seq_along(values$files$name)) {
        file <- values$files[i, ]
        df <- NA
        if(grepl(".xlsx", file$datapath) | grepl(".xls", file$datapath)){
          print("HERE")
          df <- read_excel(file$datapath)
        } else if (grepl(".csv", file$datapath)){
          df <- read.csv(file$datapath, na.strings = c("", "NA"), fill = TRUE)  # Disable automatic renaming
        }
        
        progress_updater <- function(progress_fraction, detail) {
          incProgress(progress_fraction, detail = detail)
        }
        # Use file name as the list element name
        values$df_list[[file$name]] <- df  # Store the dataframe with the selected column
        values$cols[[file$name]] <- input[[paste0("col_", file$name)]]
      }
      # Call the function and pass the updater
      print(names(values$df_list))
      print(values$cols)
      values$results_list <- get_search_objects(values$df_list, values$cols, values$flybase_gene_ids, progress_updater)
    })
  })
  
  observeEvent(values$results_list, {
    # This code runs only when results_list is updated
    req(values$results_list)
    updateNavbarPage(session, inputId = "tabs", selected = "Processing")
  }, ignoreNULL = TRUE)  # Set ignoreNULL to TRUE to avoid triggering when results_list is NULL
  
  output$progressUI <- renderUI({
    withProgress(message = 'Initializing...', value = 0, {
      # This is just a placeholder. Actual progress updates happen elsewhere.
    })
  })
  
  # Display results
  output$nextButtonUI <- renderUI({
    # Check if the results list is not null
    req(values$results_list)
    if (length(values$results_list) == 4*length(values$df_list)) {
      actionButton("nextButtonUI", "Download Results =>")
    }
  })
  
  observeEvent(input$nextButtonUI, {
    # Code to execute when the Next button is clicked, such as navigating to another tab
    updateNavbarPage(session, inputId = "tabs",selected = "Results")
  })
  
  output$resultsTable <- renderTable({
    req(values$results_list)
    data <- lapply(names(values$results_list), function(name) {
      df <- values$results_list[[name]]
      nrows = nrow(df)
      Size = format(object.size(df), units = "auto")
      list(
        Index = as.integer(which(names(values$results_list) == name)),
        Name = paste0(name, ".csv"),
        Size = Size,
        Rows = as.integer(nrows)
      )
    })
    results_df <- data.frame(do.call(rbind, data))
    results_df
  }, server = TRUE)
  
  
  #_______________________________________________________________________________
  
  # Download all results
  output$download_all <- downloadHandler(
    filename = function() {
      version <- format(Sys.Date(), "%m%d%y")
      paste0("BDSC-RNAi-Stock_", version, ".zip")
    },
    content = function(file) {
      # Define a unique directory for each session based on session token
      parent_folder_name <- paste0("BDSC_stock_ids_", Sys.Date(), "_", session$token)
      temp_dir <- file.path(tempdir(), parent_folder_name)
      
      # Ensure the parent directory exists
      if (!dir.exists(temp_dir)) {
        dir.create(temp_dir, recursive = TRUE)
      }
      
      files <- lapply(names(values$results_list), function(df_name) {
        csv_file_path <- file.path(temp_dir, paste0(df_name, ".csv"))
        #values$results_list[[df_name]][values$results_list[[df_name]] == "NA"] <- NA
        write.csv(values$results_list[[df_name]], csv_file_path, row.names = FALSE)
        csv_file_path  # Return the path for each CSV
      })
      
      # Flatten the list if nested
      files <- unlist(files)
      
      # Check if all files exist
      if (all(file.exists(files))) {
        print("All files exist, proceeding to zip.")
        # Change the current directory to the temporary directory
        setwd(tempdir())
        # Create a ZIP file with all CSV files in the parent directory
        zip(zipfile = file, files = parent_folder_name, extras = "-r")
        # Reset the working directory
        setwd(normalizePath("~"))
        
        # Clean up the temporary directory
        unlink(temp_dir, recursive = TRUE)
      } else {
        print("One or more files do not exist.")
      }
    }
  )
  
  
  # Restart application
  observeEvent(input$restart, {
    values$files <- NULL
    values$results_list <- NULL
    values$df_list <- NULL
    values$cols <- NULL
    deleteTempFiles()
    updateNavbarPage(session, inputId = "tabs",selected = "Upload")
    output$fileTable <- renderTable({
      data.frame(Name = character(), Size = numeric(), Rows = numeric())
    })
    shinyjs::reset("files")
    bdsc_db_list <- read_all_csvs_as_dt_recursive("data")
  })
  
  session$onSessionEnded(function() {
    # Delete input files
    deleteTempFiles()
    stopApp()
  })
  
  
}

# run ---------------------------------------------------------------------

# Run the application:

shinyApp(ui = ui, server = server)

