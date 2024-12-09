
# created by Ange Tchuisseu 
# date of creation: 12 June 2024
# Last update : 12 June 2024

  
    pheweb_function <- function(dataframe, chr, pos,  allele1, allele2, hg="38" ,
                                pheweb_website = "https://pheweb.org/UKB-Neale/variant/", file_name="heatmap") { 
      # This function takes a dataframe with SNP and chr_pos columns, queries the specified PheWeb website for each variant, 
      # and returns significant traits associated with each variant.
      # Reference for PheWeb: https://www.nature.com/articles/s41588-020-0622-5#MOESM2
      # 
      # Arguments:
      # - dataframe: A dataframe containing at least two columns: 'SNP' and 'chr_pos'. 'chr_pos' must be in the format: chr:pos-allele1-allele2.
      # - pheweb_website: The URL of the PheWeb website. Default is "https://pheweb.org/UKB-Neale/variant/". ( 2400 traits for 11 million imputed variants in 337,000 unrelated white British individuals.)
      #                   Other options:  
      #                     - https://pheweb.org/UKB-SAIGE/variant/ (for 1400 EHR-derived broad PheWAS codes for 20 million imputed variants in 400,000 white British individuals.)
      #                     - https://pheweb.org/UKB-TOPMed/variant/  (for 1400 EHR-derived broad PheWAS codes for 57 million TOPMed-imputed variants in 400,000 white British individuals.)
      # 
      # Returns:
      # A list containing significant traits associated with each variant, variants with no significant traits, and variants not found.
      
      # Load required libraries
      library(httr)
      library(jsonlite)
      library(dplyr)
      #set time
      time1 = Sys.time()
      
      # Check if 'SNP' column exists
      rename_rsid_to_snp <- function(df) {
        # Check if 'SNP' column exists (case-sensitive)
        if (!"SNP" %in% colnames(df)) {
          
          # Check for 'rsid' column (case-insensitive)
          rsid_col <- grep("^rsid$", colnames(df), ignore.case = TRUE, value = TRUE)
          
          # If 'rsid' column is found, rename it to 'SNP'
          if (length(rsid_col) > 0) {
            colnames(df)[colnames(df) == rsid_col] <- "SNP"
            message(paste("Column", rsid_col, "renamed to 'SNP'."))
          } else {
            message("No 'SNP' or 'rsid' column found.")
          }
          
        } else {
          message("'SNP' column already exists.")
        }
        
        return(df)
      }
      dataframe <- rename_rsid_to_snp(dataframe)
      
      if (hg=="38") { message ("Using hg38 by default, please check the panel again")} 
      else  { cat ("Using", hg)
        
        }
      # Initialize dataframes to store results
      data <- data.frame()
      no_trait_row <- data.frame()
      variant_not_found <- data.frame()
      # Loop through each variant and query the website
      if (nrow(dataframe)<2) { 
        cat("dataframe empty")
        data= data.frame(SNP  = NA, variant= NA, chr= NA, 
                         pos  = NA,
                         allele1  = NA, 
                         allele2= NA,
                         beta=NA, ac= NA, category=NA,
                         pval = NA, 
                         phenostring = NA, 
                         phenocode = NA, num_samples=NA)
        allele_switched= data.frame()
        no_trait_row = data.frame()
        variant_not_found = data.frame()
        keep= data.frame() 
        return(list(data = data, allele_switched= allele_switched, no_trait_row = no_trait_row, variant_not_found = variant_not_found, keep= keep ))
        
      }
      else 
      {
      for (i in 1:nrow(dataframe))  {
        variant_not_found_= data.frame()
        no_trait_row_ <- data.frame()
        # Extract the current variant identifier
        dataframe$chr_pos <- paste(chr , pos , sep = ":")
        dataframe$chr_pos <- paste(dataframe$chr_pos, allele1 , sep = "-") 
        dataframe$chr_pos <- paste(dataframe$chr_pos, allele2 , sep = "-")      
        variant <- dataframe$chr_pos[i]
        print(i)
        print(variant)
        # Print the current variant being processed
       # cat("Processing variant:", variant, "\n")
        
        # Create the query URL for the current variant
        query_url <- paste0(pheweb_website, variant, sep="")
        
        # Send GET request to the website
        response <- GET(query_url)
        
        # Check if the request was successful
        if (status_code(response) == 200) {
          # Try to parse the JSON response
          tryCatch({
            # Extract the relevant part of the HTML content containing the JSON data
            html_content <- sub("^.*?<script", "<script", content(response, as = "text"), perl=TRUE)
            html_content2 <- sub(".*UKBiobank PheWeb</title>", "<title>UKBiobank PheWeb</title>", html_content)
            html_content3 <- sub(".* type=\"text/javascript\"></script>\n\n\n\n\n\n<script type=\"text/javascript\">\n", 
                                 " type=\"text/javascript\"></script>\n\n\n\n\n\n<script type=\"text/javascript\">\n",
                                 html_content2)
            
            html_content4 <- sub(".*  window.variant = ", 
                                 "  window.variant = ",
                                 html_content3)
            html_content5 <- sub(";\n.*", "", html_content4)
            char_value <- gsub("\t", ",", html_content5)
            
            # Extract JSON portion from the character value
            json_start <- regexpr("\\{", char_value)
            json_string <- substr(char_value, json_start, nchar(char_value))
            
            # Parse JSON content
            parsed_json <- fromJSON(json_string)
            
            # Convert to dataframe
            df <- as.data.frame(parsed_json$phenos)
            
            # Filter the dataframe for pval < 0.00005 as given in opentargets meaning that -log (0.00005) = 4.3
            result <- subset(df, df$pval < 0.00005)
            print(colnames(df))
            # If there are significant traits, append them to the results dataframe
            if (nrow(result) > 0) {
              print("TRAIT FOUND ")
              result$variant <- variant
              result$SNP <- dataframe$SNP[i]
              result$chr= dataframe$chr[i]
              result$pos  = dataframe$pos[i]
             # result$allele1  = dataframe$allele1[i]
             # result$allele2  = dataframe$allele2[i]
              data <- rbind(data, result)
            } else {
              print("NO TRAIT FOUND ")
              
              # If no significant traits found, append a row with SNP set to variant and other columns as NA
              no_trait_row_ <- data.frame(SNP = dataframe$SNP[i], variant= variant, chr= chr[i], 
                                          pos  = pos[i], 
                                          allele1  = allele1[i], 
                                          allele2  = allele2[i],
                                          ac    =NA,        af=NA,           beta   =NA,      category=NA,
                                          num_cases =NA,   num_controls =NA,phenocode  =NA, 
                                          phenostring =NA, pval  =NA,       sebeta =NA,
                                          trait_is_bad =NA,tstat =NA,       url_sum=NA )
              
              no_trait_row <- rbind(no_trait_row, no_trait_row_)
            }
          }, error = function(e) {
            message("Error parsing JSON for variant ", variant, ": ", e)
          })
        } else {
          print("VARIANT NOT FOUND ")
          
          variant_not_found_ <- data.frame(SNP = dataframe$SNP[i], variant= variant, chr= chr[i], 
                                           pos  = pos[i], 
                                           allele1  = allele1[i], 
                                           allele2  = allele2[i], category= NA, 
                                           ac    =NA,        af=NA,           beta   =NA,     
                                           num_cases =NA,   num_controls =NA,phenocode  =NA, 
                                           phenostring =NA, pval  =NA,       sebeta =NA,
                                           trait_is_bad =NA,tstat =NA,       url_sum=NA       
                                         
                                           )
          variant_not_found <- rbind(variant_not_found, variant_not_found_)
        }
      }
      cat("Sorry, no  significant trait was found for variant :",  "\n", no_trait_row$SNP, "\n", "in the Pheweb manhattan plot. Check the website:", pheweb_website, "\n")
      
      message("         ", "\n", "\n", "\n")
      cat("Request failed for variant : Sorry, I couldn't find the variant:", "\n",  variant_not_found$SNP,"\n",  "in the Pheweb. Check the website:", pheweb_website, "\n")
      
      message("         ", "\n", "\n", "\n")
      cat(" There are ", nrow(data), "  phenotypes found")
      message("         ", "\n", "\n", "\n")
      
      data2= (variant_not_found)
      message ("switching the alleles for the variants that were not found")
      cat(" There are ", nrow(data2), "  variants to scan through after changing allele position. ")
      
      
      
      storage=data.frame()
      no_trait_row<- data.frame(SNP = NA, variant= NA, chr= NA, 
                                pos  = NA, 
                                
                                beta   =NA, ac=NA ,  category=NA,  
                                num_samples =NA,   phenocode  =NA, 
                                phenostring =NA, pval  =NA) 
      
      variant_not_found <- data.frame()
      keep= data.frame()
     print(data2)
      #repeat this for the same code by changing the alleles
      for (i in 1:nrow(data2)) {
        variant_not_found_= data.frame()
        no_trait_row_ <- data.frame()
        # retrieve the information from the dataframe 
        # Extract the current variant identifier
        data2$chr_pos <- paste(data2$chr , data2$pos , sep = ":")
        data2$chr_pos <- paste(data2$chr_pos, data2$allele2 , sep = "-") 
        data2$chr_pos <- paste(data2$chr_pos, data2$allele1 , sep = "-")      
        
        
        variant <- data2$chr_pos[i]
        print(i)
        # Print the current variant being processed
        # cat("Processing variant:", variant, "\n")
        
        # Create the query URL for the current variant
        query_url <- paste0(pheweb_website, variant, sep="")
        
        # Send GET request to the website
        response <- GET(query_url)
        
        # Check if the request was successful
        if (status_code(response) == 200) {
          # Try to parse the JSON response
          tryCatch({
            # Extract the relevant part of the HTML content containing the JSON data
            html_content <- sub("^.*?<script", "<script", content(response, as = "text"), perl=TRUE)
            html_content2 <- sub(".*UKBiobank PheWeb</title>", "<title>UKBiobank PheWeb</title>", html_content)
            html_content3 <- sub(".* type=\"text/javascript\"></script>\n\n\n\n\n\n<script type=\"text/javascript\">\n", 
                                 " type=\"text/javascript\"></script>\n\n\n\n\n\n<script type=\"text/javascript\">\n",
                                 html_content2)
            
            html_content4 <- sub(".*  window.variant = ", 
                                 "  window.variant = ",
                                 html_content3)
            html_content5 <- sub(";\n.*", "", html_content4)
            char_value <- gsub("\t", ",", html_content5)
            
            # Extract JSON portion from the character value
            json_start <- regexpr("\\{", char_value)
            json_string <- substr(char_value, json_start, nchar(char_value))
            
            # Parse JSON content
            parsed_json <- fromJSON(json_string)
            
            # Convert to dataframe
            df <- as.data.frame(parsed_json$phenos)
            
            # Filter the dataframe for pval < 5e-8
            result <- subset(df, df$pval < 5e-8)
           
            # If there are significant traits, append them to the results dataframe
            if (nrow(result) > 0) {
              result$variant <- variant
              result$SNP <- data2$SNP[i]
              result$chr= data2$chr[i]
              result$pos  = data2$pos[i]
              #result$allele1  = data2$allele1[i]
             # result$allele2  = data2$allele2[i]
              cat("finished")
              storage <- rbind(storage, result)
              
              # # Remove rows in df2 where SNPs are found in df1
              # variant_not_found <- variant_not_found %>% 
              #   filter(!SNP %in% storage$SNP)
              
              
             # print(nrow(data2))
            } else {
              # If no significant traits found, append a row with SNP set to variant and other columns as NA
              no_trait_row_ <- data.frame(SNP = dataframe$SNP[i], variant= variant, chr= chr[i], 
                                          pos  = pos[i], 
                                          category=NA,  
                                                    beta   =NA, ac=NA   ,  
                                          num_samples =NA,   phenocode  =NA, 
                                          phenostring =NA, pval  =NA
                                    
                                          )
              no_trait_row <- rbind(no_trait_row, no_trait_row_)
              print(no_trait_row)
              # variant_not_found <- variant_not_found %>% 
              #   filter(!SNP %in% no_trait_row$SNP)
            }
          }, error = function(e) {
            message("Error parsing JSON for variant ", variant, ": ", e)
            keep= rbind(keep, variant)
          })
        } else {
          message("Request failed for variant - not found ", variant)
          variant_not_found_ <- data.frame(SNP = dataframe$SNP[i], variant= variant, chr= chr[i], 
                                           pos  = pos[i],category=NA,  
                                           
                                           beta   =NA, ac=NA ,    
                                           num_samples =NA,   phenocode  =NA, 
                                           phenostring =NA, pval  =NA
                                           
          )
          
          variant_not_found <- rbind(variant_not_found, variant_not_found_)
        }
      }
      message("                 Final Results  ", "\n")
      cat("Sorry, no  significant trait was found for variant :",  "\n", no_trait_row$SNP, "\n", "in the Pheweb manhattan plot. Check the website:", pheweb_website, "\n")
      
      message("         ", "\n", "\n", "\n")
      cat("Request failed for variant : Sorry, I couldn't find the variant:", "\n",  variant_not_found$SNP,"\n",  "in the Pheweb. Check the website:", pheweb_website, "\n")
      
      # Print the final results
      time2=Sys.time()
      duration= time2-time1
      message("This took ", duration , "  minutes to run.",  "\n") 
      message("                  END  ", "\n") 
      #no_trait_row= no_trait_row%>%dplyr::select(-c(allele1, allele2))
      
      #plot the graph
      
      library(ggplot2)
      library(reshape2)
    
      
      print((data))
      if (pheweb_website %in% c("https://pheweb.org/UKB-SAIGE/variant/", "https://pheweb.org/UKB-TOPMed/variant/")) {
        library(dplyr)
       # print(colnames(storage))
        
        no_trait_row <- no_trait_row %>%
          mutate(num_cases = NA,
                 num_controls = NA,
                 category = NA,
                 sebeta = NA,
                 trait_is_bad = NA,
                 tstat = NA,
                 af = NA,
                 url_sum = NA)
        if ("num_samples" %in% colnames(no_trait_row)) {
          no_trait_row$num_samples <- NULL
        }
        cat ("data merged") 
        print((data))

        print((storage))
        
        print((no_trait_row))
        
        # print(colnames(no_trait_row))
        data= bind_rows(data, storage, no_trait_row)
        
       
        dat= data%>% dplyr::select(phenostring, SNP, pval, beta, category)
      # Reshape data for heatmap
      dat$log_pval <- -log10(dat$pval)
      
      print(dat)
      # Filter out rows where 'phenostring' and other relevant columns are NA
      # Filter data
      dat_filtered <- dat %>%
        group_by(SNP) %>%
        filter(
          # Keep rows where phenostring is not NA
          any(!is.na(phenostring)) |
            # Keep one row if all values are NA
            (all(is.na(phenostring) & is.na(pval) & is.na(beta) & is.na(category) & is.na(log_pval)) & row_number() == 1)
        ) %>%
        ungroup() %>%
        distinct()
      # Identify SNPs that have at least one non-NA value
      non_na_snps <- dat_filtered$SNP[!is.na(dat_filtered$log_pval)]
      
      # Remove rows where log_pval is NA for SNPs that have other non-NA rows
      df_filtered <- dat_filtered[!(is.na(dat_filtered$log_pval) & dat_filtered$SNP %in% non_na_snps), ]
      message(" Filter")
      # Removing rows where any column contains 'diabetes' in its name
      # Removing rows where 'method' contains 'diabetes' or 'Diabetes'
      df_filtered <- df_filtered[!grepl("diabetes|Diabetes", df_filtered$phenostring, ignore.case = TRUE), ]
      
      df_filtered$phenostring <- ifelse(is.na(df_filtered$phenostring), "No phenotype", df_filtered$phenostring)
      df_filtered$category <- ifelse(is.na(df_filtered$category), "No phenotype", df_filtered$category)
      df_filtered$log_pval <- ifelse(is.na(df_filtered$log_pval), 0, df_filtered$log_pval)
      df_filtered$pval <- ifelse(is.na(df_filtered$pval), 0, df_filtered$pval)
      
      heatmap_data <- dcast(df_filtered, SNP ~ phenostring, value.var = "log_pval")
      
      # Melt the data for ggplot
      heatmap_data_melted <- melt(heatmap_data, id.vars = "SNP")
      message(" category_phenostring")
      # Extract unique phenostring for each category
      category_phenostring <- dat %>%
        group_by(category) %>%
        filter(phenostring != "No phenotype") %>%
        select(category, phenostring) %>%
        distinct() %>%
        dplyr::rename(variable=phenostring)
      # Identify which rows belong to "No phenotype"
      #heatmap_data_melted$no_phenotype <- ifelse(heatmap_data_melted$variable == "No phenotype", TRUE, FALSE)
      
      # Melt data for ggplot2
      # <- melt(heatmap_data, id.vars = "SNP")
      
      # Plot
      library(ggplot2)
      # Plot with distinct colors for different p-value ranges and a cross for 'No phenotype' only if it’s the only phenotype
      library(ggplot2)
      message(" only_no_phenotype_snps ")
      # Identify SNPs that only have "No phenotype"
       only_no_phenotype_snps <- aggregate(value ~ SNP, data = heatmap_data_melted, function(x) all(is.na(x) | heatmap_data_melted$variable == "No phenotype"))
       only_no_phenotype_snps <- only_no_phenotype_snps$SNP[only_no_phenotype_snps$value == TRUE]
      # 
      
      # Get unique categories for x-axis text color
    
      library(RColorBrewer)
      # Create a named vector for category colors
      # Get unique categories for x-axis text color
      categories <- unique(dat$category)
      print(length(categories))
      # Create a color palette for categories
      # Create a color palette with enough colors for all categories
      color_palette <- colorRampPalette(brewer.pal(12, "Set3"))(length(categories))
      
      # Create a named vector for category colors
      category_colors <- setNames(color_palette, categories)
      
      # category_colors <- brewer.pal(n = length(categories), name = "Set3") # You can choose other palettes if needed
      # 
      # Create a named vector for category colors
      category_colors <- setNames(category_colors, categories)
      
      # Ensure that every unique category gets a color
      # Ensure that every unique category gets a color
      
      # Merge category information with heatmap_data_melted
      # Merge category information with heatmap_data_melted using left_join
      heatmap_data_melted <- left_join(heatmap_data_melted, category_phenostring, by = "variable")
      
      heatmap_data_melted <- heatmap_data_melted %>%
        mutate(category_color = ifelse(!is.na(category), category, "No phenotype")) %>%
        mutate(category_color = factor(category_color, levels = categories))
      
      message(" Add a flag for cells with No phenotype ")
      heatmap_data_melted$no_phenotype <- heatmap_data_melted$variable == "No phenotype"
      heatmap_data_melted$only_no_phenotype <- heatmap_data_melted$SNP %in% only_no_phenotype_snps
      
      # Determine the order of phenostring based on category
      heatmap_data_melted$variable <- factor(heatmap_data_melted$variable, 
                                             levels = unique(heatmap_data_melted$variable[order(heatmap_data_melted$category)]))
      
      # Plot with distinct colors for different p-value ranges and a cross for 'No phenotype' only if it’s the only phenotype
      library(ggplot2)
      library(ggplot2)
     p= ggplot(heatmap_data_melted, aes(y = variable, x = SNP, fill = ifelse(no_phenotype, "No phenotype", cut(value, breaks = c(0, 5, 10, 15, 20))))) +
        geom_tile(color = "white") +  # Ensure NAs are displayed as empty cells
        
        scale_fill_manual(values = c("lightblue", "skyblue", "blue", "darkblue", "white"), 
                          name = "-log10(p-value)", 
                          na.value = "white",  # Set NA cells to be white
                          labels = c("<5", "5-10", "10-15", ">15", "No phenotype")) +
        
        # Add a cross "X" for cells with "No phenotype" but only if it’s the only phenotype for that SNP
        # geom_text(data = heatmap_data_melted[heatmap_data_melted$no_phenotype & heatmap_data_melted$only_no_phenotype, ], 
        #           aes(label = "X"), color = "black") +
       geom_text(aes(label = ifelse(no_phenotype & only_no_phenotype, "X", NA)), color = "black") +
       # geom_text(data = unique(heatmap_data_melted %>% filter(category != "No phenotype")), 
       #           aes(label = variable, color = category_color), size = 3, vjust = -0.5) +
       # scale_color_manual(values = category_colors) +
       
        
        theme_minimal() +
        labs(x = "SNP", y = "Disease ") +
       theme(
         axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Increase size of x-axis text
         axis.text.y = element_text(size = 12#, color = as.factor(heatmap_data_melted$category_color)
                                    )# Increase size of y-axis text
       ) +
     coord_fixed(ratio = 1)  # Ensure cells are square
     
      print(p)
     ggsave(paste0("heatmap_", file_name), width = 12, height = 12)
      # Return the list of dataframes
      return(list(final= dat_filtered, data = data, allele_switched= storage,
                  no_trait_row = no_trait_row, variant_not_found = variant_not_found, keep= keep ))
      } 
      
      else {   
        data= rbind(data, storage, no_trait_row)
      dat= data%>% dplyr::select(phenostring, SNP, pval, beta)
      # Return the list of dataframes
        return(list( data = data, allele_switched= storage,
                    no_trait_row = no_trait_row, variant_not_found = variant_not_found, keep= keep ))}
      }
    }
    
    create_heatmap <- function(df_subset,  x_label, y_label) {
      
      # Subset the dataset to only include variables in var_test
      #df_subset <- df_clean2[, var_test]
      
      # Step 1: Compute the correlation matrix and p-values for the subset
      cor_results <- rcorr(as.matrix(df_subset))  # Hmisc package function
      cor_matrix <- cor_results$r
      p_matrix <- cor_results$P
      
      # Step 2: Melt the correlation and p-value matrices to long format for ggplot
      cor_melt <- melt(cor_matrix)
      p_melt <- melt(p_matrix)
      
      # Rename the columns after melting to reflect the actual variable names
      colnames(cor_melt) <- c("Variable1", "Variable2", "Correlation")
      colnames(p_melt) <- c("Variable1", "Variable2", "p_value")
      
      # Combine the correlation and p-value data into one data frame
      heatmap_data <- merge(cor_melt, p_melt, by = c("Variable1", "Variable2"))
      
      # Filter heatmap_data to only include rows where Variable2 starts with "prs"
      heatmap_data <- heatmap_data %>%
        filter(grepl("^(prs|PRS)", Variable1)) %>%
        filter(!grepl("^(prs|PRS)", Variable2))
      # Add a new column for significance labels
      # Add a new column for significance labels
      heatmap_data$significance <- cut(heatmap_data$p_value, 
                                       breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                                       labels = c("***", "**", "*", ""), 
                                       right = FALSE)
      
      # Create a combined label of correlation and significance for display
      heatmap_data$label <- ifelse(heatmap_data$p_value < 0.05, 
                                   paste0(round(heatmap_data$Correlation, 2), 
                                          " ", heatmap_data$significance), 
                                   "X")
      heatmap_data$Correlation <- ifelse(heatmap_data$p_value < 0.05, 
                                         heatmap_data$Correlation, 
                                         
                                         NA)
      
      # Step 3: Plot the heatmap with enhanced color differentiation and circles
      p = ggplot(heatmap_data, aes(Variable1, Variable2, fill = Correlation), size=5) +
        geom_tile(color = "white") +
        scale_fill_gradient2(low = "red", high = "blue", mid = "yellow", 
                             midpoint = 0.05, limit = c(-1, 1), 
                             name = "Correlation", na.value = "white")+
        geom_text(aes(label = label), color = "black", size = 3) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),  # Increase x-axis text size
              axis.text.y = element_text(size = 9))+
        labs(title = "Correlation Matrix ",
             x = paste0(x_label),
             y = paste0(y_label))
      print(p)
    }
    

