library(shiny)
library(reshape2)
library(stringi)
library(ggplot2)
library(magrittr)
library(ShortRead)
library(ggseqlogo)

# Define server logic required to plot various variables against mpg
shinyServer(function(input, output){
  
# read in fastq data
  reads <- eventReactive(input$goButton,{
    withProgress(message='processing',min=0,max=1,value=0,{
    inFile <- input$file
    if (is.null(inFile))
      return(NULL)
    if (input$file_type == 'PE') {
    fileNames <- sapply(inFile$name, function(x){
      stringi::stri_split(x,regex = '.fq')[[1]]
    })[1,] %>% unique
    
    sequences <- lapply(inFile$datapath, function(x){
      reads <- readFastq(x) # read reads
      sequences <- sread(reads) # extract reads
    })
    names(sequences) <- fileNames
    
    fileName_u <- sapply(fileNames, function(x){
        stringi::stri_split(x,regex = '_')[[1]]
      })[1,] %>% unique

    sequences_ap <- lapply(fileName_u, function(x){
      append(sequences[[paste0(x,'_1')]], sequences[[paste0(x,'_2')]] %>% reverse %>% complement)
    })
      names(sequences_ap) <- fileName_u
    }
    else {
      fileNames <- sapply(inFile$name, function(x){
        stringi::stri_split(x,regex = '.fq')[[1]]
      })[1,] %>% unique
      
      sequences_ap <- lapply(inFile$datapath, function(x){
        reads <- readFastq(x) # read reads
        sequences <- sread(reads) # extract reads
      })
      names(sequences_ap) <- fileNames
    }
    
    incProgress(0.1, detail = paste("Reading file"))
    Sys.sleep(0.1)
    })
    sequences_ap
  })
  
# Load example data
  example_data <- eventReactive(input$data_input,{
    withProgress(message='processing',min=0,max=1,value=0,{
      fls <- dir("example_data", "*fq*", full=TRUE)
      fileNames <- lapply(fls, function(x){
        rev_split <- stringi::stri_split(x,regex = '/')[[1]]  %>% rev %>% stringi::stri_split(regex = '.fq')
        fileNames <- rev_split[[1]][1]
      })
      sequences <- lapply(fls, function(x){
        reads <- readFastq(x) # read reads
        sequences <- sread(reads) # extract reads
      })
      names(sequences) <- fileNames
      
      fileName_u <- sapply(fileNames, function(x){
        stringi::stri_split(x,regex = '_')[[1]]
      })[1,] %>% unique
      
      sequences_ap <- lapply(fileName_u, function(x){
        append(sequences[[paste0(x,'_1')]], sequences[[paste0(x,'_2')]] %>% reverse %>% complement)
      })
      
      names(sequences_ap) <- fileName_u
      
      incProgress(0.1, detail = paste("Reading file"))
      Sys.sleep(0.1)
    })
    sequences_ap
  })
  
# RenderUI for sample select.
  output$sampleID <- renderUI({
    if (input$data_input == 'local_file') {
      inFile <- input$file
      if (is.null(inFile))
        return(NULL)
      if (input$file_type == 'PE') {
        fileNames <- sapply(inFile$name, function(x){
          stringi::stri_split(x,regex = '.fq')[[1]]
        })[1,] %>% unique
        
        fileName_u <- sapply(fileNames, function(x){
          stringi::stri_split(x,regex = '_')[[1]]
        })[1,] %>% unique
      }
      else {
        fileName_u <- sapply(inFile$name, function(x){
          stringi::stri_split(x,regex = '.fq')[[1]]
        })[1,] %>% unique
      }
    }
    else {
      fls <- dir("example_data", "*fq*", full=TRUE)
      fileNames <- lapply(fls, function(x){
        rev_split <- stringi::stri_split(x,regex = '/')[[1]]  %>% rev %>% stringi::stri_split(regex = '.fq')
        fileNames <- rev_split[[1]][1]
      })
      fileName_u <- sapply(fileNames, function(x){
        stringi::stri_split(x,regex = '_')[[1]]
      })[1,] %>% unique
    }
    names(fileNames) = fileNames
    selectInput('sample_id', 'Which sample to show:', fileName_u)
  })

  
# show some sequences  
  output$Seq <- renderDataTable({
    if (input$sample_id %>% is.null) 
      return(NULL)
    if (input$data_input == 'local_file'){
      data <- reads()[[input$sample_id]]
    }
    else {
      data <- example_data()[[input$sample_id]]
    }
    reads <- data[1:1000] %>% as.data.frame
    colnames(reads) <- 'Sequences'
    reads
  },options = list(pageLength = 10, autoWidth = TRUE, scrollX=TRUE))
  
# Select wantted sequences  
  flt_seq <- eventReactive(input$anButton,{
    if (input$data_input == 'local_file') {
      data <- reads()
    }
    else {
      data <- example_data()
    }
    withProgress(message='processing',min=0,max=1,value=0,{
    
    fseq <- lapply(data, function(x){
      idx <- grepl(paste0(input$Seq5, '.*', input$Seq3), x)
      x[idx]
    })
    incProgress(0.1, detail = paste("Finding target reads"))
    # cut the sequences
    sub_seq <- lapply(fseq, function(x){
      font_seq <- gsub(paste0('.*', input$Seq5), replacement = input$Seq5, x) # cut 5' start
      beh_seq <- gsub(paste0(input$Seq3, '.*'), replacement = input$Seq3, font_seq) # cut 3' end
      sub_seq <- substr(beh_seq, nchar(input$Seq5) + 1, nchar(beh_seq) - nchar(input$Seq3)) # remove the sart and end base pair
    })
    incProgress(0.1, detail = paste("Extractting target reads"))
    Sys.sleep(0.1)
    })
    sub_seq
  })
  
  # RenderUI for sample select.
  output$Flt_sampleID <- renderUI({
    if (input$data_input == 'local_file') {
      inFile <- input$file
      if (is.null(inFile))
        return(NULL)
      if (input$file_type == 'PE') {
        fileNames <- sapply(inFile$name, function(x){
          stringi::stri_split(x,regex = '.fq')[[1]]
        })[1,] %>% unique
        
        fileName_u <- sapply(fileNames, function(x){
          stringi::stri_split(x,regex = '_')[[1]]
        })[1,] %>% unique
      }
      else {
        fileName_u <- sapply(inFile$name, function(x){
          stringi::stri_split(x,regex = '.fq')[[1]]
        })[1,] %>% unique
      }
    }
    else {
      fls <- dir("example_data", "*fq*", full=TRUE)
      fileNames <- lapply(fls, function(x){
        rev_split <- stringi::stri_split(x,regex = '/')[[1]]  %>% rev %>% stringi::stri_split(regex = '.fq')
        fileNames <- rev_split[[1]][1]
      })
      fileName_u <- sapply(fileNames, function(x){
        stringi::stri_split(x,regex = '_')[[1]]
      })[1,] %>% unique
    }
    names(fileNames) = fileNames
    selectInput('Flt_sample_id', 'Which sample to show:', fileName_u)
  })
  
# Show filted sequences
  output$Flt_seq <- renderDataTable({
    reads <- flt_seq()[[input$Flt_sample_id]][1:1000] %>% as.data.frame
    colnames(reads) <- 'Sequences'
    reads
  },options = list(pageLength = 15, autoWidth = TRUE))
  
# Split WT SNP and Indel
  # RenderUI for sample select.
  output$Var_sampleID <- renderUI({
    if (input$data_input == 'local_file') {
      inFile <- input$file
      if (is.null(inFile))
        return(NULL)
      if (input$file_type == 'PE') {
        fileNames <- sapply(inFile$name, function(x){
          stringi::stri_split(x,regex = '.fq')[[1]]
        })[1,] %>% unique
        
        fileName_u <- sapply(fileNames, function(x){
          stringi::stri_split(x,regex = '_')[[1]]
        })[1,] %>% unique
      }
      else {
        fileName_u <- sapply(inFile$name, function(x){
          stringi::stri_split(x,regex = '.fq')[[1]]
        })[1,] %>% unique
      }
    }
    else {
      fls <- dir("example_data", "*fq*", full=TRUE)
      fileNames <- lapply(fls, function(x){
        rev_split <- stringi::stri_split(x,regex = '/')[[1]]  %>% rev %>% stringi::stri_split(regex = '.fq')
        fileNames <- rev_split[[1]][1]
      })
      fileName_u <- sapply(fileNames, function(x){
        stringi::stri_split(x,regex = '_')[[1]]
      })[1,] %>% unique
    }
    names(fileNames) = fileNames
    selectInput('Var_sample_id', 'Which sample to show:', fileName_u)
  })
  
  WT <- eventReactive(input$SpButton,{
    WT <- lapply(flt_seq(), function(x){
      x[x == input$WT]
    })
  })
  
  output$WT <- renderDataTable({
    reads <- WT()[[input$Var_sample_id]][1:100] %>% as.data.frame
    colnames(reads) <- 'WT'
    reads
  },options = list(pageLength = 10, autoWidth = TRUE))
  
  SNP <- eventReactive(input$SpButton,{
    SNP <- lapply(flt_seq(), function(x){
      x[x != input$WT & nchar(x) == nchar(input$WT)]
    })
  })
  
  output$SNP <- renderDataTable({
    reads <- SNP()[[input$Var_sample_id]][1:100] %>% as.data.frame
    colnames(reads) <- 'SNP'
    reads
  },options = list(pageLength = 10, autoWidth = TRUE))
  
  Indel <- eventReactive(input$SpButton,{
    Indel <- lapply(flt_seq(), function(x){
      x[nchar(x) != nchar(input$WT)]
    })
  })
  
  output$Indel <- renderDataTable({
    reads <- Indel()[[input$Var_sample_id]][1:100] %>% as.data.frame
    colnames(reads) <- 'Indel'
    reads
  },options = list(pageLength = 10, autoWidth = TRUE))
  
# Plot total var
  Total <- reactive({
    WT_table <- sapply(WT(), length) %>% data.frame
    colnames(WT_table) <- 'WT'
    SNP_table <- sapply(SNP(), length) %>% data.frame
    colnames(SNP_table) <- 'SNP'
    Indel_table <- sapply(Indel(), length) %>% data.frame
    colnames(Indel_table) <- 'Indel'
    Total <- cbind(WT_table, SNP_table, Indel_table)
  })
  
  Total_t <- reactive({
    Total <- Total()
    for (i in Total %>% colnames) {
      New_n <- paste(i, 'Proportion', sep = ' ')
      Total[New_n] <- Total[i] / Total %>% rowSums
    }
    Total
  })
  
  output$Total_t <- renderDataTable({
    Total_t()
  },options = list(pageLength = 10, autoWidth = TRUE))
  
  # Download total table
  output$downloadTotal_t <- downloadHandler(
    filename = function() {
      paste0(input$Total_n,".csv")
    },
    content = function(file) {
      write.csv(Total_t(), file, row.names=T)
    }
  )
  
  Total_p <- reactive({
    if (dim(Total())[1] != 1) {
      if (input$tp == 'Number') {
        Total_n <- data.frame(Total(), samples = rownames(Total()))
        Total_n <- reshape2::melt(Total_n, id.vars='samples')
        ggplot(Total_n, aes(x=samples, y=value, fill=variable))+
          geom_bar(stat = "identity", width = 0.6)+
          labs(x = 'Samples', y = 'Numbers')+
          theme_bw()+
          theme(aspect.ratio = input$T_width / input$T_height)
      }
      else{
        Total_r <- data.frame(Total() / rowSums(Total()), samples = rownames(Total()))
        Total_r <- reshape2::melt(Total_r, id.vars='samples')
        ggplot(Total_r, aes(x=samples, y=value, fill=variable))+
          geom_bar(stat = "identity", width = 0.6)+
          labs(x = 'Samples', y = 'Proportion')+
          theme_bw()+
          theme(aspect.ratio = input$T_width / input$T_height)
      }
    }
    else {
      if (input$tp == 'Number') {
        Total_n <- Total() %>% reshape2::melt()
        Total_n$variable <- Total_n$variable %>% as.character
        ggplot(Total_n, aes(x=variable, y=value, fill=variable))+
          geom_bar(stat = "identity", width = 0.6)+
          labs(x = 'Samples', y = 'Numbers')+
          theme_bw()+
          theme(aspect.ratio = input$T_width / input$T_height)
      }
      else{
        Total_r <- (Total() / rowSums(Total())) %>% reshape2::melt()
        Total_r$variable <- Total_r$variable %>% as.character
        ggplot(Total_r, aes(x=variable, y=value, fill=variable))+
          geom_bar(stat = "identity", width = 0.6)+
          labs(x = 'Samples', y = 'Proportion')+
          theme_bw()+
          theme(aspect.ratio = input$T_width / input$T_height)
      }
    }
  })
  
  output$Total_p <- renderPlot({
    Total_p()
  })
  
# Dowload total plot
  output$downloadTotalPng <- downloadHandler(
    filename = function()  {paste0(input$Total_n,".png")},
    content = function(file) {
      p <- Total_p()
      ggsave(file, p, width = input$T_width, height = input$T_height)
    }
  )
  
  output$downloadTotalPdf <- downloadHandler(
    filename = function() {
      paste0(input$Total_n,".pdf")
    },
    content = function(file) {
      p <- Total_p()
      ggsave(file, p, width = input$T_width, height = input$T_height)
    }
  )
  
# SNP table
  output$SNP_sampleID <- renderUI({
    if (input$data_input == 'local_file') {
      inFile <- input$file
      if (is.null(inFile))
        return(NULL)
      if (input$file_type == 'PE') {
        fileNames <- sapply(inFile$name, function(x){
          stringi::stri_split(x,regex = '.fq')[[1]]
        })[1,] %>% unique
        
        fileName_u <- sapply(fileNames, function(x){
          stringi::stri_split(x,regex = '_')[[1]]
        })[1,] %>% unique
      }
      else {
        fileName_u <- sapply(inFile$name, function(x){
          stringi::stri_split(x,regex = '.fq')[[1]]
        })[1,] %>% unique
      }
    }
    else {
      fls <- dir("example_data", "*fq*", full=TRUE)
      fileNames <- lapply(fls, function(x){
        rev_split <- stringi::stri_split(x,regex = '/')[[1]]  %>% rev %>% stringi::stri_split(regex = '.fq')
        fileNames <- rev_split[[1]][1]
      })
      fileName_u <- sapply(fileNames, function(x){
        stringi::stri_split(x,regex = '_')[[1]]
      })[1,] %>% unique
    }
    names(fileNames) = fileNames
    selectInput('SNP_sample_id', 'Which sample to show:', fileName_u)
  })
  
  SNP_t <- reactive({
    SNP_uniq <- lapply(SNP(), unique)
    
    split_table <- lapply(SNP_uniq, function(x){
      len <- x %>% nchar %>% unique
      
      split_table <- sapply(x, function(y){
        substring(y, 1:len, 1:len)
      }) %>% t %>% as.data.frame
    })
    
    SNP_t <- lapply(split_table, function(x){
      len <- dim(split_table[[1]])[2]
      table_counts <- apply(x, 2, table %>% as.data.frame)
      table_rb <- do.call(rbind, table_counts)
      table_rb['col_n'] <- rownames(table_rb) %>% stringi::stri_replace(replacement = '', regex = '\\..*')
      tb_count <- reshape2::dcast(table_rb,  value.Var1 ~ col_n, value.var = 'value.Freq', fill = 0) %>% as.data.frame
      tb_count <- data.frame(tb_count[,-1], row.names = tb_count[,1])
      colnames(tb_count) <- tb_count %>% colnames %>% stringi::stri_replace(replacement = '', regex = 'V')
      snp_t <- tb_count[,seq(1, len) %>% as.character]
    })
    SNP_t
  })
  
  output$SNP_t <- renderDataTable({
    if (input$SNP_sample_id %>% is.null) 
      return(NULL)
    
    SNP_t()[[input$SNP_sample_id]]
  }, options = list(autoWidth = TRUE, scrollX=TRUE))
  
  output$downloadSnp_t <- downloadHandler(
    filename = function() {
      paste0(input$SNP_n,".csv")
    },
    content = function(file) {
      write.csv(SNP_t()[[input$SNP_sample_id]], file, row.names=T)
    }
  )
  
# plot SNP var
  SNP_p <- reactive({
    SNP_uniq <- lapply(SNP(), unique)
    if (input$Sp == 'prob') {
      ggseqlog_plot <- lapply(SNP_uniq, function(x){
        ggplot()+
          geom_logo(x, method=input$Sp, seq_type="dna", col_scheme = 'nucleotide')+
          theme_logo()+
          theme(aspect.ratio = input$S_width / input$S_height)
      })
    }
    else {
      ggseqlog_plot <- lapply(SNP_uniq, function(x){
        ggplot()+
          geom_logo(x, method=input$Sp, seq_type="dna", col_scheme = 'nucleotide')+
          theme_logo()+
          theme(aspect.ratio = input$S_width / input$S_height)
      })
    }
    
  })
  
  output$SNP_p <- renderPlot({
    if (input$SNP_sample_id %>% is.null) 
      return(NULL)
    
    SNP_p()[[input$SNP_sample_id]]
  })

# Download SNP plot  
  output$downloadSnpPng <- downloadHandler(
    filename = function()  {paste0(input$SNP_n,".png")},
    content = function(file) {
      p <- SNP_p()[[input$SNP_sample_id]]
      ggsave(file, p, width = input$S_width, height = input$S_height)
    }
  )

  output$downloadSnpPdf <- downloadHandler(
    filename = function() {
      paste0(input$SNP_n,".pdf")
    },
    content = function(file) {
      p <- SNP_p()[[input$SNP_sample_id]]
      ggsave(file, p, width = input$S_width, height = input$S_height)
    }
  )
  
# Plot Indel var
  Indels <- reactive({
    Indel_uniq <- lapply(Indel(), unique)
    Indel <- lapply(Indel_uniq, function(x){
      Indel_type <- nchar(x) %% 3 %>% table %>% data.frame
      Indel_type <- Indel_type$Freq %>% t %>% as.data.frame
      colnames(Indel_type) <- c('3N', '3N+1', '3N+2')
      Indel_type %>% t
    })
    Indel <- data.frame(Indel) %>% t %>% as.data.frame
  })
  
  Indel_t <- reactive({
    Indels <- Indels()
    for (i in Indels %>% colnames) {
      New_n <- paste(i, 'Proportion', sep = ' ')
      Indels[New_n] <- Indels[i] / Indels %>% rowSums
    }
    Indels
  })
  
  output$Indel_t <- renderDataTable({
    Indel_t()
  },options = list(pageLength = 10, autoWidth = TRUE))
  
  output$downloadIndel_t <- downloadHandler(
    filename = function() {
      paste0(input$Indel_n,".csv")
    },
    content = function(file) {
      write.csv(Indel_t(), file, row.names=T)
    }
  )
  
  Indel_p <- reactive({
    if (dim(Indels())[1] != 1) {
      if (input$Ip == 'Number') {
        Indel_n <- data.frame(Indels(), samples = rownames(Indels()))
        colnames(Indel_n) <- c('3N', '3N+1', '3N+2', 'samples')
        Indel_n <- reshape2::melt(Indel_n, id.vars='samples')
        ggplot(Indel_n, aes(x=samples, y=value, fill=variable))+
          geom_bar(stat = "identity", width = 0.6)+
          labs(x = 'Samples', y = 'Numbers')+
          theme_bw()+
          theme(aspect.ratio = input$I_width / input$I_height)
      }
      else{
        Indel_r <- data.frame(Indels() / rowSums(Indels()), samples = rownames(Indels()))
        colnames(Indel_r) <- c('3N', '3N+1', '3N+2', 'samples')
        Indel_r <- reshape2::melt(Indel_r, id.vars='samples')
        ggplot(Indel_r, aes(x=samples, y=value, fill=variable))+
          geom_bar(stat = "identity", width = 0.6)+
          labs(x = 'Samples', y = 'Proportion')+
          theme_bw()+
          theme(aspect.ratio = input$I_width / input$I_height)
      }
    }
    else {
      if (input$Ip == 'Number') {
        Indel_n <- Indels() %>% reshape2::melt()
        ggplot(Indel_n, aes(x=variable, y=value, fill=variable))+
          geom_bar(stat = "identity", width = 0.6)+
          labs(x = 'Samples', y = 'Numbers')+
          theme_bw()+
          theme(aspect.ratio = input$I_width / input$I_height)
      }
      else{
        Indel_r <- (Indels() / rowSums(Indels())) %>% reshape2::melt()
        ggplot(Indel_r, aes(x=variable, y=value, fill=variable))+
          geom_bar(stat = "identity", width = 0.6)+
          labs(x = 'Samples', y = 'Proportion')+
          theme_bw()+
          theme(aspect.ratio = input$I_width / input$I_height)
      }
    }
  })
  
  
  output$Indel_p <- renderPlot({
    Indel_p()
  })

# Download indel plot  
  output$downloadIndelPng <- downloadHandler(
    filename = function()  {paste0(input$Indel_n,".png")},
    content = function(file) {
      p <- Indel_p()
      ggsave(file, p, width = input$I_width, height = input$I_height)
    }
  )
  
  output$downloadIndelPdf <- downloadHandler(
    filename = function() {
      paste0(input$Indel_n,".pdf")
    },
    content = function(file) {
      p <- Indel_p()
      ggsave(file, p, width = input$I_width, height = input$I_height)
    }
  )

})
