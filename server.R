library(shiny)
library(shinydashboard)
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
    reads <- readFastq(inFile$datapath) # read reads
    incProgress(0.1, detail = paste("Reading file"))
    sequences <- sread(reads) # extract reads
    incProgress(0.1, detail = paste("Extractting reads"))
    Sys.sleep(0.1)
    })
    sequences
  })
  
# Load example data
  example_data <- eventReactive(input$data_input,{
    withProgress(message='processing',min=0,max=1,value=0,{
      reads <- readFastq('example_data/example.fq.gz')
      incProgress(0.1, detail = paste("Reading file"))
      sequences <- sread(reads) # extract reads
      incProgress(0.1, detail = paste("Extractting reads"))
      Sys.sleep(0.1)
    })
    sequences
  })
  
  
# show some sequences  
  output$Seq <- renderDataTable({
    if (input$data_input == 'local_file'){
      data <- reads()
    }
    else {
      data <- example_data()
    }
    reads <- as.data.frame(data[1:1000])
    colnames(reads) <- 'Sequences'
    reads
  }, options = list(pageLength = 15))
  
# Select wantted sequences  
  flt_seq <- eventReactive(input$anButton,{
    if (input$data_input == 'local_file') {
      data <- reads()
    }
    else {
      data <- example_data()
    }
    withProgress(message='processing',min=0,max=1,value=0,{
    gseq <- paste0(input$Seq5, '.*', input$Seq3)
    idx <- grepl(gseq, data)
    fseq <- data[idx]
    incProgress(0.1, detail = paste("Finding target reads"))
    font_seq <- gsub(paste0('.*', input$Seq5), replacement = input$Seq5, fseq)
    incProgress(0.1, detail = paste("Delete 5' reads"))
    beh_seq <- gsub(paste0(input$Seq3, '.*'), replacement = input$Seq3, font_seq)
    incProgress(0.1, detail = paste("Delete 3' reads"))
    sub_seq <- substr(beh_seq, nchar(input$Seq5) + 1, nchar(beh_seq) - nchar(input$Seq3))
    incProgress(0.1, detail = paste("Extractting target reads"))
    Sys.sleep(0.1)
    })
    sub_seq
  })
  
# Show filted sequences
  output$Flt_seq <- renderDataTable({
    reads <- as.data.frame(flt_seq()[1:1000])
    colnames(reads) <- 'Sequences'
    reads
  }, options = list(pageLength = 15))
  
# Split WT SNP and Indel
  WT <- eventReactive(input$SpButton,{
    flt_seq()[flt_seq() == input$WT]
  })
  
  output$WT <- renderDataTable({
    reads <- as.data.frame(WT()[1:100])
    colnames(reads) <- 'WT'
    reads
  }, options = list(pageLength = 10))
  
  SNP <- eventReactive(input$SpButton,{
    flt_seq()[flt_seq() != input$WT & nchar(flt_seq()) == nchar(input$WT)]
  })
  
  output$SNP <- renderDataTable({
    reads <- as.data.frame(SNP()[1:100])
    colnames(reads) <- 'SNP'
    reads
  }, options = list(pageLength = 10))
  
  Indel <- eventReactive(input$SpButton,{
    flt_seq()[nchar(flt_seq()) != nchar(input$WT)]
  })
  
  output$Indel <- renderDataTable({
    reads <- as.data.frame(Indel()[1:100])
    colnames(reads) <- 'Indel'
    reads
  }, options = list(pageLength = 10))
  
# Plot total var
  Total <- reactive({
    Total <- data.frame(WT=length(WT()), SNP=length(SNP()), Indel=length(Indel())) %>% t %>% as.data.frame
    colnames(Total) <- 'Numbers'
    Total$Freq <- Total$Numbers / sum(Total$Numbers)
    Total$Type <- rownames(Total)
    Total
  })
  
  output$Total_t <- renderDataTable({
    Total()
  })
  
  output$Total_p <- renderPlot({
    ggplot(Total(), aes(x=Type, y=Freq, fill=Type))+
      geom_bar(stat = "identity")+
      theme_bw()
  })
# Dowload total plot
  output$downloadTotalPng <- downloadHandler(
    filename = function()  {paste0(input$Total_n,".png")},
    content = function(file) {
      png(file)
      p <- ggplot(Total(), aes(x=Type, y=Freq, fill=Type))+
              geom_bar(stat = "identity")+
              theme_bw()
      print(p)
      dev.off()
    }
  )
  
  output$downloadTotalPdf <- downloadHandler(
    filename = function() {
      paste0(input$Total_n,".pdf")
    },
    content = function(file) {
      pdf(file)
      p <- ggplot(Total(), aes(x=Type, y=Freq, fill=Type))+
              geom_bar(stat = "identity")+
              theme_bw()
      print(p)
      dev.off()
    }
  )
  
# plot SNP var
  SNP_u <- reactive({
    SNP_uniq <- unique(SNP())
  })
  
  output$SNP_p <- renderPlot({
    ggplot()+geom_logo(SNP_u(), method="prob", seq_type="dna", col_scheme = 'nucleotide')+theme_logo()
  })

# Download SNP plot  
  output$downloadSnpPng <- downloadHandler(
    filename = function()  {paste0(input$SNP_n,".png")},
    content = function(file) {
      png(file)
      p <- ggplot(Total(), aes(x=Type, y=Freq, fill=Type))+
        geom_bar(stat = "identity")+
        theme_bw()
      print(p)
      dev.off()
    }
  )

  output$downloadSnpPdf <- downloadHandler(
    filename = function() {
      paste0(input$SNP_n,".pdf")
    },
    content = function(file) {
      pdf(file)
      p <- ggplot(Total(), aes(x=Type, y=Freq, fill=Type))+
        geom_bar(stat = "identity")+
        theme_bw()
      print(p)
      dev.off()
    }
  )
  
# Plot Indel var
  Indels <- reactive({
    Indel_uniq <- unique(Indel())
    Indel_type <- nchar(Indel_uniq) %% 3 %>% table %>% data.frame
    colnames(Indel_type) <- c('Type', 'Numbers')
    Indel_type$Type <- c('3N', '3N+1', '3N+2')
    Indel_type$Freq <- Indel_type$Numbers / sum(Indel_type$Numbers)
    Indel_type
  })
  
  output$Indel_t <- renderDataTable({
    Indels()
  })
  
  output$Indel_p <- renderPlot({
    ggplot(Indels(), aes(x=Type, y=Freq, fill=Type))+
      geom_bar(stat = "identity")+
      theme_bw()
  })

# Download indel plot  
  output$downloadIndelPng <- downloadHandler(
    filename = function()  {paste0(input$Indel_n,".png")},
    content = function(file) {
      png(file)
      p <- ggplot(Total(), aes(x=Type, y=Freq, fill=Type))+
        geom_bar(stat = "identity")+
        theme_bw()
      print(p)
      dev.off()
    }
  )
  
  output$downloadIndelPdf <- downloadHandler(
    filename = function() {
      paste0(input$Indel_n,".pdf")
    },
    content = function(file) {
      pdf(file)
      p <- ggplot(Total(), aes(x=Type, y=Freq, fill=Type))+
        geom_bar(stat = "identity")+
        theme_bw()
      print(p)
      dev.off()
    }
  )

})
