library(DT)
library(shinydashboard)
library(shinycssloaders)

options(shiny.sanitize.errors = FALSE)
options(shiny.maxRequestSize=1000*1024^2) 

dashboardPage(skin = "black",
              # start of Header
              dashboardHeader(
                title = "HT-PCR Tools"
              ),
              # start of Sidebar
              dashboardSidebar(
                sidebarMenu(
                  menuItem('1. Introduction', tabName = 'intro', icon = icon('info-circle')
                  ),
                  menuItem('2. Upload data', tabName = 'load', icon = icon('cloud-upload-alt')
                  ),
                  menuItem('3. Data processing', tabName = 'dp', icon = icon('hourglass-start')
                  ),
                  menuItem('4. Plotting results', tabName = 'plot', icon = icon('chart-area')
                  )
                )
              ),
              # start of Body
              dashboardBody(
                tabItems(
                  tabItem(tabName = 'intro',
                    fluidRow(
                      box(title = "User Guide", width = 12, solidHeader = T, status = "primary",
                        column(12,
                          includeMarkdown("intro.Rmd")
                        )
                      )
                    )
                  ),
                  tabItem(tabName = 'load',
                      tabPanel('Upload you data',
                        fluidRow(
                          column(3,
                            box(title = 'Load in data',solidHeader = T, 
                                status = "primary", collapsible = T,width = 12,
                                h4(strong('Select fastq.gz/fq.gz')),
                                hr(),
                                radioButtons('data_input', 'Use example file or upload your own data',
                                  c('Upload local file'='local_file',
                                    'Example Data'='example_file'),
                                  selected = 'local_file'
                                ),
                                conditionalPanel(condition="input.data_input=='local_file'",
                                  radioButtons('file_type', 'Single or Pair-end sample',
                                               c('Single end'='SE',
                                                 'Pair end'='PE'),
                                               selected = 'PE'),
                                  p('* Please wait for the upload complete!'),
                                  fileInput('file', '', accept = c('text/gz',
                                    'text/comma-separated-values,text/plain', 
                                    '.gz'), multiple = TRUE
                                  ),
                                  actionButton("goButton", "Upload data"),
                                  tags$li("(This will take a few minutes.)")
                                ),
                                conditionalPanel(condition="input.data_input=='example_file'"
                                )
                            )
                          ),
                          column(9,
                            box(title = 'Raw sequences',solidHeader = T, 
                              status = "primary", collapsible = T,width = 12,
                              uiOutput('sampleID'),
                              withSpinner(dataTableOutput("Seq"))
                            )
                          )
                        )
                      )
                  ),
                  tabItem(tabName = 'dp',
                    tabsetPanel(
                      tabPanel(title = '3.1 Total variant',
                        fluidRow(
                          column(3,
                            box(title = 'Extract target reads', solidHeader = T, 
                              status = "primary", collapsible = T,width = 12,
                              textInput("Seq5","input 5' Seq",
                                value='GTAACTAAGA'),
                              textInput("Seq3","input 3' Seq",
                                value='TACTTTATCT'),
                              actionButton("anButton", "Extract !"),
                              tags$li("(This will take a few minutes.)")
                            )
                          ),
                          column(9, 
                            box(title = 'Target sequences',solidHeader = T, 
                              status = "primary", collapsible = T,width = 12,
                              uiOutput('Flt_sampleID'),
                              withSpinner(dataTableOutput("Flt_seq"))
                            )
                          )
                        )
                      ),
                      tabPanel(title = '3.2 Split varivattions',
                        fluidRow(
                          column(3,
                            box(title = 'Input WT sequences',solidHeader = T, 
                                status = "primary", collapsible = T,width = 12,
                              textInput('WT', 'Please enter WT sequence',
                                value = 'TCTGATTCCTGAAGACATTCCTCTGGAGCCTCAGTAAATT'),
                              actionButton("SpButton", "Start to Split"),
                              hr(),
                              uiOutput('Var_sampleID')
                            )
                          ),
                          
                          column(9,
                            box(title = 'WT sequences',solidHeader = T, 
                                status = "primary", collapsible = T,width = 12,
                                dataTableOutput('WT')
                                # withSpinner(dataTableOutput('WT'))
                            ),
                            box(title = 'SNP sequences',solidHeader = T, 
                                status = "primary", collapsible = T,width = 12,
                                withSpinner(dataTableOutput('SNP'))
                            ),
                            box(title = 'Indel sequences',solidHeader = T, 
                                status = "primary", collapsible = T,width = 12,
                                withSpinner(dataTableOutput('Indel'))
                            )
                          )
                        )
                      )
                    )
                  ),
                  tabItem(tabName = 'plot',
                    tabsetPanel(
                      tabPanel(title = 'Barplot',
                        fluidRow(
                          column(12,
                            box(title = 'Total variations table',solidHeader = T, 
                              status = "primary", collapsible = T,width = 6,
                              withSpinner(dataTableOutput('Total_t')),
                              hr(),
                              downloadButton('downloadTotal_t','Download .csv', class = "btn btn-warning")
                            ),
                            box(title = 'Indel variations table',solidHeader = T, 
                                status = "primary", collapsible = T,width = 6,
                                withSpinner(dataTableOutput('Indel_t')),
                                hr(),
                                downloadButton('downloadIndel_t','Download .csv', class = "btn btn-warning")
                            )
                          ),
                          column(12,
                            box(title = 'Total variations plot',solidHeader = T, 
                                status = "primary", collapsible = T,width = 6,
                                column(4, selectInput('tp', 'Chose value to plot:',
                                            choices = c(
                                              'Number'='Number',
                                              'Proportion' = 'Proportion'
                                            ), selected = 'Proportion')
                                ),
                                column(4, numericInput('T_height', 'Figure Width:', min = 1, max = 20, value = 6)),
                                column(4, numericInput('T_width', 'Figure Height:', min = 1, max = 20, value = 4)),
                                hr(),
                                withSpinner(plotOutput('Total_p')),
                                hr(),
                                textInput('Total_n', 'Output file name:','Total_variations_plot'),
                                downloadButton('downloadTotalPng','Download .png', class = "btn btn-warning"),
                                downloadButton('downloadTotalPdf','Download .pdf', class = "btn btn-warning")
                            ),
                            box(title = 'Indel variations plot',solidHeader = T, 
                                status = "primary", collapsible = T,width = 6,
                                column(4, selectInput('Ip', 'Chose value to plot:',
                                            choices = c(
                                              'Number'='Number',
                                              'Proportion' = 'Proportion'
                                            ), selected = 'Proportion')
                                ),
                                column(4, numericInput('I_height', 'Figure Width:', min = 1, max = 20, value = 6)),
                                column(4, numericInput('I_width', 'Figure Height:', min = 1, max = 20, value = 4)),
                                hr(),
                                withSpinner(plotOutput('Indel_p')),
                                hr(),
                                textInput('Indel_n', 'Output file name:','Indel_variations_plot'),
                                downloadButton('downloadIndelPng','Download .png', class = "btn btn-warning"),
                                downloadButton('downloadIndelPdf','Download .pdf', class = "btn btn-warning")
                            )
                          )
                        )
                      ),
                      tabPanel(title = 'Seqlogo Plot',
                        fluidRow(
                          column(12,
                                 box(title = 'SNP variations table',solidHeader = T, 
                                     status = "primary", collapsible = T,width = 12,
                                     uiOutput('SNP_sampleID'),
                                     withSpinner(dataTableOutput('SNP_t')),
                                     textInput('SNP_n', 'Output file name:','SNP_variations_plot'),
                                     downloadButton('downloadSnp_t','Download .csv', class = "btn btn-warning")
                                 ),
                                 box(title = 'SNP variations plot',solidHeader = T, 
                                     status = "primary", collapsible = T,width = 12,
                                     column(4, selectInput('Sp', 'Chose method to plot:',
                                                           choices = c(
                                                             'bits'='bits',
                                                             'prob' = 'prob'
                                                           ), selected = 'prob')
                                     ),
                                     column(4, numericInput('S_height', 'Figure Width:', min = 1, max = 30, value = 15)),
                                     column(4, numericInput('S_width', 'Figure Height:', min = 1, max = 10, value = 3)),
                                     column(12, withSpinner(plotOutput('SNP_p'))),
                                     downloadButton('downloadSnpPng','Download .png', class = "btn btn-warning"),
                                     downloadButton('downloadSnpPdf','Download .pdf', class = "btn btn-warning")
                                 )
                          )
                        )
                      )
                     )
                    )
                  )
                )      # end of body
)
