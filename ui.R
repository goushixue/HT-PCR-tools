library(shiny)
library(shinydashboard)
library(shinycssloaders)

options(shiny.maxRequestSize=300*1024^2) 

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
                                  p('* Please wait for the upload complete!'),
                                  fileInput('file', '', accept = c('text/gz',
                                    'text/comma-separated-values,text/plain', 
                                    '.gz')
                                  ),
                                  actionButton("goButton", "Upload data"),
                                  tags$li("(This will take a few minutes.)")
                                ),
                                conditionalPanel(condition="input.data_input=='example_file'",
                                  fileInput('file', '', accept = c('text/gz',
                                    'text/comma-separated-values,text/plain', 
                                    '.gz')
                                  )
                                )
                                # fileInput('file', 'Choose input File:'),
                                # actionButton("goButton", "Upload data"),
                                # tags$li("(This will take a few minutes.)")
                            )
                          ),
                          column(9,
                            box(title = 'Raw sequences',solidHeader = T, 
                              status = "primary", collapsible = T,width = 12,
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
                                value='CCAGCAGCTC'),
                              textInput("Seq3","input 3' Seq",
                                value='TCCAGTGTCA'),
                              actionButton("anButton", "Extract !"),
                              tags$li("(This will take a few minutes.)")
                            )
                          ),
                          column(9, 
                            box(title = 'Target sequences',solidHeader = T, 
                              status = "primary", collapsible = T,width = 12,
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
                                value = 'GGGAGCCCAGGTGGGCGGATCCATCTCCTCTGGCTCCTCCGCC'),
                              actionButton("SpButton", "Start to Split")
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
                      tabPanel(title = 'Plot figure',
                        fluidRow(
                          column(6,
                            box(title = 'Total variations table',solidHeader = T, 
                              status = "primary", collapsible = T,width = 12,
                              withSpinner(dataTableOutput('Total_t'))
                            ),
                            box(title = 'Total variations plot',solidHeader = T, 
                                status = "primary", collapsible = T,width = 12,
                                withSpinner(plotOutput('Total_p')),
                                hr(),
                                textInput('Total_n', 'Output file name:','Total_variations_plot'),
                                downloadButton('downloadTotalPng','Download .png', class = "btn btn-warning"),
                                downloadButton('downloadTotalPdf','Download .pdf', class = "btn btn-warning")
                            )
                          ),
                          column(6,
                            box(title = 'Indel variations table',solidHeader = T, 
                              status = "primary", collapsible = T,width = 12,
                              withSpinner(dataTableOutput('Indel_t'))
                              ),
                            box(title = 'Indel variations plot',solidHeader = T, 
                                status = "primary", collapsible = T,width = 12,
                                withSpinner(plotOutput('Indel_p')),
                                hr(),
                                textInput('SNP_n', 'Output file name:','SNP_variations_plot'),
                                downloadButton('downloadIndelPng','Download .png', class = "btn btn-warning"),
                                downloadButton('downloadIndelPdf','Download .pdf', class = "btn btn-warning")
                            )
                          ),
                          column(12,
                            box(title = 'SNP variations plot',solidHeader = T, 
                              status = "primary", collapsible = T,width = 12,
                              withSpinner(plotOutput('SNP_p')),
                              hr(),
                              textInput('Indel_n', 'Output file name:','Indel_variations_plot'),
                              downloadButton('downloadSnpPng','Download .png', class = "btn btn-warning"),
                              downloadButton('downloadSnpPdf','Download .pdf', class = "btn btn-warning")
                            )
                          )
                        )
                      )
                    )
                  )
                )
               # end of body
)
