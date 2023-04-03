### Shiny : TCGA-GSEA website
## Arata Hayashi
## 15 June 2022



library(shiny)
library(shinyFiles)
library(shinycssloaders)
library(shinycustomloader)
library(shinythemes)
library(shinydashboard)
library(rintrojs)
library(DT)
library(shinyWidgets)
library(shinyhelper)
library(shinyalert)
library(shinyjs)



library(msigdbr)
library(clusterProfiler)
library(ggplot2)
#library(enrichplot)
library(fgsea)
#need for enrichment plot
library(RColorBrewer)
library(gridExtra)
library(cowplot)
#need to download exel file
library(writexl)
#library(dplyr)
library(reshape2)
library(stringr)

#define icon for intermediate result
okIcon <- '<span class="glyphicon glyphicon-check"></span>'
noIcon <- '<span class="glyphicon glyphicon-unchecked"></span>'

# Define code for opening a new window
js_code <- "
shinyjs.browseURL = function(url) {
  window.open(url,'_blank');
}
"

jscode <- "shinyjs.toTop = function() {document.body.scrollTop = 0;}"

# load pics from open source
geni_png <- "https://www.linkpicture.com/q/geni_logo_layers_2023_1_1.png"
lab_pic2_png  <- "https://www.linkpicture.com/q/lab_pic2.png"

# load tables
mytable <- as.data.frame(readRDS("dat/studies.RDS"))
rownames(mytable) <- unlist(mytable$name)
genesetList <- readRDS("dat/gene.sets.choise.RDS")
conversion <- readRDS("dat/conversion_updated.RDS")
help_text <- readRDS("dat/helptext.RDS")


ui  <- fluidPage(
  
  titlePanel(tags$head(tags$link(rel = "icon", type = "image/png", href = "https://www.linkpicture.com/q/GENI_only_GENI.png"),
                       tags$title("GENI - geni - gene enrichment identifier"))
                       
  ),
  
  theme = shinytheme("flatly"), 
  collapsible = TRUE, 
  color = "#000000",
  
  
  useShinyjs(),
  introjsUI(),
  
  extendShinyjs(text = js_code, functions = 'browseURL'),
  extendShinyjs(text = jscode, functions = 'toTop'),
  
  
  navbarPage(
    title = actionLink("reload", div(img(src = "https://www.linkpicture.com/q/GENI_only_GENI.png",
                                         width="45", height="80",
                                         filetype = "image/png",
                                         style = "margin-top:20px;
                                         padding-right: -5px;
                                         padding-bottom: 30px"),
                                     "GENI",
                                     style="margin-top: -35px;
                                          font-size: 38px;")),
    
    id = "chosentab",
    
    tabPanel( # first top tab for TCGA GSEA analysis ----
              
              title = "Start GENI", 
              
              icon = span(class = "glyphicon glyphicon-search"),
              
              
              
              sidebarLayout( ## side bar layout for TCGA
                
                
                sidebarPanel( # Sidebar panel for inputs 
                  
                  width = 3, 
                  
                  br(),
                  
                  
                  
                  div(textInput( ### gene name insert ----
                                 inputId = "GeneID",
                                 label = "Find your gene",
                                 width = '100%'
                  ),style="font-size:14px;text-align: center;font-family:Arial,sans-serif;") %>%
                    helper(type = "inline",
                           icon = "question-circle",
                           title = "",
                           fade = TRUE,
                           content = c("<br><font size='+1'>Enter the name of your gene of interest or NCBI gene ID.",
                                       "<font color=\"#000000\">Example: EXT1(ext1), DPYSL2(dpysl2)",

                                       "<br>If not recognized, please look for another gene symbol.",
                                       "<font color=\"#000000\">Example: P53 = <b>TP53</b>, HER2 = <b>ERBB2<b/></font>"),
                           size = "m"),
                  div(
                    htmlOutput("check_gene_id"),
                    align = "center"
                  ),
                  
                  br(),
                  br(),
                  
                  div(strong("Search tissue and Select your study"), font = "Arial", style = "font-size:14px;text-align: center;font-family:Arial,sans-serif;"),
                  
                  selectInput( ### gene set name select ----
                               "tissue_select",
                               "",
                               width = '100%',
                               c("All","Adrenal Gland","Bile Duct","Bladder","Bowel",
                                 "Breast","Cervix","CNS/Brain","Esophagus","Eye",
                                 "Head and Neck","Kidney", "Liver", "Lung","Lymphoid",
                                 "Myeloid","Ovary", "Pancreas","Pleura","Prostate",
                                 "Skin", "Soft tissue", "Stomach","Testis","Thymus",
                                 "Thyroid","Uterus"),
                               selected = "All"
                  ) %>% 
                    helper(type = "inline",
                           icon = "question-circle",
                           title = "",
                           fade = TRUE,
                           content = c("<br><font size='+1'>Here you choose one of the <b>TCGA studies</b> for GSEA.",
                                       "You can look for studies according to the tissue type.",
                                       "<br>If you want to see all the studies,",
                                       "please select <b>All</b>, then click <b>Show studies</b>.",
                                       "<br><font size='+3'>Selected column will be colored dark blue.</font>"),
                           size = "m"),
                  
                  div(DTOutput("chosen_study"),style = "font-size:10px; "),
                  
                  
                  br(),
                  htmlOutput("check_study"),
                  br(),
                  div(strong("Select gene set"), font = "Arial", style = "font-size:14px;text-align: center;font-family:Arial,sans-serif;"),
                  selectInput( ### gene set name select ----
                                   inputId = "geneset",
                                   "",
                                   width = '100%',
                                   genesetList$name,
                                   selected = "hallmark gene sets (H)"
                  ) %>%  
                    helper(type = "inline",
                           icon = "question-circle",
                           title = "",
                           fade = TRUE,
                           content = c("<br><font size='+1'>Here you choose one of <b>gene sets</b> for GSEA.",
                                       "<br>For more details about the gene set, please check <a href = https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp target = _blank ><b>MSigDB web site</b></a>."),
                           size = "m"),
                  br(),
                  
                  
                  # advance option
                  
                  div(dropdownButton(
                    div( 
                      tags$h3("Advance settings"),
                      
                      numericInput(
                      "nPerm",
                      "Change permutations",
                      width = '100%',
                      min = 1000,
                      max = 100000,
                      value = 10000,
                      step = 100
                    ),align = "center"),
                    
                    
                    div(numericInput(
                      "minGSSize",
                      "Change minimum number of genes in set",
                      width = '80%',
                      min = 1,
                      max = 50,
                      value = 15,
                      step = 1
                    ),align = "center"),
                    
                    
                    div(numericInput(
                      "maxGSSize",
                      "Change maximum number of genes in set",
                      width = '80%',
                      min = 100,
                      max = 1000,
                      value = 500,
                      step = 100
                    ),align = "center"),
                    
                    
                    div(selectInput( ### add option for Normalization
                      "normalization",
                      "Select normalization option",
                      width = '80%',
                      c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                      selected = "BH"
                    ),align = "center"),
                    
                    div(selectInput( ### add option for Normalization
                      "exponent",
                      "Select exponent option",
                      width = '80%',
                      c("classic", "weighted", "weighted_p2", "weighted_p3"),
                      selected = "weighted"
                    ),align = "center"),
                    
                    
                    div(numericInput(
                      "pcutoff",
                      "Change p value cut off",
                      width = '80%',
                      min = 0,
                      max = 1,
                      value = 1,
                      step = 0.01
                    ),align = "center"),
                    
                    div(actionButton("reset_advance",
                                     "Reset"),align = "right"),
                    
                    
                    
                    
                    
                    inputId = "mydropdown",
                    label = "Advanced settings",
                    icon = icon("sliders-h", verify_fa = FALSE),
                    circle = FALSE, 
                    status = "danger", 
                    size = "sm",
                    up = TRUE,
                    #icon = icon("cog"), 
                    width = "300px",
                    tooltip = FALSE
                    
                  ),align = "center") %>% 
                    helper(type = "inline",
                           icon = "question-circle",
                           title = "",
                           fade = TRUE,
                           content = c("<br><font size='+1'>Here you can change details for GSEA.",
                                       "Default:",
                                       "<img src = 'advance_default.png'>"),
                           size = "m"),
                  
                  br(),
                  
                  div(actionButton( ### action button for GSEA ----
                                    "make_intermediatesheet", 
                                    "Apply GENI",
                                    width = "80%"
                  ),align = "center"),
                  
                  
                ),# end of sidebar panels for input
                
                mainPanel(# Main panel for  outputs ----
                          
                          width = 9,
                          
                          tabsetPanel( # tab set panel for main panel
                            type = "hidden",
                            id = "GSEAanalyses",
                            
                            tabPanelBody( ## first main tab ----
                                          "first_tab",
                                          
                                         
                                          
                                          column(9,
                                                 div(img(src = geni_png, height = 200), align = "center"),
                                                 
                                                 br(),
                                                 h1("Welcome to",tags$b("GENI", style = 'color:#3C5271'),align = "center"),
                                                 h4("~ Gene ENrichment Identifier ~",align = "center"),
                                                 br(),
                                                 h3("An open web tool for easy analysis of Gene Set Enrichment Analysis of clinical data.", align = "center"),
                                                 br(),
                                                 br(),
                                                 h2("What can GENI do for you?", align = "center"),
                                                 column(8, offset = 2,
                                                        br(),
                                                        h4(tags$ul(
                                                   tags$li("Find the gene sets that are enriched for your gene of interest from cancer patient data (TCGA)."),
                                                   br(),
                                                   tags$li("Apply Gene Set Enrichment Analysis (GSEA) easily on your rank file."),
                                                   br(),
                                                   tags$li("Find the gene sets that are enriched for your gene of interest from cancer cell line data (CCLE)",tags$b("(Upcoming)",style = "color: #93DEE5;")),
                                                   
                                                  br(),
                                                  br(),
                                                 ))),
                                                 
                                                 br(),
                                                 div(actionButton(inputId = "intro", label = "Introduction tour with example", icon = icon("info-circle",verify_fa = FALSE)),align = "center"),
                                                 
                                                 
                                                 ),
                                          column(3,
                                                 br(),
                                                 br(),
                                                 h4("Updates:",align = "center"),
                                                 #h4(paste("Data modified on: ",as.Date(file.info("dat/studies.RDS")$mtime),sep = "")), # I added data modified date
                                                 h5("2022.06.15 GENI beta version uploaded"),
                                                 h5("2022.09.18 GENI beta2 version uploaded"),
                                                 h5("2022.10.19 GENI v1.0 uploaded"),
                                                 h5("2022.10.22 Bug in select study fixed"),
                                                 h5("2022.11.29 Search by NCBI Gene ID has been added"),
                                                 h5("2022.12.10 Introduction with example has been added"),
                                                 h5("2022.12.18 GA Cookies are disabled"),
                                                 h5("2023.01.23 Bug in download fixed"),
                                                 h5("2023.02.02 Gene search method is modified"),
                                                 
                                                 br(),
                                                 br(),
                                                 br(),
                                                 br(),
                                                 br(),
                                                 h4("Info:",align = "center"),
                                                 h4("Arata Hayashi,",align = "center",style = "color: #000000;"),
                                                 h4("Department of Biochemistry and Molecular Biology, Faculty of Medicine",align = "center"),
                                                 h4("The Hebrew University-Hadassah Medical School",align = "center"),
                                                 h4("Ein Karem, 9112001, Jerusalem, Israel.",align = "center"),
                                                 div(a(href = "mailto:arata.hayashi@mail.huji.ac.il","arata.hayashi@mail.huji.ac.il",target="_blank",style = "color: #000000;"),align = "center"),
                                                 br(),
                                                 
                                                 br(),
                                                 div(img(src = "https://www.linkpicture.com/q/lab_pic.png",  width="70%",height="70%"), align = "center"),
                                                 div(img(src = "https://www.linkpicture.com/q/uni_pic.png",  width="80%",height="80%"), align = "center"),
                                                 br(),
                                                 br(),
                                                 style = "background-color:powderblue;"),
                                          
                                         
                                          
                                          
                            ),
                            
                            
                            tabPanelBody( ## 2nd main tab ----
                                          "2nd_tab",
                                          
                                          h1("Choose the Gene Set to get details",align = "center"),
                                          
                                          br(),
                                          column(width = 12, style = "overflow-x: scroll;", ## allow to scroll inside column
                                                 
                                                 ### output: GSEA summary result ----
                                                 withLoader(DT::dataTableOutput('intermediate_result'),type = "html",loader = "dnaspin") %>% 
                                                   helper(type = "inline",
                                                          icon = "question-circle",
                                                          title = "",
                                                          fade = TRUE,
                                                          content = c("<br><font size='+3'>Here you can see summary table of GSEA.",
                                                                      "<br>To see details, click the row.",
                                                                      "<br>To download, select right column."),
                                                          size = "l")
                                          ),
                                          br(),
                                          br(),
                                          
                                          
                                          column(
                                            12,
                                            plotOutput("summary_dotplot_up") %>% 
                                              helper(id = "summary_dotplot_up_helper",
                                                     type = "inline",
                                                     icon = "question-circle",
                                                     title = "",
                                                     fade = TRUE,
                                                     content = c("<br><font size='+3'>You can hide/show summary plots by clicking the button below."),
                                                     size = "l")
                                          ),
                                          column(
                                            12,
                                            plotOutput("summary_dotplot_down")
                                          ),
                                          
                                          column(#separation of main panel to 2 columns ->1
                                            8,
                                            column(
                                              5,
                                              div(
                                              actionButton( ### action button for GSEA ----
                                                            "show", 
                                                            "Show/Hide summary plots"
                                              ),
                                              br(),
                                              br(),
                                              actionBttn( ### action button for gene set details
                                                
                                                "detail",
                                                HTML("Click for gene set details<br/>(jump to MSigDB page)"),
                                                
                                                
                                              ),
                                              
                                              htmlOutput("link"),
                                              br(),
                                              br(),
                                              tableOutput("details_table"), ### output: small table result ----
                                              br(),
                                              br(),
                                              h3("Download selected gene sets"),
                                              downloadButton("zipdownloadData", label = "Download",align = "center"), #### download all result ----
                                              h5("It may takes up to few minuites", style = "color: #F5918F"),
                                              align = "center"
                                              ),
                                              
                                              
                                            ),
                                            column(
                                              7,
                                              plotOutput("plotgraph1"), ### output: GSEA result plot ----
                                              checkboxInput("pvalue_table",
                                                            label = "put p value table") #### select p value ----
                                            )
                                          ),
                                          column(#separation of main panel to 2 columns ->1
                                            4,
                                            DTOutput("rank_table"), ### output: Rank and Core Enrichment table ----
                                            br()
                                          )
                                          
                                          
                            )# end of 2nd main tab
                            
                          )# end of tab set panel for main panel
                          
                )# end of main panel for outputs
                
              )# side bar layout for TCGA
              
    ), ## end first top tab for TCGA GSEA analysis
    
    
    
    tabPanel( # second top tab for TCGA GSEA analysis ----
              
              title = "Apply GSEA on your own data", 
              
              icon = span(class = "glyphicon glyphicon-search"),
              
              sidebarLayout( ## side bar layout for TCGA
                
                
                sidebarPanel( # Sidebar panel for inputs 
                  width = 3, 
                  
                  div(fileInput( ### file upload ----
                                 "rnkfile_input", 
                                 "Upload your rank file",
                                 multiple = FALSE,
                                 accept = c(".csv",".txt",".xlsx")
                  ),style="font-size:14px;text-align: center;font-family:Arial,sans-serif;") %>% 
                    helper(type = "inline",
                           icon = "question-circle",
                           title = "",
                           fade = TRUE,
                           content = c("<br><font size='+3'>Please upload your pre-ranked file.",
                                       "You can upload <b>txt</b>, <b>csv</b> or <b>xlsx</b> file,",
                                       "with header without spaces between words.",
                                       "Example: 'Gene_Name' and 'Fold_Change'"),
                           size = "m"),
                  
                  
                  div(htmlOutput("check_input_file"),align = "center"),
                  
                  
                  br(),
                  div(selectInput( ### gene set name select ----
                                   "geneset2",
                                   "Gene sets for analysis",
                                   width = '100%',
                                   genesetList$name,
                                   selected = "hallmark gene sets (H)"
                  ),style="font-size:14px;text-align: center;font-family:Arial,sans-serif;") %>% 
                    helper(type = "inline",
                           icon = "question-circle",
                           title = "",
                           fade = TRUE,
                           content = c("<br><font size='+3'>Here you choose one of the <b>gene sets</b> for GSEA.",
                                       "<br>For more details about the gene set, please check the <a href = https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp target = _blank ><b>MSigDB web site</b></a>."),
                           size = "m"),
                  br(),
                  
                  
                  ## advance options
                  
                  div(dropdownButton(
                    tags$h3("Advance options"),
                    
                    
                    
                    div(numericInput(
                      "nPerm_uploaded",
                      "Change permutations",
                      width = '80%',
                      min = 1000,
                      max = 100000,
                      value = 10000,
                      step = 100
                    ),align = "center"),
                    
                    
                    div(numericInput(
                      "minGSSize_uploaded",
                      "Change minimum number of genes in set",
                      width = '80%',
                      min = 1,
                      max = 50,
                      value = 15,
                      step = 1
                    ),align = "center"),
                    
                    
                    div(numericInput(
                      "maxGSSize_uploaded",
                      "Change maximum number of genes in set",
                      width = '80%',
                      min = 100,
                      max = 1000,
                      value = 500,
                      step = 100
                    ),align = "center"),
                    
                    
                    div(selectInput( ### add option for Normalization
                      "normalization_uploaded",
                      "Select normalization option",
                      width = '80%',
                      c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                      selected = "BH"
                    ),align = "center"),
                    
                    div(selectInput( ### add option for Normalization
                      "exponent_uploaded",
                      "Select exponent option",
                      width = '80%',
                      c("classic", "weighted", "weighted_p2", "weighted_p3"),
                      selected = "weighted"
                    ),align = "center"),
                    
                    div(numericInput(
                      "pcutoff_uploaded",
                      "Change p value cut off",
                      width = '80%',
                      min = 0,
                      max = 1,
                      value = 1,
                      step = 0.01
                    ),align = "center"),
                    
                    div(actionButton("reset_advance_uploaded",
                                     "Reset"),align = "right"),
                    br(),
                    br(),
                    
                    
                    inputId = "mydropdown_uploaded",
                    label = "Advanced settings",
                    icon = icon("sliders-h", verify_fa = FALSE),
                    circle = FALSE, 
                    status = "danger", 
                    size = "sm",
                    up = FALSE,
                    #icon = icon("cog"), 
                    width = "300px",
                    tooltip = FALSE
                    
                  ),align = "center") %>% 
                    helper(type = "inline",
                           icon = "question-circle",
                           title = "",
                           fade = TRUE,
                           content = c("<br><font size='+3'>Here you can change details for GSEA.",
                                       "Default:",
                                       "<img src = 'advance_default.jpg'>"),
                           size = "m"),
                  br(),
                  
                  div(actionButton( ### action button for GSEA: uploaded file ----
                                    "make_intermediatesheet_uploaded",
                                    "Apply GSEA",
                                    width = "56%"
                  ),align = "center")
                  
                ),# end of sidebar panels for input
                
                mainPanel(# Main panel for  outputs ----
                          
                          width = 9,
                          
                          tabsetPanel( # tab set panel for main panel
                            type = "hidden",
                            id = "GSEAanalyses_uploaded",
                            
                            tabPanelBody( ## first main tab ----
                                          "first_tab_uploaded",
                                          column(9, 
                                          div(img(src = geni_png, height = 200), align = "center"),
                                          br(),
                                          h2("You are now in 'Apply GSEA on your own data' Tab",align = "center"),
                                          br(),
                                          br(),
                                          h3("Here you can easily apply GSEA on your pre-ranked data.",align = "center"),
                                          
                                          
                                          h3("Please upload your rank data, and choose the study of interests.", align = "center"),
                                          
                                          br(),
                                          br(),
                                          br(),
                                          column(4,offset = 4,
                                                 
                                                 
                                                 div(
                                                   h3("Download example file from here"),
                                                   selectInput( ### select download format ----
                                                                "format_example",
                                                                "Download format",
                                                                c("xlsx",
                                                                  "csv",
                                                                  "txt"),
                                                                width = "30%"
                                                                
                                                   ),
                                                   downloadButton("example",
                                                                  "Download"), align = "center"))
                            )),
                            
                            
                            tabPanelBody( ## 2nd main tab ----
                                          "2nd_tab_uploaded",
                                          
                                          h1("Choose Gene group to get details",align = "center"),
                                          
                                          br(),
                                          column(width = 12, style = "height:500px; overflow-y: scroll;overflow-x: scroll;",
                                                 withLoader(DTOutput('intermediate_result_uploaded'),type = "html",loader = "dnaspin") %>% 
                                                   helper(type = "inline",
                                                          icon = "question-circle",
                                                          title = "",
                                                          fade = TRUE,
                                                          content = c("<br><font size='+3'>Here you can see summary table of GSEA.",
                                                                      "<br>To see details, click the row.",
                                                                      "<br>To download, select right column."),
                                                          size = "m")
                                          ),
                                          
                                          
                                          br(),
                                          
                                          column(
                                            12,
                                            plotOutput("summary_dotplot_up_uploaded") %>% 
                                              helper(id = "summary_dotplot_up_uploaded_helper",
                                                     type = "inline",
                                                     icon = "question-circle",
                                                     title = "",
                                                     fade = TRUE,
                                                     content = c("<br><font size='+3'>You can hide/show summary plots by clicking the button below."),
                                                     size = "m")
                                          ),
                                          column(
                                            12,
                                            plotOutput("summary_dotplot_down_uploaded")
                                          ),
                                          
                                          column(#separation of main panel to 2 columns ->1
                                            8,
                                            
                                            
                                            column(
                                              5,
                                              actionButton( ### action button for GSEA ----
                                                            "show_upload", 
                                                            "Show/Hide summary plots"
                                              ),
                                              br(),
                                              br(),
                                              actionBttn( ### action button for gene set details
                                                
                                                "detail_uploaded",
                                                HTML("Click for gene set detail<br/>(jump to MSigDB page)")
                                                
                                              ),
                                              
                                              h4("Download selected gene sets"),
                                              
                                              downloadButton("zipdownloadData_uploaded", label = "Download",align = "center"), #### download all result ----
                                              br(),
                                              br(),
                                              br(),
                                              tableOutput("details_table_uploaded") ### output: small table result ----
                                            ),
                                            column(
                                              7,
                                              plotOutput("plotgraph1_uploaded"), ### output: GSEA result plot ----
                                              checkboxInput("pvalue_table_uploaded",
                                                            label = "put p value table") #### select p value ----
                                            )
                                          ),
                                          column(#separation of main panel to 2 columns ->1
                                            4,
                                            DTOutput("rank_table_uploaded"), ### output: Rank and Core Enrichment table ----
                                            br()
                                          )
                                          
                                          
                            )# end of 2nd main tab
                            
                          )# end of tab set panel for main panel
                          
                )# end of main panel for outputs
                
              )# side bar layout for TCGA
    ), ## end first top tab for TCGA GSEA analysis
    
    
    tabPanel(# -----second top tab------
             title = "Download gene expression data",
             icon = span(class = "glyphicon glyphicon-download-alt"),
             #

             div(
               img(src = geni_png, height = 200),
               h1("You can download gene expression data here"),
               selectInput( ### select study to download ----
                            "Study2",
                            "Chose your interest study",
                            width = '30%',
                            readRDS("dat/studies.RDS")$name,
                            selected = "Adrenocortical Carcinoma (TCGA, Firehose Legacy)"

               ),

               br(),
               br(),

               selectInput( ### select download format ----
                            "format3",
                            "Chose the download format",
                            c("csv",
                              "txt"),
                            width = "10%"
               ),

               downloadButton('downloadData3', 'Download'),align = "center"), ### download original expression data button ----

    ), # end of third top panel
    
    
    tabPanel(# -----4th top tab------
             
             title = "HELP",
             icon = span(class = "glyphicon glyphicon-question-sign"),
             fluidPage(
               column(width = 10,offset = 1,style = "background-color:powderblue;",
                      
                       htmlOutput("help")
               )
             )
             
    ),
    navbarMenu("About",
               tabPanel(# -----About 1st tab------
                        
                        title = "Contact",
                        icon = span(class = "glyphicon glyphicon-user"),
                        div(
                          img(src = geni_png, height = 200),
                          h1("Yoav Shaul's lab"),
                          br(),
                          br(),
                          h4("We are working on metabolism and cancer."),
                          br(),
                          
                          br(),
                          a(href = "https://www.shaullab.com/home","Visit our lab web page",target="_blank"),
                          br(),
                          a(href = "https://www.shaullab.com/home", img(height = 120, 
                                                                        width = 90,
                                                                        src = lab_pic2_png),target="_blank"),
                          
                          align = "center")
               ),
               tabPanel(# -----About 2nd tab------
                        
                        title = "Terms of Use",
                        column(width = 8, offset = 2,
                               #icon = span(class = "glyphicon glyphicon-user"),
                               div(
                                 img(src = geni_png, height = 200),
                                 h1("Terms of Use"),
                                 br(),
                                 h4(),
                                 br(),
                                 h4("The use of GENI web site tools and services are free and open to all users. For technical support, please feel free to contact us."),
                                 h4("GENI also contains links to other websites, GENI does not have any knowledge of or control over the information in such other websites."),
                                 
                                 h4(""),
                                 br(),
                                 
                                 align = "center"))
               ),
               tabPanel(
                 title = "Acknowledgements",
                 column(width = 8, offset = 2,
                        #icon = span(class = "glyphicon glyphicon-user"),
                        div(
                          img(src = geni_png, height = 200),
                          h1("Acknowledgements"),
                          br(),
                          br(),
                          h4("Several of GENI's key functionalities were enabled because of the shiny, ggplot2, enrichplot and clusterProfiler R packages as well as msigdb."),
                          h4("I would like to thank my colleagues and friends, expecially Emilia Malachi for the lovely logo."),
                          br(),
                          
                          align = "center"))
               )
    )# end if navbarMenu
  )#end of navber page
)
























server <- function(input, output, session){
  
  proxy_select_study <- dataTableProxy("chosen_study")
  conversion <- readRDS("dat/conversion_updated.RDS")
 

  
  observeEvent(input$reload, {
    refresh()
  })

  #session$onSessionEnded(stopApp)
  
  # load observe_helpers to add help widget 
  observe_helpers(withMathJax = FALSE)
  

  ## set value of demo run
  observeEvent(input$intro,{
    
    demo_gene <- "DPYD"
    updateTextInput(session,input = "GeneID",value = demo_gene)
    
    demo_tissue <- "Breast"
    updateSelectInput(session, input = "tissue_select", selected = demo_tissue)
    
    
    })


  
  ## move to next tab
  observeEvent(input$select_study, { # change tab with click: ----
    # tab name to move
    newvalue <- "chose_study"
    # name of id for tabs
    updateTabItems(session, "GSEAanalyses", newvalue)
  })
  
  

   output$help <- renderUI(includeHTML("dat/help_tab.html"))
  
  ## introduction tour (not working)
  observeEvent(input$intro,{
    introjs(session, options = list("nextLabel" = "Continue",
                                    "prevLabel" = "Previous",
                                    "doneLabel" = "Alright. Let's go",
                                    steps=help_text)
    )
  })
  

  observeEvent(input$make_intermediatesheet, { # change tab with click: ----
    #move to top of the page
    js$toTop();
    # Very important! We do not the gene set selected for download to persist once we select a new analysis.
    selected(c())
    
    # show error message with click "Apply GSEA" if there is luck in information
    if( length(Interested_gene_id()) >= 2 & is.null(input$chosen_study_rows_selected)){ # no gene no study
      
      shinyalert("Oops!", "Please input gene and select study.", type = "error",closeOnClickOutside = TRUE)
      
    }else if (length(Interested_gene_id()) == 1 & is.null(input$chosen_study_rows_selected)){ # no study
      
      shinyalert("Oops!", "Please select study.", type = "error",closeOnClickOutside = TRUE)
      
    }else if (length(Interested_gene_id())>=2  & is.null(input$chosen_study_rows_selected) == FALSE){ # no gene
      
      shinyalert("Oops!", "Please select gene.", type = "error",closeOnClickOutside = TRUE)
      
    } else if (length(Interested_gene_id()) == 1 & is.null(input$chosen_study_rows_selected) == FALSE){
      
      # tab name to move
      newvalue <- "2nd_tab"
      # name of id for tabs
      updateTabItems(session, "GSEAanalyses", newvalue)
    }
    
  })
  
  

  
  
  observeEvent(input$select_study, { # change tab with click: ----
    # tab name to move
    newvalue <- "chose_study_uploaded"
    # name of id for tabs
    updateTabItems(session, "GSEAanalyses_uploaded", newvalue)
  })
  
  observeEvent( input$make_intermediatesheet_uploaded, { # change tab with click: ----
    # move back to the top of the page
    js$toTop();
    
    # Very important! We do not the gene set selected for download to persist once we select a new analysis.
    selected(c())
    infile <- input$rnkfile_input
    
    # show error message according to input file with click "Apply GSEA"
    if (is.null(infile)) {
      shinyalert("Oops!", "Please upload your rank file.", type = "error")
      
    } else if (endsWith(infile$name, '.xlsx') | endsWith(infile$name, '.csv') | endsWith(infile$name, '.txt')){
      
      # tab name to move
      newvalue <- "2nd_tab_uploaded"
      # name of id for tabs
      updateTabItems(session, "GSEAanalyses_uploaded", newvalue)
      
    }  else { 
      # if none of read option is not fit, return message in red.
      
      shinyalert("Oops!", "Please upload one of these formats: txt  csv  xlsx", type = "error")
    }
  })
  
  observeEvent(input$select_study, { # change tab with click: ----
    # Very important! We do not the gene set selected for download to persist once we select a new analysis.
    selected(c())
  })
  
  
  
  # reset advanced settings to the default
  observeEvent(input$reset_advance,{
    updateNumericInput(session, "nPerm", value = 10000)
    updateNumericInput(session, "minGSSize", value = 15)
    updateNumericInput(session, "maxGSSize", value = 500)
    updateSelectInput(session, "normalization",selected = "BH")
    updateNumericInput(session, "pcutoff", value = 1)
    
  })
  # reset advanced settings to the default
  observeEvent(input$reset_advance_uploaded,{
    updateNumericInput(session, "nPerm_uploaded", value = 10000)
    updateNumericInput(session, "minGSSize_uploaded", value = 15)
    updateNumericInput(session, "maxGSSize_uploaded", value = 500)
    updateSelectInput(session, "normalization_uploaded",selected = "BH")
    updateNumericInput(session, "pcutoff_uploaded", value = 1)
    
  })
  
  # show/hide plot output with click
  observeEvent(input$show,{
    toggle("summary_dotplot_up")
    toggle("summary_dotplot_down")
    toggle("summary_dotplot_up_helper")
  })
  observeEvent(input$show_upload,{
    toggle("summary_dotplot_up_uploaded")
    toggle("summary_dotplot_down_uploaded")
    toggle("summary_dotplot_up_uploaded_helper")
  })
  
  
  
  # cehck inputed file type and show error message
  output$check_input_file <- renderText({
    
    infile <- input$rnkfile_input
    if (is.null(infile)) {
      return(NULL)
      # change read function according to file name
    } else if (endsWith(infile$name, '.xlsx')){
      
      paste("You uploaded xlsx file")
      
      
    } else if (endsWith(infile$name, '.csv')){
      
      paste("You uploaded csv file")
      
      
    } else if (endsWith(infile$name, '.txt')){
      
      paste("You uploaded txt file")
      
    } else { 
      # if none of read option is not fit, return message in red.
      paste("Please upload one of these formats: txt  csv  xlsx")
    }
  })
  
  
  
  
  
  
  # look for gene id if inputted gene name
  Interested_gene_id <- reactive({
    #req(input$GeneID)
    if(is.na(as.numeric(input$GeneID))){# if input is not numeric
      if(is.na(conversion[match(toupper(input$GeneID),conversion$Hugo_Symbol),"Entrez_Gene_Id"]) == FALSE){
        conversion[match(toupper(input$GeneID),conversion$Hugo_Symbol),"Entrez_Gene_Id"]


      } else {
         conversion[grep(toupper(input$GeneID),conversion$Hugo_Symbol),"Entrez_Gene_Id"]


      }
    } else {# if not ( if input is numeric)
      if(is.na(conversion[match(input$GeneID,conversion$Entrez_Gene_Id),"Entrez_Gene_Id"]) == FALSE){
        conversion[match(input$GeneID,conversion$Entrez_Gene_Id),"Entrez_Gene_Id"]
        
      } else{
         conversion[grep(input$GeneID,conversion$Entrez_Gene_Id),"Entrez_Gene_Id"]
        
      }
      
    }
    
    
  })
  
  
  
  output$check_gene_id <- renderText({ # output: check inputted GeneID ----
    
    # if gene name is found in conversion table show ID in green
    if( length(Interested_gene_id()) > 20 ){
      
      paste("Please enter your gene of interest.Example: EXT1, IL6, 2131, 1806")
      
      
    }else if (length(Interested_gene_id()) >= 2 ){
      if(is.na(as.numeric(input$GeneID))){
      paste("# ",conversion[grep(toupper(input$GeneID),conversion$Hugo_Symbol),"Hugo_Symbol"])
      } else {
        paste("# ",conversion[grep(input$GeneID,conversion$Entrez_Gene_Id),"Hugo_Symbol"])
      }
      #print(conversion[grep(toupper(input$GeneID),conversion$Hugo_Symbol),1],sep="\n")
      #print.data.frame(data.frame(conversion[grep(toupper(input$GeneID),conversion$Hugo_Symbol),1]))
    }else if(length(Interested_gene_id()) == 1 ){
      
      paste("<font color=\"#32CD32\"><b>","Gene ID is valid.")
      
    } else {
      
      # if no gene id is found, show message in red
      paste("<font color=\"#FF0000\"><b>Please try another gene name or NCBI gene ID.")
      
    }
    
  }) ## end of output: check inputted GeneID 
  
  
  
  
  select_study <- reactive({
    # ask input
    req(input$tissue_select)
    if (input$tissue_select != "All"){
      table <- filter(mytable, Tissue == input$tissue_select)
      table[,c("name","Number of patients")]
    } else{
      mytable[,c("name","Number of patients")]
    }
    
  })
  
  
  
  output$chosen_study <- renderDT(
    
    datatable(select_study(),selection = 'single', escape = FALSE,rownames=F, 
              options = list(paging = TRUE))
    
  )
  
  # set row name 
  selected_study <- reactive({
    req(select_study())
      
      select_study()[input$chosen_study_rows_selected,"name"]
    
    
  })
  
  
  

  
  
  output$check_study <- renderText({
    # ask input file
    req(selected_study())
    
    
    paste("<font color=\"#32CD32\"><b>","You have selected ",selected_study(),".")
  })
  
  
  
  
  
  datainput <- eventReactive(input$make_intermediatesheet_uploaded | input$make_intermediatesheet ,{ # read data from upload file / TCGA file ----
    # condition = tab name
    if(input$chosentab =="Apply GSEA on your own data") {## read data from upload file ----
      # ask input file
      req(input$rnkfile_input)  
      
      infile <- input$rnkfile_input
      # check inside file, if it's empty, return null
      if (endsWith(infile$name, '.xlsx')){
        library(readxl)
        UploadFile <- read_excel(infile$datapath)
        
      } else if (endsWith(infile$name, '.csv')){
        
        UploadFile <- read.csv(infile$datapath)
        
        
      } else if (endsWith(infile$name, '.txt')){
        
        UploadFile <- read.table(infile$datapath, header = T)
        
        
      }
      
      # check if Gene Name is inputed or Gene ID is inputed
      if(length(grep("G", UploadFile[,1])) != 0){ # if input contains Gene Name
        # set colmn names as gene name
        colnames(UploadFile) <- c("GeneName", "FoldChange")
        # add new column for gene id but fill with NA
        UploadFile$GeneID <- as.character(NA)
        
        n <- 1
        # look for gene name from conversion and add gene id in GeneID column,
        # if no ID is found, add NA
        try(for(i in 1:nrow(UploadFile)){ # convert Gene Name to Gene ID
          
          UploadFile[i,3] <- conversion[match(UploadFile[i,1],conversion$Hugo_Symbol),]$Entrez_Gene_Id
          
          n <- n + 1
          
        })
        
        
        # chose only GeneID and FoldChange column
        UploadFile <- data.frame(UploadFile[, c("GeneID","FoldChange")])
        
      } else if(length(grep("G",UploadFile[,1])) == 0){ # if input contains Gene ID
        # set column names as gene id
        colnames(UploadFile) <- c("GeneID", "FoldChange")
        
      }
      
      UploadFile <- UploadFile[complete.cases(UploadFile),]
      
      ## end of read data from upload file
    } else if(input$chosentab == "Start GENI") {## read data from TCGA file ----
      # ask input for read
      req(selected_study())
      # according to inputted study id, read expression data from dat file
      data <- readRDS(paste0( "dat/", unlist(mytable[selected_study(),"RDS"])))
      
      # delete patiant id column
      data <- data[,-1]
      # extract expression data of interested gene 
      genes <- colnames(data)
      
      
      if(is.na(match(Interested_gene_id(),genes)) == TRUE){
        
        shinyalert("Oops!", "Selected gene is not appeared in this study", type = "error",closeOnClickOutside = TRUE)
        stop()
      }
      gene1 <- data[,match(Interested_gene_id(),genes)]
      # delete expression data of interested gene from file
      data <- data[,-match(Interested_gene_id(),genes)]
      
      # calculate Spearman correlation
      corr <- cor(gene1, data, method="spearman")
      # change x-y axis (row name will be gene and one column of correlation)
      corr <- data.frame(t(corr))
      # put new column gene name from row name
      corr$name <- rownames(corr)
      # omit na row
      corr <- na.omit(corr)
      # set colmun name
      colnames(corr) <- c("Spearman","GeneID")
      # change order
      corr <- corr[,c("GeneID","Spearman")]
      
    }## end of read data from TCGA file
    
  })## end of read data from upload file / TCGA file 
  
  term2g <- reactive({
    if(input$chosentab =="Apply GSEA on your own data") {## read data from upload file ----
      # ask input file ( gene set name)
      req(input$geneset2)
      
      geneset <- input$geneset2
    } else if(input$chosentab =="Start GENI") {
      
      # ask input file ( gene set name)
      req(input$geneset)
      
      geneset <- input$geneset
    }
    
    
    # define category
    category <- genesetList[match(geneset,genesetList$name),]$category
    # define sub category
    if(is.na(genesetList[match(geneset,genesetList$name),]$subcategory) == TRUE){
      
      subcategory <- NA
      
    } else {
      
      subcategory <- genesetList[match(geneset,genesetList$name),]$subcategory
      
    }
    # pull down gene set data from msigdbr api
    msigdbr(species = "Homo sapiens", category = category, subcategory = subcategory) %>% select(gs_name, entrez_gene)
    
  })
  
  
  GSEA_result <- reactive({
    req(datainput())
    if(input$chosentab =="Apply GSEA on your own data") {## make list from uploaded file and apply GSEA ----
      # ask input
      req(input$make_intermediatesheet_uploaded)
      # call correlation data 
      UploadFile <- datainput()
      # restrict digits for GSEA and make list
      original_gene_list <- round(as.numeric(UploadFile$FoldChange), digits = 6)
      # add gene id for same list
      names(original_gene_list) <- UploadFile$GeneID
      
      # sort by fold change
      gene_list <- sort(original_gene_list, decreasing = TRUE)
      # delete duplicate gene names 
      gene_list <- gene_list[!duplicated(names(gene_list))]
      expo_uploaded <- switch(input$exponent_uploaded,
                              "classic" = 0,
                              "weighted" = 1,
                              "weighted_p2" = 2,
                              "weighted_p3" = 3)
      # GSEA analysis Hallmark
      gsea <-GSEA(gene_list,
                  TERM2GENE     = term2g(),
                  minGSSize     = input$minGSSize_uploaded,
                  maxGSSize     = input$maxGSSize_uploaded,
                  pAdjustMethod = input$normalization_uploaded,
                  pvalueCutoff  = input$pcutoff_uploaded,
                  seed          = TRUE,
                  eps           = 1e-100,
                  exponent      = expo_uploaded,
                  nPermSimple   = input$nPerm_uploaded,
                  by            = "fgsea")
      
    } else if(input$chosentab =="Start GENI") {## make list from TCGA file and apply GSEA ----
      req(input$make_intermediatesheet)
      # call correlation data 
      Spearman <- datainput()
      # restrict digits for GSEA and make list
      original_gene_list <- round(as.numeric(Spearman$Spearman), digits = 6)
      # add gene id for same list
      names(original_gene_list) <- Spearman$GeneID
      
      # sort by fold change
      gene_list <- sort(original_gene_list, decreasing = TRUE)
      # delete duplicate gene names 
      gene_list <- gene_list[!duplicated(names(gene_list))]
      
      expo <- switch(input$exponent,
                     "classic" = 0,
                     "weighted" = 1,
                     "weighted_p2" = 2,
                     "weighted_p3" = 3)
      # GSEA analysis Hallmark
      gsea <-GSEA(gene_list,
                  TERM2GENE     = term2g(),
                  minGSSize     = input$minGSSize,
                  maxGSSize     = input$maxGSSize,
                  pAdjustMethod = input$normalization,
                  pvalueCutoff  = input$pcutoff,
                  seed          = TRUE,
                  eps           = 1e-100,
                  nPermSimple   = input$nPerm,
                  exponent      = expo,     
                  by            = "fgsea")
    }
    
    #trace(GSEA,edit = TRUE)
    
  })
  
  
  # make GSEA result summary sheet
  gsea_result_summary <- reactive({
    # ask input
    req(GSEA_result())
    # put as data frame
    as.data.frame(GSEA_result())
    
  })
  
  output_intermediate_table <- reactive({ # make intermediate result of GSEA to select gene set
    # ask input
    req(gsea_result_summary())
    
    # select columns to present as intermediate table
    gsea_result_sum <- gsea_result_summary()[,c(1,3,5,6,7,8)]
    # change column names
    colnames(gsea_result_sum) <- c("GeneSet","Size","NES","NOM p-val","FWER p-val","FDR q-val")
    
    
    
    # remove row names(IDs)
    rownames(gsea_result_sum) <- NULL
    # add Download selection column
    cbind(Download=noIcon, gsea_result_sum)
    
  })
  
  
  
  
  # output of intermediate result
  
  output$intermediate_result_uploaded <- renderDT(selection = 'single',# output of intermediate result ----
                                                  escape = FALSE, 
                                                  output_intermediate_table(),
                                                  rownames=F,
                                                  options = list(paging = FALSE,
                                                                 columnDefs = list(list(className = 'dt-center',
                                                                                        targets = 0)))
  ) 
  
  output$intermediate_result <- renderDT(selection = 'single',# output of intermediate result ----
                                         escape = FALSE, 
                                         output_intermediate_table(),
                                         rownames=F,
                                         options = list(paging = TRUE,
                                                        columnDefs = list(list(className = 'dt-center', 
                                                                               targets = 0))
                                         )
  ) 
  
  
  # set row name 
  selected_geneset <- reactive({ # specify selected row name (ID) ----
    req(output_intermediate_table())
    
    if(input$chosentab =="Apply GSEA on your own data") {
      req(input$intermediate_result_uploaded_rows_selected)
      
      output_intermediate_table()[input$intermediate_result_uploaded_rows_selected, "GeneSet"]
    } else if(input$chosentab =="Start GENI") {
      req(input$intermediate_result_rows_selected)
      
      output_intermediate_table()[input$intermediate_result_rows_selected, "GeneSet"]
    }
  })
  
  
  # link to details
  
  
  
  
  observeEvent(input$detail, {
    req(selected_geneset())
    js$browseURL(paste0("https://www.gsea-msigdb.org/gsea/msigdb/cards/",selected_geneset(),".html"))
  })
  
  
  observeEvent(input$detail_uploaded, {
    req(selected_geneset())
    js$browseURL(paste0("https://www.gsea-msigdb.org/gsea/msigdb/cards/",selected_geneset(),".html"))
  })
  
  
  
  
  
  # modify GSEA plot function ----
  nrank <- reactive(nrow(datainput()))
  nhigh <- reactive(max(datainput()[,2])/2)
  
  # modify GSEA plot function ----
  
  

  
  plot1 <- reactive({ # make plot for selected gene set of GSEA
    # ask input
    req(GSEA_result())
    req(selected_geneset())
    # print(nrank())
    # print(nhigh())
    source("make_enrichment_plot.R")
    # make GSEA plot for selected gene set
    if(input$chosentab =="Apply GSEA on your own data") {
      
      plot <- make_enrichment_plot(GSEA_result(),
                        title = selected_geneset(),
                        geneSetID = selected_geneset(),
                        pvalue_table = input$pvalue_table_uploaded,
                        nrank = nrank(),
                        nhigh = nhigh())
      
    } else if(input$chosentab =="Start GENI") {
      
      plot <- make_enrichment_plot(GSEA_result(),
                        title = selected_geneset(),
                        geneSetID = selected_geneset(),
                        pvalue_table = input$pvalue_table,
                        nrank = nrank(),
                        nhigh = nhigh())
    }
    plot
  }) # make plot for selected gene set of GSEA
  
  output$plotgraph1_uploaded <- renderPlot({ # output: plot of GSEA result
    # ask input
    req(plot1())
    
    plot1()
    
  },width = "auto",
  height = "auto")# end of output: plot of GSEA result
  
  
  output$plotgraph1 <- renderPlot({ # output: plot of GSEA result
    # ask input
    req(plot1())
    
    plot1()
    
  },width = "auto",
  height = "auto")# end of output: plot of GSEA result
  
  
  
  
  
  res_table <- reactive({ # make result summary with values ----
    # ask inputs
    req(GSEA_result())
    req(selected_geneset())
    
    
    # make vector contain names for 1st column
    Name <- c("Enrichment Score (ES)","Normalized Enrichment Score (NES)","Nominal p-value","FWER p-value","FDR q-value")
    # extract values from main result table
    Value <- c(paste(round(GSEA_result()[selected_geneset(), "enrichmentScore"], digits = 4)),
                paste(round(GSEA_result()[selected_geneset(), "NES"], digits = 4)),
                paste(formatC(GSEA_result()[selected_geneset(), "pvalue"], format = "e", digits = 3)),
                paste(formatC(GSEA_result()[selected_geneset(), "p.adjust"], format = "e", digits = 3)),
                paste(formatC(GSEA_result()[selected_geneset(), "qvalue"], format = "e", digits = 3)))
    # make data frame from vectors
    data.frame(Name,Value)
    
  }) # end of make result summary with values
  
  
  
  output$details_table_uploaded <- renderTable({ # output: result summary with values ----
    # ask input
    req(res_table())
    
    res_table()
    
  })
  
  
  output$details_table <- renderTable({ # output: result summary with values ----
    # ask input
    req(res_table())
    
    res_table()
    
  })
  
  
  
  
  
  
  rank_table <- reactive({ # make result table with rank, GeneID, core enrichment ----
    # ask inputs
    req(GSEA_result())
    req(selected_geneset())
    req(term2g())
    # load input data (uploaded / TCGA)
    data <- datainput()
    colnames(data) <- c("GeneID", "diff")
    
    
    term2g <- term2g()
    
    
    
    # extract relevant gene set and delete duplicated gene ids
    gene_set <- unique(filter(term2g, gs_name == selected_geneset()))
    # add gene name column and fill with NA
    gene_set$GeneName <- as.character(NA)
    # add rank column and fill with NA
    gene_set$rank <- as.double(NA)
    # add core enrichment column and fill with "No"
    gene_set$CoreEnrichment <- "No"
    
    # look for gene name and correlation / fold change
    for( i in 1:nrow(gene_set)){
      # look for gene name
      gene_set[i,3] <- conversion[match(gene_set[i,2], conversion$Entrez_Gene_Id), ]$Hugo_Symbol
      # look for correlation or fold change
      gene_set[i,4] <- data[match(gene_set[i,2], data$GeneID), ]$diff
      
    }
    
    # load gsea result summary sheet
    gsea_result_summary <- gsea_result_summary()
    # select core enrichment column and convert to vector of gene ids
    coreID <- str_split(gsea_result_summary[selected_geneset(),"core_enrichment"], pattern = "/")
    # loop for making Yes row
    for( i in 1:nrow(gene_set)){
      
      if(is.na(match(gene_set[i,2],coreID[[1]])) == FALSE){
        # if gene id is found, replace No to Yes
        gene_set[i,5]<- "Yes"
        
      }
      
    }
    # change column order
    gene_set <- gene_set[,c(3,2,4,5)]
    # change row order according to NES
    if(gsea_result_summary[selected_geneset(),"NES"] > 0){
      # if NES positive, order will be from high to low 
      gene_set <- gene_set[order(gene_set$CoreEnrichment, gene_set$rank,decreasing=T),]
      
    } else if (gsea_result_summary[selected_geneset(),"NES"] < 0){
      # if NES negative, order will be from low to high
      gene_set <- gene_set[order(gene_set$CoreEnrichment, gene_set$rank,decreasing=F),]
      
    }
    # order as Yes rows first
    gene_set <- gene_set[order(gene_set$CoreEnrichment, gene_set$rank,decreasing=T),]
    colnames(gene_set) <- c("Symbol","Entrez gene ID","Rank Metric Score", "Core Enrichment")
    
    na.omit(gene_set)
    
    
  }) # end of make result table with rank, GeneID, core enrichment 
  
  
  output$rank_table_uploaded <- renderDT({ # output: result table with rank, GeneID, core enrichment ----
    # ask inputs
    req(selected_geneset())
    req(rank_table())
    # add color to core enrichment column
    datatable(rank_table(),rownames=F, selection = "none",options = list(lengthMenu = c(10,25,50,75,100,150),
                                                                         pageLength = 10))%>%
      formatStyle(
        
        columns = "Core Enrichment",
        background = styleEqual(c("Yes","No"),c("lightgreen","lightpink"))
        
        
      )
  })
  
  output$rank_table <- renderDT({ # output: result table with rank, GeneID, core enrichment ----
    # ask inputs
    req(selected_geneset())
    req(rank_table())
    # add color to core enrichment column
    datatable(rank_table(),rownames=F, selection = "none", options = list(lengthMenu = c(10,25,50,75,100,150),
                                                                          pageLength = 10
    ))%>%
      formatStyle(
        
        columns = "Core Enrichment",
        background = styleEqual(c("Yes","No"),c("lightgreen","lightpink"))
        
        
      )
  })
  
  
  
  
  
  
  # Stores the row indexes to be downloaded
  selected <- reactiveVal(c())
  
  proxy_intermediate_result_uploaded <- dataTableProxy("intermediate_result_uploaded")
  
  proxy_intermediate_result <- dataTableProxy("intermediate_result") 
  
  # Technical variable for updating the table without re-rendering it from scratch.
  
  
  
  
  observeEvent( input$intermediate_result_uploaded_cell_clicked, {
    colC <- input$intermediate_result_uploaded_cell_clicked$col 
    rowC <- input$intermediate_result_uploaded_cell_clicked$row
    
    
    
    # If did not click on the index of the Download column, do nothing.
    req( colC  == 0 )
    
    # If the row is newly added, add it
    if(!rowC %in% selected())
      selected(c(selected(), rowC))
    # else, remove it from the downloads (toggle state)
    else 
      selected( setdiff(selected(), rowC))
    
    
    # Create a static version of output table
    tableToReplace <- output_intermediate_table()
    # Add a download icon to the selected rows
    tableToReplace[selected(),"Download"] <- okIcon
    
    
    
    # Efficiently replace the data.
    proxy_intermediate_result_uploaded %>% replaceData(tableToReplace, resetPaging=F, rownames=F)
    # Select the row clicked upon 
    #(this is a minor technical point. If this were omitted, then the row would be selected for download, but the graph would not be displayed automatically)
    proxy_intermediate_result_uploaded %>% selectRows(rowC)
    
    
  })
  
  observeEvent(input$intermediate_result_cell_clicked , {
    
    # Store the clicked cell coordinates in variables with simple names  
    colC <- input$intermediate_result_cell_clicked$col 
    rowC <- input$intermediate_result_cell_clicked$row
    
    
    
    # If did not click on the index of the Download column, do nothing.
    req( colC  == 0 )
    
    # If the row is newly added, add it
    if(!rowC %in% selected())
      selected(c(selected(), rowC))
    # else, remove it from the downloads (toggle state)
    else 
      selected( setdiff(selected(), rowC))
    
    
    # Create a static version of output table
    tableToReplace <- output_intermediate_table()
    # Add a download icon to the selected rows
    tableToReplace[selected(),"Download"] <- okIcon
    
    
    
    # Efficiently replace the data.
    proxy_intermediate_result %>% replaceData(tableToReplace, resetPaging=F, rownames=F)
    # Select the row clicked upon 
    #(this is a minor technical point. If this were omitted, then the row would be selected for download, but the graph would not be displayed automatically)
    proxy_intermediate_result %>% selectRows(rowC)
    
    
  })
  
  
  
  
  summary_plot_up <- reactive({ # make summary plot for NES positive gene sets ----
    # ask input
    req(gsea_result_summary())
    # extract needed columns to make summary dot plot
    ggdata <- gsea_result_summary()[,c("Description","NES","p.adjust","qvalue","setSize")]
    # set column names
    colnames(ggdata) <- c("Description","NES","FWER p-val","FDR q-val","Size")
    # chose only up regulated gene sets (NES positive)
    ggdata_up <- ggdata[gsea_result_summary()$NES > 0,]
    # change order according to adjusted p value
    ggdata_up$Description<-reorder(ggdata_up$Description, ggdata_up$NES, decreasing = F)
    
    if(nrow(ggdata_up)>15){
      
      ggdata_up <- ggdata_up[1:15,]
      
    }
    ggplot(ggdata_up,aes(x = Description, y = NES, size = Size, color = `FDR q-val`, fill = `FDR q-val`))+
      geom_point(shape = 21)+
      scale_size(range = c(3,8))+
      scale_color_continuous(low = 'red', high = 'blue')+
      scale_fill_continuous(low = 'red', high = 'blue')+
      xlab('Gene set (ordered by FWER p-val)')+ 
      ylab('Enrichment score')+
      labs(title = 'Summary dot plot (upregulated)')+
      theme_bw(base_size = 15)+
      theme(
        legend.position = 'right',
        plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1, hjust = 0.5, color = "black"),
        axis.text.x = element_text(angle = 0, size = 11, face = 'bold', hjust = 0.5, color = "black"),
        axis.text.y = element_text(angle = 0, size = 11, face = 'bold', vjust = 0.5, color = "black"),
        axis.title = element_text(size = 11, face = 'bold'),
        axis.title.y = element_text(size = 11, face = 'bold'),
        legend.key.size = unit(0.5, "cm"), # Sets overall area/size of the legend
        legend.text = element_text(size = 12, face = "bold"), # Text size
        title = element_text(size = 14, face = "bold"))+
      guides(size = guide_legend(order = 1))+
      coord_flip()
    
  })
  
  summary_plot_down <- reactive({ # make summary plot for NES positive gene sets ----
    # ask input
    req(gsea_result_summary())
    
    
    # extract needed columns to make summary dot plot
    ggdata <- gsea_result_summary()[,c("Description","NES","p.adjust","qvalue","setSize")]
    # set column names
    colnames(ggdata) <- c("Description","NES","FWER p-val","FDR q-val","Size")
    # chose only down regulated gene sets (NES negative)
    ggdata_down <- ggdata[gsea_result_summary()$NES < 0,]
    # change order according to adjusted p value
    ggdata_down$Description<-reorder(ggdata_down$Description, ggdata_down$NES, decreasing = T)
    if(nrow(ggdata_down)>15){
      
      ggdata_down <- ggdata_down[1:15,]
      
    }
    ggplot(ggdata_down,aes(x = Description, y = NES, size = Size, color = `FDR q-val`, fill = `FDR q-val`))+
      geom_point(shape = 21)+
      scale_size(range = c(3,8))+
      scale_color_continuous(low = 'red', high = 'blue')+
      scale_fill_continuous(low = 'red', high = 'blue')+
      xlab('Gene set (ordered by FWER p-val)')+ 
      ylab('Enrichment score')+
      labs(title = 'Summary dot plot (downregulated)')+
      theme_bw(base_size = 15)+
      theme(
        legend.position = 'right',
        plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1, hjust = 0.5, color = "black"),
        axis.text.x = element_text(angle = 0, size = 11, face = 'bold', hjust = 0.5, color = "black"),
        axis.text.y = element_text(angle = 0, size = 11, face = 'bold', vjust = 0.5, color = "black"),
        axis.title = element_text(size = 11, face = 'bold'),
        axis.title.y = element_text(size = 11, face = 'bold'),
        legend.key.size = unit(0.5, "cm"), # Sets overall area/size of the legend
        legend.text = element_text(size = 12, face = "bold"), # Text size
        title = element_text(size = 14, face = "bold"))+
      guides(size = guide_legend(order = 1))+
      coord_flip()
    
  })
  
  output$summary_dotplot_up_uploaded <- renderPlot({
    req(summary_plot_up())
    summary_plot_up()
  })
  output$summary_dotplot_down_uploaded <- renderPlot({
    req(summary_plot_down())
    summary_plot_down()
  })
  
  output$summary_dotplot_up <- renderPlot({
    req(summary_plot_up())
    summary_plot_up()
  })
  output$summary_dotplot_down <- renderPlot({
    req(summary_plot_down())
    summary_plot_down()
  })
  
  
  
  output$zipdownloadData_uploaded <- downloadHandler( # zip-download all result ----
                                                      # set name for saved file
                                                      filename = function(){
                                                        if( length(selected()) == 0){
                                                          shinyalert("Oops!", "Please select gene set by clicking first colum.", type = "error",closeOnClickOutside = TRUE)
                                                          
                                                        } else {
                                                          infilename <- input$rnkfile_input$name
                                                          paste0(input$geneset2,"_",infilename, ".zip", sep="")
                                                        }
                                                        
                                                      },
                                                      # make files for the zip file
                                                      content = function(fname) {
                                                        # ask inputs
                                                        req(gsea_result_summary())
                                                        req(datainput())
                                                        req(selected())
                                                        # empty vector for file names
                                                        fs <- c()
                                                        # set directory to save file temporally (change every time)
                                                        fdir <- tempdir() %>% paste0("/", runif(1) , "/")
                                                        # create directory
                                                        dir.create(fdir) 
                                                        
                                                        # common summary result for all gene sets
                                                        fname1 <- paste(fdir,"summary_result.xlsx", sep = "")
                                                        # summary dot plot for up regulated gene sets
                                                        fname2 <- paste(fdir,"summary_dotplot_up.pdf")
                                                        # summary dot plot for down regulated gene sets
                                                        fname3 <- paste(fdir,"summary_dotplot_down.pdf")
                                                        # output of correlation / fold change used for GSEA
                                                        fname4 <- paste(fdir, "rank_file_all.xlsx", sep = "")
                                                        # make out put of Excel file for all gene sets result
                                                        gsea_result_summary <- gsea_result_summary()
                                                        colnames(gsea_result_summary) <- c("GeneSet","Description","Size","ES","NES","NOM p-val","FWER p-val","FDR q-val", "rank","leading_edge","core_enrichment")
                                                        write_xlsx(gsea_result_summary, path = fname1)
                                                        # make out put of dot plot for up regulated gene sets
                                                        ggsave(filename = fname2,
                                                               plot = summary_plot_up(),
                                                               height = 13,
                                                               width = 13,
                                                               device = 'pdf')
                                                        # make out put of dot plot for down regulated gene sets
                                                        ggsave(filename = fname3,
                                                               plot = summary_plot_down(),
                                                               height = 13,
                                                               width = 13,
                                                               device = 'pdf')
                                                        # call data frame which used for GSEA
                                                        data <- datainput()
                                                        # add new column for gene name and fill with NA
                                                        data$GeneName <- as.character(NA)
                                                        # look for gene name and add
                                                        for (i in 1:nrow(data)){
                                                          
                                                          data[i,"GeneName"] <- conversion[match(data[i,"GeneID"],
                                                                                                 conversion$Entrez_Gene_Id),
                                                          ]$Hugo_Symbol
                                                          
                                                        }
                                                        
                                                        # change column order
                                                        data <- data[,c(3,1,2)]
                                                        # make output of correlation / fold change used for GSEA
                                                        write_xlsx(data, path = fname4)
                                                        # add all common file names to fs
                                                        fs <- c(fname1,fname2,fname3,fname4)
                                                        # loop for making GSEA plot result and rank output ----
                                                        
                                                        ## make rank result ----
                                                        for (i in 1:length(selected())){
                                                          # loop for selected gene set name
                                                          geneset_name <- output_intermediate_table()[selected()[i], "GeneSet"]
                                                          # rank result of selected gene sets
                                                          path <- paste(fdir, geneset_name,"_rank.xlsx", sep = "")
                                                          # GSEA plot of selected gene sets
                                                          path2 <- paste(fdir,geneset_name,".pdf", sep = "")
                                                          # add file names to fs
                                                          fs <- c(fs, path, path2)
                                                          # call data frame which used for GSEA
                                                          Spearman <- datainput()
                                                          # set column names 
                                                          colnames(Spearman) <- c("GeneID", "Spearman")
                                                          # make condition according to the input
                                                          if(input$chosentab =="Apply GSEA on your own data") {## load gene set data for uploaded file ----
                                                            # ask inputs
                                                            req(term2g())
                                                            term2g <- term2g()
                                                            
                                                          } else if(input$chosentab =="Start GENI") {## load gene set data for TCGA file ----
                                                            # ask inputs
                                                            req(term2g())
                                                            term2g <- term2g()
                                                            
                                                          }
                                                          # extract relevant gene set and delete duplicated gene ids
                                                          gene_set <- unique(filter(term2g, gs_name == geneset_name))
                                                          # add gene name column and fill with NA
                                                          gene_set$GeneName <- as.character(NA)
                                                          # add rank column and fill with NA
                                                          gene_set$rank <- as.double(NA)
                                                          # add core enrichment column and fill with "No"
                                                          gene_set$CoreEnrichment <- "No"
                                                          
                                                          # look for gene name and correlation / fold change
                                                          for( i in 1:nrow(gene_set)){
                                                            # look for gene name
                                                            gene_set[i,3] <- conversion[match(gene_set[i,2], conversion$Entrez_Gene_Id), ]$Hugo_Symbol
                                                            # look for correlation or fold change
                                                            gene_set[i,4] <- Spearman[match(gene_set[i,2], Spearman$GeneID), ]$Spearman
                                                            
                                                          }
                                                          
                                                          # load gsea result summary sheet
                                                          gsea_result_summary <- gsea_result_summary()
                                                          # select core enrichment column and convert to vector of gene ids
                                                          coreID <- str_split(gsea_result_summary[geneset_name,"core_enrichment"], pattern = "/")
                                                          # loop for making Yes row
                                                          for( i in 1:nrow(gene_set)){
                                                            
                                                            if(is.na(match(gene_set[i,2],coreID[[1]])) == FALSE){
                                                              # if gene id is found, replace No to Yes
                                                              gene_set[i,5]<- "Yes"
                                                              
                                                            }
                                                            
                                                          }
                                                          # change column order
                                                          gene_set <- gene_set[,c(3,2,4,5)]
                                                          # change row order according to NES
                                                          if(gsea_result_summary[geneset_name,"NES"] > 0){
                                                            # if NES positive, order will be from high to low 
                                                            gene_set <- gene_set[order(gene_set$CoreEnrichment, gene_set$rank,decreasing=T),]
                                                            
                                                          } else if (gsea_result_summary[geneset_name,"NES"] < 0){
                                                            # if NES negative, order will be from low to high
                                                            gene_set <- gene_set[order(gene_set$CoreEnrichment, gene_set$rank,decreasing=F),]
                                                            
                                                          }
                                                          # order as Yes rows first
                                                          gene_set <-gene_set[order(gene_set$CoreEnrichment, gene_set$rank,decreasing=T),]
                                                          # add column names
                                                          colnames(gene_set) <- c("Symbol","Entrez gene ID","Rank Metric Score", "Core Enrichment")
                                                          # delete empty rows
                                                          gene_set <- na.omit(gene_set)
                                                          ## end of make rank result ----
                                                          
                                                          # make GSEA result plot for selected gene sets
                                                          source("make_enrichment_plot.R")
                                                          plot_to_save <- make_enrichment_plot(GSEA_result(),
                                                                                               title = geneset_name,
                                                                                               geneSetID = geneset_name,
                                                                                               pvalue_table = input$pvalue_table_uploaded,
                                                                                               nrank = nrank(),
                                                                                               nhigh = nhigh())
                                                          
                                                          
                                                          # make output of rank result 
                                                          write_xlsx(gene_set, path = path)
                                                          # make out put of GSEA result for selected gene sets
                                                          ggsave(filename = path2,
                                                                 plot = plot_to_save,
                                                                 height = 6,
                                                                 width = 6,
                                                                 device = 'pdf')
                                                          
                                                        }# end of loop
                                                        
                                                        # zip the files
                                                        zip(zipfile=fname, files=fs, extras = "-j")
                                                        
                                                      },
                                                      # set export type
                                                      contentType = "application/zip"
                                                      
  )# end of zip-download all result
  
  output$zipdownloadData <- downloadHandler( # zip-download all result ----
                                             # set name for saved file
                                             
                                             
                                             filename = function(){
                                               if( length(selected()) == 0){
                                                 shinyalert("Oops!", "Please select gene set by clicking first colum.", type = "error",closeOnClickOutside = TRUE)
                                                 
                                               }else {
                                                 paste(input$GeneID,"_",select_study()[input$chosen_study_rows_selected,"name"],"_",input$geneset, ".zip", sep="")
                                               } 
                                             },
                                             # make files for the zip file
                                             content = function(fname) {
                                               # ask inputs
                                               req(gsea_result_summary())
                                               req(datainput())
                                               req(selected())
                                               # empty vector for file names
                                               fs <- c()
                                               # set directory to save file temporally (change every time)
                                               fdir <- tempdir() %>% paste0("/", runif(1) , "/")
                                               # create directory
                                               dir.create(fdir) 
                                               
                                               # common summary result for all gene sets
                                               fname1 <- paste(fdir,"summary_result.xlsx", sep = "")
                                               # summary dot plot for up regulated gene sets
                                               fname2 <- paste(fdir,"summary_dotplot_up.pdf")
                                               # summary dot plot for down regulated gene sets
                                               fname3 <- paste(fdir,"summary_dotplot_down.pdf")
                                               # output of correlation / fold change used for GSEA
                                               fname4 <- paste(fdir, "rank_file_all.xlsx", sep = "")
                                               # make out put of Excel file for all gene sets result
                                               gsea_result_summary <- gsea_result_summary()
                                               colnames(gsea_result_summary) <- c("GeneSet","Description","Size","ES","NES","NOM p-val","FWER p-val","FDR q-val", "rank","leading_edge","core_enrichment")
                                               write_xlsx(gsea_result_summary, path = fname1)
                                               # make out put of dot plot for up regulated gene sets
                                               ggsave(filename = fname2,
                                                      plot = summary_plot_up(),
                                                      height = 13,
                                                      width = 13,
                                                      device = 'pdf')
                                               # make out put of dot plot for down regulated gene sets
                                               ggsave(filename = fname3,
                                                      plot = summary_plot_down(),
                                                      height = 13,
                                                      width = 13,
                                                      device = 'pdf')
                                               # call data frame which used for GSEA
                                               data <- datainput()
                                               # add new column for gene name and fill with NA
                                               data$GeneName <- as.character(NA)
                                               # look for gene name and add
                                               for (i in 1:nrow(data)){
                                                 
                                                 data[i,"GeneName"] <- conversion[match(data[i,"GeneID"],
                                                                                        conversion$Entrez_Gene_Id),
                                                 ]$Hugo_Symbol
                                                 
                                               }
                                               
                                               # change column order
                                               data <- data[,c(3,1,2)]
                                               # make output of correlation / fold change used for GSEA
                                               write_xlsx(data, path = fname4)
                                               # add all common file names to fs
                                               fs <- c(fname1,fname2,fname3,fname4)
                                               # loop for making GSEA plot result and rank output ----
                                               
                                                
                                               
                                               ## make rank result ----
                                               for (i in 1:length(selected())){
                                                 # loop for selected gene set name
                                                 geneset_name <- output_intermediate_table()[selected()[i], "GeneSet"]
                                                 # rank result of selected gene sets
                                                 path <- paste(fdir, geneset_name,"_rank.xlsx", sep = "")
                                                 # GSEA plot of selected gene sets
                                                 path2 <- paste(fdir,geneset_name,".pdf", sep = "")
                                                 # add file names to fs
                                                 fs <- c(fs, path, path2)
                                                 # call data frame which used for GSEA
                                                 Spearman <- datainput()
                                                 # set column names 
                                                 colnames(Spearman) <- c("GeneID", "Spearman")
                                                 # make condition according to the input
                                                 if(input$chosentab =="Apply GSEA on your own data") {## load gene set data for uploaded file ----
                                                   # ask inputs
                                                   req(term2g())
                                                   term2g <- term2g()
                                                   
                                                 } else if(input$chosentab =="Start GENI") {## load gene set data for TCGA file ----
                                                   # ask inputs
                                                   req(term2g())
                                                   term2g <- term2g()
                                                   
                                                 }
                                                 # extract relevant gene set and delete duplicated gene ids
                                                 gene_set <- unique(filter(term2g, gs_name == geneset_name))
                                                 # add gene name column and fill with NA
                                                 gene_set$GeneName <- as.character(NA)
                                                 # add rank column and fill with NA
                                                 gene_set$rank <- as.double(NA)
                                                 # add core enrichment column and fill with "No"
                                                 gene_set$CoreEnrichment <- "No"
                                                 
                                                 # look for gene name and correlation / fold change
                                                 for( i in 1:nrow(gene_set)){
                                                   # look for gene name
                                                   gene_set[i,3] <- conversion[match(gene_set[i,2], conversion$Entrez_Gene_Id), ]$Hugo_Symbol
                                                   # look for correlation or fold change
                                                   gene_set[i,4] <- Spearman[match(gene_set[i,2], Spearman$GeneID), ]$Spearman
                                                   
                                                 }
                                                 
                                                 # load gsea result summary sheet
                                                 gsea_result_summary <- gsea_result_summary()
                                                 # select core enrichment column and convert to vector of gene ids
                                                 coreID <- str_split(gsea_result_summary[geneset_name,"core_enrichment"], pattern = "/")
                                                 # loop for making Yes row
                                                 for( i in 1:nrow(gene_set)){
                                                   
                                                   if(is.na(match(gene_set[i,2],coreID[[1]])) == FALSE){
                                                     # if gene id is found, replace No to Yes
                                                     gene_set[i,5]<- "Yes"
                                                     
                                                   }
                                                   
                                                 }
                                                 # change column order
                                                 gene_set <- gene_set[,c(3,2,4,5)]
                                                 # change row order according to NES
                                                 if(gsea_result_summary[geneset_name,"NES"] > 0){
                                                   # if NES positive, order will be from high to low 
                                                   gene_set <- gene_set[order(gene_set$CoreEnrichment, gene_set$rank,decreasing=T),]
                                                   
                                                 } else if (gsea_result_summary[geneset_name,"NES"] < 0){
                                                   # if NES negative, order will be from low to high
                                                   gene_set <- gene_set[order(gene_set$CoreEnrichment, gene_set$rank,decreasing=F),]
                                                   
                                                 }
                                                 # order as Yes rows first
                                                 gene_set <-gene_set[order(gene_set$CoreEnrichment, gene_set$rank,decreasing=T),]
                                                 # add column names
                                                 colnames(gene_set) <- c("Symbol","Entrez gene ID","Rank Metric Score", "Core Enrichment")
                                                 # delete empty rows
                                                 gene_set <- na.omit(gene_set)
                                                 ## end of make rank result ----
                                                 source("make_enrichment_plot.R")
                                                 # make GSEA result plot for selected gene sets
                                                 plot_to_save <-  make_enrichment_plot(GSEA_result(),
                                                                                       title = geneset_name,
                                                                                       geneSetID = geneset_name,
                                                                                       pvalue_table = input$pvalue_table,
                                                                                       nrank = nrank(),
                                                                                       nhigh = nhigh())
                                                
                                                 # make output of rank result 
                                                 write_xlsx(gene_set, path = path)
                                                 # make out put of GSEA result for selected gene sets
                                                 ggsave(filename = path2,
                                                        plot = plot_to_save,
                                                        height = 6,
                                                        width = 6,
                                                        device = 'pdf')
                                                 
                                               }# end of loop
                                               
                                               # zip the files
                                               zip(zipfile=fname, files=fs, extras = "-j")
                                               
                                             },
                                             # set export type
                                             contentType = "application/zip"
                                             
  )# end of zip-download all result
  
  
  
  
  
  
  
  output$downloadData3 = downloadHandler(
    filename = function() {
      paste(input$Study2,".",input$format3, sep = "")
    },
    content = function(file) {
      if(input$format3 == "csv"){
        idx <- (mytable$name == input$Study2) %>% which()
        
        data <- as.data.frame(readRDS(paste0( "dat/", mytable[idx, ]$RDS)))
        write.csv(data, file, row.names = FALSE)
      }else if(input$format3 == "txt"){
        idx <- (mytable$name == input$Study2) %>% which()
        
        data <- as.data.frame(readRDS(paste0( "dat/", mytable[idx, ]$RDS)))
        write.table(data, file)
      }
    }
  )
  
  output$example <- downloadHandler(
    filename = function() {
      paste("Example.",input$format_example, sep = "")
    },
    content = function(file) {
      if(input$format_example == "xlsx"){
        data <- readRDS("dat/Example.RDS")
        write_xlsx(data, file,)
      }else if(input$format_example == "csv"){
        data <- readRDS("dat/Example.RDS")
        write.csv(data, file, row.names = F)
      }else if(input$format_example == "txt"){
        data <- readRDS("dat/Example.RDS")
        write.table(data, file, row.names = F)
      }
    }
  ) 
  
  
}## end of function



shinyApp(ui = ui, server=server)
