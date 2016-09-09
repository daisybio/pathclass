library(shiny)
library(shinyjs)
library(shinydashboard)
library(org.Hs.eg.db)

dashboardPage(
    dashboardHeader(
    ),
    dashboardSidebar(
        sidebarMenu(
            menuItem("About", tabName = "About", icon = icon("info")),
            menuItem("Selected Pathways", tabName = "Pathways", icon = icon("bookmark")),
            menuItem("TCGA BRCA Data", tabName = "TCGA_array", icon = icon("th")),
            menuItem("Gene Expression Omnibus", tabName = "GEO", icon = icon("bus")),
            menuItem("Upload Data", tabName = "CUSTOM", icon = icon("cloud-upload"))
        )
    ),
    dashboardBody(
        useShinyjs(),
        tags$header(tags$link(rel = "stylesheet", type = "text/css", href = "main.css"),
                    tags$script(HTML('
                        $(document).ready(function(){
                            for(i = 2; i < 9; i++){
                                $("#custom_nav li:nth-child(" + i + ") a").hide();
                            }
                            for(i = 2; i < 7; i++){
                                $("#geo_nav li:nth-child(" + i + ") a").hide();
                            }
                        });
                    '))),
        tabItems(
            tabItem(tabName = "About",
                HTML(aboutText)
            ),
            tabItem(tabName = "Pathways",
                    h2("Selected Pathways"),
                    tabsetPanel(
                        tabPanel("Overview", plotOutput("pathways_overview")),
                        tabPanel("Pathway Members",
                                 fluidRow(box(
                        selectInput("selp_dataset", "Select a predictor",
                                    choices = pathway.sources),
                        selectInput("selp_pathway", "Select a pathway",
                                    choices = NULL)
                                 )),
                        dataTableOutput("selected_pathway_members")
                        )
                    )
            ),

            tabItem(tabName = "TCGA_array",
                    h2("TCGA BRCA data (training data)"),
                    tabsetPanel(
                        tabPanel("Setup", box(
                            selectInput("tcga_predictors",
                                        "Select a reference / gold standard",
                                        choices = all.sources,
                                        selected = "pam50")
                        )),
                        tabPanel("Predictions Table",
                                 box(width = 2, downloadButton("download_tcga_subtypes", "Download")),
                                 dataTableOutput("gssea_result_tcga")),
                        tabPanel("Predictions Plot", value = "custom_tab2",
                                 plotOutput("gssea_result_tcga_heatmap")),
                        tabPanel("Performance (subtypes)",
                                 plotOutput("class_error_subt_plot_tcga_array")
                        ),
                        tabPanel("Performance (predictors)",
                                 plotOutput("class_error_pred_plot_tcga_array")
                        ),
                        tabPanel("Confusion matrices",
                                 plotOutput("class_error_cm_plot_tcga_array")
                        )
                    )

            ),
            tabItem(tabName = "GEO",
                    h2("Import data from the Gene Expression Omnibus for subtyping"),
                    tabsetPanel(
                        id = "geo_nav",
                        tabPanel("Setup",
                             fluidRow(
                                 column(width = 4,
                                 box(width = NULL,
                                     title = "Enter a GEO identifier",
                                     textInput("geo_id", "GSE", "GSE45827"),
                                     shinysky::shinyalert("geo_error"),
                                     actionButton("geo_dl_button", "Download", icon = icon("download"))
                                 )),
                                 column(width = 4,
                                 conditionalPanel("output.geo_file_uploaded",
                                 box(title = "Select reference / gold standard",
                                     width = NULL,
                                     selectInput("geo_predictors", "Predictor",
                                                 c(pathway.sources, "GEO annotation"),
                                                 "GEO annotation")
                                 ),
                                 box(
                                     title = "BRCA class labels",
                                     width = NULL,
                                     selectInput("geo_pheno_columns", "Select a class label column (optional)", choices = NULL),
                                     selectInput("subtype_LumA", "LumA", choices = NULL),
                                     selectInput("subtype_LumB", "LumB", choices = NULL),
                                     selectInput("subtype_Basal", "Basal", choices = NULL),
                                     selectInput("subtype_Her2", "Her2", choices = NULL)
                                 ))),
                                 column(width = 2,
                                        conditionalPanel("output.geo_file_uploaded",
                                                         box(background = "red",
                                                             width = NULL,
                                                             actionButton("geo_button", "Start", icon = icon("rocket"))
                                                         )))
                             )
                        ),
                        tabPanel("Predictions Table",
                                 box(width = 2, downloadButton("download_geo_subtypes", "Download")),
                                 dataTableOutput("gssea_result_geo")),
                        tabPanel("Predictions Plot",
                            plotOutput("gssea_result_geo_heatmap")
                        ),
                        tabPanel("Performance (subtypes)",
                            plotOutput("class_error_subt_plot_geo")
                        ),
                        tabPanel("Performance (predictors)",
                                 plotOutput("class_error_pred_plot_geo")
                        ),
                        tabPanel("Confusion matrices",
                                 plotOutput("class_error_cm_plot_geo")
                        )
                    )
            ),

            tabItem(tabName = "CUSTOM",
                    h2("Upload custom data for subtyping"),
                    tabsetPanel(
                        id = "custom_nav",
                        tabPanel("Setup",
                                 fluidRow(
                                     column(width = 2,
                                         box(
                                             width = NULL,
                                             title = "Upload a dataset",
                                             fileInput("custom_file", "File for subtyping"),
                                             selectInput("sep", "Column separator",
                                                         choices = c("tab" = "\t", "semicolon ;" = ";", "comma ," = ","),
                                                         selected = "\t"),
                                             shinysky::shinyalert("custom_error")
                                         )
                                     ),
                                     column(width = 2,
                                     conditionalPanel("output.custom_file_uploaded",
                                     box(
                                        title = "ID mapping",
                                        width = NULL,
                                        selectInput("custom_gene_id_type", "Gene ID type",
                                                    choices = AnnotationDbi::columns(org.Hs.eg.db),
                                                    selected = "SYMBOL")
                                     ))),
                                     column(width = 2,
                                     conditionalPanel("output.custom_file_uploaded",
                                     box(title = "Select reference / gold standard",
                                         width = NULL,
                                         selectInput("custom_predictors", "Predictor",
                                                     c(all.sources),
                                                     "pam50")
                                     ), box(title = "Upload custom class labels",
                                            width = NULL,
                                            fileInput("custom_class_labels", "File with expected class labels (optional)")),
                                     conditionalPanel("output.custom_labels_uploaded",
                                                      box(width = NULL,
                                                          title = "Custom BRCA class labels",
                                                          selectInput("custom_subtype_LumA", "LumA", choices = NULL),
                                                          selectInput("custom_subtype_LumB", "LumB", choices = NULL),
                                                          selectInput("custom_subtype_Basal", "Basal", choices = NULL),
                                                          selectInput("custom_subtype_Her2", "Her2", choices = NULL)
                                                      )))
                                     ),
                                     column(width = 2,
                                     conditionalPanel("output.custom_file_uploaded",
                                     box(background = "red",
                                         width = NULL,
                                         actionButton("custom_button", "Start", icon = icon("rocket"))
                                     )))
                                 )
                        ),
                        tabPanel("Data Table", dataTableOutput("raw_data_custom")),
                        tabPanel("Gene ID mapping", dataTableOutput("raw_data_id_map")),
                        tabPanel("Predictions Table",
                                 box(width = 2, downloadButton("download_custom_subtypes", "Download")),
                                 dataTableOutput("gssea_result_custom")),
                        tabPanel("Predictions Plot",
                                 plotOutput("gssea_result_custom_heatmap")
                        ),
                        tabPanel("Performance (subtypes)",
                                 plotOutput("class_error_subt_plot_custom")
                        ),
                        tabPanel("Performance (predictors)",
                                 plotOutput("class_error_pred_plot_custom")
                        ),
                        tabPanel("Confusion matrices",
                                 plotOutput("class_error_cm_plot_custom")
                        )
                    )
                )
        )
    )
)

