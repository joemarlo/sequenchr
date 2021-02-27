UI <- fluidPage(
  
  # download fonts
  HTML('<link rel="preconnect" href="https://fonts.gstatic.com">'),
  HTML('<link rel="stylesheet" href="//fonts.googleapis.com/css?family=Roboto:400,300,700,400italic">'),
  HTML('<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Libre+Barcode+128+Text&display=swap">'),

  # load custom CSS file
  includeCSS(file.path("www", "custom_css.css")),
  
  # set top left title
  titlePanel(
    title = "",
    windowTitle = "sequenchr"
  ),
  
  # main UI
  sidebarLayout(
    sidebarPanel(
      width = 3,
      HTML("<h2 id='logo' style='font-family: \"Libre Barcode 128 Text\", cursive; font-weight:500; color: #8c8c8c; font-size: 5rem;'>sequenchr</h2>"),
      br(),
      h4("Data summary"),
      htmlOutput(outputId = 'summary_table'),
      checkboxInput(inputId = 'plotting_check_cluster',
                    label = 'Cluster the data?'),
      conditionalPanel(
        condition = 'input.plotting_check_cluster == true',
        selectInput(inputId = 'clustering_select_distanceMethod',
                    label = 'Distance method',
                    choices = c('OM', 'HAM', 'DHD'),
                    selected = 'OM'),
        conditionalPanel(
          condition = 'input.clustering_select_distanceMethod == "OM" || input.clustering_select_distanceMethod == "DHD"',
          selectInput(inputId = 'clustering_select_substitutionMethod',
                      label = 'Substitution method',
                      choices = 'TRATE',
                      selected = 'TRATE')
        ),
        conditionalPanel(
          condition = 'input.clustering_select_distanceMethod == "OM"',
          selectInput(inputId = 'clustering_slider_indel',
                      label = 'Insertion/deletion cost',
                      choices = c('auto', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10'))
          ),
        selectInput(inputId = 'clustering_select_clustering_method',
                    label = 'Clustering method',
                    choices = c('ward.D2', 'complete', 'average'),
                    selected = 'ward.D2'),
        actionButton(inputId = 'clustering_button_cluster',
                     label = 'Cluster the data'),
        uiOutput(outputId = 'clustering_UI')
      ),
      br(), br(),
      HTML('<details><summary>Additional settings</summary>'),
      br(),
      sliderInput(inputId = 'clustering_slider_facet_ncol',
                  label = 'Number of columns for plot facets',
                  min = 1,
                  max = 20,
                  value = 1,
                  ticks = FALSE),
      sliderInput(inputId = 'clustering_slider_dendrogram_depth',
                  label = 'Depth to plot dendrogram. Deeper = slower peformance',
                  min = 0,
                  max = 500,
                  step = 5,
                  value = 50,
                  ticks = FALSE),
      sliderInput(inputId = 'clustering_slider_separation_range',
                  label = 'n clusters to calculate silhouette width at',
                  min = 2,
                  max = 20,
                  step = 1,
                  value = c(2, 10),
                  ticks = FALSE),
      br(),
      uiOutput(outputId = 'clustering_button_UI'),
      br(),
      h4("Citations"),
      p("....lorem ipsum..."),
      HTML('</details><br>')
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        id = 'plotting_tabs',
        type = 'tabs',
        tabPanel(
          title = 'Sequence index plot',
          br(),
          plotOutput(outputId = 'plotting_plot_sequence',
                     height = 750)
        ),
        # tabPanel(
        #   title = 'Top 10',
        #   br(),
        #   plotOutput(outputId = 'plotting_plot_common',
        #              height = 750)
        # ),
        tabPanel(
          title = 'State distribution',
          br(),
          plotOutput(outputId = 'plotting_plot_state',
                     height = 750)
        ),
        tabPanel(
          title = 'Modal states',
          br(),
          plotOutput(outputId = 'plotting_plot_modal',
                     height = 750)
        ),
        tabPanel(
          title = 'Covariates',
          br(),
          plotOutput(outputId = 'plotting_plot_covariates',
                     height = 750)
        ),
        tabPanel(
          title = "Transitions",
          br(),
          h3("The below two plots visualize the frequency of transitions between states"),
          br(),
          sliderInput(inputId = 'plotting_slider_chord',
                      label = 'Periods to include in transition rate calculation:',
                      min = 1,
                      max = 10,
                      value = c(1, 10),
                      step = 1,
                      width = '50%',
                      animate = TRUE),
          tabsetPanel(
            id = 'explore_tabs',
            type = 'tabs',
            tabPanel(
              title = "Chord plot",
              br(),
              # h3("Transition between states"),
              h5("The outside arc (node) represents the frequency of the states"),
              h5("The connection between the nodes is the transition rates between the two states"),
              br(), 
              chorddiag::chorddiagOutput(
                outputId = 'explore_plot_chord',
                height = 700
                )
            ),
            tabPanel(
              title = 'Matrix',
              br(),
              plotOutput(outputId = 'explore_plot_matrix',
                         height = 700)
            )
          )
        ),
        tabPanel(
          title = 'Legend',
          br(),
          plotOutput(outputId = 'plotting_plot_legend',
                     height = 400,
                     width = 320)
        )
      )
    )
  )
)
