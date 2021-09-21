shinyServer(function(input, output, session) {

  # pre-processing ----------------------------------------------------------

  # set ggplot theme
  ggplot2::theme_set(ggplot2::theme_minimal(base_size = 16))

  # retrieve pass-through variables
  sequence_data <- shiny::getShinyOption("sequence_data")
  covariates_data <- shiny::getShinyOption("covariates_data")

  # establish color mapping
  color_mapping <- viridis::viridis_pal()(length(alphabet(sequence_data)))
  names(color_mapping) <- TraMineR::alphabet(sequence_data)

  # tidy the data
  tidy_data <- tidy_sequence_data(sequence_data)
  if (!is.null(covariates_data)){
    tidy_cov_data <- covariates_data %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(sequenchr_seq_id = dplyr::row_number()) %>%
      tidyr::pivot_longer(cols = -sequenchr_seq_id)
  } else tidy_cov_data <- NULL

  # initialize store for reactive variables
  store <- reactiveValues()


  # summary table -----------------------------------------------------------

  # summary table
  output$summary_table <- renderText({

    row_names <- c("n sequences", "n unique sequences", "n periods")
    row_values <- c(nrow(sequence_data),
                    nrow(dplyr::distinct(as.data.frame(sequence_data))),
                    ncol(sequence_data))
    row_values <- round(row_values, 2)
    html_table <- bootstrap_table(row_names, row_values)

    return(html_table)
  })


  # plotting ----------------------------------------------------------------

  # render the sequence index plot
  output$plotting_plot_sequence <- renderPlot({

    if (isFALSE(input$plotting_check_cluster)){
      p <- plot_sequence_index(seq_def_tidy = tidy_data,
                               color_mapping = color_mapping)
    } else {
      p <- plot_sequence_index(
        seq_def_tidy = tidy_data,
        color_mapping = color_mapping,
        cluster_labels = label_clusters_reactive(),
        n_col_facets = input$clustering_slider_facet_ncol
      )
    }

    return(p)
  })

  # plot of just the legend colors
  output$plotting_plot_legend <- renderPlot({

    # plot it
    p <- plot_legend(color_mapping = color_mapping)

    return(p)
  })

  # state distribution plot
  output$plotting_plot_state <- renderPlot({

    if (isFALSE(input$plotting_check_cluster)){
      p <- plot_state(seq_def_tidy = tidy_data,
                      color_mapping = color_mapping)
    } else {
      p <- plot_state(
        seq_def_tidy = tidy_data,
        color_mapping = color_mapping,
        cluster_labels = label_clusters_reactive(),
        n_col_facets = input$clustering_slider_facet_ncol
      )
    }

    return(p)
  })

  # modal plot
  output$plotting_plot_modal <- renderPlot({

    if (isFALSE(input$plotting_check_cluster)){
      p <- plot_modal(seq_def_tidy = tidy_data,
                      color_mapping = color_mapping)
    } else {
      p <- plot_modal(
        seq_def_tidy = tidy_data,
        color_mapping = color_mapping,
        cluster_labels = label_clusters_reactive(),
        n_col_facets = input$clustering_slider_facet_ncol
      )
    }

    return(p)
  })

  # covariates plot
  output$plotting_plot_covariates <- renderPlot({

    validate(need(is(tidy_cov_data, 'data.frame'),
                  'Covariates data not provided'))

    if (isFALSE(input$plotting_check_cluster)){

      # plot without clustering
      p <- tidy_cov_data %>%
        ggplot2::ggplot(ggplot2::aes(x = value)) +
        ggplot2::geom_density() +
        ggplot2::facet_wrap(~name, scales = 'free')

    } else {
      # plot with clustering
      p <- dplyr::tibble(cluster = label_clusters_reactive(),
                         sequenchr_seq_id = 1:length(label_clusters_reactive())) %>%
        dplyr::right_join(tidy_cov_data, by = 'sequenchr_seq_id') %>%
        ggplot2::ggplot(ggplot2::aes(x = value, group = cluster, color = cluster)) +
        ggplot2::geom_density() +
        ggplot2::facet_wrap(~name, scales = 'free')
    }

    return(p)
  })


  # clustering --------------------------------------------------------------

  # render correct options for substituion cost based on distance method input
  observeEvent(input$clustering_select_distanceMethod, {

    if (input$clustering_select_distanceMethod == 'DHD'){
      updateSelectInput(session = session,
                        inputId = 'clustering_select_substitutionMethod',
                        choices = 'TRATE',
                        selected = 'TRATE')
    } else if (input$clustering_select_distanceMethod == 'OM') {
      updateSelectInput(session = session,
                        inputId = 'clustering_select_substitutionMethod',
                        choices = c('TRATE', 'CONSTANT', 'INDELS', 'INDELSLOG'),
                        selected = 'TRATE')
    }
  })

  # cluster the data
  observeEvent(input$clustering_button_cluster, {

    # message to console
    message(paste0("1/3 Computing distance matrix with method ",
                   input$clustering_select_distanceMethod))

    # compute the distance matrix
    if (input$clustering_select_distanceMethod == 'DHD'){
      store$dist_matrix <- TraMineR::seqdist(
        seqdata = sequence_data,
        method = "DHD"
      )
    } else if (input$clustering_select_distanceMethod == 'HAM'){
      store$dist_matrix <- TraMineR::seqdist(
        seqdata = sequence_data,
        method = "HAM"
      )
    } else if (input$clustering_select_distanceMethod == 'OM'){

      indel <- input$clustering_slider_indel
      if (indel != 'auto') indel <- as.numeric(indel)

      store$dist_matrix <- TraMineR::seqdist(
        seqdata = sequence_data,
        method = "OM",
        indel = indel,
        sm = input$clustering_select_substitutionMethod
      )
    }

    # message to console
    message("2/3 Clustering the data")

    # cluster the data
    store$cluster_model <- hclust(
      d = as.dist(store$dist_matrix),
      method = input$clustering_select_clustering_method
    )

    # message to console
    message("3/3 Clustering finished")

    # render the clustering UI
    output$clustering_UI <- renderUI({
      tagList(
        br(),
        sliderInput(inputId = 'clustering_slider_n_clusters',
                    label = 'Number of clusters',
                    min = 2,
                    max = 20,
                    step = 1,
                    value = 1,
                    ticks = FALSE),
        br(),
        actionButton(inputId = 'clustering_button_separation',
                     label = 'Calculate validity statistics')
      )
    })

    # add the download button
    output$clustering_button_UI <- renderUI({
      downloadButton(outputId = 'clustering_button_download',
                     label = 'Download cluster assignments')
    })

    # remove and add dendrogram tab
    removeTab(inputId = 'plotting_tabs',
              target = 'Dendrogram')
    insertTab(
      inputId = 'plotting_tabs',
      target = 'Legend',
      position = 'after',
      select = TRUE,
      tab = tabPanel(
        title = 'Dendrogram',
        br(),
        plotOutput(outputId = 'clustering_plot_dendrogram',
                   height = 600)
      )
    )
  })

  # returns the current cluster assignments
  label_clusters_reactive <- reactive({

    # stop here if clustering hasn't been run yet
    validate(need(is(store$cluster_model, 'hclust'),
                  'Cluster the data first'))

    # set k number of clusters
    k <- input$clustering_slider_n_clusters
    if (is.null(k)) k <- 2

    # get cluster cluster assignments
    cluster_labels <- label_clusters(.model = store$cluster_model, k = k)

    return(cluster_labels)
  })

  # plot the dendrogram
  output$clustering_plot_dendrogram <- renderPlot({

    # stop here if clustering hasn't been run yet
    validate(need(is(store$cluster_model, 'hclust'),
                  'Cluster the data first'))

    # default hcl_k to 2
    hcl_k <- input$clustering_slider_n_clusters
    if (is.null(hcl_k)) hcl_k <- 2

    # plot it
    p <- plot_dendrogram(cluster_model = store$cluster_model,
                         k = hcl_k,
                         h = input$clustering_slider_dendrogram_depth)

    return(p)
  })

  # compute and plot silhouette width
  observeEvent(input$clustering_button_separation, {

    # message to console
    message('1/2 Calculating cluster validity statistics')

    # compute cluster stats
    store$separation_metrics <- cluster_stats(
      dist_matrix = as.dist(store$dist_matrix),
      cluster_model = store$cluster_model,
      k_min = input$clustering_slider_separation_range[1],
      k_max = input$clustering_slider_separation_range[2]
    )

    # remove and add silhouette plot tab
    removeTab(inputId = 'plotting_tabs',
              target = 'Cluster validity')
    insertTab(
      inputId = 'plotting_tabs',
      target = 'Dendrogram',
      position = 'after',
      select = TRUE,
      tab = tabPanel(
        title = 'Cluster validity',
        br(),
        plotOutput(outputId = 'clustering_plot_separation',
                   height = 600)
      )
    )

    # message to console
    message('2/2 Cluster validity finished')
  })

  # plot the silhouette width
  output$clustering_plot_separation <- renderPlot({

    # stop here if silhouette width hasn't been run yet
    validate(need(is(store$separation_metrics, 'data.frame'),
                  'Calculate the silhouette width first'))

    # set breaks
    x_brks <- seq(isolate(input$clustering_slider_separation_range[[1]]),
                isolate(input$clustering_slider_separation_range[[2]]),
                by = 1)

    # plot it
    p <- store$separation_metrics %>%
      dplyr::rename(`Calinski and Harabasz index` = ch_norm,
                    `Silhouette width` = silhouette_norm) %>%
      tidyr::pivot_longer(cols = c("Calinski and Harabasz index", "Silhouette width")) %>%
      ggplot2::ggplot(ggplot2::aes(x = k, y = value, group = name, color = name, fill = name)) +
      ggplot2::geom_line(size = 1.2) +
      ggplot2::geom_point(size = 5, pch = 21, color = 'white', stroke = 1.3) +
      ggplot2::scale_x_continuous(breaks = x_brks) +
      ggplot2::scale_fill_discrete(name = NULL) +
      ggplot2::labs(title = "Cluster validity statistics",
                    subtitle = 'Maximum values == optimal number of clusters',
                    x = 'n clusters',
                    y = 'Normalized index',
                    color = NULL) +
      ggplot2::theme(legend.position = 'bottom')

    return(p)
  })

  # download the clusters
  output$clustering_button_download <- downloadHandler(

    # use plot title as file name but only retain alpha-numeric characters
    filename <- function() {
      time <- gsub("-|:| ", "", Sys.time())
      paste0(time, '_cluster_assignments.csv')
    },

    # dataframe of clusters to download
    content <- function(file) {
      label_clusters_reactive() %>%
        as.data.frame() %>%
        dplyr::mutate(row = dplyr::row_number(),
                      cluster = sub("  \\|.*", "", `.`)) %>%
        dplyr::select(-`.`) %>%
        write.csv(., file, row.names = FALSE)
    }
  )

  # update the max number of rows to plot based on the number of clusters
  observeEvent(input$clustering_slider_n_clusters, {
    updateSliderInput(session = session,
                      inputId = 'clustering_slider_facet_ncol',
                      max = input$clustering_slider_n_clusters)
  })


  # chord plot --------------------------------------------------------------

  # transition matrix
  transition_matrix_reactive <- reactive({

    if (isFALSE(input$plotting_check_cluster)){
      TRATE_tidy <- transition_matrix(
        tidy_data,
        period_min = input$plotting_slider_chord[1],
        period_max = input$plotting_slider_chord[2],
        cluster_labels = NULL)
    } else {
      TRATE_tidy <- transition_matrix(
        tidy_data,
        period_min = input$plotting_slider_chord[1],
        period_max = input$plotting_slider_chord[2],
        cluster_labels = label_clusters_reactive())
    }

    return(TRATE_tidy)
  })

  # render the transition plot
  output$explore_plot_matrix <- renderPlot({

    # get the transition matrix
    TRATE_mat <- transition_matrix_reactive()

    # plot it
    p <- plot_transition_matrix(
      TRATE_mat,
      n_col_facets = input$clustering_slider_facet_ncol
      )

    return(p)
  })

  # on load, update the slider with the correct number of periods based on the data
  updateSliderInput(session = session,
                    inputId = 'plotting_slider_chord',
                    max = ncol(sequence_data),
                    value = c(1, ncol(sequence_data)))

  # only render the chord plot if package is installed; otherwise render a message
  if ('chorddiag' %in% rownames(installed.packages())){

    # render the chord plot
    output$explore_plot_chord <- chorddiag::renderChorddiag({

      validate(need(isFALSE(input$plotting_check_cluster),
                    "Chord plot currently only supported without clustering"))

      # get the transition matrix
      trans_tidy <- transition_matrix_reactive()

      # convert from tidy df to matrix
      TRATE_wide <- dplyr::select(trans_tidy, -cluster)
      TRATE_wide <- tidyr::pivot_wider(TRATE_wide, values_from = 'n', names_from = 'previous')
      TRATE_mat <- as.matrix(TRATE_wide[,-1])
      rownames(TRATE_mat) <- TRATE_wide$current

      # create the color vector
      states_included <- base::intersect(names(color_mapping), rownames(TRATE_mat))
      colors_chord <- as.vector(color_mapping[states_included])

      # plot the chord diagram
      p <- chorddiag::chorddiag(
        data = TRATE_mat,
        groupColors = colors_chord,
        groupnamePadding = 20,
        groupnameFontsize = 12,
        precision = 4
      )

      return(p)
    })

    output$explore_chord_UI <- renderUI({
      tagList(
        h5("The outside arc (node) represents the frequency of the states"),
        h5("The connection between the nodes is the transition rates between the two states"),
        br(),
        chorddiag::chorddiagOutput(
          outputId = 'explore_plot_chord',
          height = 700
        )
      )
    })
  } else {
    output$explore_chord_UI <- renderUI({
      HTML("chorddiag package not detected. Install via devtools::install_github('mattflor/chorddiag')")
    })
  }

})
