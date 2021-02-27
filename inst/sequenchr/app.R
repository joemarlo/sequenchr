# set ggplot theme
ggplot2::theme_set(ggplot2::theme_minimal(base_size = 16))

# load UI
source('main_ui.R', local = TRUE)

# retrieve pass-through variables
sequence_data <- shiny::getShinyOption("sequence_data")
covariates_data <- shiny::getShinyOption("covariates_data")


# pre-processing ----------------------------------------------------------

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


# app ---------------------------------------------------------------------

shinyApp(
    ui = UI,
    server = function(input, output, session){

        # initialize store for reactive variables
        store <- reactiveValues()

        # summary table
        output$summary_table <- renderText({
            data.frame(
                n_sequences = nrow(sequence_data),
                n_unique_sequences = nrow(dplyr::distinct(as.data.frame(sequence_data))),
                n_periods = ncol(sequence_data)
            ) %>%
                t() %>%
                `rownames<-`(c('n sequences', 'n unique sequences', 'n periods')) %>%
                knitr::kable(digits = 2, format = 'html') %>%
                kableExtra::kable_styling(bootstrap_options = c("hover", "condensed"))
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
              cluster_assignments = cluster_assignments(),
              n_col_facets = input$clustering_slider_facet_ncol
            )
          }

          return(p)
        })

        # render the top 10 most common sequences
        # TODO: should show frequency of sequence somehow
        # output$plotting_plot_common <- renderPlot({
        #
        #     if (isFALSE(input$plotting_check_cluster)){
        #
        #         # plot without clustering
        #         p <- tidy_data%>%
        #             group_by(sequenchr_seq_id) %>%
        #             summarize(seq_collapsed = paste0(value, collapse = 'SE3P'),
        #                       .groups = 'drop') %>%
        #             count(seq_collapsed) %>%
        #             arrange(desc(n)) %>%
        #             slice_head(n = 10) %>%
        #             separate(seq_collapsed, into = paste0('p', 1:ncol(sequence_data)), sep = "SE3P") %>%
        #             mutate(sequenchr_seq_id = row_number()) %>%
        #             pivot_longer(cols = setdiff(colnames(.), c('n', "sequenchr_seq_id"))) %>%
        #             mutate(name = as.numeric(gsub('p', '', name))) %>%
        #             rename(period = name) %>%
        #             ggplot(aes(x = period, y = sequenchr_seq_id, fill = value)) +
        #             geom_tile() +
        #             scale_fill_manual(values = color_mapping) +
        #             scale_y_discrete(labels = NULL, breaks = NULL) +
        #             labs(title = "Top 10 most common sequences",
        #                  x = 'Period',
        #                  y = 'Sequence (ranked by count)',
        #                  fill = NULL)
        #     } else {
        #
        #         # plot with clustering
        #         p <- tidy_data%>%
        #             left_join(data.frame(cluster = factor(sub("  \\|.*", "", cluster_assignments()),
        #                                                   levels = paste0('Cluster ', 1:length(cluster_assignments()))),
        #                                  sequenchr_seq_id = 1:length(cluster_assignments())),
        #                       by = "sequenchr_seq_id") %>%
        #             group_by(sequenchr_seq_id, cluster) %>%
        #             summarize(seq_collapsed = paste0(value, collapse = 'SE3P'),
        #                       .groups = 'drop') %>%
        #             count(cluster, seq_collapsed) %>%
        #             group_by(cluster) %>%
        #             arrange(desc(n)) %>%
        #             slice_head(n = 10) %>%
        #             ungroup() %>%
        #             separate(seq_collapsed, into = paste0('p', 1:ncol(sequence_data)), sep = "SE3P") %>%
        #             mutate(id = row_number()) %>%
        #             pivot_longer(cols = setdiff(colnames(.), c('n', "id", 'cluster'))) %>%
        #             mutate(name = as.numeric(gsub('p', '', name))) %>%
        #             rename(period = name) %>%
        #             ggplot(aes(x = period, y = id, fill = value)) +
        #             geom_tile() +
        #             scale_fill_manual(values = color_mapping) +
        #             scale_y_discrete(labels = NULL, breaks = NULL) +
        #             facet_wrap(~cluster, scales = 'free_y', ncol = input$clustering_slider_facet_ncol) +
        #             labs(title = "Top 10 most common sequences by cluster",
        #                  x = 'Period',
        #                  y = 'Sequence (ranked by count)',
        #                  fill = NULL)
        #
        #     }
        #
        #     return(p)
        # })

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
              cluster_assignments = cluster_assignments(),
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
              cluster_assignments = cluster_assignments(),
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
              ggplot2::ggplot(ggplot2::aes(x = value, group = name)) +
              ggplot2::geom_density()

          } else {
            # plot with clustering
            p <- dplyr::tibble(cluster = cluster_assignments(),
                               sequenchr_seq_id = 1:length(cluster_assignments())) %>%
              dplyr::right_join(tidy_cov_data, by = 'sequenchr_seq_id') %>%
              ggplot2::ggplot(ggplot2::aes(x = value, group = name)) +
              ggplot2::geom_density() +
              ggplot2::facet_wrap( ~ cluster,
                          scales = 'free_y',
                          ncol = input$clustering_slider_facet_ncol)
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

          # cluster the data
          store$cluster_model <- fastcluster::hclust(
            d = as.dist(store$dist_matrix),
            method = input$clustering_select_clustering_method
            )

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
                                 label = 'Calculate separation metrics')
                )
            })

          # add the download button
          output$clustering_button_UI <- renderUI({
            downloadButton(outputId = 'clustering_button_download',
                           label = 'Download cluster assignments')
            })
        })

        # returns the current cluster assignments
        cluster_assignments <- reactive({

            # stop here if clustering hasn't been run yet
            validate(need(is(store$cluster_model, 'hclust'),
                          'Cluster the data first'))

            # get the cluster assignments
            hcl_k <- input$clustering_slider_n_clusters
            cluster_assignments <- stats::cutree(store$cluster_model, k = hcl_k)

            # reorder clusters to match dendrogram left to right
            cluster_to_dend_mapping <- dplyr::tibble(cluster = cluster_assignments[store$cluster_model$order]) %>%
              tidyr::nest(-cluster) %>%
              dplyr::mutate(cluster_dend = dplyr::row_number()) %>%
              tidyr::unnest(data) %>%
              dplyr::distinct()
            cluster_sorted <- dplyr::tibble(cluster = cluster_assignments) %>%
              dplyr::left_join(cluster_to_dend_mapping, by = 'cluster') %>%
              dplyr::pull(cluster_dend)

            # add label
            cluster_ns <- base::table(cluster_sorted)
            cluster_names <- factor(
                cluster_sorted,
                labels = paste("Cluster", 1:hcl_k, " | n = ", cluster_ns)
            )

            return(cluster_names)
        })

        # plot the dendrogram
        output$clustering_plot_dendrogram <- renderPlot({

            # stop here if clustering hasn't been run yet
            validate(need(is(store$cluster_model, 'hclust'),
                          'Cluster the data first'))

            # retrieve the current cluster model and k cluster value
            cluster <- store$cluster_model
            k <- input$clustering_slider_n_clusters
            h <- input$clustering_slider_dendrogram_depth

            # plot it
            p <- plot_dendrogram(cluster, k, h)

            return(p)
        })

        # compute and plot silhouette width
        observeEvent(input$clustering_button_separation, {
            # get optimal cluster sizes by calculating silhouette width
            # store$s_width <- NbClust::NbClust(
            #     data = NULL,
            #     diss = as.dist(store$dist_matrix),
            #     distance = NULL,
            #     method = 'ward.D2',
            #     max.nc = input$clustering_slider_separation_range[2],
            #     min.nc = input$clustering_slider_separation_range[1],
            #     index = 'silhouette'
            # )

            store$separation_metrics <- cluster_stats(
                dist_matrix = as.dist(store$dist_matrix),
                cluster_model = store$cluster_model,
                k_min = input$clustering_slider_separation_range[1],
                k_max = input$clustering_slider_separation_range[2]
            )

            # update slider with best k value
            # updateSelectInput(session = session,
            #                   inputId = 'clustering_slider_n_clusters',
            #                   selected = store$s_width$Best.nc[['Number_clusters']])

            # remove and add silhouette plot tab
            removeTab(inputId = 'plotting_tabs',
                      target = 'Separation plot')
            insertTab(
                inputId = 'plotting_tabs',
                target = 'Dendrogram',
                position = 'after',
                select = TRUE,
                tab = tabPanel(
                    title = 'Separation plot',
                    br(),
                    plotOutput(outputId = 'clustering_plot_separation',
                               height = 600)
                )
            )
        })

        # plot the silhouette width
        output$clustering_plot_separation <- renderPlot({

          # stop here if silhouette width hasn't been run yet
          validate(need(is(store$separation_metrics, 'data.frame'),
                        'Calculate the silhouette width first'))

          # plot it
          p <- store$separation_metrics %>%
            dplyr::rename(`CH index` = ch_norm,
                   `Silhouette width` = silhouette_norm) %>%
            tidyr::pivot_longer(cols = c("CH index", "Silhouette width")) %>%
            ggplot2::ggplot(ggplot2::aes(x = k, y = value, group = name, color = name)) +
            ggplot2::geom_line(size = 1.2) +
            ggplot2::labs(title = "Cluster seperation measured by Calinski-Harabasz index and silhouette width",
                 subtitle = 'Optimal clusters: minimum CH, maximum silhouette width',
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
            cluster_assignments() %>%
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

        # # render the d3 plot
        # output$explore_d3_chord <- renderD3({
        #     # r2d3(data = runif(5, 0, input$bar_max),
        #     #      script = file.path('d3_plots', 'chord.js')
        #     r2d3(data = matrix(round(runif(input$bar_max, 1, 10000)),
        #                        ncol = 4, nrow = 4),
        #          script = file.path('d3_plots', 'chord.js'))
        # })

        # transition matrix
        transition_matrix <- reactive({

            # filter the data to the periods specfiied by the input slider
            # add NA filler rows after each group before calculating transition matrix
            # this prevents end of day looping back to beginning of day for next group
            freq_data <- tidy_data%>%
              dplyr::filter(period >= input$plotting_slider_chord[1],
                            period <= input$plotting_slider_chord[2]) %>%
              dplyr::mutate(value = as.character(value)) %>%
              dplyr::group_by(sequenchr_seq_id) %>%
              dplyr::group_split() %>%
              lapply(X = ., FUN = function(df){
                df %>% dplyr::add_row(sequenchr_seq_id = NA, value = NA, period = NA)
              }) %>%
              dplyr::bind_rows()

            # calculate transition matrix
            n <- nrow(freq_data)
            TRATE_mat <- base::table(data.frame(previous = freq_data$value[1:(n-1)],
                                          current = freq_data$value[2:n]))
            TRATE_mat <- TRATE_mat / sum(TRATE_mat)

            # ensure matrix contains all the states (b/c above filters may remove some)
            unique_states <- unique(tidy_data$value) %>% as.vector()
            TRATE_filled <- tidyr::crossing(previous = unique_states, current = unique_states) %>%
                dplyr::left_join(dplyr::as_tibble(TRATE_mat),
                                 by = c('previous', 'current')) %>%
                tidyr::replace_na(list(n = 0)) %>%
                tidyr::pivot_wider(names_from = previous, values_from = n)
            TRATE_filled_mat <- as.matrix(TRATE_filled[, -1])
            rownames(TRATE_filled_mat) <- TRATE_filled[[1]]

            return(list(TRATE_filled, TRATE_filled_mat))
        })

        # render the chord plot
        output$explore_plot_chord <- chorddiag::renderChorddiag({

            # get the transition matrix
            trans_mat <- transition_matrix()
            freq_data <- trans_mat[[1]]
            TRATE_mat <- trans_mat[[2]]

            # create the color vector
            states_included <- base::intersect(names(color_mapping), rownames(TRATE_mat)) #unique(freq_data$value))
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

        # render the transition plot
        output$explore_plot_matrix <- renderPlot({

            # get the transition matrix
            TRATE_mat <- transition_matrix()[[1]]

            # plot it
            p <- plot_transition_matrix(TRATE_mat)

            return(p)
        })

        # on load, update the slider with the correct number of periods based on the data
        updateSliderInput(session = session,
                          inputId = 'plotting_slider_chord',
                          max = ncol(sequence_data),
                          value = c(1, ncol(sequence_data)))

    }
)
