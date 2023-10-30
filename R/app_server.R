#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny plotly data.table ggplot2
#' @noRd
app_server <- function(input, output, session) {
  # Your application server logic

  # reactive values; set default
  v <- reactiveValues(data=NULL)
                        #list(DataStore(points = data.table::fread("/Users/xx20081/git/SlopeExplorer/inst/extdata/example_1.csv"))))

  # listen for changing data points things (threshold change or source change)
  reassess_data_points <- reactive({
    list(input$xp_thresh, input$data_source)
  })

  # process changing data points things
  observeEvent(reassess_data_points(), {
    cat("Calling observeEvent(input$data_source)\n")
    if(input$data_source=="file") {
      # get from file
    } else if (input$data_source=="example_1") {
      point_data <- data.table::fread("/Users/xx20081/git/SlopeExplorer/inst/extdata/example_1.csv")
    } else if (input$data_source=="composite_1") {
      point_data <- data.table::fread("/Users/xx20081/git/SlopeExplorer/inst/extdata/composite_1_autosomes.gz")
    }
    # filter on threshold, dont keep all of the points
    xp_thresh <- input$xp_thresh
    point_data <- point_data[i.P <= xp_thresh, ]
    # add/update to store
    v$data <- list(DataStore(points = point_data))
    updateSliderInput(session, "iter", value=0, min=0, max=length(v$data)-1, step=1)
  })


  surface_grid <- reactive({
    cat("Calling surface_grid <- reactive\n")
    points <- v$data[[1]][["points"]]
    lim = c(max(abs(c(points[["i.BETA"]], points[["p.BETA"]])), na.rm=TRUE))
    n_points <- 500
    surface_grid <- expand.grid(x = seq(-lim, lim, length.out=n_points), y = seq(-lim, lim, length.out=n_points))
    surface_grid <- data.table::as.data.table(surface_grid)
    return(surface_grid)
  })


  # do the E-step of the EM algorithm
  # function to work out the normalised z0 and z1 values given a set of distribution parameters
  point_data <- reactive({
    cat("Calling point_data <- reactive\n")
    # distribution 1 - Gi / SNPs associated with incidence
    ux0  <- input$ux0
    uy0  <- input$uy0
    sx0  <- input$sx0
    sy0  <- input$sy0
    sxy0 <- input$sxy0
    dir0 <- input$dir0
    pi0  <- input$pi0
    # distribution 2 - Gip / SNPs associated with incidence and progression
    ux1  <- input$ux1
    uy1  <- input$uy1
    sx1  <- input$sx1
    sy1  <- input$sy1
    sxy1 <- input$sxy1
    dir1 <- input$dir1
    # iter step
    iter <- input$iter
    # the points
    point_data <- v$data[[iter+1]][["points"]]
    # mean vectors
    u0 <- c(ux0, uy0)
    u1 <- c(ux1, uy1)
    # cov matrices
    if(iter==0) {
      sxy0 = sxy0 # allow playing with the controls if on the first iter; thereafter display SlopeHunter method
    } else {
      sxy0 = sx0*sy0*dir0*0.95 # the x & y perfectly correlated under 1st component  #===========
    }
    sigma0 = matrix(c(sx0^2,sxy0,sxy0,sy0^2), 2, 2)
    sigma1 = matrix(c(sx1^2,sxy1,sxy1,sy1^2), 2, 2)

    point_data <- gen_point_z_data(point_data,u0,u1,sigma0,sigma1,pi0)
  })

  gen_point_z_data <- function(point_data, mean_vec1, mean_vec2, sigma_mat1, sigma_mat2, pi0) {
    cat("Calling gen_point_z_data()\n")
    # 1st component
    point_data[, f0 := mclust::dmvnorm(point_data[, c("i.BETA", "p.BETA")], mean=mean_vec1, sigma=sigma_mat1)]
    point_data[, f0 := ifelse(f0<1e-300, 1e-300, f0)]
    # 2nd component
    point_data[, f1 := mclust::dmvnorm(point_data[, c("i.BETA", "p.BETA")], mean=mean_vec2, sigma=sigma_mat2)]
    point_data[, f1 := ifelse(f1<1e-300, 1e-300, f1)]
    ## proportional contribution of density of f0 (for every point) to the total mixture
    point_data[, pt := pi0*f0/(pi0*f0+(1-pi0)*f1)]
    point_data[, pt := ifelse(pt>0.9999999, 0.9999999, pt)]
    return(point_data)
  }

  # do the E-step of the EM algorithm
  # function to work out the normalised z0 and z1 values given a set of distribution parameters
  surface_data <- reactive({
    cat("Calling surface_data <- reactive\n")
    # distribution 1 - Gi / SNPs associated with incidence
    ux0  <- input$ux0
    uy0  <- input$uy0
    sx0  <- input$sx0
    sy0  <- input$sy0
    sxy0 <- input$sxy0
    dir0 <- input$dir0
    pi0  <- input$pi0
    # distribution 2 - Gip / SNPs associated with incidence and progression
    ux1  <- input$ux1
    uy1  <- input$uy1
    sx1  <- input$sx1
    sy1  <- input$sy1
    sxy1 <- input$sxy1
    dir1 <- input$dir1
    # surface grid data
    point_data <- surface_grid()
    # mean vectors
    u0 <- c(ux0, uy0)
    u1 <- c(ux1, uy1)
    # cov matrices
    sigma0 <- matrix(c(sx0^2,sxy0,sxy0,sy0^2), 2, 2)
    sigma1 <- matrix(c(sx1^2,sxy1,sxy1,sy1^2), 2, 2)
    # create the Z data
    point_data[, z0 := mclust::dmvnorm(point_data[, c("x","y")], mean=u0, sigma=sigma0)]
    point_data[, z1 := mclust::dmvnorm(point_data[, c("x","y")], mean=u1, sigma=sigma1)]
    # normalise 0-1
    point_data[, z0 := (z0-min(z0)) / (max(z0)-min(z0))]
    point_data[, z1 := (z1-min(z1)) / (max(z1)-min(z1))]

    # matrices
    dim_x <- length(unique(point_data[["x"]]))
    z0_mat <- matrix(point_data[["z0"]], nrow=dim_x, byrow=TRUE)
    z1_mat <- matrix(point_data[["z1"]], nrow=dim_x, byrow=TRUE)
    rownames(z0_mat) = rownames(z1_mat) = unique(point_data[["x"]])
    colnames(z0_mat) = colnames(z1_mat) = unique(point_data[["y"]])

    # return a list of matrices
    return(list(z0 = z0_mat, z1 = z1_mat))
  })


  # the main display output
  output$main_plot <- renderPlotly({
    cat("Calling main_plot()\n")

    # the iter value
    iter <- input$iter
    cat("Plotting:", iter, "\n")

    # the iter data
    iter_data <- v$data[[iter+1]]

    # update the GUI if going through the algo iterations
    if(iter>0) {
      isolate(updateNumericInput(session, "ux0",  value = 0))
      isolate(updateNumericInput(session, "uy0",  value = 0))
      isolate(updateNumericInput(session, "ux1",  value = 0))
      isolate(updateNumericInput(session, "uy1",  value = 0))
      isolate(updateNumericInput(session, "sx0",  value = iter_data[["sx0"]]))
      isolate(updateNumericInput(session, "sy0",  value = iter_data[["sy0"]]))
      isolate(updateNumericInput(session, "sx1",  value = iter_data[["sx1"]]))
      isolate(updateNumericInput(session, "sy1",  value = iter_data[["sy1"]]))
      isolate(updateNumericInput(session, "sxy0", value = iter_data[["sxy0"]]))
      isolate(updateNumericInput(session, "sxy1", value = iter_data[["sxy1"]]))
      isolate(updateNumericInput(session, "dir0", value = iter_data[["dir0"]]))
      isolate(updateNumericInput(session, "pi0",  value = iter_data[["pi0"]]))
    }

    # the points data
    point_data <- v$data[[1]][["points"]] # only keep one copy of the points coordinates
    if(is.null(iter_data[["pt_data"]])) { # if null pt data (i.e. not run SH yet), then generate the pt data from the GUI controls
      d <- point_data()
      point_data[, pt := d[["pt"]]]
    } else { # else use the pt data stored when we ran the SH algo.
      point_data[, pt := iter_data[["pt_data"]]]
    }

    # create the Hunted factor
    point_data[, clusters := factor(ifelse(pt>=0.5, "Hunted", "Pleiotropic"), levels=c("Hunted", "Pleiotropic"))]

    # get the surface data; will use the updated GUI values
    surface_data <- surface_data()

    # downsample Pleiotropic SNPs for plotting
    set.seed(123)
    point_data <- point_data[clusters == "Hunted" | runif(.N) < 1000 / .N, ]

    # the b line
    if(all(!is.null(iter_data[["dir0"]]) & !is.null(iter_data[["sy0"]]) & !is.null(iter_data[["sx0"]]))) {
      b = iter_data[["dir0"]]*iter_data[["sy0"]]/iter_data[["sx0"]]
      x_line <- seq(min(point_data[["i.BETA"]]), max(point_data[["i.BETA"]]), length.out = 100)
      line_data <- data.table::data.table(x = x_line, y = b * x_line)
    } else {
      line_data <- NULL
    }

    # create the plot
    p <- plot_base(point_data, surface_data, line_data)

    if(input$surface_type=="3d") {
      p <- plot_3d(p, point_data, surface_data, line_data)
    } else {
      p <- plot_2d(p, point_data, surface_data, line_data)

      # # ggplot doesnt want matrices
      # g <- expand.grid(x = as.numeric(rownames(surface_data[["z0"]])), y = as.numeric(colnames(surface_data[["z0"]])))
      # g[["z0"]] <- as.vector(t(surface_data[["z0"]]))
      # g[["z1"]] <- as.vector(t(surface_data[["z1"]]))
      # p <- plot_2d(point_data, g, line_data)
    }

    # show the plot
    p
  })

  plot_base <- function(point_data, surface_data=NULL, line_data=NULL) {
    cat("Calling plot_base()\n")

    gg <- plot_ly() |>

      # plot layout options
      layout(
        plot_bgcolor = "white",
        scene = list(
          xaxis = list(title = "\u03B2 incidence"),
          yaxis = list(title = "\u03B2 progression"),
          zaxis = list(title = "Normalised density"))
      )



    return(gg)
  }

  # produce the 3D plot
  plot_3d <- function(gg, point_data, surface_data=NULL, line_data=NULL) {
    cat("Calling plot_3d()\n")

    # the points
    gg <- gg |>
      add_markers(data       = point_data,
                  x          = ~i.BETA,
                  y          = ~p.BETA,
                  z          = ~pt,
                  mode       = "markers",
                  color      = ~clusters,
                  colors     = c("red", "grey"),
                  marker     = list(size = 5))

    # the line
    if(!is.null(line_data)) {
      gg <- gg |>
        add_trace(type = "scatter3d",
                  data = line_data,
                  x = ~x,
                  y = ~y,
                  z = rep(0, nrow(line_data)),
                  mode = "lines",
                  line = list(color = "black", width = 2),
                  name = "SlopeHunter adjustment factor")
    }

    # surface 1
    if(!is.null(surface_data)) {
      gg <- gg |>

        add_surface(x          = as.numeric(rownames(surface_data[["z0"]])),
                    y          = as.numeric(colnames(surface_data[["z0"]])),
                    z          = surface_data[["z0"]],
                    opacity    = 0.5,
                    colorscale = list(c(0, 1), c("lightblue", "blue")),
                    showscale  = FALSE,
                    showlegend = TRUE,
                    name       = "Incidence only\nSNP distribution")
    }

    # surface 2
    if(!is.null(surface_data)) {
      gg <- gg |>
        add_surface(x          = as.numeric(rownames(surface_data[["z1"]])),
                    y          = as.numeric(colnames(surface_data[["z1"]])),
                    z          = surface_data[["z1"]],
                    opacity    = 0.5,
                    colorscale = list(c(0, 1), c("lightgrey", "darkgrey")),
                    showscale  = FALSE,
                    showlegend = TRUE,
                    name       = "Incidence & progression\n (Pleiotropic) SNP distribution")
    }

    # return plot
    return(gg)
  }


  # produce the 2D plot
  plot_2d <- function(gg, point_data, surface_data=NULL, line_data=NULL) {
    cat("Calling plot_2d()\n")

    # the points
    gg <- gg |>
      add_markers(data       = point_data,
                  x          = ~i.BETA,
                  y          = ~p.BETA,
                  type       = "scatter",
                  mode       = "markers",
                  color      = ~clusters,
                  colors     = c("red", "grey"),
                  marker     = list(size = 5))

    # the line
    if(!is.null(line_data)) {
      gg <- gg |>
        add_trace(type = "scatter",
                  data = line_data,
                  x = ~x,
                  y = ~y,
                  mode = "lines",
                  line = list(color = "black", width = 2),
                  name = "SlopeHunter adjustment factor")
    }

    # contour 1
    if(!is.null(surface_data)) {
      gg <- gg |>
        add_contour(x          = as.numeric(rownames(surface_data[["z0"]])),
                    y          = as.numeric(colnames(surface_data[["z0"]])),
                    z          = surface_data[["z0"]],
                    line       = list(width = 2, color = "blue"),
                    opacity    = 0.5,
                    colorscale = list(c(0, 1), c("rgba(0,0,0,0.0)", "rgba(0,0,0,0.0)")),
                    showscale  = FALSE,
                    showlegend = TRUE,
                    name       = "Incidence only\nSNP distribution")

      # contour 2
      gg <- gg |>
        add_contour(x          = as.numeric(rownames(surface_data[["z1"]])),
                    y          = as.numeric(colnames(surface_data[["z1"]])),
                    z          = surface_data[["z1"]],
                    line       = list(width = 2, color = "lightgrey"),
                    opacity    = 0.5,
                    colorscale = list(c(0, 1), c("rgba(0,0,0,0.0)", "rgba(0,0,0,0.0)")),
                    showscale  = FALSE,
                    showlegend = TRUE,
                    name       = "Incidence & progression\n (Pleiotropic) SNP distribution")
    }

    # return plot
    return(gg)
  }




  observeEvent(input$run_button, {
    cat("Running SlopeHunter...\n")

    withProgress(message = 'Running SlopeHunter', value = 0, {

      # get the initialising parameters
      xp_thresh <- input$xp_thresh
      pi0       <- input$pi0
      sxy1      <- input$sxy1_slope_init
      n         <- 50 # guess number of iterations needed

      # get the point data
      point_data <- v$data[[1]][["points"]]

      # filter by the xp_threshold
      point_data <- point_data[i.P <= xp_thresh, ]

      # if there are significant points in the incidence GWAS, continue
      if(nrow(point_data)==0) {
        cat("No significant incidence SNPs at xp_threshold", xp_thresh, "\n")
        return()
      }

      # set the initial sdev and cov
      sx0 = sx1 = sd( point_data[["i.BETA"]] )
      sy0 = sy1 = sd( point_data[["p.BETA"]] )
      dir0 = sign(  cov(point_data[["i.BETA"]], point_data[["p.BETA"]]) )
      if (dir0==0) stop("All associations with at least either x or y are constant")

      # convergence criterion
      loglkl_ck = 0

      # data store
      data_store_iters = list()

      # create the baseline data situation and add as iter=0
      iter0 = DataStore(
        points = point_data,
        pt_data= rep(0, nrow(point_data)),
        sx0    = sx0,
        sy0    = sy0,
        dir0   = dir0,
        sxy0   = sx0*sy0*dir0*0.95,
        sx1    = sx1,
        sy1    = sy1,
        sxy1   = sxy1,
        pi0    = pi0
      )
      data_store_iters <- append(data_store_iters, list(iter0))

      ### EM algorithm
      for(iter in 1:50000){
        cat("Iter step:", iter, "\n")
        incProgress(1/n, detail = paste("Iteration ", iter))

        #### The E step:
        # covariance matrix for the target component (f0)
        sxy0 = sx0*sy0*dir0*0.95       # the x & y perfectly correlated under 1st component  #===========
        sigma0 = matrix(c(sx0^2,sxy0,sxy0,sy0^2), 2, 2)
        # covariance matrix for the component (f1)
        sigma1 = matrix(c(sx1^2,sxy1,sxy1,sy1^2), 2, 2)

        # get the proportional contribution of density of f0 (for every point) to the total mixture
        point_data <- gen_point_z_data(point_data,c(0,0),c(0,0),sigma0,sigma1,pi0)

        # loglik of the mixture model: pi0 * f0 + (1-p0) * f1
        loglkl = sum(log(pi0*point_data[["f0"]]+(1-pi0)*point_data[["f1"]]))

        # add some iteractions to the store, not all as potentially loads of memory
        if(iter %in% 1:20 | iter%%20==0 | (loglkl - loglkl_ck)/loglkl < 1e-10) {
          iter_n = DataStore(
            pt_data= point_data[['pt']],
            sx0    = sx0,
            sy0    = sy0,
            dir0   = dir0,
            sxy0   = sxy0,
            sx1    = sx1,
            sy1    = sy1,
            sxy1   = sxy1,
            pi0    = pi0
          )
          data_store_iters <- append(data_store_iters, list(iter_n))
        }

        #### The M step:
        # update pi0
        pi0 = mean(point_data[['pt']])
        if (pi0<0.0001) pi0 = 0.0001
        if (pi0>0.9999) pi0 = 0.9999

        # update sx0 & sy0
        sx0  = sqrt(sum(point_data[['pt']]*(point_data[["i.BETA"]]^2))/sum(point_data[['pt']]))
        sy0  = sqrt(sum(point_data[['pt']]*(point_data[["p.BETA"]]^2))/sum(point_data[['pt']]))
        dir0 = sign(sum(point_data[['pt']]*point_data[["p.BETA"]]*point_data[["i.BETA"]])/sum(point_data[['pt']]) )
        if (dir0==0) dir0=sample(c(1,-1), 1)   # avoid slope = 0 (horizontal line)

        # update sx1, sy1 & sxy1
        sx1  = sqrt(sum((1-point_data[['pt']])*(point_data[["i.BETA"]]^2))/(length(point_data[["i.BETA"]])-sum(point_data[['pt']])))
        sy1  = sqrt(sum((1-point_data[['pt']])*(point_data[["p.BETA"]]^2))/(length(point_data[["p.BETA"]])-sum(point_data[['pt']])))
        sxy1 = sum((1-point_data[['pt']])*point_data[["i.BETA"]]*point_data[["p.BETA"]])/(length(point_data[["i.BETA"]])-sum(point_data[['pt']]))
        if (abs(sxy1) > 0.75*sx1*sy1)  sxy1 = sign(sxy1)*0.75*sx1*sy1     #===========

        ## Check convergence
        if (iter%%10==0){
          if ((loglkl - loglkl_ck)/loglkl < 1e-10){
            break
          } else {
            loglkl_ck = loglkl
          }
        }
      }

      # Diagnosis
      if (iter == 50000) warning("The algorithm may not have converged.\n")

      # add list to reactiveValues data
      v$data <- data_store_iters

      # set the max iter of the slider
      updateSliderInput(session, "iter", value=1, min=0, max=length(data_store_iters)-1, step=1)

      cat("Finished SlopeHunter run\n")

    })
  })




  # else {
  #
  #   # the iter value
  #   iter <- input$iter
  #   cat("Plotting:", iter, "\n")
  #
  #   # the points data - only stored in the first iter
  #   points_data <- v$data[[1]][["points"]]
  #
  #   # the iter data4345
  #   # structure:
  #   # List of 10
  #   # $ points :Classes ‘data.table’ and 'data.frame':	100 obs. of  7 variables:
  #   #   ..$ id    : int [1:100] 1 2 3 4 5 6 7 8 9 10 ...
  #   # ..$ i.BETA: num [1:100] -0.9164 -1.0962 -0.1146 -0.0183 0.5002 ...
  #   # ..$ i.SE  : num [1:100] 0.329 -0.337 0.348 -0.411 -0.308 ...
  #   # ..$ i.P   : num [1:100] 0.04316 0.11304 0.10182 0.17469 0.00513 ...
  #   # ..$ p.BETA: num [1:100] -0.348 0.748 0.7 2.267 0.337 ...
  #   # ..$ p.SE  : num [1:100] 0.0569 -0.251 0.7399 0.0752 0.528 ...
  #   # ..$ p.P   : num [1:100] 0.10342 0.00765 0.06166 0.17929 0.00187 ...
  #   # ..- attr(*, ".internal.selfref")=<externalptr>
  #   # ...
  #   # $ pt_data: num 0.111
  #   # $ sx0    : num 135
  #   # $ sy0    : num 1.09
  #   # $ sxy0   : num -1.04
  #   # $ dir0   : num -1
  #   # $ sx1    : num 1
  #   # $ sy1    : num 1.09
  #   # $ sxy1   : num 1e-05
  #   # $ pi0    : num 0.6
  #   # - attr(*, "class")= chr "data_store"
  #   iter_data <- v$data[[iter+1]]
  #
  #   # update the distribution parameters
  #   # means are always zero
  #   isolate(updateNumericInput(session, "Gi_i_mean",  value = 0))
  #   isolate(updateNumericInput(session, "Gi_p_mean",  value = 0))
  #   isolate(updateNumericInput(session, "Gip_i_mean", value = 0))
  #   isolate(updateNumericInput(session, "Gip_p_mean", value = 0))
  #   # std deviations
  #   isolate(updateNumericInput(session, "Gi_i_sdev",  value = iter_data[["sx0"]]))
  #   isolate(updateNumericInput(session, "Gi_p_sdev",  value = iter_data[["sy0"]]))
  #   isolate(updateNumericInput(session, "Gip_i_sdev", value = iter_data[["sx1"]]))
  #   isolate(updateNumericInput(session, "Gip_p_sdev", value = iter_data[["sx1"]]))
  #   # covariances
  #   isolate(updateNumericInput(session, "Gi_sxy0",    value = iter_data[["sxy0"]]))
  #   isolate(updateNumericInput(session, "Gi_dir0",    value = iter_data[["dir0"]]))
  #   isolate(updateNumericInput(session, "Gip_sxy1",   value = iter_data[["sxy1"]]))
  #   # pi
  #   isolate(updateNumericInput(session, "pi0",        value = iter_data[["pi0"]]))
  #
  #
  #   # the grid for the surfaces
  #   grid_axes_data <- surface_grid_axes()
  #
  #   # the matrix of z (normalised distribution density) for the surfaces
  #   i_surface_data  <- matrix(surface_z_data()[["z_Gi" ]], nrow=length(unique(grid_axes_data[["x"]])), byrow=TRUE)
  #   ip_surface_data <- matrix(surface_z_data()[["z_Gip"]], nrow=length(unique(grid_axes_data[["x"]])), byrow=TRUE)
  #
  #   # the b line
  #   b = iter_data[["dir0"]]*iter_data[["sy0"]]/iter_data[["sx0"]]
  #   isolate(updateNumericInput(session, "b_value", value = b))
  #   x_line <- seq(min(points_data[["i.BETA"]]), max(points_data[["i.BETA"]]), length.out = 100)
  #   y_line <- b * x_line
  #
  #   # create the basic plot
  #   gg <- plot_ly() |>
  #     add_markers(data       = points_data,
  #                 x          = ~i.BETA,
  #                 y          = ~p.BETA,
  #                 z          = iter_data[["pt_data"]],
  #                 mode       = "markers",
  #                 color      = factor(ifelse(iter_data[["pt_data"]]>=0.5,"Hunted", "Pleiotropic"), levels=c("Hunted", "Pleiotropic")),
  #                 colors     = c("red", "grey"),
  #                 marker     = list(size = 5)) |>
  #     add_trace(type = "scatter3d",
  #               x = x_line,
  #               y = y_line,
  #               z = rep(0, length(x_line)),
  #               mode = "lines",
  #               line = list(color = "black", width = 2),
  #               text = paste0("b = ", round(b, digits=2)),
  #               name = "SlopeHunter adjustment factor"
  #     ) |>
  #     add_surface(data       = grid_axes_data,
  #                 x          = ~x,
  #                 y          = ~y,
  #                 z          = i_surface_data,
  #                 opacity    = 0.3,
  #                 colorscale = list(c(0, 1), c("lightblue", "blue")),
  #                 showscale  = FALSE,
  #                 showlegend = TRUE,
  #                 name       = "Incidence only\nSNP distribution") |>
  #     add_surface(data       = grid_axes_data,
  #                 x          = ~x,
  #                 y          = ~y,
  #                 z          = ip_surface_data,
  #                 opacity    = 0.3,
  #                 colorscale = list(c(0, 1), c("lightgrey", "darkgrey")),
  #                 showscale  = FALSE,
  #                 showlegend = TRUE,
  #                 name       = "Incidence & progression\n (Pleiotropic) SNP distribution") |>
  #     layout(scene = list(
  #       xaxis = list(title = "\u03B2 incidence"),
  #       yaxis = list(title = "\u03B2 progression"),
  #       zaxis = list(title = "Normalised density")
  #     ))
  #
  #   # return / show
  #   gg
  #
  # }
  #
  #   })



}
