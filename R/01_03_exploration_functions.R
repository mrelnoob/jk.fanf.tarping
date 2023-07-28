### Sourcing previous scripts:
source(file = here::here("R/01_02_cleaning_functions.R")) # I need to do that because one of the functions
# below uses one of the functions from the sourced R file.



# _________________________________________
### Creation of a function that draws univariate boxplots for all numeric variables in a given dataset:

#' Univariate boxplots
#'
#' @description The `uni.boxplots` function draws, within a single panel, an independent boxplot for each
#' numeric (continuous or discrete) variable in a given dataset. It is particularly useful for data
#' exploration (e.g. Zuur \emph{et al.}, 2010). For instance, to simultaneously observe the
#' distributions of all numeric variables or \strong{to detect their univariate outliers}.
#'
#' @details The `uni.boxplots` function only modifies the graphical parameters of the
#' \code{\link[graphics:boxplot]{boxplot}} function in the `graphics` package to match some predefined
#' preferences. Therefore, default values of `uni.boxplots` create nice looking boxplots but retain
#' default \emph{heuristic} aspects of `boxplot` (such as the length of whiskers or the plotting of
#' outliers). These aspects can however be changed as in \code{\link[graphics:boxplot]{boxplot}}. \cr
#' On the other hand, panel parameters are internally controlled using `par`. However, to avoid unforeseen
#' conflicts with other internal parameters, it is not possible to tune panel parameters as we would
#' do with `par`. Instead, parametrization is only possible with the given subset of parameters.
#'
#' @note To avoid \emph{recursive argument errors}, internal arguments should be called using upper case
#' letters (e.g. CEX.LAB = 0.9) whereas other arguments from the `boxplot` function should be called with
#' their normal case writing (e.g. outline = FALSE)!
#'
#' @param dataset The input dataset containing all variables to be plotted. It may contain all kinds of
#' variables, the `uni.boxplots` function will automatically detect and plot numeric variables (columns).
#' @param ... Any other parameter that can be incorporated in \code{\link[graphics:boxplot]{boxplot}}.
#' @param MAR A numerical vector of the form `c(bottom, left, top, right)` which gives the number of lines
#' of margin to be specified on the four sides of the plot. The default is `c(0.5,4.1,1.1,1.5)`.
#' @param CEX.LAB The magnification to be used for x and y labels relative to the current setting of `cex`.
#' @param FONT.LAB The font to be used for x and y labels.
#' @param BTY A character string which determined the type of box which is drawn about plots. If `BTY` is
#' one of "o", "l", "7", "c", "u", or "]" the resulting box resembles the corresponding upper
#' case letter. A value of "n" suppresses the box (the default).
#' @param FG The color to be used for the foreground of plots. This is the default color used for things
#' like axes and boxes around plots (defaults to "gray35").
#' @param COL.AXIS The color to be used for axis annotation. Defaults to "gray35".
#' @param COL.LAB The color to be used for x and y labels. Defaults to "gray20".
#' @param CEX.PAR A numerical value giving the amount by which plotting text and symbols should be
#' magnified relative to the default (for `par`, the panel manager). This starts as 1 when a device
#' is opened, and is reset when the layout is changed, e.g. by setting `mfrow`. Defaults to 0.8.
#' @param TCL The length of tick marks as a fraction of the height of a line of text. The default
#' value is -0.3.
#' @param MGP The margin line (in `mex` units) for the axis title, axis labels and axis line.
#' Note that `mgp[1]` affects title whereas `mgp[2:3]` affect axis. The default is c(2.4, 0.6, 0).
#' @param OMA A vector of the form `c(bottom, left, top, right)` giving the size of the outer margins
#' in lines of text.
#' @param TYPE The type of boxplot to draw. Default is "n".
#' @param BORDER An optional vector of colors for the outlines of the boxplots. The values in border
#' are recycled if the length of border is less than the number of plots. Default is "lightcoral".
#' @param COL If col is non-null it is assumed to contain colors to be used to colour the bodies of
#' the boxplots. Default is "moccasin".
#' @param LTY The line type. Line types can either be specified as an integer (0=blank, 1=solid (default),
#' 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash) or as one of the character strings "blank",
#' "solid", "dashed", "dotted", "dotdash", "longdash", or "twodash", where "blank" uses ‘invisible lines’
#' (i.e., does not draw them).
#' @param STAPLEWEX Staple line width expansion, proportional to box width. Default is 0.
#' @param WHISKLWD Whisker line width expansion. Default is 2.
#' @param BOXWEX A scale factor to be applied to all boxes. When there are only a few groups, the
#' appearance of the plot can be improved by making the boxes narrower. Default is 0.7.
#' @param BOXLWD Width of boxplot outer lines. Default is 0.1.
#' @param MEDLWD Width of the median line. Default is 2.6.
#' @param PCH The type of points to be drawn for outliers. Default is 19. See \code{\link[graphics:points]{points}}
#' for possible values and their interpretation.
#'
#' @return A panel of univariate boxplots.
#' @export
#' @import graphics
#'
#' @examples
#' data("mtcars")
#' uni.boxplots(dataset = mtcars)
uni.boxplots <- function(dataset, MAR=c(0.5,4.1,1.1,1.5), CEX.LAB=1, FONT.LAB=2, BTY = "n", FG = "gray35",
                         COL.AXIS = "gray35", COL.LAB = "gray20", CEX.PAR = 0.8, TCL = -0.3,
                         MGP = c(2.4, 0.6, 0), OMA = c(1, 0, 0, 0),
                         TYPE = "n", BORDER = "lightcoral", COL = "moccasin",  LTY = 1, STAPLEWEX = 0,
                         WHISKLWD = 2, BOXWEX = 0.7, BOXLWD = 0.1, MEDLWD = 2.6, PCH = 19, ...){
  num.data <- dataset[, sapply(dataset, is.numeric)]
  nam <- names(num.data)
  ncol.data <- ncol(num.data)
  ncol.adjust <- ceiling(x = ncol.data/4) # Round to the next integer (e.g. ceiling(x = 7.12) returns 8)!

  graphics::par(mfrow= c(ncol.adjust,4), mar=MAR, cex.lab=CEX.LAB, font.lab=FONT.LAB, bty=BTY, fg=FG,
                col.axis=COL.AXIS, col.lab=COL.LAB, cex=CEX.PAR, tcl=TCL, mgp=MGP, oma=OMA)
  for (i in c(1:ncol.data)) {
    graphics::boxplot(num.data[,i],ylab =(nam[i]), type=TYPE, border=BORDER, col = COL,
                      lty=LTY, staplewex=STAPLEWEX, whisklwd=WHISKLWD, boxwex=BOXWEX, boxlwd=BOXLWD,
                      medlwd=MEDLWD, pch=PCH, cex=0.7, ...) }
}





# _________________________________________
### Creation of a function that draws univariate Cleveland dotplots for all numeric variables in a
### given dataset:

#' Univariate Cleveland dotplots
#'
#' @description The `uni.dotplots` function draws, within a single panel, an independent Cleveland dotplot
#' (i.e. plotting value against rank) for each numeric (continuous or discrete) variable in a given
#' dataset. It is particularly useful for data exploration (e.g. Zuur \emph{et al.}, 2010). For instance,
#' to simultaneously observe the distributions of all numeric variables or \strong{to detect their
#' univariate outliers}.
#'
#' @details The `uni.dotplots` function only modifies the graphical parameters of the
#' \code{\link[graphics:plot.default]{plot}} function in the `graphics` package to match some predefined
#' preferences. Therefore, default values of `uni.dotplots` create nice looking dotplots but retain
#' some aspects of the original `plot` function. These aspects can however be changed as in
#' \code{\link[graphics:plot.default]{plot}}. \cr
#' On the other hand, panel parameters are internally controlled using `par`. However, to avoid unforeseen
#' conflicts with other internal parameters, it is not possible to tune panel parameters as we would
#' do with `par`. Instead, parametrization is only possible with the given subset of parameters.
#'
#' @note To avoid \emph{recursive argument errors}, internal arguments should be called using upper case
#' letters (e.g. CEX.LAB = 0.9) whereas other arguments from the `plot` function should be called with
#' their normal case writing (e.g. sub = "My subtitle")!
#'
#' @param dataset The input dataset containing all variables to be plotted (must be a `data.frame` with
#' at least 2 variables). It may contain all kinds of columns, the `uni.dotplots` function will
#' automatically detect and plot numeric variables (columns).
#' @param MAR A numerical vector of the form `c(bottom, left, top, right)` which gives the number of lines
#' of margin to be specified on the four sides of the plot. The default is `c(0.5,4.1,1.1,1.5)`.
#' @param CEX.LAB The magnification to be used for x and y labels relative to the current setting
#' of `CEX.PAR`.
#' @param FONT.LAB The font to be used for x and y labels.
#' @param BTY A character string which determined the type of box which is drawn about plots. If `BTY` is
#' one of "o", "l", "7", "c", "u", or "]" the resulting box resembles the corresponding upper
#' case letter. A value of "n" suppresses the box (the default).
#' @param FG The color to be used for the foreground of plots. This is the default color used for things
#' like axes and boxes around plots (defaults to "gray35").
#' @param COL.AXIS The color to be used for axis annotation. Defaults to "gray35".
#' @param COL.LAB The color to be used for x and y labels. Defaults to "gray20".
#' @param CEX.PAR A numerical value giving the amount by which plotting text and symbols should be
#' magnified relative to the default (for `par`, the panel manager). This starts as 1 when a device
#' is opened, and is reset when the layout is changed, e.g. by setting `mfrow`. Defaults to 0.8.
#' @param TCL The length of tick marks as a fraction of the height of a line of text. The default
#' value is -0.3.
#' @param MGP The margin line (in `mex` units) for the axis title, axis labels and axis line.
#' Note that `mgp[1]` affects title whereas `mgp[2:3]` affect axis. The default is c(2.4, 0.6, 0).
#' @param OMA A vector of the form `c(bottom, left, top, right)` giving the size of the outer margins
#' in lines of text.
#' @param LAB A numerical vector of the form `c(x, y, len)` which modifies the default way that axes
#' are annotated. The values of `x` and `y` give the (approximate) number of tickmarks on the x and y
#' axes and len specifies the label length. The default is `c(5, 5, 7)`. Note that this only affects
#' the way the parameters xaxp and yaxp are set when the user coordinate system is set up, and is
#' not consulted when axes are drawn. `len` is unimplemented in `R`.
#' @param COL.PCH The color to be used for points. Default is "lightcoral".
#' @param PCH The type of points to be drawn. Default is 19. See \code{\link[graphics:points]{points}}
#' for possible values and their interpretation.
#' @param COL.GRID The color of the background grid. Default is "lavender".
#' @param NX The number of lines for the grid in x. Default value is 5.
#' @param NY The number of lines for the grid in Y. Default value is 9.
#' @param LTY The type of lines to be drawn in the background grid. Default value is 6.
#' @param ... Any other parameter that can be incorporated in \code{\link[graphics:plot.default]{plot}}.
#'
#' @return A panel of univariate dotplots.
#' @export
#' @import graphics
#'
#' @examples
#' data("mtcars")
#' uni.dotplots(dataset = mtcars, COL.GRID = "lightblue", LTY = 1, NX = 10, NY = 20)
uni.dotplots <- function(dataset, MAR=c(3,2,0.5,1.5), CEX.LAB = 1.2, FONT.LAB = 2, BTY = "n",
                         FG = "gray35", COL.AXIS = "gray35", COL.LAB = "gray20", CEX.PAR = 0.6,
                         TCL = -0.3, MGP = c(1.7, 0.6, 0.1), OMA = c(1, 0, 1, 0), LAB = c(5, 10, 7),
                         COL.PCH = "lightcoral", PCH = 19, COL.GRID = "lavender", NX = 5, NY = 9, LTY = 6,
                         ...){
  num.data <- dataset[, sapply(dataset, is.numeric)]
  nam <- names(num.data)
  ncol.data <- ncol(num.data)
  ncol.adjust <- ceiling(x = ncol.data/4) # Round to the next integer (e.g. ceiling(x = 7.12) returns 8)!
  num.data <- as.matrix(num.data)

  graphics::par(mfrow= c (ncol.adjust,4), mar=MAR, cex.lab = CEX.LAB, font.lab=FONT.LAB, bty = BTY, fg = FG,
                col.axis = COL.AXIS, col.lab = COL.LAB, cex = CEX.PAR, tcl = TCL,
                mgp = MGP, oma = OMA, lab = LAB)
  for (i in c(1:ncol(num.data))) {
    graphics::plot(x = num.data[,i], y = 1:length(num.data[,i]), type = "p", xlab = nam[i], ylab = "",
                   col = COL.PCH, pch = PCH, panel.first = {
                     grid(col=COL.GRID,nx = NX,ny = NY, lty = LTY)
                   }, ...) }
  # Here, the argument panel.first={} is used to draw the grid first, so behind the points!
}





# _________________________________________
### Creation of a function that generate random samples from various distributions based on my variables
### parameters:

#' Random sample simulation
#'
#' @description The `uni.simudistrib` function **automatically generates 5 Cleveland dotplots** of random
#' samples from different distributions (either `normal`, `log-normal`, or `poisson`) based on the
#' parameters of the variables in the input `data.frame` or `matrix` (see \emph{Details}).  \cr
#' This function is useful to see whether variables' extreme values are actual outliers or whether they
#' lie in a range of values possible for a random sample drawn from a `normal`, `log-normal`, or a `poisson`
#' distribution. \emph{In fine}, it may help determine if the original variable can be approximated by
#' these distribution with or without a transformation.
#'
#' @details The `uni.simudistrib` function extracts some key parameters from the input variables (sample
#' size, mean and standard deviation) and generates random samples based on these parameters. For instance,
#' if `simu.var` contains \emph{i} variables `X1`, `X2`, ... `Xi` and if `distribution = "normal"`, the
#' function will return a panel of \emph{i}x5 plots:
#'
#' * The 1st row will contain five dotplots for five random samples with \emph{n} = `length(X1)` and drawn
#' from a Normal distribution with the same mean and standard deviation as `X1`.
#' * The 2nd row will contain five dotplots for five random samples with \emph{n} = `length(X2)` and drawn
#' from a Normal distribution with the same mean and standard deviation as `X2`.
#' * Etc.
#'
#' \strong{Warning}: the function may fail for `log-normal` and `poisson` distributions if input
#' variables contain negative values (because these distributions are by definition positive). Additionally,
#' if `distribution = "poisson"`, the resulting plots will return \strong{integer} values as Poisson
#' is a discrete probability distribution.
#'
#' @param simu.var A `data.frame` or a `matrix`. For obvious layout and readability reasons, `simu.var`
#' should not include too many variables (`p` < 12 is advised) if the plot is to be printed in a page
#' or window of limited dimensions. For any use in HTML documents (e.g. with `RMarkdown`), the number of
#' input variables should not be a problem.
#' @param distribution Either `"normal"`, `"log-normal"` or `"poisson"` (case sensitive).
#'
#' @return A panel of `p`*5 plots, where `p` is the number of variables in `simu.var`.
#' @export
#' @import graphics
#' @importFrom stats na.omit sd rnorm rlnorm rpois
#'
#' @examples
#' uni.simudistrib(simu.var = iris[,1:4], distribution = "normal")
uni.simudistrib <- function(simu.var, distribution){
  if (!is.data.frame(simu.var) && !is.matrix(simu.var)) {
    warning("The input dataset (i.e. simu.var) must either be a data.frame or a matrix!")
  }
  if(missing(distribution)){
    stop("Please specify the desired distribution. See ?uni.simudistrib for details.")
  }

  ncol.data <- ncol(simu.var)
  nam <- names(simu.var)
  num.data <- as.matrix(simu.var)

  graphics::par(mfrow = c(ncol.data,5), mar = c(3,2,0.1,1.5), cex.lab = 1, font.lab=2, bty = "n",
                fg = "gray35", col.axis = "gray35", col.lab = "gray20", cex = 0.6, tcl = -0.3,
                mgp = c(1.8, 0.6, 0.1), oma = c(0.2, 0.2, 0, 0), lab = c(5, 10, 7))
  for (i in c(1:ncol(num.data))) {
    mu <- mean(stats::na.omit(num.data[,i]))
    std <- stats::sd(stats::na.omit(num.data[,i]))
    n <- length(num.data[,i])
    rrr <- NULL

    if (distribution == "normal") {
      for (j in 1:5) {
        rrr[[j]] <- stats::rnorm(n = n, mean = mu, sd = std)
        graphics::plot(x = rrr[[j]], y = 1:length(rrr[[j]]), type = "p",
                       xlab = paste("rnorm based on", nam[i], sep = " "), ylab = "", col = "hotpink3",
                       pch = 19, panel.first = {grid(col="lavender",nx = 9,ny = 3, lty = 6)})}
    }

    if (distribution == "log-normal") {
      for (j in 1:5) {
        rrr[[j]] <- stats::rlnorm(n = n,
                                  meanlog = log(mu^2 / sqrt(std^2 + mu^2)),
                                  sdlog = sqrt(log(1 + (std^2 / mu^2)))) # That's how you use rlnorm!!!!!!
        graphics::plot(x = rrr[[j]], y = 1:length(rrr[[j]]), type = "p",
                       xlab = paste("rlnorm based on", nam[i], sep = " "), ylab = "", col = "hotpink3",
                       pch = 19, panel.first = {grid(col="lavender",nx = 9,ny = 3, lty = 6)})}
    }

    if (distribution == "poisson") {
      for (j in 1:5) {
        rrr[[j]] <- stats::rpois(n = n, lambda = mu)
        graphics::plot(x = rrr[[j]], y = 1:length(rrr[[j]]), type = "p",
                       xlab = paste("rpois based on", nam[i], sep = " "), ylab = "", col = "hotpink3",
                       pch = 19, panel.first = {grid(col="lavender",nx = 9,ny = 3, lty = 6)})}
    }
  }
}





# _________________________________________
### Creation of a function that generates histograms for each numeric variable:

#' Histogram Panel
#'
#' @description The `uni.histograms` function draws, within a single panel, an independent histogram
#' or each numeric (continuous or discrete) variable in a given dataset. It is particularly useful
#' for data exploration (e.g. Zuur \emph{et al.}, 2010). For instance, to simultaneously observe
#' the distributions of all numeric variables and determine which one will require transformation.
#'
#' @param dataset The input dataset containing all variables to be plotted (must be a `data.frame` with
#' at least 2 variables). It may contain all kinds of columns, the `uni.histograms` function will
#' automatically detect and plot numeric variables (columns).
#' @param MAR A numerical vector of the form `c(bottom, left, top, right)` which gives the number of lines
#' of margin to be specified on the four sides of the plot. The default is `c(0.5,4.1,1.1,1.5)`.
#' @param CEX.LAB The magnification to be used for x and y labels relative to the current setting
#' of `CEX.PAR`.
#' @param FONT.LAB The font to be used for x and y labels.
#' @param BTY A character string which determined the type of box which is drawn about plots. If `BTY` is
#' one of "o", "l", "7", "c", "u", or "]" the resulting box resembles the corresponding upper
#' case letter. A value of "n" suppresses the box (the default).
#' @param FG The color to be used for the foreground of plots. This is the default color used for things
#' like axes and boxes around plots (defaults to "gray35").
#' @param COL.AXIS The color to be used for axis annotation. Defaults to "gray35".
#' @param COL.LAB The color to be used for x and y labels. Defaults to "gray20".
#' @param CEX.PAR A numerical value giving the amount by which plotting text and symbols should be
#' magnified relative to the default (for `par`, the panel manager). This starts as 1 when a device
#' is opened, and is reset when the layout is changed, e.g. by setting `mfrow`. Defaults to 0.8.
#' @param TCL The length of tick marks as a fraction of the height of a line of text. The default
#' value is -0.3.
#' @param MGP The margin line (in `mex` units) for the axis title, axis labels and axis line.
#' Note that `mgp[1]` affects title whereas `mgp[2:3]` affect axis. The default is c(2.4, 0.6, 0).
#' @param OMA A vector of the form `c(bottom, left, top, right)` giving the size of the outer margins
#' in lines of text.
#' @param LAB A numerical vector of the form `c(x, y, len)` which modifies the default way that axes
#' are annotated. The values of `x` and `y` give the (approximate) number of tickmarks on the x and y
#' axes and len specifies the label length. The default is `c(5, 5, 7)`. Note that this only affects
#' the way the parameters xaxp and yaxp are set when the user coordinate system is set up, and is
#' not consulted when axes are drawn. `len` is unimplemented in `R`.
#' @param BREAKS One of:
#' * a vector giving the breakpoints between histogram cells,
#' * a function to compute the vector of breakpoints,
#' * a single number giving the number of cells for the histogram,
#' * a character string naming an algorithm to compute the number of cells (see
#' \code{\link[graphics:hist]{hist}}),
#' * a function to compute the number of cells.
#'
#' In the last three cases the number is a suggestion only; as the breakpoints will be set to pretty
#' values, the number is limited to 1e6 (with a warning if it was larger). If breaks is a function,
#' the x vector is supplied to it as the only argument (and the number of breaks is only limited by
#' the amount of available memory).
#' @param COL The color of the bar of the histograms (bins).
#' @param BORDER The color of the border of the bars.
#'
#' @return A panel of histograms.
#' @import graphics
#' @export
#'
#' @examples
#' uni.histograms(dataset = iris[,1:4])
uni.histograms <- function(dataset, MAR=c(3,2,0.5,1.5), CEX.LAB = 1.2, FONT.LAB = 2, BTY = "n",
                           FG = "gray35", COL.AXIS = "gray35", COL.LAB = "gray20", CEX.PAR = 0.6,
                           TCL = -0.3, MGP = c(1.7, 0.6, 0.1), OMA = c(1, 0, 1, 0), LAB = c(5, 10, 7),
                           BREAKS = 10, COL = "moccasin", BORDER = "white"){
  num.data <- dataset[, sapply(dataset, is.numeric)]
  nam <- names(num.data)
  ncol.data <- ncol(num.data)
  ncol.adjust <- ceiling(x = ncol.data/4) # Round to the next integer (e.g. ceiling(x = 7.12) returns 8)!
  num.data <- as.matrix(num.data)

  graphics::par(mfrow= c (ncol.adjust,4), mar=MAR, cex.lab = CEX.LAB, font.lab=FONT.LAB, bty = BTY, fg = FG,
                col.axis = COL.AXIS, col.lab = COL.LAB, cex = CEX.PAR, tcl = TCL,
                mgp = MGP, oma = OMA, lab = LAB)
  for (i in c(1:ncol(num.data))) {
    graphics::hist(num.data[,i], breaks = BREAKS, col = COL, border = BORDER,
                   main = "", xlab = nam[i], ylab = "")
  }
}





# ___________________________________________________________
### Creation of a function that generates conditional boxplots between the given variables:

#' Bivariate Conditional Boxplots
#'
#' @description The `cond.boxplots` function draws, within a single panel, conditional boxplots (i.e.
#' a plot of a continuous Y variable against one or several categorical X(s) variable(s) (factor),
#' or \strong{discretized} numeric variable(s)) in order to explore bivariate relationships and
#' assess \strong{heteroscedasticity}. \cr
#' If some of the X variables are quantitative (numeric), they will consequently be transformed into
#' factors and discretized into a given number of classes, set approximately by the `N` parameter. This
#' factorisation is required because boxplots are not meant to plot two quantitative variables (without
#' factorisation, the function would plot as many boxplots as there are values in X).
#'
#' @param dataset The input `data.frame` containing the variables to plot (both the Y and the X(s)). Must
#' only contain numeric or factor variables!
#' @param Y The number of the column of `dataset` containing the variable to be plotted as Y. This
#' parameter should be specified as an integer. For instance, if the variable to be used as Y is called
#' "blip" and is the second column of the `dataset`, the `Y` parameter should be `2` (and NOT "two",
#' `dataset$blip`, or `dataset[,2]`).
#' @param Xs The numbers of the columns to be plotted as Xs (e.g. `4:10`, if all the columns from the
#' fourth to the tenth are to be plotted against `Y`).
#' @param outlier Logical (`TRUE` or `FALSE`). Should outliers be plotted or not?
#' @param N Integer. The number of bins on X. Default is 6.
#'
#' @param MAR A numerical vector of the form `c(bottom, left, top, right)` which gives the number of lines
#' of margin to be specified on the four sides of the plot. The default is `c(2.2,2.1,0.5,1.7)`.
#' @param CEX.LAB The magnification to be used for x and y labels relative to the current setting
#' of `CEX.PAR`.
#' @param FONT.LAB The font to be used for x and y labels.
#' @param BTY A character string which determined the type of box which is drawn about plots. If `BTY` is
#' one of "o", "l", "7", "c", "u", or "]" the resulting box resembles the corresponding upper
#' case letter. A value of "n" suppresses the box (the default).
#' @param FG The color to be used for the foreground of plots. This is the default color used for things
#' like axes and boxes around plots (defaults to "gray35").
#' @param COL.AXIS The color to be used for axis annotation. Defaults to "gray35".
#' @param COL.LAB The color to be used for x and y labels. Defaults to "gray20".
#' @param CEX.PAR A numerical value giving the amount by which plotting text and symbols should be
#' magnified relative to the default (for `par`, the panel manager). This starts as 1 when a device
#' is opened, and is reset when the layout is changed, e.g. by setting `mfrow`. Defaults to 0.8.
#' @param TCL The length of tick marks as a fraction of the height of a line of text. The default
#' value is -0.3.
#' @param MGP The margin line (in `mex` units) for the axis title, axis labels and axis line.
#' Note that `mgp[1]` affects title whereas `mgp[2:3]` affect axis. The default is c(1.2, 0.4, 0.2).
#' @param OMA A vector of the form `c(bottom, left, top, right)` giving the size of the outer margins
#' in lines of text.
#' @param TYPE The type of boxplot to draw. Default is "n".
#' @param BORDER An optional vector of colors for the outlines of the boxplots. The values in border
#' are recycled if the length of border is less than the number of plots. Default is "lightcoral".
#' @param COL If col is non-null it is assumed to contain colors to be used to colour the bodies of
#' the boxplots. Default is "moccasin".
#' @param LTY The line type. Line types can either be specified as an integer (0=blank, 1=solid (default),
#' 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash) or as one of the character strings "blank",
#' "solid", "dashed", "dotted", "dotdash", "longdash", or "twodash", where "blank" uses ‘invisible lines’
#' (i.e., does not draw them).
#' @param STAPLEWEX Staple line width expansion, proportional to box width. Default is 0.
#' @param WHISKLWD Whisker line width expansion. Default is 2.
#' @param BOXWEX A scale factor to be applied to all boxes. When there are only a few groups, the
#' appearance of the plot can be improved by making the boxes narrower. Default is 0.7.
#' @param BOXLWD Width of boxplot outer lines. Default is 0.1.
#' @param MEDLWD Width of the median line. Default is 2.6.
#' @param PCH The type of points to be drawn for outliers. Default is 19. See \code{\link[graphics:points]{points}}
#' for possible values and their interpretation.
#' @param ... Any other graphical parameter of `boxplot()`.
#'
#' @return A panel of conditional boxplots.
#' @export
#' @import graphics
#' @importFrom fields bplot.xy
#'
#' @examples
#' cond.boxplots(dataset = iris, Y = 1, Xs = 2:5, outlier = TRUE)
cond.boxplots <- function(dataset, Y, Xs, outlier, N = 6,
                          MAR = c(2.2,2.1,0.5,1.7), CEX.LAB = 0.9, FONT.LAB = 2, BTY = "n", FG = "gray35",
                          COL.AXIS = "gray35", COL.LAB = "gray20", CEX.PAR = 0.8, TCL = -0.3,
                          MGP = c(1.2, 0.4, 0.2), OMA = c(1, 0, 0, 0),
                          TYPE = "n", BORDER = "moccasin", COL = "gray50", LTY = 1, STAPLEWEX = 0,
                          WHISKLWD = 2, BOXWEX = 0.5, BOXLWD = 0.1, MEDLWD = 2.6, PCH = 19, ...){
  mydata <- as.data.frame(dataset)

  if(!is.numeric(mydata[,Y])){
    stop("The Y variable must be numeric! Please check the class of your columns before proceeding.")
  }
  if(!any(sapply(mydata, is.numeric)) && !any(sapply(mydata, is.factor))){
    stop("Input dataset should only contain numeric or factor variables!  Please check the class of your columns before proceeding.")
  }
  if(missing(Y)){
    stop("Please specify the Y and X(s) variables to be used and respect the required format. See ?cond.boxplots for details.")
  }

  nam <- names(mydata)
  ncol.data <- ncol(mydata[,Xs])
  ncol.adjust <- ceiling(x = ncol.data/4) # Round to the next integer (e.g. ceiling(x = 7.12) returns 8)!
  minX <- min(Xs)
  maxX <- max(Xs)
  yyy <- as.matrix(mydata[,Y])

  graphics::par(mfrow= c(ncol.adjust,4), mar=MAR, cex.lab=CEX.LAB, font.lab=FONT.LAB, bty=BTY, fg=FG,
                col.axis=COL.AXIS, col.lab=COL.LAB, cex=CEX.PAR, tcl=TCL, mgp=MGP, oma=OMA)
  for (i in c(minX:maxX)) {

    if (is.numeric(mydata[,i])) {
      fields::bplot.xy(y = mydata[,Y], x = mydata[,i], outlier=outlier, N=N,
                       xlab=(nam[i]), type=TYPE, border=BORDER, col=COL,lty=LTY, staplewex=STAPLEWEX,
                       whisklwd=WHISKLWD, boxwex=BOXWEX, boxlwd=BOXLWD, medlwd=MEDLWD, pch=PCH, cex=0.7, ...)
    }else if (is.factor(mydata[,i])) {
      xxx <- as.matrix(mydata[,i])
      graphics::boxplot(yyy~xxx, outline=outlier, ylab="",
                        xlab=(nam[i]), type=TYPE, border=BORDER, col=COL, lty=LTY, staplewex=STAPLEWEX,
                        whisklwd=WHISKLWD, boxwex=BOXWEX, boxlwd=BOXLWD,medlwd=MEDLWD, pch=PCH, cex=0.7, ...)
    }
  }
}
