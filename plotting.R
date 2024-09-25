# Colors
col1 = "darkmagenta"
col2 = "black"

col_disp1 = "deepskyblue2"
col_disp2 = "deepskyblue4"
col_shift1 = "orange1"
col_shift1_sat = "orange3"
col_shift2 = "orangered1"
col_shift2_sat = "orangered3"

grey = transp_grey <- rgb(0, 0, 0, 0.2)
darkgrey <- rgb(0.6, 0.6, 0.6)

# Labels
quantile = "Quantile"
qlevel = "Quantile level"
coverage = "Coverage"
densities = "Densities"

################################################################################
# Helper function to add individual letters in circles to plots
library(plotrix)

text_in_circle <- function(x, y, txt, col, cex = 0.75){
  draw.circle(x, y, 0.75*strwidth("M", cex = 0.75), border = col, col = "white", lwd = 1)
  text(x, y, txt, col = col, cex = 0.75)
}

################################################################################
# Plot AVM decomposition in quantile spread plot with densities
plot.avm_decomp = function(qF,qG,colF = col1,colG = col2,
                           col_dispG = col_disp2,col_dispF = col_disp1,
                           col_shiftG = col_shift2,col_shiftG_double = col_shift2_sat,
                           col_shiftF = col_shift1,col_shiftF_double = col_shift1_sat,
                           lab_F = "F",lab_G = "G",use_legend = TRUE,legend_loc = "bottomleft",ylim = NULL,...){
  use_sat_col = TRUE

  tol = 10^-10
  x = seq(0.001,0.999,0.001) # 1 - coverage
  
  disp1 = pmax(0,qF(1-x/2) - qF(x/2) - qG(1-x/2) + qG(x/2))
  disp2 = pmax(0,qG(1-x/2) - qG(x/2) - qF(1-x/2) + qF(x/2))
  shift1 = pmax(0,pmin(qF(x/2) - qG(x/2),qF(1-x/2) - qG(1-x/2)))
  shift2 = pmax(0,pmin(qG(x/2) - qF(x/2),qG(1-x/2) - qF(1-x/2)))
  shift = (shift1 + shift2)
  
  outerQ = function(p) ifelse(p <= 0.5, pmin(qF(p),qG(p)), pmax(qF(p),qG(p)))
  innerQ = function(p) ifelse(p <= 0.5, pmax(qF(p),qG(p)), pmin(qF(p),qG(p)))
  
  p_u = 1-x/2
  p_l = x/2
  outer_u = outerQ(p_u)
  inner_u = innerQ(p_u)
  outer_l = outerQ(p_l)
  inner_l = innerQ(p_l)
  
  if(is.null(ylim)){
    ylim = c(min(outer_l),max(outer_u))
    print(ylim)
  }
  
  plot(NULL, type = "l", xlim = 0:1, ylim = ylim,
       xlab = coverage, ylab = quantile ,...)
  
  # positive Disp
  i = which(disp1 > tol & outer_u > inner_u + shift + tol)
  li = split(i,i-1:length(i))
  lapply(li,function(i) polygon(c(1-x[i],rev(1-x[i])),c(outer_u[i],rev(inner_u[i] + shift[i])),
          border = NA,col = col_dispF))
  i = which(disp1 > tol & outer_l < inner_l - shift - tol)
  li = split(i,i-1:length(i))
  lapply(li,function(i) polygon(c(1-x[i],rev(1-x[i])),c(outer_l[i],rev(inner_l[i] - shift[i])),
          border = NA,col = col_dispF))
  
  # negative Disp
  i = which(disp2 > tol & outer_u > inner_u + shift + tol)
  li = split(i,i-1:length(i))
  lapply(li,function(i) polygon(c(1-x[i],rev(1-x[i])),c(outer_u[i],rev(inner_u[i] + shift[i])),
          border = NA,col = col_dispG))
  i = which(disp2 > tol & outer_l < inner_l - shift - tol)
  li = split(i,i-1:length(i))
  lapply(li,function(i) polygon(c(1-x[i],rev(1-x[i])),c(outer_l[i],rev(inner_l[i] - shift[i])),
          border = NA,col = col_dispG))
  
  # positive Shift
  i = shift1 > tol
  polygon(c(1-x[i],rev(1-x[i])),c(inner_u[i],rev(inner_u[i] + shift[i])),
          border = NA,col = col_shiftF)
  polygon(c(1-x[i],rev(1-x[i])),c(inner_l[i],rev(inner_l[i] - shift[i])),
          border = NA,col = col_shiftF)
  if(use_sat_col){
    i = i & inner_u < inner_l
    polygon(c(1-x[i],rev(1-x[i])),c(inner_l[i],rev(inner_u[i])),
            border = NA,col = col_shiftF_double)
  }
  
  # negative Shift
  i = shift2 > tol
  polygon(c(1-x[i],rev(1-x[i])),c(inner_u[i],rev(inner_u[i] + shift[i])),
          border = NA,col = col_shiftG)
  polygon(c(1-x[i],rev(1-x[i])),c(inner_l[i],rev(inner_l[i] - shift[i])),
          border = NA,col = col_shiftG)
  if(use_sat_col){
    i = i & inner_u < inner_l
    polygon(c(1-x[i],rev(1-x[i])),c(inner_l[i],rev(inner_u[i])),
            border = NA,col = col_shiftG_double)
  }
  
  lines(1-x, qF(p_l), col = colF, lwd = 2)
  lines(1-x, qF(p_u), col = colF, lwd = 2)
  lines(1-x, qG(p_l), col = colG, lwd = 2)
  lines(1-x, qG(p_u), col = colG, lwd = 2)
  
  atF = 0.9
  atG = 0.95
  text_in_circle(atF, qF((1+atF)/2), lab_F, col = colF)
  text_in_circle(atF, qF((1-atF)/2), lab_F, col = colF)
  text_in_circle(atG, qG((1+atG)/2), lab_G, col = colG)
  text_in_circle(atG, qG((1-atG)/2), lab_G, col = colG)
  
  if(use_legend){
    legend(legend_loc, col = c(col_dispF, col_shiftF, col_dispG, col_shiftG),
           legend = c(expression(Disp["+"]),expression(Shift["+"]),
                      expression(Disp["-"]),expression(Shift["-"])), 
           pch = 15, bty = "n", ncol = 2)
  }
}

plot.densities = function(dF,dG,pmF = NULL,pmG = NULL,colF = col1,colG = col2,ylim = NULL,...){
  x.grid = seq(ylim[1],ylim[2],0.01)
  xlim = c(0,max(c(dF(x.grid),dG(x.grid)),pmF$p,pmG$p))*1.1
  if(is.null(ylim)) ylim = c(-3,3)
  plot(NULL,xlim = xlim,ylim = ylim,xlab = densities,...)

  x.grid = seq(ylim[1],ylim[2],0.01)
  polygon(c(0,dF(x.grid),0),c(ylim[1],x.grid,ylim[2]),border = colF,col = adjustcolor(colF, alpha.f=0.5))
  polygon(c(0,dG(x.grid),0),c(ylim[1],x.grid,ylim[2]),border = colG,col = adjustcolor(colG, alpha.f=0.5))
  
  # Plot point masses
  if(!is.null(pmF)){
    segments(rep(0,length(pmF$p)),pmF$q,pmF$p,pmF$q,col = colF)
    points(pmF$p,pmF$q,pch = 16,col = colF)
  }
  if(!is.null(pmG)){
    segments(rep(0,length(pmG$p)),pmG$q,pmG$p,pmG$q,col = colG)
    points(pmG$p,pmG$q,pch = 16,col = colG)
  }
}

plot.densities.decomp = function(dF,dG,qF,qG,pmF = NULL,pmG = NULL,ylim = NULL,colF = col1,colG = col2,
                                 lab_F = "F",lab_G = "G",use_legend = TRUE,legend_loc = "bottomleft",...){
  layout(matrix(c(1,2),nrow = 1),widths = c(1,2.5))
  
  par(mar = c(2.2,2.2,0.2,0),mgp = c(1.2,0.4,0))
  plot.densities(dF,dG,pmF,pmG,colF = colF,colG = colG,ylim = ylim,ylab = quantile,xaxt = "n")
  
  par(mar = c(2.2,0,0.2,0),mgp = c(1.2,0.4,0))
  plot.avm_decomp(qF,qG,colF = colF,colG = colG,ylim = ylim,yaxt = "n",lab_F = lab_F,lab_G = lab_G,use_legend = use_legend,legend_loc = legend_loc,...)
}
################################################################################
# Plot CD decomposition slice in quantile spread plot with densities
plot.CD_decomp = function(qF,qG,beta,colF = col1,colG = col2,
                           col_dispG = col_disp2,col_dispF = col_disp1,
                           col_shiftG = col_shift2,col_shiftG_double = col_shift2_sat,
                           col_shiftF = col_shift1,col_shiftF_double = col_shift1_sat,
                           lab_F = "F",lab_G = "G",use_legend = TRUE, ylim = NULL,...){
  use_sat_col = TRUE

  tol = 10^-10
  x = seq(0.001,0.999,0.001) # 1 - coverage
  
  disp1 = ifelse(x >= 1-beta,
                 pmax(0,qF(1-x/2) - qF(x/2) - qG(1-(1-beta)/2) + qG((1-beta)/2)),
                 0)
  disp2 = ifelse(x <= 1-beta, pmax(0,qG(1-(1-beta)/2) - qG((1-beta)/2) - qF(1-x/2) + qF(x/2)),0)
  shift1 = pmax(0,pmin(qF(x/2) - qG((1-beta)/2),qF(1-x/2) - qG(1-(1-beta)/2)))
  shift2 = pmax(0,pmin(qG((1-beta)/2) - qF(x/2),qG(1-(1-beta)/2) - qF(1-x/2)))
  shift = (shift1 + shift2)

  outerQ = function(p) ifelse(p <= 0.5, pmin(qF(p),qG((1-beta)/2)), pmax(qF(p),qG(1-(1-beta)/2)))
  innerQ = function(p) ifelse(p <= 0.5, pmax(qF(p),qG((1-beta)/2)), pmin(qF(p),qG(1-(1-beta)/2)))
  
  p_u = 1-x/2
  p_l = x/2
  outer_u = outerQ(p_u)
  inner_u = innerQ(p_u)
  outer_l = outerQ(p_l)
  inner_l = innerQ(p_l)
  
  if(is.null(ylim)){
    ylim = c(min(outer_l),max(outer_u))
    print(ylim)
  }
  
  plot(NULL, type = "l", xlim = 0:1, ylim = ylim,
       xlab = coverage, ylab = quantile ,...)
  
  # positive Disp
  i = disp1 > tol & outer_u > inner_u + shift + tol
  polygon(c(1-x[i],rev(1-x[i])),c(outer_u[i],rev(inner_u[i] + shift[i])),
          border = NA,col = col_dispF)
  i = disp1 > 10^-10 & outer_l < inner_l - shift - tol
  polygon(c(1-x[i],rev(1-x[i])),c(outer_l[i],rev(inner_l[i] - shift[i])),
          border = NA,col = col_dispF)
  
  # negative Disp
  i = disp2 > tol & outer_u > inner_u + shift + tol
  polygon(c(1-x[i],rev(1-x[i])),c(outer_u[i],rev(inner_u[i] + shift[i])),
          border = NA,col = col_dispG)
  i = disp2 > 10^-10 & outer_l < inner_l - shift - tol
  polygon(c(1-x[i],rev(1-x[i])),c(outer_l[i],rev(inner_l[i] - shift[i])),
          border = NA,col = col_dispG)
  
  # positive Shift
  i = shift1 > tol & 1-x <= beta
  polygon(c(1-x[i],rev(1-x[i])),c(inner_u[i],rev(inner_u[i] + shift[i])),
          border = NA,col = col_shiftF)  
  if(use_sat_col){
    i = i & inner_u < inner_l
    polygon(c(1-x[i],rev(1-x[i])),c(inner_l[i],rev(inner_u[i])),
            border = NA,col = col_shiftF_double)
  }
  i = shift1 > tol & 1-x >= beta
  polygon(c(1-x[i],rev(1-x[i])),c(inner_l[i],rev(inner_l[i] - shift[i])),
          border = NA,col = col_shiftF)
  if(use_sat_col){
    i = i & inner_u < inner_l
    polygon(c(1-x[i],rev(1-x[i])),c(inner_l[i],rev(inner_u[i])),
            border = NA,col = col_shiftF_double)
  }

  # negative Shift
  i = shift2 > tol & 1-x >= beta
  polygon(c(1-x[i],rev(1-x[i])),c(inner_u[i],rev(inner_u[i] + shift[i])),
          border = NA,col = col_shiftG)
  if(use_sat_col){
    i = i & inner_u < inner_l
    polygon(c(1-x[i],rev(1-x[i])),c(inner_l[i],rev(inner_u[i])),
            border = NA,col = col_shiftG_double)
  }
  i = shift2 > tol & 1-x <= beta
  polygon(c(1-x[i],rev(1-x[i])),c(inner_l[i],rev(inner_l[i] - shift[i])),
          border = NA,col = col_shiftG)
  if(use_sat_col){
    i = i & inner_u < inner_l
    polygon(c(1-x[i],rev(1-x[i])),c(inner_l[i],rev(inner_u[i])),
            border = NA,col = col_shiftG_double)
  }
  
  lines(1-x, qF(p_l), col = colF, lwd = 2)
  lines(1-x, qF(p_u), col = colF, lwd = 2)
  lines(1-x, qG(p_l), col = adjustcolor(colG,alpha = 0.3), lwd = 2,lty = 2)
  lines(1-x, qG(p_u), col = adjustcolor(colG,alpha = 0.3), lwd = 2,lty = 2)
  abline(h = qG((1-beta)/2), col = colG, lwd = 2)
  abline(h = qG(1-(1-beta)/2), col = colG, lwd = 2)
  abline(v = beta,col = "grey",lwd = 1)
  
  atF = 0.9
  atG = 0.95
  text_in_circle(atF, qF((1+atF)/2), lab_F, col = colF)
  text_in_circle(atF, qF((1-atF)/2), lab_F, col = colF)
  text_in_circle(atG, qG((1+atG)/2), lab_G, col = adjustcolor(col2,alpha = 0.3))
  text_in_circle(atG, qG((1-atG)/2), lab_G, col = adjustcolor(col2,alpha = 0.3))
  
  if(use_legend){
    legend("bottomleft", col = c(col_dispF, col_shiftF, col_dispG, col_shiftG),
           legend = c(expression(Disp["+"]),expression(Shift["+"]),
                      expression(Disp["-"]),expression(Shift["-"])), 
           pch = 15, bty = "n", ncol = 2)
  }
}

plot.densities.decompCD = function(dF,dG,qF,qG,beta,use_legend = TRUE,ylim = NULL){
  layout(matrix(c(1,2),nrow = 1),widths = c(1,2.5))
  
  par(mar = c(2,2,0,-0.2) + 0.2,mgp = c(1.2,0.4,0))
  plot.densities(dF,dG,ylim = ylim,ylab = quantile,xaxt = "n")
  
  par(mar = c(2,-0.2,0,0) + 0.2,mgp = c(1.2,0.4,0))
  plot.CD_decomp(qF,qG,beta,ylim = ylim,yaxt = "n",use_legend = use_legend)
}

################################################################################
# Illustrations at the interval level [Remove in public replication?]

# Function to add lines representing decomposition for a single pair of intervals
avm_lines <- function(l_F, u_F, l_G, u_G, at_l, at_u, # interval ends and an x coordinate (at)
                      col_F_disp = col_disp1,
                      col_G_disp = col_disp2,
                      col_F_larger = col_shift1,
                      col_G_larger = col_shift2, ...){

  # compute interval widths
  d1 <- u_F - l_F
  d2 <- u_G - l_G

  # initialize vectors where line segments wil be stored
  F_disp_l <- F_disp_u <- G_disp_l <- G_disp_u <-
    F_larger_l <- F_larger_u <- G_larger_l <- G_larger_u <- NA

  # case where interval 1 is wider
  if(d1 >= d2){
    d <- d1 - d2

    # run through cases
    if(u_F < l_G){
      F_disp_l <- c(l_F, l_F + d)
      G_larger_u <- c(u_G, u_F)
      G_larger_l <- c(l_G, l_F + d)
    }

    if(l_F < l_G & l_G < u_F & u_F < u_G){
      F_disp_l <- c(l_F, l_F + d)
      G_larger_l <- c(l_F + d, l_G)
      G_larger_u <- c(u_F, u_G)
    }

    if(l_F < l_G & u_G < u_F){
      F_disp_l <- c(l_F, l_G)
      F_disp_u <- c(u_G, u_F)
    }

    if(l_G < l_F & l_F < u_G & u_G < u_F){
      F_disp_u <- c(u_F, u_F - d)
      F_larger_u <- c(u_F - d, u_G)
      F_larger_l <- c(l_F, l_G)
    }

    if(u_G < l_F){
      F_disp_u <- c(u_F, u_F - d)
      F_larger_u <- c(u_F - d, u_G)
      F_larger_l <- c(l_F, l_G)
    }
  }

  # case where interval 2 is wider
  if(d2 >= d1){
    d <- d2 - d1

    # run through cases
    if(u_G < l_F){
      G_disp_l <- c(l_G, l_G + d)
      F_larger_u <- c(u_F, u_G)
      F_larger_l <- c(l_F, l_G + d)
    }

    if(l_G < l_F & l_F < u_G & u_G < u_F){
      G_disp_l <- c(l_G, l_G + d)
      F_larger_l <- c(l_G + d, l_F)
      F_larger_u <- c(u_G, u_F)
    }

    if(l_G < l_F & u_F < u_G){
      G_disp_l <- c(l_G, l_F)
      G_disp_u <- c(u_F, u_G)
    }

    if(l_F < l_G & l_G < u_F & u_F < u_G){
      G_disp_u <- c(u_G, u_G - d)
      G_larger_u <- c(u_G - d, u_F)
      G_larger_l <- c(l_G, l_F)
    }

    if(u_F < l_G){
      dispu <- c(u_G, u_G - d)
      G_larger_u <- c(u_G - d, u_F)
      G_larger_l <- c(l_G, l_F)
    }
  }

  # add lines
  lines(rep(at_l, length(F_disp_l)), F_disp_l, col = col_F_disp, ...)
  lines(rep(at_l, length(G_disp_l)), G_disp_l, col = col_G_disp, ...)
  lines(rep(at_u, length(F_disp_u)), F_disp_u, col = col_F_disp, ...)
  lines(rep(at_u, length(G_disp_u)), G_disp_u, col = col_G_disp, ...)
  lines(rep(at_l, length(F_larger_l)), F_larger_l, col = col_F_larger, ...)
  lines(rep(at_l, length(G_larger_l)), G_larger_l, col = col_G_larger, ...)
  lines(rep(at_u, length(F_larger_u)), F_larger_u, col = col_F_larger, ...)
  lines(rep(at_u, length(G_larger_u)), G_larger_u, col = col_G_larger, ...)
}

# function to plot an interval (vertically)
add_interval <- function(l, u, at, wdt = 0.1, col  ="black"){
  lines(rep(at, 2), c(l, u), col = col)
  lines(at + c(wdt, -wdt)/2, rep(l, 2), col = col)
  lines(at + c(wdt, -wdt)/2, rep(u, 2), col = col)
}

# function wrapping up plotting of intervals and lines
add_intervals_plus_avm_lines <- function(l_F, u_F, l_G, u_G, # interval ends
                                         at1, at2, # x-coordinates
                                         col_F = col1, col_G = col2, # colours
                                         label_F = NULL, label_G = NULL, # labels
                                         wdt = 0.1){ # width
  add_interval(l_F, u_F, at1, wdt = wdt, col = col_F)
  if(!is.null(label_F)) text_in_circle(at1, (l_F + u_F)/2, label_F, col = col_F)
  add_interval(l_G, u_G, at2, wdt = wdt, col = col_G)
  if(!is.null(label_G)) text_in_circle(at2, (l_G + u_G)/2, label_G, col = col_G)
  avm_lines(l_F, u_F, l_G, u_G, at_l = at1/3 + 2*at2/3, at_u = 2*at1/3 + at2/3, lwd = 3)
}

# For CD

# function wrapping up plotting of intervals and lines
# function to add lines representing decomposition for a single pair of intervals
cd_lines <- function(l_F, u_F, l_G, u_G, at_l, at_u, # interval ends and an x coordinate (at)
                     col_F_disp = col_disp1, # colours
                     col_G_disp = col_disp2,
                     col_F_larger = col_shift1,
                     col_G_larger = col_shift2, ...){

  # compute interval widths
  d1 <- u_F - l_F
  d2 <- u_G - l_G

  # initialize vectors where line segments wil be stored
  F_disp_l <- F_disp_u <- G_disp_l <- G_disp_u <-
    F_larger_l <- F_larger_u <- G_larger_l <- G_larger_u <- NA

  # case where interval 1 is wider
  if(d1 >= d2){
    d <- d1 - d2

    # run through cases
    if(u_F < l_G){
      F_disp_l <- c(l_F, l_F + d)
      G_larger_u <- c(u_G, u_F)
      G_larger_l <- c(l_G, u_F)
    }

    if(l_F < l_G & l_G < u_F & u_F < u_G){
      F_disp_l <- c(l_F, l_F + d)
      G_larger_l <- c(l_F + d, l_G)
      G_larger_u <- c(u_F, u_G)
    }

    if(l_F < l_G & u_G < u_F){
      F_disp_l <- c(l_F, l_G)
      F_disp_u <- c(u_G, u_F)
    }

    if(l_G < l_F & l_F < u_G & u_G < u_F){
      F_disp_u <- c(u_F, u_F - d)
      F_larger_u <- c(u_F - d, u_G)
      F_larger_l <- c(l_F, l_G)
    }

    if(u_G < l_F){
      F_disp_u <- c(u_F, u_F - d)
      F_larger_u <- c(u_F - d, u_G)
      F_larger_l <- c(l_F, u_G)
    }
  }

  # case where interval 2 is wider
  if(d2 >= d1){
    d <- d2 - d1

    # run through cases
    if(u_G < l_F){
      G_disp_l <- c(l_G, l_G + d)
      F_larger_u <- c(u_F, u_G)
      F_larger_l <- c(l_F, u_G)
    }

    if(l_G < l_F & l_F < u_G & u_G < u_F){
      G_disp_l <- c(l_G, l_G + d)
      F_larger_l <- c(l_G + d, l_F)
      F_larger_u <- c(u_G, u_F)
    }

    if(l_G < l_F & u_F < u_G){
      G_disp_l <- c(l_G, l_F)
      G_disp_u <- c(u_F, u_G)
    }

    if(l_F < l_G & l_G < u_F & u_F < u_G){
      G_disp_u <- c(u_G, u_G - d)
      G_larger_u <- c(u_G - d, u_F)
      G_larger_l <- c(l_G, l_F)
    }

    if(u_F < l_G){
      dispu <- c(u_G, u_G - d)
      G_larger_u <- c(u_G - d, u_F)
      G_larger_l <- c(l_G, u_F)
    }
  }

  # add lines
  lines(rep(at_l, length(F_disp_l)), F_disp_l, col = col_F_disp, ...)
  lines(rep(at_l, length(G_disp_l)), G_disp_l, col = col_G_disp, ...)
  lines(rep(at_u, length(F_disp_u)), F_disp_u, col = col_F_disp, ...)
  lines(rep(at_u, length(G_disp_u)), G_disp_u, col = col_G_disp, ...)
  if(u_F < l_G | u_G < l_F){
    lines(rep(at_l, length(F_larger_l)), F_larger_l, col = col_F_larger, ...)
    lines(rep(at_l, length(G_larger_l)), G_larger_l, col = col_G_larger, ...)
  }
  lines(rep(at_u, length(F_larger_u)), F_larger_u, col = col_F_larger, ...)
  lines(rep(at_u, length(G_larger_u)), G_larger_u, col = col_G_larger, ...)
}

add_intervals_plus_cd_lines <- function(l_F, u_F, l_G, u_G, # interval ends
                                        at1, at2, # x-coordinates
                                        col_F = col1, col_G = col2, # colours
                                        label_F = NULL, label_G = NULL, # labels
                                        wdt = 0.1){ # width
  add_interval(l_F, u_F, at1, wdt = wdt, col = col_F)
  if(!is.null(label_F)) text_in_circle(at1, (l_F + u_F)/2, label_F, col = col_F)
  add_interval(l_G, u_G, at2, wdt = wdt, col = col_G)
  if(!is.null(label_G)) text_in_circle(at2, (l_G + u_G)/2, label_G, col = col_G)
  cd_lines(l_F, u_F, l_G, u_G, at_l = at1/3 + 2*at2/3, at_u = 2*at1/3 + at2/3, lwd = 3)
}

################################################################################
# Quantile plot with densities
plot.QFs = function(qF,qG,ylim = c(-2,2),colF = col1,colG = col2,
                    lab_F = "F",lab_G = "G",fold = FALSE,...){
  xlim = c(0.001,0.999)
  x.grid = seq(xlim[1],xlim[2],0.001)
  x.grid.folded = ifelse(x.grid >= 0.5,x.grid,1-x.grid)
  
  plot(NULL,xlim = c(0,1),ylim = ylim,xlab = qlevel,...)
  
  if(fold){
    polygon(c(1-xlim[1],x.grid.folded,xlim[2],xlim[2],rev(x.grid.folded),1-xlim[1]),
            c(qF(xlim[1]),qF(x.grid),qF(xlim[2]),qG(xlim[2]),qG(sort(x.grid,decreasing = TRUE)),qG(xlim[1])),
            border = NA,col = transp_grey)
    x.grid.dark = x.grid[qF(1-x.grid) > qG(x.grid) & x.grid >= 0.5]
    polygon(c(x.grid.dark,rev(x.grid.dark)),c(qF(1-x.grid.dark),qG(rev(x.grid.dark))),border = NA,col = darkgrey)
    plot(qF,xlim = c(0.5,1),col = col1,add =TRUE,n = 1000*(xlim[2] - xlim[1])+1, lwd = 2)
    plot(qG,xlim = c(0.5,1),col = col2,add =TRUE,n = 1000*(xlim[2] - xlim[1])+1, lwd = 2)
    plot(function(x) qF(1-x),xlim = c(0.5,1),col = col1,add =TRUE,n = 1000*(xlim[2] - xlim[1])+1, lwd = 2)
    plot(function(x) qG(1-x),xlim = c(0.5,1),col = col2,add =TRUE,n = 1000*(xlim[2] - xlim[1])+1, lwd = 2)
  }
  else{
    polygon(c(xlim[1],x.grid,xlim[2],xlim[2],sort(x.grid,decreasing = TRUE),xlim[1]),
            c(qF(xlim[1]),qF(x.grid),qF(xlim[2]),qG(xlim[2]),qG(sort(x.grid,decreasing = TRUE)),qG(xlim[1])),
            border = NA,col = transp_grey)
    plot(qF,xlim = c(0,1),col = col1,add =TRUE,n = 1000*(xlim[2] - xlim[1])+1, lwd = 2)
    plot(qG,xlim = c(0,1),col = col2,add =TRUE,n = 1000*(xlim[2] - xlim[1])+1, lwd = 2)
  }
  
  atF = 0.9
  atG = 0.95
  text_in_circle((1+atF)/2, qF((1+atF)/2), lab_F, col = colF)
  text_in_circle((1+atG)/2, qG((1+atG)/2), lab_G, col = colG)
}

plot.densities.QF = function(dF,dG,qF,qG,ylim = NULL,colF = col1,colG = col2,
                             lab_F = "F",lab_G = "G",fold = FALSE){
  layout(matrix(c(1,2),nrow = 1),widths = c(1,2.5))
  
  par(mar = c(2.2,2.2,0.2,0),mgp = c(1.2,0.4,0))
  plot.densities(dF,dG,colF = colF,colG = colG,ylim = ylim,ylab = quantile,xaxt = "n")
  
  par(mar = c(2.2,0,0.2,0),mgp = c(1.2,0.4,0))
  plot.QFs(qF,qG,colF = colF,colG = colG,ylim = ylim,yaxt = "n",lab_F = lab_F,lab_G = lab_G,fold = fold)
}



