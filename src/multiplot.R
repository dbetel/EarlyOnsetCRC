# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
# Copied from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/

multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL, ttl=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
      if(is.null(ttl)) {
          layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                           ncol = cols, nrow = ceiling(numPlots/cols))
      } else {
        
          layout <- matrix(c(rep(1,cols), seq(2, 1 + (cols * ceiling(numPlots/cols)))),
                           ncol = cols, nrow = 1 + ceiling(numPlots/cols), byrow=TRUE)
      }
  }

  if (numPlots==1) {
      print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    if(!is.null(ttl)){
      pushViewport(viewport(layout=grid.layout(nrow(layout), ncol(layout),
                              heights=c(0.1, rep(0.9/(nrow(layout)-1),nrow(layout)-1 )))))
    } else {
      pushViewport(viewport(layout=grid.layout(nrow(layout), ncol(layout))))
    }
    
    if(!is.null(ttl)){
        grid.text(ttl, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:cols, gp=gpar(fontsize=20, fontface='bold')))
    }

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      
      ## Get the i,j matrix positions of the regions that contain this subplot
      if(!is.null(ttl)){
        matchidx <- as.data.frame(which(layout == 1+ i, arr.ind = TRUE))
      }else {
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      }

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



Title_multiplot <- function(..., plotlist=NULL, cols=1, ttl) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
  layout <- matrix(c(rep(1,cols), seq(2, 1 + (2 *(cols * ceiling(numPlots/cols))))),
                           ncol = cols, nrow = 1 + 2*(ceiling(numPlots/cols)), byrow=TRUE)
      
  if (numPlots==1) {
      print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    display.rows <-  (nrow(layout) -1)/2
    pushViewport(viewport(layout=grid.layout(nrow(layout), ncol(layout),
                              heights=c(0.03, rep(c((0.97 - (display.rows*0.02) )/display.rows, 0.02) ,
                                display.rows)))))
    
    if(!is.null(ttl)){
        grid.text(ttl, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:cols, gp=gpar(fontsize=28, fontface='bold')),vjust=5)
    }

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      
      ## Get the i,j matrix positions of the regions that contain this subplot
      plotnumber <- i + (as.integer((i-1)/cols) * cols) +1
      matchidx <- as.data.frame(which(layout == plotnumber, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
      grid.text(paste0('[',LETTERS[i], ']'), vp= viewport(layout.pos.row = matchidx$row +1,
                        layout.pos.col = matchidx$col, gp=gpar(fontsize=20, fontface='bold')),hjust=10.5, vjust=-2)
    }
  }
}
