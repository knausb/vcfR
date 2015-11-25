

#' @rdname chromo_plot
#' @title Plot chromR object
#' @name chromo_plot
#' @export
#' @aliases chromo
#' 
#' @param x an object of class chrom
#' @param nsum logical for whether nsum will be displayed
#' @param DP logical for whether cumulative depth will be displayed
#' @param QUAL logical for whether variant quality will be displayed
#' @param MQ logical for whether mapping quality will be displayed
#' @param HWE logical for whether Hardy-Weinberg equilibrium should be displayed
#' @param NE logical for whether effective size will be displayed
#' @param TPI logical for whether Theta sub pi will be displayed
#' @param TAJD logical for whether Tajima's D will be displayed
#' @param FWH logical for whether Fay and Wu's H will be displayed
#' @param SNPDEN logical for whether variant density will be displayed
#' @param NUC logical for whether nucleotide content will be displayed
#' @param ANN logical for whether annotation will be displayed
#' @param x1 numeric vector for custom tracks
#' @param y1 numeric vector for custom tracks
#' @param label1 optional string label for dot track one.
#' @param x2 numeric vector for custom tracks
#' @param y2 numeric vector for custom tracks
#' @param label2 optional string label for dot track two.
#' @param dot.alpha hexadecimal [00, ff] indicating transparency of dots.
#' @param verbose logical stating whether to produce verbose output
#' @param ... arguments
#' 
chromo <- function(x = x,
                   verbose=TRUE,
                   nsum=TRUE,
                   DP=FALSE,
                   QUAL=FALSE,
                   MQ=FALSE,
                   HWE=FALSE,
                   NE=FALSE,
                   TPI=FALSE,
                   TAJD=FALSE,
                   FWH=FALSE,
                   SNPDEN=FALSE,
                   NUC=FALSE,
                   ANN=FALSE,
                   x1=NULL,
                   y1=NULL,
                   label1=NULL,
                   x2=NULL,
                   y2=NULL,
                   label2=NULL,
                   dot.alpha=22,
                   ...){
  brows <- 0
  srows <- 0
  #
  if(length( x@var.info$DP[x@var.info$mask])>0 & DP  ){brows <- brows+1} # dp
  if(length( x@var.info$MQ[x@var.info$mask])>0 & MQ  ){brows <- brows+1} # mq
  if(length(x@vcf@fix[,'QUAL'][x@var.info$mask])>0 & QUAL){brows <- brows+1} # qual
  #
  if(length( x@var.info$hwe.Da[x@var.info$mask])>0 & HWE ){brows <- brows+1} # HWE
  #
  if(length(       x@var.info$Ne[x@var.info$mask])>0 & NE  ){brows <- brows+1} # Ne
  if(length( x@var.info$theta_pi[x@var.info$mask])>0 & TPI ){brows <- brows+1} # Theta_pi
  if(length(x@var.info$tajimas_d[x@var.info$mask])>0 & TAJD){brows <- brows+1} # Tajima's D
  if(length(  x@var.info$faywu_h[x@var.info$mask])>0 & FWH ){brows <- brows+1} # Fay and Wu's
  #
  if( length(x@win.info$variants)>0 & SNPDEN){brows <- brows+1}
  if(length(x@win.info$A)>0 & NUC   ){brows <- brows+1}
  #
  #  if( sum(x1 == FALSE) > 0 & sum(y1 == FALSE) > 0 ){brows <- brows+1}
  if(is.null(x1) == FALSE & is.null(y1) == FALSE){brows <- brows+1}
  if(is.null(x2) == FALSE & is.null(y2) == FALSE){brows <- brows+1}
  #  if( sum(x2 == FALSE) > 0 & sum(y2 == FALSE) > 0 ){brows <- brows+1}
  #  if(x2 != FALSE & y2 != FALSE){brows <- brows+1}
  #
  if(length(x@ann)>0 & ANN){srows <- srows+1}
  if(nrow(x@seq.info$nuc.win)>0   ){srows <- srows+1}
  #
  if(verbose){
    cat('  Chromo\n')
    cat(paste("  Wide rows   (brows): ", brows, "\n"))
    cat(paste("  Narrow rows (srows): ", srows, "\n"))
  }
  #
  layout(matrix(c(1:((brows+srows)*2)), ncol=2, byrow=T),
         widths=rep(c(1,0.1), times=c(brows+srows)),
         heights=c(rep(1,times=brows),rep(0.4, times=srows)))
  par(mar=c(0,0,0,0))
  par(oma=c(4,4,3,1))
  #
  #  if(sum(x1 == FALSE) > 0 & sum(y1 == FALSE) > 0 ){
  if(is.null(x1) == FALSE & is.null(y1) == FALSE){
    if(is.null(label1)){label1 <- "Custom track 1"}
    plot(x = c(1, x@len),
         y = c(min(y1, na.rm=TRUE), max(y1, na.rm=TRUE)),
         axes = F, 
         frame.plot = TRUE,
         ylab = "",
         type = 'n',
         ...,
    )
    points(x = x1,
           y = y1,
           pch = 20,
           col = paste("#FF8000", dot.alpha, sep = ""),
    )
    title(main = label1, line = -1)
    axis(side = 2, las = 2)
    #
    boxplot(y1, axes = FALSE, frame.plot = TRUE, col = "#FF8000")
  }
  #
  if(is.null(x2) == FALSE & is.null(y2) == FALSE){
    if(is.null(label2)){label2 <- "Custom track 2"}
    plot(x = c(1, x@len),
         y = c(min(y2, na.rm=TRUE), max(y2, na.rm=TRUE)),
         axes = FALSE,
         frame.plot = TRUE,
         ylab="",
         type = 'n',
         ...,
    )
    points(x = x2,
           y = y2,
           pch = 20,
           col = paste("#228B22", dot.alpha, sep=""),
    )
    title(main = label2, line = -1)
    axis(side = 2, las = 2)
    #
    boxplot(y2, axes = FALSE, frame.plot = T, col = "#228B22")
  }
  #
  if(length( x@var.info$hwe.Da[x@var.info$mask])>0 & HWE ){ # HWE
    colv <- rep('#008000', times=length(x@var.info$hwe.p))
    colv[x@var.info$hwe.p < 0.05] <- "#ff0000"
    #    plot(x@vcf.fix$POS[x@var.info$mask], x@var.info$hwe.Da[x@var.info$mask], pch=20,
    plot(x = c(1, x@len),
         y = c(0, max(x@var.info$hwe.chisq[x@var.info$mask], na.rm=TRUE)),
         type = "n",
         axes = FALSE,
         frame.plot = TRUE,
         ylab = "",
         ...
    )
    #    plot(x = x@vcf.fix$POS[x@var.info$mask],
    points(x = x@var.info$POS[x@var.info$mask],
           y = x@var.info$hwe.chisq[x@var.info$mask],
           pch = 20,
           col = paste(colv[x@var.info$mask], dot.alpha, sep="")
    )
    #         col=paste("#800080", dot.alpha, sep=""),
    #         ylim=c(0, max(x@var.info$hwe.chisq[x@var.info$mask], na.rm = TRUE)),
    axis(side = 2, las = 2)
    title(main = "H-W (dis)equilibrium (Chi-square)", line = -1)
    #    boxplot(as.numeric(x@vcf.fix[x@mask,6]), axes=FALSE, frame.plot=T, col="#800080")
    #    boxplot(as.numeric(x@var.info$hwe.Da[x@var.info$mask]), axes=FALSE, frame.plot=T, col="#008000")    
    boxplot(as.numeric(x@var.info$hwe.chisq[x@var.info$mask]), axes = FALSE, frame.plot = TRUE, col="#008000")    
  }
  #  
  if(length(x@var.info$DP[x@var.info$mask])>0 & DP){ # dp
    #    plot(x@vcf.fix[x@mask,2], x@vcf.info[x@mask,1], pch=20, col="#0080ff22", axes=F, frame.plot=T, ylab="", ...)
    #    plot(x = x@vcf.fix$POS[x@var.info$mask],
    plot(x=c(1,x@len),
         y=c(min(x@var.info$DP[x@var.info$mask], na.rm=TRUE), 
             max(x@var.info$DP[x@var.info$mask], na.rm=TRUE)),
         type='n',
         axes = FALSE, 
         frame.plot = TRUE, 
         ylab = "",
         ...
    )
    points(x = x@var.info$POS[x@var.info$mask],
           y = x@var.info$DP[x@var.info$mask],
           pch = 20,
           #         xlim = c(1,x@len),
           col = paste("#0080ff", dot.alpha, sep="")
    )
    axis(side = 2, las = 2)
    title(main = "Read depth (DP)", line = -1)
    #    boxplot(x@vcf.info[x@mask,1], axes=FALSE, frame.plot=T, col="#0080ff")
    boxplot(x@var.info$DP[x@var.info$mask], axes = FALSE, frame.plot = TRUE, col = "#0080ff")
  }
  #
  if(length(x@var.info$MQ[x@var.info$mask])>0 & MQ){ # MQ
    #    plot(x@vcf.fix[x@mask,2], x@vcf.info[x@mask,2], pch=20, col="#3CB37122", axes=F, frame.plot=T, ylab="", ...)
    if(sum(is.na(x@var.info$MQ[x@var.info$mask])) < length(x@var.info$MQ[x@var.info$mask])){
      plot(x = c(1, x@len),
           y = c(min(x@var.info$MQ[x@var.info$mask], na.rm=TRUE), 
                 max(x@var.info$MQ[x@var.info$mask], na.rm=TRUE)),
           axes = FALSE,
           frame.plot = TRUE,
           ylab="",
           type='n',
           ...
      )
      points(x = x@var.info$POS[x@var.info$mask],
             y = x@var.info$MQ[x@var.info$mask],
             pch = 20, 
             col = paste("#3CB371", dot.alpha, sep="")
      )
      axis(side = 2, las = 2)
      title(main = "Mapping quality (MQ)", line = -1)
      #    boxplot(x@vcf.info[x@mask,2], axes=FALSE, frame.plot=T, col="#3CB371")
      boxplot(x@var.info$MQ[x@var.info$mask], axes = FALSE, frame.plot = TRUE, col = "#3CB371")
    } else {
      plot(1,1, type = 'n')
      text(1,1,"No mapping qualities found")
      plot(1,1, type = 'n', axes = FALSE, frame.plot = FALSE)
    }
  }
  #
  if(length(x@vcf@fix[x@var.info$mask, 'QUAL'])>0 & QUAL){ # qual
    #    plot(x@vcf.fix[x@mask,2], x@vcf.fix[x@mask,6], pch=20, col="#80008022", axes=F, frame.plot=T, ylab="", ...)
    plot(x = c(0, x@len),
         y = c(min(x@vcf@fix[x@var.info$mask, 'QUAL'], na.rm=TRUE), 
               max(x@vcf@fix[x@var.info$mask, 'QUAL'], na.rm=TRUE)),
         axes = FALSE,
         frame.plot = TRUE,
         ylab = "",
         type = 'n',
         ...
    )
    points(x = x@var.info$POS[x@var.info$mask],
           y = x@vcf@fix[x@var.info$mask, 'QUAL'],
           pch = 20,
           col = paste("#800080", dot.alpha, sep="")
    )
    axis(side = 2, las = 2)
    title(main = "QUAL", line = -1)
    #    boxplot(as.numeric(x@vcf.fix[x@mask,6]), axes=FALSE, frame.plot=T, col="#800080")
    boxplot(as.numeric(x@vcf@fix[x@var.info$mask, 'QUAL']), axes = FALSE, frame.plot = TRUE, col = "#800080")
  }
  #
  #
  if(length(x@var.info$Ne[x@var.info$mask])>0 & NE){ # Ne
    #    plot(x@vcf.fix[x@mask,2], x@vcf.stat[x@mask,6], pch=20, col="#00008B22", axes=F, frame.plot=T, ylab="", ...)
    plot(x = c(1, x@len),
         y = c(min(x@var.info$Ne[x@var.info$mask], na.rm=TRUE),
               max(x@var.info$Ne[x@var.info$mask], na.rm=TRUE)),
         axes = FALSE,
         frame.plot = TRUE,
         ylab = "",
         type = 'n',
         ...
    )
    points(x = x@var.info$POS[x@var.info$mask],
           y = x@var.info$Ne[x@var.info$mask],
           pch = 20, 
           col = paste("#00008B", dot.alpha, sep="")
    )
    title(main = "Ne", line = -1)
    axis(side = 2, las = 2)
    #    boxplot(as.numeric(x@vcf.stat[x@mask,6]), axes=FALSE, frame.plot=T, col="#00008B")
    boxplot(as.numeric(x@var.info$Ne[x@var.info$mask]), axes = FALSE, frame.plot = TRUE, col = "#00008B")
  }
  if(length(x@var.info$theta_pi[x@var.info$mask])>0 & TPI){ # Theta_pi
    #    plot(x@vcf.fix[x@mask,2], x@vcf.stat[x@mask,7], pch=20, col="#FF8C0022", axes=F, frame.plot=T, ylab="", ...)
    plot(x = c(1, x@len),
         y = c(min(x@var.info$theta_pi[x@var.info$mask], na.rm=TRUE),
               max(x@var.info$theta_pi[x@var.info$mask], na.rm=TRUE)),
         axes = FALSE,
         frame.plot = TRUE,
         ylab="",
         type = 'n',
         ...
    )
    points(x = x@var.info$POS[x@var.info$mask],
           y = x@var.info$theta_pi[x@var.info$mask],
           pch = 20,
           col = paste("#FF8C00", dot.alpha, sep="")
    )
    #    title(main=expression(paste(theta[pi], pi, "Theta_pi")), line=-1)
    title(main = "Theta_pi", line = -1)
    axis(side = 2, las = 2)
    #    boxplot(as.numeric(x@vcf.stat[x@mask,7]), axes=FALSE, frame.plot=T, col="#FF8C00")
    boxplot(as.numeric(x@var.info$theta_pi[x@var.info$mask]), axes = FALSE, frame.plot = TRUE, col = "#FF8C00")
  }
  #
  if(length(x@var.info$tajimas_d[x@var.info$mask])>0 & TAJD){ # Tajima's D
    #    plot(x@vcf.fix[x@mask,2], x@vcf.stat[x@mask,10], pch=20, col="#00640022", axes=F, frame.plot=T, ylab="", ...)
    plot(x = c(1, x@len),
         y = c(min(x@var.info$tajimas_d[x@var.info$mask], na.rm=TRUE),
               max(x@var.info$tajimas_d[x@var.info$mask], na.rm=TRUE)),
         axes = FALSE,
         frame.plot = TRUE,
         ylab="",
         type = 'n',
         ...
    )
    points(x = x@var.info$POS[x@var.info$mask],
           y = x@var.info$tajimas_d[x@var.info$mask],
           pch = 20,
           col = paste("#006400", dot.alpha, sep="")
    )
    abline(0, 0, lty = 2)
    title(main = "Tajima's D", line = -1)
    axis(side = 2, las = 2)
    #    boxplot(as.numeric(x@vcf.stat[x@mask,10]), axes=FALSE, frame.plot=T, col="#006400")
    boxplot(as.numeric(x@var.info$tajimas_d[x@var.info$mask]), axes = FALSE, frame.plot = TRUE, col = "#006400")
  }
  #
  if(length(x@var.info$faywu_h[x@var.info$mask])>0 & FWH){ # Fay and Wu's H
    #    plot(x@vcf.fix[x@mask,2], x@vcf.stat[x@mask,11], pch=20, col="#8B008B22", axes=F, frame.plot=T, ylab="", ...)
    plot(x = c(1, x@len),
         y = c(min(x@var.info$faywu_h[x@var.info$mask], na.rm=TRUE),
               max(x@var.info$faywu_h[x@var.info$mask], na.rm=TRUE)),
         axes = FALSE,
         frame.plot = TRUE,
         ylab="",
         type = 'n',
         ...
    )
    points(x = x@var.info$POS[x@var.info$mask],
           y = x@var.info$faywu_h[x@var.info$mask],
           pch = 20, 
           col = paste("#8B008B", dot.alpha, sep="")
    )
    abline(0, 0, lty = 2)
    title(main = "Fay and Wu's H", line = -1)
    axis(side = 2, las = 2)
    #    boxplot(as.numeric(x@vcf.stat[x@mask,11]), axes=FALSE, frame.plot=T, col="#8B008B")
    boxplot(as.numeric(x@var.info$faywu_h[x@var.info$mask]), axes = FALSE, frame.plot = TRUE, col = "#8B008B")
  }
  if(length(x@win.info$variants)>0 & SNPDEN){
    # SNP density.
    snpden <- x@win.info$variants/x@win.info$length
    plot(x = c(0,x@len),
         y = c(0,max(snpden)),
         type = 'n',
         xlab = "",
         ylab = "",
         axes = FALSE,
         frame.plot = TRUE,
         ...)
    abline(h = seq(0.1, 1, by = 0.1), col = "#a0a0a0")
    abline(h = seq(0.02, 0.08, by = 0.02), lty = 3, col = "#a0a0a0")
    rect(xleft = x@win.info$start,
         ybottom = 0,
         xright = x@win.info$end,
         ytop = snpden,
         col = "#cc0000",
         border=NA)
    axis(side = 2, las = 2)
    title(main = "Variants per site", line = -1)
#    title(main = paste(sum(x@win.info$variants), "total variants"), line = -2)
    title(main = paste(sum(x@var.info$mask), "total variants"), line = -2)
    boxplot(snpden, axes = FALSE, frame.plot = TRUE, col = "#cc0000", ylim = c(0,max(snpden)))
  }
  if(length(x@win.info$A)>0 & NUC){
    # GC and AT content.
    AT <- rowSums(cbind(x@win.info$A, x@win.info$T))/x@win.info$length
    GC <- rowSums(cbind(x@win.info$G, x@win.info$C))/x@win.info$length
    plot(x = c(0,x@len),
         y = c(0,1),
         type = 'n',
         xlab = "",
         ylab = "",
         axes = FALSE,
         frame.plot = TRUE,
         ...)
    rect(xleft = x@win.info$start,
         ybottom = 0,
         xright = x@win.info$end,
         ytop = GC,
         col = "#0000cc",
         border = NA)
    rect(xleft = x@win.info$start,
         ybottom = GC,
         xright = x@win.info$end,
         ytop = GC+AT,
         col = "#ffd700",
         border = NA)
    segments(x0 = x@win.info$start,
             y0 = AT+GC,
             x1 = x@win.info$end,
             y1 = GC+AT,
             col = "#000000")
    axis(side = 2, las = 2)
    title(main = "GC and AT content", line = -1)
    boxplot(x = cbind(GC, AT),
            axes=FALSE,
            frame.plot = TRUE,
            col=c("#0000cc", "#ffd700"),
            ylim=c(0,1),
            border=c("#0000cc", "#ffd700")
    )
    text(1, 0.1, "G/C")
    text(2, 0.1, "A/T")
  }
  #
  if(length(x@ann)>0 & ANN){
    # Annotations.
    plot(x = c(0,x@len),
         y = c(-1,1),
         type = 'n',
         xlab = "",
         ylab = "",
         las = 1,
         axes = FALSE,
         frame.plot = TRUE,
         ...)
    lines(x = c(0,x@len),
          y = c(0, 0),
          lwd=2)
    if(nrow(x@ann) == 0){
      title(main = "No annotations found", line = -1)
      boxplot(x = 1,
              axes = FALSE,
              frame.plot = TRUE)
    } else {
      rect(xleft = as.numeric(as.character(x@ann[,4])),
           ybottom = -1,
           xright = as.numeric(as.character(x@ann[,5])),
           ytop = 1,
           col = "#b22222",
           border = NA)
      title(main = "Annotations", line = -1)
      #
      plot(x = 1:10,
           y = 1:10,
           type = 'n',
           axes = FALSE,
           xlab = "",
           ylab = "",
           ...)
      if(nsum){
        text(5,10,"Genic bases:")
        text(5,9,sum(as.numeric(as.character(x@ann[,5]))-as.numeric(as.character(x@ann[,4]))))
        text(5,8,format(sum(as.numeric(as.character(x@ann[,5]))-as.numeric(as.character(x@ann[,4])))/x@len, digits=3))
      }
    }
  }
  #
  #  if(length(x@acgt.w)>0){
  if(nrow(x@seq.info$nuc.win)>0){    
    # Chromosome.
    plot(x = c(0,x@len),
         y = c(-1,1),
         type = 'n',
         xlab = "",
         ylab = "",
         las = 1,
         axes = FALSE,
         frame.plot = TRUE,
         ...)
    lines(c(0, x@len),c(0, 0), lwd=2)
    #    rect(x@acgt.w[,1], -0.7, x@acgt.w[,2], 0.7, col="#00cc00", border=NA)
    #    rect(x@n.w[,1], -0.4, x@n.w[,2], 0.4, col="#ff6666", border=NA)
    rect(xleft = x@seq.info$nuc.win[,1],
         ybottom = -0.7,
         xright = x@seq.info$nuc.win[,2],
         ytop = 0.7,
         col = "#00cc00",
         border = NA)
#    if(!is.na(x@seq.info$N.win[1,1])){
    if(nrow(x@seq.info$N.win) > 0){
      rect(xleft = x@seq.info$N.win[,1],
           ybottom = -0.4,
           xright = x@seq.info$N.win[,2],
           ytop = 0.4,
           col = "#ff6666",
           border = NA)
    }
    axis(side = 1)
    title(xlab = "Base pairs", line = 2, outer = T)
    title(main = "Green = called bases; Red = n", line = -1)
    #
    plot(x = 1:10,
         y = 1:10,
         type = 'n',
         axes = FALSE,
         xlab = "",
         ylab = "")
    if(nsum){
      text(5, 9, "Length (bp)")
      text(5, 8, x@len)
      text(5, 6, "GC (0.25; 0.5; 0.75)")
      text(5, 5, paste(format(quantile(x@nuccomp.w$gcf, probs=c(0.25,0.5,0.75)), digits=3), collapse="; "))
      text(5, 3, "SNPs (0.25; 0.5; 0.75)")
      text(5, 2, paste(format(quantile(x@snpden.w$density, probs=c(0.25,0.5,0.75)), digits=3), collapse="; "))
    }
  }
  title(main = x@name, outer = TRUE, line = 1)
  #
  # Return to defaults
  layout(matrix(c(1), ncol = 1, byrow = TRUE))
  par(mar=c(5, 4, 4, 2))
  par(oma=c(0, 0, 0, 0))
}


#' @rdname chromo_plot
#' @export
#' @aliases chromoqc
#'
chromoqc <- function(x, nsum = FALSE, ...){
  chromo(x = x,
#         verbose = TRUE,
         nsum = FALSE,
         DP = TRUE,
         QUAL = TRUE, 
         MQ = TRUE,
         SNPDEN = TRUE, 
         NUC = TRUE, 
         ANN = TRUE,
         #         x1=FALSE, y1=FALSE, x2=FALSE, y2=FALSE,
         ...)
}


#' @rdname chromo_plot
#' @export
#' @aliases chromoqc
#'
chromohwe <- function(x, nsum = FALSE, ...){
  chromo(x, 
         verbose = TRUE, 
         nsum = FALSE, 
         DP = TRUE,
         #         QUAL=TRUE, MQ=TRUE,
         HWE = TRUE,
         SNPDEN = TRUE, 
         NUC = TRUE, 
         ANN = TRUE,
         #         x1=FALSE, y1=FALSE, x2=FALSE, y2=FALSE,
         ...)
}


#' @rdname chromo_plot
#' @export
#' @aliases chromodot
#'
chromodot <- function(x, nsum = FALSE, x1 = NULL, y1 = NULL, x2 = NULL, y2 = NULL, ...){
  chromo(x = x,
         verbose = TRUE,
         nsum = FALSE,
         DP = TRUE,
         #         QUAL=FALSE, MQ=FALSE, 
         SNPDEN = TRUE,
         NUC = TRUE,
         ANN = TRUE,
         x1 = x1,
         y1 = y1,
         x2 = x2,
         y2 = y2,
         #         label1=NULL, label2=NULL,
         ...)
}


#' @rdname chromo_plot
#' @export
#' @aliases chromopop
#'
chromopop <- function(x, ...){
  chromo(x = x,
         ANN=TRUE,
         nsum=FALSE, 
         NE=TRUE, 
         TPI=TRUE,
         TAJD=TRUE, 
         FWH=TRUE,
         SNPDEN=TRUE,
         ...)
}

chromoall <- function(x, ...){
  chromo(x = x,
         ANN=TRUE,
         nsum=FALSE,
         DP=TRUE,
         QUAL=TRUE,
         MQ=TRUE,
         NE=TRUE,
         TPI=TRUE,
         TAJD=TRUE,
         FWH=TRUE,
         SNPDEN=TRUE,
         NUC=TRUE, ...)
}







