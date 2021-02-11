################################################################################################
# Prepared for the textbook:
# Data Analysis for Business, Economics, and Policy
# by Gabor BEKES (Central Europen University) and  Gabor KEZDI (University of Michigan)
# Cambridge University Press 2020

# License: Free to share, modify and use for educational purposes. 
# Not to be used for business purposes
# 
#
###############################################################################################x

twoClassSummaryExtended <- function (data, lev = NULL, model = NULL)
{
    lvls <- levels(data$obs)
    rmse <- sqrt(mean((data[, lvls[1]] - ifelse(data$obs == lev[2], 0, 1))^2))
    c(defaultSummary(data, lev, model), "RMSE" = rmse)
}


createRocPlot <- function(r, file_name,  myheight_small = 5.625, mywidth_small = 7.5) {
    all_coords <- coords(r, x="all", ret="all", transpose = FALSE)
    
    roc_plot <- ggplot(data = all_coords, aes(x = fpr, y = tpr)) +
        geom_line(color='purple3', size = 0.7) +
        geom_area(aes(fill = 'navyblue', alpha=0.4), alpha = 0.3, position = 'identity', color = 'purple3') +
        scale_fill_viridis(discrete = TRUE, begin=0.6, alpha=0.5, guide = FALSE) +
        xlab("False Positive Rate (1-Specifity)") +
        ylab("True Positive Rate (Sensitivity)") +
        geom_abline(intercept = 0, slope = 1,  linetype = "dotted", col = "black") +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .1), expand = c(0, 0.01)) +
        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .1), expand = c(0.01, 0)) +
        theme_minimal()
    #+    theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size=13),
    #        axis.title.x = element_text(size=13), axis.title.y = element_text(size=13))
    #save_fig(file_name, output, "small")
    
    #ggsave(plot = roc_plot, paste0(file_name, ".png"),      width=mywidth_small, height=myheight_small, dpi=1200)
    #cairo_ps(filename = paste0(file_name, ".eps"),    #        width = mywidth_small, height = myheight_small, pointsize = 12,    #       fallback_resolution = 1200)
    #print(roc_plot)
    #dev.off()
    
    roc_plot
}

create_calibration_plot <- function(data, prob_var, actual_var, y_lab = "Actual event probability" , n_bins = 10, breaks = NULL) {
    
    if (is.null(breaks)) {
        breaks <- seq(0,1,length.out = n_bins + 1)
    }
    
    binned_data <- data %>%
        mutate(
            prob_bin = cut(!!as.name(prob_var), 
                           breaks = breaks,
                           include.lowest = TRUE)
        ) %>%
        group_by(prob_bin, .drop=FALSE) %>%
        summarise(mean_prob = mean(!!as.name(prob_var)), mean_actual = mean(!!as.name(actual_var)), n = n())
    
    p <- ggplot(data = binned_data) +
        geom_line(aes(mean_prob, mean_actual), color='navyblue', size=0.6, show.legend = TRUE) +
        geom_point(aes(mean_prob,mean_actual), color = 'navyblue', size = 1, shape = 16, alpha = 0.7, show.legend=F, na.rm = TRUE) +
        geom_segment(x=min(breaks), xend=max(breaks), y=min(breaks), yend=max(breaks), color='purple3', size=0.3) +
        theme_minimal() +
        labs(x= "Predicted event probability",
             y= y_lab) +
        coord_cartesian(xlim=c(0,1), ylim=c(0,1))+
        expand_limits(x = 0.01, y = 0.01) +
        scale_y_continuous(expand=c(0.01,0.01),breaks=c(seq(0,1,0.1))) +
        scale_x_continuous(expand=c(0.01,0.01),breaks=c(seq(0,1,0.1))) 
    
    #save_fig(file_name, output, "small")
    p
}

createLossPlot <- function(r, best_coords, file_name,  myheight_small = 5.625, mywidth_small = 7.5) {
    t <- best_coords$threshold[1]
    sp <- best_coords$specificity[1]
    se <- best_coords$sensitivity[1]
    n <- rowSums(best_coords[c("tn", "tp", "fn", "fp")])[1]
    
    all_coords <- coords(r, x="all", ret="all", transpose = FALSE)
    all_coords <- all_coords %>%
        mutate(loss = (fp*FP + fn*FN)/n)
    l <- all_coords[all_coords$threshold == t, "loss"]
    
    loss_plot <- ggplot(data = all_coords, aes(x = threshold, y = loss)) +
        geom_line(color='navyblue', size=0.7) +
        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
        geom_vline(xintercept = t , color = 'purple3' ) +
        annotate(geom = "text", x = t, y= min(all_coords$loss),
                 label=paste0("best threshold: ", round(t,2)),
                 colour='purple3', angle=90, vjust = -1, hjust = -0.5, size = 7) +
        annotate(geom = "text", x = t, y= l,
                 label= round(l, 2), hjust = -0.3, size = 7) +
        theme_minimal()
    
    
    #  ggsave(plot = loss_plot, paste0(file_name,".png"), width=mywidth_small, height=myheight_small, dpi=1200)
    #  cairo_ps(filename = paste0(file_name,".eps"), width = mywidth_small, height = myheight_small, pointsize = 12, fallback_resolution = 1200)
    #  print(loss_plot)
    #  dev.off()
    
    loss_plot
}

#createRocPlotWithOptimal <- function(r, best_coords, file_name,  mywidth_large=12, myheight_large = 9) {
createRocPlotWithOptimal <- function(r, best_coords, file_name,  myheight_small = 5.625, mywidth_small = 7.5) {
    
    all_coords <- coords(r, x="all", ret="all", transpose = FALSE)
    t <- best_coords$threshold[1]
    sp <- best_coords$specificity[1]
    se <- best_coords$sensitivity[1]
    
    roc_plot <- ggplot(data = all_coords, aes(x = specificity, y = sensitivity)) +
        geom_line(color='navyblue', size=0.7) +
        scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
        scale_x_reverse(breaks = seq(0, 1, by = 0.1)) +
        geom_point(aes(x = sp, y = se)) +
        annotate(geom = "text", x = sp, y = se,
                 label = paste(round(sp, 2),round(se, 2),sep = ", "),
                 hjust = 1, vjust = -1, size = 7) +
        theme_minimal()
    #  + theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20),
    #          axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))
    
    
    #  ggsave(plot = roc_plot, paste0(file_name, ".png"),         width=mywidth_small, height=myheight_small, dpi=1200)
    # cairo_ps(filename = paste0(file_name, ".eps"),           width = mywidth_small, height = myheight_small, pointsize = 12,           fallback_resolution = 1200)
    #print(roc_plot)
    #dev.off()
    
    roc_plot
}


draw_confusion_matrix <- function(cm) {
    
    layout(matrix(c(1,1,2)))
    par(mar=c(2,2,2,2))
    plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
    title('CONFUSION MATRIX', cex.main=2)
    
    # create the matrix 
    rect(150, 430, 240, 370, col='#3F97D0')
    text(195, 435, 'Class1', cex=1.2)
    rect(250, 430, 340, 370, col='#F7AD50')
    text(295, 435, 'Class2', cex=1.2)
    text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
    text(245, 450, 'Actual', cex=1.3, font=2)
    rect(150, 305, 240, 365, col='#F7AD50')
    rect(250, 305, 340, 365, col='#3F97D0')
    text(140, 400, 'Class1', cex=1.2, srt=90)
    text(140, 335, 'Class2', cex=1.2, srt=90)
    
    # add in the cm results 
    res <- as.numeric(cm$table)
    text(195, 400, res[1], cex=1.6, font=2, col='white')
    text(195, 335, res[2], cex=1.6, font=2, col='white')
    text(295, 400, res[3], cex=1.6, font=2, col='white')
    text(295, 335, res[4], cex=1.6, font=2, col='white')
    
    # add in the specifics 
    plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
    text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
    text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
    text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
    text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
    text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
    text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
    text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
    text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
    text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
    text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
    
    # add in the accuracy information 
    text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
    text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
    text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
    text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}  
