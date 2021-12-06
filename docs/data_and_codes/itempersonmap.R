

# Function to draw and item-person map
# mod is a mirt model returned from the mirt function
# For meaningful results, Partial Credit Model must be
# used in estimating item thresholds.

itempersonmap <- function(mod) {
  require("mirt")
  require("dplyr")
  require("reshape2")
  require("ggplot2")
  require("cowplot")
  
  if(unique(mod@Model$itemtype) != "Rasch") {
    stop('You must select itemtype = "Rasch" for all items')
  } else {
    
    pars <- as.data.frame(coef(mod, IRTpars = TRUE, simplify = TRUE)$items)
    
    pars <- pars %>%
      select(-a) %>%
      mutate(item = row.names(pars)) %>%
      melt(data = ., 
           id.vars = "item",
           variable.name = "threshold",
           value.name = "parameter")
    
    pars_mean <- pars %>%
      group_by(item) %>%
      summarise(mean_threshold = mean(parameter))
    
    theta <- as.data.frame(fscores(mod)) %>%
      rename(theta = F1)
    
    # Histogram of latent trait distribution
    p1 <- ggplot(data = theta,
                 aes(x = theta)) +
      geom_histogram(bins = 30, fill = "royalblue2", colour = "lightgray") +
      theme_bw(base_size = 13) +
      theme(axis.text.x = element_text(angle = 45)) +
      labs(x = "Latent Trait", y = "") +
      scale_x_continuous(limits = c(-5, 5), breaks = seq(-5, 5, by = 1)) +
      coord_flip() + scale_y_reverse()
    
    # Dot and line plot of item thresholds
    p2 <- ggplot(data = pars, 
                 aes(x = item, y = parameter)) +
      geom_line() +
      geom_point(aes(shape = threshold), size = 3, colour = "indianred1") +
      geom_point(data = pars_mean, aes(x = item, y = mean_threshold), size = 3, 
                 colour = "black", shape = 8) +
      theme_bw(base_size = 13) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      labs(x = "", y = "Item Thresholds", shape = "Threshold") +
      scale_y_continuous(position = "right", limits = c(-5, 5), breaks = seq(-5, 5, by = 1))
    
    # Combine the plots together
    cowplot::plot_grid(p1, p2, align = "h")
    
  }
  
}
