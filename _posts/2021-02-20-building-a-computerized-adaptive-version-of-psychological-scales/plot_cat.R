
# model is an object of randomCAT in the catR package

plot_cat <- function(model) {
  require("ggplot2")
  
  cat_summary <- data.frame(
    items = factor(model$testItems, levels = model$testItems),
    thetas = model$thetaProv,
    se = model$seProv,
    theta_lb = model$thetaProv - 1.96 * model$seProv,
    theta_ub = model$thetaProv + 1.96 * model$seProv,
    wr = factor(model$pattern, levels=c(0, 1), labels=c("Wrong", "Right")))
  
  p <- ggplot(data=cat_summary, aes(x=items, y=thetas, color=wr)) + 
    geom_point(aes(size=se)) +
    geom_linerange(aes(ymin=theta_lb, ymax=theta_ub), linetype=3, size = 1) +
    geom_point(aes(x=tail(items, n = 1), y=tail(thetas, n = 1)), color="black", pch=4, size=4) +
    geom_hline(aes(yintercept = model$trueTheta), color="#808080", linetype="dashed") +
    coord_cartesian(ylim=c(-4, 4)) + 
    scale_size_continuous(range=c(1, 3)) +
    labs(x = "Items", y = expression(paste("Estimated ", theta)), color = "Response") +
    guides(size=F, alpha=F) + 
    theme_bw() + 
    theme(legend.key=element_blank())
  
  return(p)
}

# Assuming your model is called "model 1"
# plot_cat(model1)
