plot_sfr <- function(type, fit, save = FALSE){
  caption = glue("PH Test For  {type}")
  p <- ggcoxzph(caption = caption, 
                fit = fit, 
                ggtheme = theme_bw()
                )
  
  
  if (save == TRUE){
    
    filename = glue("schoenfeld_residuals_{type}.png")
    path = "figures/schoenfeld_residual_plot/"
    
    print(glue("Saving {filename} into {path}"))

    png(file.path(path, filename), width = 5*300, height = 7*300, res = 300)
    print(p$plot)   # ensures grid graphics render
    dev.off()
  }
  
  p
}
