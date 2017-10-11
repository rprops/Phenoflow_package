## ----contrasts, fig.width=7,fig.height=6,fig.cap="Contrast in density between start-up and non-operational phase of cooling water microbial community"----
suppressPackageStartupMessages(library(Phenoflow))
data("CoolingTower")

### Check which parameters are in the fingerprint
CoolingTower@param

### Lets run with the standard c("FL1-H","FL3-H") and lets evaluate how the 
### microbial community compares between the start-up and control phase of the
### cooling water system.
comp <- fp_contrasts(CoolingTower, comp1 = c(5:9), comp2 = c(1:4), 
                     param=c("FL1-H","FL3-H"), thresh=0.01)

### Plot
v <- ggplot2::ggplot(comp, ggplot2::aes(`FL1-H`, `FL3-H`, z = Density))+
  ggplot2::geom_tile(ggplot2::aes(fill=Density)) + 
  ggplot2::geom_point(colour="gray", alpha=0.4)+
  ggplot2::scale_fill_distiller(palette="RdBu", na.value="white") + 
  ggplot2::stat_contour(ggplot2::aes(fill=..level..), geom="polygon", binwidth=0.1)+
  ggplot2::theme_bw()+
  ggplot2::geom_contour(color = "white", alpha = 1)

### Red/positive values indicate higher density for the reactor start-up.
### blue/negative values indicate lower density during start-up.
print(v)

