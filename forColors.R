library(getopt)
 
Opts <- matrix(
  c("colorsRegion",  "a", 1, "integer"),
  byrow=TRUE, ncol=4)
 opt <- getopt(spec=Opts)
print(opt$colorsRegion)
rwColors <- rainbow(opt$colorsRegion)
palette(rwColors)
colors=palette()

rgb.palette <- colorRampPalette(c("red", "green", "blue"), space = "rgb") 
colors<-rgb.palette(opt$colorsRegion)
write.table(colors,file='colors.txt',sep=' ') 
