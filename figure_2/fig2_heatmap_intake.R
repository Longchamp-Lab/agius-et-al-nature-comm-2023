## Create an heatmap 
library(pheatmap)
library(readxl)
library(dplyr)
library(scales)
library(viridis)
library(superheat)

setwd("D:/OneDrive - Université de Lausanne/Documents - LongchampLab/Thomas_Kidney_IRI/20210118_Baseline50_12SucroseLP_Martine/Kcal intake")

# Create a simulated dataset (e.g., 20 rows and 40 columns to match the aspect ratio)
set.seed(123) # for reproducibility
data <- as.data.frame(read_xlsx("intake_data.xlsx", sheet = 1))
rownames(data) <- data$group

data <- data[,-1]

# Assuming your data is in the format where columns are samples and rows are features/genes/observations
# Normalize the data to keep the distribution but scale between 0 and 1
normalized_data <- as.data.frame(lapply(data, rescale))

# Convert the normalized dataframe to a matrix for the heatmap
data_t <- as.data.frame(t(as.matrix(normalized_data)))

rownames(data_t) <- colnames(data)

# identify all scaled values that fall below -0.3
x.col <- data_t < 0.3
# set all values that satisfy the condition to "white"
x.col <- gsub("TRUE", "white", x.col)
# set all values that do not satisfy the condition to "black"
x.col <- gsub("FALSE", "black", x.col)
# convert to matrix
x.col <- matrix(x.col, ncol = ncol(data_t))


## 
png("intake_heatmap.png", height = 6, width = 3.5, units = "in", res = 300)
p <- superheat(data_t,
          scale = F,
          heat.pal = viridis(10),
          heat.lim = c(0, 1),
          heat.pal.values = c(0, 0.5, 1),
          
          X.text = round(t(as.matrix(data)),1),
          X.text.col = x.col,
          
          n.clusters.rows = 1,
          n.clusters.cols = 1,
          grid.hline.size = 1,
          grid.vline.size = 1,
          
          legend.height = 0.2,
          legend.width = 1,
          legend.text.size = 14,
          
          #bottom.label.text.angle = 90,
          #left.label.col = "white",
          bottom.label.col = "white",
          left.label = "none",
          bottom.label = "variable"
          )
dev.off()



### other

png("heatmap.png", width = 10, height = 10, units = "in", res = 300)

# Generate the heatmap without column annotations
heatmap_plot <- pheatmap(data_matrix, 
                         color = viridis(255), # Use viridis color palette
                         cluster_rows = FALSE, 
                         cluster_cols = FALSE,
                         show_rownames = TRUE, 
                         show_colnames = TRUE,
                         scale = "none",
                         #annotation_col = col_annotation,
                         border_color = NA) # No border around cells for now

# Add a black border around the entire plot. You'll have to adjust the dimensions according to your plot size.
#library(grid)
#grid.rect(x = 0.5, y = 0.5, width = 0.99, height = 0.99, just = "center", gp = gpar(col="black", lwd=2))

dev.off() # Close the plotting device

