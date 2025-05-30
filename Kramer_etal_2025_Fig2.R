---
  title: "Upset Plot for EXPORTS 18S samples"
# Sasha Kramer
# guide: https://cran.r-project.org/web/packages/UpSetR/readme/README.html
# and: https://upset.app/implementations/#python-libraries
---

## Clear workspace:
rm(list=ls())

## Set working directory
setwd("~/18S/upset/")

## Load library
library(UpSetR)
library(ggplot2)

## Load data
atl <- read.csv("upset_atlantic_nodeep.csv",header = T,sep = ",",check.names = "FALSE")
pac <- read.csv("pacific_upset_nodeep.csv",header = T,sep = ",",check.names = "FALSE")

## Color plot
upset(atl,sets = c("Surface water","Bulk particles","Ind. particles"),matrix.color="black",main.bar.color="black",sets.bar.color = c("#dcaae8","#ac395b","#e09cc3"),queries = list(list(query = intersects,params = list("Ind. particles"), color = "#e09cc3", active = T),list(query = intersects,params = list("Bulk particles"), color = "#ac395b", active = T),list(query = intersects,params = list("Surface water"), color = "#dcaae8", active = T)),point.size = 6, line.size = 1, mainbar.y.label = "ASVs",sets.x.label = "ASVs",order.by = "freq",text.scale = c(4,4,4,4,4,4),keep.order = TRUE)
upset(pac,sets = c("Surface water","Bulk particles","Ind. particles"),matrix.color="black",main.bar.color="black",sets.bar.color = c("#b7e2d6","#50738a","#97b2c9"),queries = list(list(query = intersects,params = list("Ind. particles"), color = "#97b2c9", active = T),list(query = intersects,params = list("Bulk particles"), color = "#50738a", active = T),list(query = intersects,params = list("Surface water"), color = "#b7e2d6", active = T)),point.size = 3.5, line.size = 1, mainbar.y.label = "ASVs",sets.x.label = "ASVs",order.by = "freq",text.scale = c(4,4,4,4,4,4),keep.order = TRUE)
