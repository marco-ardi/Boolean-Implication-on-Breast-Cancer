library(stringr)
library(dplyr)
library(purrr)
library(igraph)
options(warn=-1)

#uploading datasets and splitting by cancer subtype
brca_clinical <- read.table("Data/brca_clinical.txt", header = TRUE, sep = "\t", dec = ".")
brca_clinical <- na.omit(brca_clinical)
brca_clinical.LumA <- brca_clinical[brca_clinical$SUBTYPE=="LumA",]
brca_clinical.LumB <- brca_clinical[brca_clinical$SUBTYPE=="LumB",]
brca_clinical.Basal <- brca_clinical[brca_clinical$SUBTYPE=="Basal",]
#check.names=FALSE is MANDATORY, otherwise R transforms "-" in "." in IDs
brca_expressions <- read.table("Data/brca_expressions.txt", header = TRUE, sep = "\t", check.names=FALSE)
brca_proteomics <- read.table("Data/brca_proteomics.txt", header = TRUE, sep = "\t", dec = ".", check.names=FALSE)
#intersection in order to split

#mRNA
brca_expressions.LumA <- brca_expressions[1]
brca_expressions.LumB <- brca_expressions[1]
brca_expressions.Basal <- brca_expressions[1]

for(i in 1:length(brca_clinical.LumA$ID)){
  brca_expressions.LumA[i] <- brca_expressions[colnames(brca_expressions) %in% brca_clinical.LumA$ID[i]]
}
for(i in 1:length(brca_clinical.LumB$ID)){
brca_expressions.LumB[i] <- brca_expressions[colnames(brca_expressions) %in% brca_clinical.LumB$ID[i]]
}
for(i in 1:length(brca_clinical.Basal$ID)){
brca_expressions.Basal[i] <- brca_expressions[colnames(brca_expressions) %in% brca_clinical.Basal$ID[i]]
}
#proteins
brca_proteomics.LumA <- brca_proteomics[1]
brca_proteomics.LumB <- brca_proteomics[1]
brca_proteomics.Basal <- brca_proteomics[1]

for(i in 1:length(brca_clinical.LumA$ID)){
  brca_proteomics.LumA[i] <- brca_proteomics[colnames(brca_proteomics) %in% brca_clinical.LumA$ID[i]]
}
for(i in 1:length(brca_clinical.LumB$ID)){
  brca_proteomics.LumB[i] <- brca_proteomics[colnames(brca_proteomics) %in% brca_clinical.LumB$ID[i]]
}
for(i in 1:length(brca_clinical.Basal$ID)){
  brca_proteomics.Basal[i] <- brca_proteomics[colnames(brca_proteomics) %in% brca_clinical.Basal$ID[i]]
}
#creating files
#quote=FALSE is MANDATORY in order to not have names between ""
#mRNA
write.table(brca_expressions.LumA, "Data/brca_expressions.LumA.txt", sep="\t", dec=".", quote=FALSE)
write.table(brca_expressions.LumB, "Data/brca_expressions.LumB.txt", sep="\t", dec=".", quote=FALSE)
write.table(brca_expressions.Basal, "Data/brca_expressions.Basal.txt", sep="\t", dec=".", quote=FALSE)
#proteins
write.table(brca_proteomics.LumA, "Data/brca_proteomics.LumA.txt", sep="\t", dec=".", quote=FALSE)
write.table(brca_proteomics.LumB, "Data/brca_proteomics.LumB.txt", sep="\t", dec=".", quote=FALSE)
write.table(brca_proteomics.Basal, "Data/brca_proteomics.Basal.txt", sep="\t", dec=".", quote=FALSE)
#now we run java's script
#shell() for windows, system() for unix
shell(cmd="javac -cp . *.java")
#StepMiner on mRNA
shell(cmd="java StepMiner -i Data/brca_expressions.LumA.txt -o MyResults/brca_expressions.LumA.discr.txt")
shell(cmd="java StepMiner -i Data/brca_expressions.LumB.txt -o MyResults/brca_expressions.LumB.discr.txt")
shell(cmd="java StepMiner -i Data/brca_expressions.Basal.txt -o MyResults/brca_expressions.Basal.discr.txt")
#StepMiner on proteins
shell(cmd="java StepMiner -i Data/brca_proteomics.LumA.txt -o MyResults/brca_proteomics.LumA.discr.txt")
shell(cmd="java StepMiner -i Data/brca_proteomics.LumB.txt -o MyResults/brca_proteomics.LumB.discr.txt")
shell(cmd="java StepMiner -i Data/brca_proteomics.Basal.txt -o MyResults/brca_proteomics.Basal.discr.txt")
#BooleanNet on discretized mRNA matrixes
shell(cmd="java BooleanNet -i MyResults/brca_expressions.LumA.discr.txt -o MyResults/brca_expressions.LumA.BoolNet.txt -s 2.0")
shell(cmd="java BooleanNet -i MyResults/brca_expressions.LumB.discr.txt -o MyResults/brca_expressions.LumB.BoolNet.txt -s 2.0")
shell(cmd="java BooleanNet -i MyResults/brca_expressions.Basal.discr.txt -o MyResults/brca_expressions.Basal.BoolNet.txt -s 2.0")
#BooleanNet on discretized proteins matrixes
shell(cmd="java BooleanNet -i MyResults/brca_proteomics.LumA.discr.txt -o MyResults/brca_proteomics.LumA.BoolNet.txt -s 2.0")
shell(cmd="java BooleanNet -i MyResults/brca_proteomics.LumB.discr.txt -o MyResults/brca_proteomics.LumB.BoolNet.txt -s 2.0")
shell(cmd="java BooleanNet -i MyResults/brca_proteomics.Basal.discr.txt -o MyResults/brca_proteomics.Basal.BoolNet.txt -s 2.0")

#some helpful functions

#from txt booleanet file to df divided by genes
from_txt_to_df <- function(file, header=TRUE, sep = "\t", dec = ".", check.names=FALSE){
  prova <- read.table(file=file, header = header, sep = sep, dec = dec, check.names=check.names)
  prova <- prova[,1]
  prova <- str_split_fixed(prova, "AND",n=2)
  for(j in 1:length(prova)){
    if(prova[j] == "")
      prova[j] <- prova[-j]
    else
      prova[length(prova)+1] <- prova[j]
    
  }
  prova <- unique(prova) #here it gives a lot of warnings, but since we apply unique() it is ok
  df <- data.frame("Gene1"=prova, "Gene2"="")
  df <- str_split_fixed(df$Gene1, "=>",n=2)
  df <- data.frame(df)
  colnames(df)<- c('Gene1','Gene2')
  tmp_gene1 <- df$Gene1
  tmp_gene2 <- df$Gene2
  Gene1 <- str_replace_all(df$Gene1, " ([a-z])*", replace="")
  Gene2 <- str_replace_all(df$Gene2, " ([a-z])*", replace="")
  df <- data.frame(Gene1, Gene2)
  df
}
#takes genes values from dataframe
from_df_to_values <- function(dataframe){
  Gene1 <- str_replace_all(tmp_gene1, "([A-Z0-9])* ", replace="")
  Gene2 <- str_replace_all(tmp_gene2, "([A-Z0-9])* ", replace="")
  values <- data.frame(Gene1, Gene2)
  results <-rep("-1", length(values$Gene1))
  for(j in 1:length(values$Gene1)){
    if(values[j,1]=="low" && values[j,2]=="low")
      results[j] <- 1
    else if(values[j,1]=="low" && values[j,2]=="high")
      results[j] <- 2
    else if(values[j,1]=="high" && values[j,2]=="low")
      results[j] <- 3
    else if(values[j,1]=="high" && values[j,2]=="high")
      results[j] <- 4
  }
  results
}

#given 2 dataframes, genes and their kind of edge, displays the relative graph
display_network <-function(df, values, title){
  net <- graph.data.frame(df, directed=T)
  net <- set_edge_attr(net, "Legame", index=E(net), value=values)
  
  E(net)$color[E(net)$Legame == 1] <- 'green'
  E(net)$color[E(net)$Legame == 2] <- 'red'
  E(net)$color[E(net)$Legame == 3] <- 'blue'
  E(net)$color[E(net)$Legame == 4] <- 'orange'
  
  png(file=paste("Plots/",title, ".png"), width=1440, height=1440)
  plot(net, edge.arrow.size=1, edge.curved=.3, edge.width=2.5, vertex.size=6.5, main=title)
  legend("bottomright",
         c("Low=>Low", "Low=>High", "High=>Low", "High=>High"),
         col = c("green", "red", "blue", "orange"),
         cex = 2.0,
         lwd = 1, lty = 1)
  dev.off() #gives a warning but everything works
  net
}
#intersect and maintain colors
#works with 2 or 3 graphs
intersect_colors <- function(graph1, graph2, graph3=NULL, title){
  if(is.null(graph3)){
  intersezione <- intersection(graph1, graph2, keep.all.vertices=FALSE)
  E(intersezione)$color <- ifelse(is.na(E(intersezione)$color_1),
                                  E(intersezione)$color_2,E(intersezione)$color_1)
  }
  else{
    intersezione <- intersection(graph1, graph2, graph3, keep.all.vertices=FALSE)
    E(intersezione)$color <- ifelse(is.na(E(intersezione)$color_1),
                                    E(intersezione)$color_2,E(intersezione)$color_1)
  }
  png(file=paste("Plots/",title, ".png"), width=1440, height=1440)
  plot(intersezione, edge.arrow.size=1, edge.curved=.3, edge.width=2.5, vertex.size=6.5, 
       main=title)
  legend("bottomright",
         c("Low=>Low", "Low=>High", "High=>Low", "High=>High"),
         col = c("green", "red", "blue", "orange"),
         cex = 2.0,
         lwd = 1, lty = 1)
  dev.off()
  intersezione
  
}

#now let's test out function passing boolnet files
#mRNA expressions
df_BoolNet_expressions.LumA <- from_txt_to_df("MyResults/brca_expressions.LumA.BoolNet.txt")
values_BoolNet_expressions.LumA <- from_df_to_values(df_BoolNet_expressions.LumA)
df_BoolNet_expressions.LumB <- from_txt_to_df("MyResults/brca_expressions.LumB.BoolNet.txt")
values_BoolNet_expressions.LumB <- from_df_to_values(df_BoolNet_expressions.LumB)
df_BoolNet_expressions.Basal <- from_txt_to_df("MyResults/brca_expressions.Basal.BoolNet.txt")
values_BoolNet_expressions.Basal <- from_df_to_values(df_BoolNet_expressions.Basal)
#proteins values
df_BoolNet_proteomics.LumA <- from_txt_to_df("MyResults/brca_proteomics.LumA.BoolNet.txt")
values_BoolNet_proteomics.LumA <- from_df_to_values(df_BoolNet_proteomics.LumA)
df_BoolNet_proteomics.LumB <- from_txt_to_df("MyResults/brca_proteomics.LumB.BoolNet.txt")
values_BoolNet_proteomics.LumB <- from_df_to_values(df_BoolNet_proteomics.LumB)
df_BoolNet_proteomics.Basal <- from_txt_to_df("MyResults/brca_proteomics.Basal.BoolNet.txt")
values_BoolNet_proteomics.Basal <- from_df_to_values(df_BoolNet_proteomics.Basal)
#display expressions
net_expressions.LumA <- display_network(df_BoolNet_expressions.LumA, values_BoolNet_expressions.LumA, "LumA - Expressions")
net_expressions.LumB <- display_network(df_BoolNet_expressions.LumB, values_BoolNet_expressions.LumB, "LumB - Expressions")
net_expressions.Basal <- display_network(df_BoolNet_expressions.Basal, values_BoolNet_expressions.Basal, "Basal - Expressions")
#display proteins
net_proteomics.LumA <- display_network(df_BoolNet_proteomics.LumA, values_BoolNet_proteomics.LumA, "LumA - Proteomics")
net_proteomics.LumB <- display_network(df_BoolNet_proteomics.LumB, values_BoolNet_proteomics.LumB, "LumB - Proteomics")
net_proteomics.Basal <- display_network(df_BoolNet_proteomics.Basal, values_BoolNet_proteomics.Basal, "Basal - Proteomics")
#intersection between graphs
inter_LumA <- intersect_colors(net_expressions.LumA, net_proteomics.LumA, title="Intersection - LumA")
inter_LumB <- intersect_colors(net_expressions.LumB, net_proteomics.LumB, title="Intersection - LumB")
inter_Basal <- intersect_colors(net_expressions.Basal, net_proteomics.Basal, title="Intersection - Basal")

inter_expressions <- intersect_colors(net_expressions.LumA, net_expressions.LumB, net_expressions.Basal, title="Intersection - Expressions")
inter_proteomics <- intersect_colors(net_proteomics.LumA, net_proteomics.LumB, net_proteomics.Basal, title="Intersection - Proteins")
