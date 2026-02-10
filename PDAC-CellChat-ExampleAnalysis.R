library(CellChat)
library(Seurat)

## CellChat Object
data <- subset(data.combined, subset = nclust == "niche1")
cellchat <- createCellChat(object = data,
                           group.by = "labels", 
                           datatype = "spatial")

cellchat@DB <- subsetDB(CellChatDB.human)

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)

cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = TRUE, interaction.range = 250, scale.distance = 2,
                              contact.dependent = TRUE, contact.range = 100)
cellchat <- aggregateNet(cellchat)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Setup Colors
cellchat.colors <- ccolss
names(cellchat.colors) <- paste0("cluster", names(cellchat.colors))

# Interaction Network
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), 
                   weight.scale = T, label.edge= F, color.use = cellchat.colors, top = 0.05,
                   title.name = "Interaction Weights/Strength")

# Cluster x Cluster Heatmap
netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues", color.use = cellchat.colors)


# Pathway Heatmap - All
netAnalysis_signalingRole_heatmap(cellchat, pattern = "all", color.use = cellchat.colors)

# Pathway Heatmap - Incoming
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", color.use = cellchat.colors)

# Pathway Heatmap - Outgoing
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", color.use = cellchat.colors)

# Network Interactions per Cluster (qgraph)
library(qgraph)
target_cluster = "cluster4"

df <- as.data.frame(cellchat@net$weight) %>%
    tibble::rownames_to_column("source") %>%
    pivot_longer(-source, names_to = "target", values_to = "weight")

cluster.df <- df %>% 
    filter(source == cl.name | target == cl.name) %>%
    as.matrix()

qgraph(cluster.df,
           edgelist = T, directed = T, arrows = T,
           labels = unique(df$source), loopRotation = pi,
           layout = "spring", layoutScale = c(0.75, 0.75)




