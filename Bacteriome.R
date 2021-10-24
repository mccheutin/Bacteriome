### Phylo-tree------
setwd("/Users/marie-charlottecheutin/Drive/Thesis/Presentations/Bacteriome")
load("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/Physeq_objects/gut_core.RData")
fish_tree <- read.tree(file = '/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/phylogeny/fish_tree_modif.tre')

# Species by colored diet
db <- sample_data(gut_core)
write.table(db, file ="./db.txt", sep ="\t", row.names = F)
db <- read.table( file ="./db.txt", sep ="\t", header = T)
library(ggtree)
grp <- list(FC = levels(factor(db[db$diet1 == "FC",]$tax2)),
            MI = levels(factor(db[db$diet1 == "MI",]$tax2)),
            SI = levels(factor(db[db$diet1 == "SI",]$tax2)),
            OM = levels(factor(db[db$diet1 == "OM",]$tax2)),
            PK =  levels(factor(db[db$diet1 == "PK",]$tax2)),
            HD =  levels(factor(db[db$diet1 == "HD",]$tax2)),
            H =  levels(factor(db[db$diet1 == "H",]$tax2)))


figtree <- ggtree(fish_tree, layout = 'fan',open.angle = 180, branch.length='none')
pdf(file = "./tree_fish_diet.pdf" , he= 15, wi= 13)
groupOTU(figtree, grp, 'Diet') + 
  geom_tiplab(size= 4, aes(angle=angle), hjust = -0.05)+
  geom_tippoint(aes(color = Diet) ,size= 3, alpha=.75)+
  scale_color_manual(values = c("darkblue",
                                "darkred",
                                "darkkhaki",
                                "darkgreen",
                                "darkorange", 
                                "darkblue",
                                "gray",
                                "darkorchid")) +
  theme_tree2(legend.position="right")
dev.off()  



## Constitution ----
core_physeq_order <- gut_core %>%
  tax_glom(taxrank = "Order") %>%                     # agglomerate at order level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Order)

core_physeq_order$Abundance <- core_physeq_order$Abundance/378
core_physeq_order %>%  filter(Order == "Desulfovibrionales" ) %>% select(Abundance) %>% sum()

core_physeq_order$Phylum = as.character(core_physeq_order$Phylum) # Avoid error message with factor for next step
sort(table(factor(core_physeq_order$Phylum)), T)
core_physeq_order[!core_physeq_order$Phylum %in% 
                    c("Proteobacteria","Bacteroidetes","Cyanobacteria","Firmicutes",
                      "Planctomycetes", "Spirochaetes", "Verrucomicrobia","Fusobacteria","Tenericutes"),
                  which(names(core_physeq_order) == "Phylum", T)] <- "Other" #Change phylum to other for those not included in the list

group <-  core_physeq_order$Phylum
subgroup <- core_physeq_order$Order
value <- core_physeq_order$Abundance

gut_core_treemap_data=data.frame(group,subgroup,value)


## SAmples ----
db <- read.table( file ="./db.txt", sep ="\t", header = T)
freq <- as.data.frame(table(db$tax1))
colnames(freq) <- c("tax1" , "N")
db_frq <- db %>% select(tax1 , family, diet2, diet3, region)
db_frq$region <- str_replace_all(db_frq$region , 
                                 c("Europa" = "Indian ocean", 
                                   "Juan_de_nova" = "Indian ocean",
                                   "Seychelles" = "Indian ocean",
                                   "Martinique" = "Caribbean sea"))

db_frq <- db_frq %>% unique()
db_frq <- left_join(freq, db_frq, by = "tax1")
write.table(db_frq , file = "db_frq.txt", sep ="\t", row.names = F)


## Core ----
load("/Users/marie-charlottecheutin/Drive/Thesis/Presentations/Bacteriome/global_core.RData")
load("/Users/marie-charlottecheutin/Drive/Thesis/Presentations/Bacteriome/sw_core.RData")
load("/Users/marie-charlottecheutin/Drive/Thesis/Presentations/Bacteriome/sediment_core.RData")

source("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram3.r")

merged = merge_samples(global_core, "type") #merge samples for herbivores, carnivores and the macroalgae
rownames(merged@otu_table@.Data)
min(rowSums(merged@otu_table@.Data)) #how many reads per sample
set.seed(10000)
venn = rarefy_even_depth(merged, sample.size = 290704)

gut <- otu_table(subset_samples(global_core, type == "gut"))
sand <- otu_table(subset_samples(global_core, type == "sand"))
sw <- otu_table(subset_samples(global_core, type == "seawater"))

venn_comp <- venn_diagram3( gut,sand,sw, "gut", "sand", "sw", colors= c("darkred","gray","darkblue"), euler=F)



