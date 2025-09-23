base_dir <- switch(Sys.info()[["nodename"]],
                   "DESKTOP-6HPT8FH" = "C:/Abbie/research/seurat/prepostnatal",
                   "gauss" = "/home/abbiew/single_cell/velmeshev",
                   "mac.lan"     = "/Users/abbiesiru/Desktop/research/shen lab/autism/files/allen brain institute",
                   "."
)

meta <- read.csv(file.path(base_dir, "metadata_human_multiple_cortical_areas.csv"), sep = ",", header = TRUE)

length(unique(meta$sample_name))
### samples: 49417

length(unique(meta$external_donor_name_label))
unique(meta$external_donor_name_label)
### donors: 3
# "H200.1025" "H200.1030" "H200.1023"

length(unique(meta$cell_type_alias_label))
unique(meta$cell_type_alias_label)
### cell types: 20
# [1] ""                "VIP"             "LAMP5"           "IT"             
# [5] "PAX6"            "Oligodendrocyte" "Astrocyte"       "L5/6 IT Car3"   
# [9] "L5/6 NP"         "SST"             "L6 CT"           "OPC"            
# [13] "PVALB"           "L6b"             "Microglia"       "L5 ET"          
# [17] "Pericyte"        "Endothelial"     "L4 IT"           "VLMC"    

### cell subtypes: 121
# [1] ""                           "Inh L2-5 VIP TOX2"         
# [3] "Inh L1 LAMP5 GGT8P"         "Inh L1 LAMP5 NDNF"         
# [5] "Inh L1-3 VIP ZNF322P1"      "Inh L3 VIP CBLN1"          
# [7] "Inh L1-4 LAMP5 DUSP4"       "Exc L2-3 LINC00507 RPL9P17"
# [9] "Inh L1 SST CXCL14"          "Inh L1 PAX6 GRIP2"         
# [11] "Inh L1-2 VIP PPAPDC1A"      "Oligo L4-6 OPALIN"         
# [13] "Inh L1 PAX6 CA4"            "Inh L1 ADARB2 ADAM33"      
# [15] "Inh L1-4 VIP CHRNA2"        "Astro L1-6 FGFR3 ETNPPL"   
# [17] "Inh L2-6 VIP VIP"           "Inh L1-6 LAMP5 CA13"       
# [19] "Exc L5-6 THEMIS GPR21"      "Exc L5-6 FEZF2 MYBPHL"     
# [21] "Exc L4-5 RORB RPL31P31"     "Exc L4-5 RORB LCN15"       
# [23] "Inh L4-6 SST MTHFD2P6"      "Exc L6 THEMIS LINC00343"   
# [25] "Exc L6 FEZF2 FAM95C"        "Exc L4-5 RORB LINC01474"   
# [27] "OPC L1-6 MYT1"              "Inh L5-6 PVALB FAM150B"    
# [29] "Exc L6 FEZF2 KRT17"         "Inh L5 PVALB CNTNAP3P2"    
# [31] "Inh L5-6 LAMP5 SFTA3"       "Exc L5 RORB SNHG7"         
# [33] "Exc L3-4 RORB SEMA6D"       "Inh L1-5 VIP KCNJ2"        
# [35] "Inh L1-3 VIP SSTR1"         "Inh L1-3 PVALB WFDC2"      
# [37] "Astro L1 FGFR3 MT1G"        "Inh L1 VIP PRSS8"          
# [39] "Exc L3-4 RORB PRSS12"       "Inh L5-6 SST TH"           
# [41] "Inh L3-5 SST MAFB"          "Exc L5 RORB LINC01202"     
# [43] "Inh L5-6 SST ISOC1"         "Oligo L4-6 MOBP COL18A1"   
# [45] "Inh L5-6 SST KLHL14"        "Exc L5-6 FEZF2 CYP26B1"    
# [47] "Exc L3-5 RORB CMAHP"        "Micro L1-6 C1QC"           
# [49] "Exc L6 THEMIS C6orf48"      "Exc L5-6 THEMIS TMEM233"   
# [51] "Exc L5-6 RORB LINC00320"    "Exc L3-4 RORB FOLH1B"      
# [53] "Exc L6 FEZF2 TBC1D26"       "Exc L5 FEZF2 SCN7A"        
# [55] "Exc L6 FEZF2 SLITRK6"       "Inh L1-2 PAX6 SCGN"        
# [57] "Exc L6 FEZF2 P4HA3"         "Exc L5-6 FEZF2 ANKRD20A1"  
# [59] "Exc L4-5 RORB HNRNPA1P46"   "Inh L5-6 PVALB STON2"      
# [61] "Inh L6 LAMP5 C1QL2"         "Astro L1 FGFR3 FOS"        
# [63] "Exc L6 FEZF2 CPZ"           "Inh L6 SST NPY"            
# [65] "Inh L1-3 VIP GGH"           "Inh L2-4 PVALB C8orf4"     
# [67] "Exc L3 RORB CARTPT"         "Exc L6 FEZF2 ETV4"         
# [69] "Exc L5-6 FEZF2 CABP7"       "Inh L6 LHX6 GLP1R"         
# [71] "Inh L1-6 PVALB SCUBE3"      "Inh L2-4 VIP DSEL"         
# [73] "Inh L6 LAMP5 ANKRD20A11P"   "Inh L1 VIP SOX11"          
# [75] "Inh L1 ADARB2 DISP2"        "Inh L1 VIP PCDH20"         
# [77] "Inh L3-6 VIP KCTD13"        "Exc L5-6 THEMIS THTPA"     
# [79] "Inh L2-4 SST AHR"           "Inh L3-6 PVALB MFI2"       
# [81] "Peri L1-6 MUSTN1"           "Inh L3-4 PVALB HOMER3"     
# [83] "Inh L2-4 VIP LGI2"          "Exc L5-6 THEMIS OR1J1"     
# [85] "Endo L2-5 CLDN5"            "Inh L1-2 VIP RPL41P3"      
# [87] "Inh L1-6 VIP RGS16"         "Inh L1-3 VIP ACHE"         
# [89] "Inh L1 VIP TNFAIP8L3"       "Exc L3-5 RORB HSPB3"       
# [91] "Exc L3-5 THEMIS ELOF1"      "Exc L2-4 RORB GRIK1"       
# [93] "Inh L1-6 VIP PENK"          "Inh L1-2 PVALB TAC1"       
# [95] "Inh L1-3 PAX6 NABP1"        "Inh L1-3 VIP CCDC184"      
# [97] "Inh L1-6 VIP RCN1"          "Exc L3-5 FEZF2 ONECUT1"    
# [99] "Exc L6 FEZF2 TBCC"          "Exc L3-5 RORB CD24"        
# [101] "VLMC L1-3 CYP1B1"           "Exc L5-6 THEMIS IL7R"      
# [103] "Exc L4 RORB BHLHE22"        "Exc L4 RORB CACNG5"        
# [105] "Exc L4-5 RORB AIM2"         "Exc L4-6 RORB HPCA"        
# [107] "Exc L6 THEMIS EGR3"         "Exc L6 FEZF2 VWA2"         
# [109] "Inh L4-5 PVALB TRIM67"      "Exc L4-5 RORB ASCL1"       
# [111] "Exc L5-6 FEZF2 RSAD2"       "Exc L3 LINC00507 PSRC1"    
# [113] "Exc L3-4 RORB RPS3P6"       "Exc L5 FEZF2 MORN2"        
# [115] "Exc L3-5 LINC00507 SLN"     "Exc L3-5 THEMIS UBE2F"     
# [117] "Exc L3-5 FEZF2 DCN"         "Exc L4 RORB CCDC168"       
# [119] "Exc L3 LINC00507 CTXN3"     "Exc L3 THEMIS PLA2G7"      
# [121] "Exc L5 FEZF2 DYRK2"  

length(unique(meta$region_label))
unique(meta$region_label)
### regions: 8
# "MTG"  "V1C"  "CgG"  "M1lm" "S1ul" "S1lm" "M1ul" "A1C" 

length(unique(meta$cortical_layer_label))
unique(meta$cortical_layer_label)
### cortical layers: 13
# [1] "L1"   "L5"   "L4"   "L2"   "L3"   "L6"   "L4ab" "L4c"  "L6a"  "L6b"  "L5a"  "L5b" 
# [13] "WM" 
