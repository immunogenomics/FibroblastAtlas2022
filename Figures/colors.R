palette_global <- c(
    Mean = 'black',
    
    ## Greens
    lung = '#31a354',
    Lung = '#31a354',
    LungDisease = '#31a354',
    `Lung Disease` = '#31a354',
    LungControl = '#a1d99b',
    `Early-stage ILD` = '#31a354',
    `End-stage IPF/RA-ILD` = '#69be77',

    
    ## Reds
    synovium = '#de2d26',
    Synovium = '#de2d26',
    RheumatoidArthritis = '#de2d26',
    Osteoarthritis = '#fc9272',

#     ## Purples
#     gut = '#756bb1',
#     Gut = '#756bb1',
#     GutInflamed = '#756bb1',
#     GutNonisnflamed = '#bcbddc',
#     GutControl = '#efedf5',
#     Inflamed = '#756bb1',
#     `Non-inflamed` = '#bcbddc',

    ## Yellows 
    gut = '#edc948',
    Gut = '#edc948',
    GutInflamed = '#edc948',
    GutNonisnflamed = '#F8DB89',
    GutControl = '#FEEDC4',
    Inflamed = '#edc948',
    `Non-inflamed` = '#F8DB89',
    
    
    ## Lip
    salivarygland = '#3182bd',
    SalivaryGland = '#3182bd',
    lip = '#3182bd',
    Lip = '#3182bd',
    PSS = '#3182bd',
    pSS = '#3182bd',
    SICCA = '#9ecae1',
    Sicca = '#9ecae1',

    ## Greys
    Dermal = '#2c7fb8',
    Lesional = '#2c7fb8',
    `Non-lesional` = '#7fcdbb',
#     Healthy = '#edf8b1',
    Healthy = '#feb24c',
    

    ## Cell type lineages
    Endothelial = '#4E79A7',
    endothelial = '#4E79A7',
    vascular_endothelial = '#4E79A7',
    Leukocyte = '#A0CBE8',
    Immune = '#A0CBE8',
    hematopoietic = '#A0CBE8',
    Mural = '#F28E2B',
    mural = '#F28E2B',
    Fibroblast = '#FFBE7D',
    fibroblast = '#FFBE7D',
    Stromal = '#FFBE7D',
    lymphatic_endothelial = '#499894', 
    Epithelial = '#59A14F',
    epithelial = '#59A14F',
    Plasma = '#8CD17D',
    Glial = '#499894',
    proliferating = 'black',
    
    
    ## Tissue specific fibroblast subtypes 
    `Lining (SC-F4)` = '#4E79A7',
    `Sublining` = '#A0CBE8',
    `CD34+ sublining (SC-F1)` = '#F28E2B',
    `DKK3+ sublining (SC-F3)` = '#FFBE7D',
    `HLA-DRAhi sublining (SC-F2)` = '#59A14F',
    `Fibroblast` = '#8CD17D',
    `Myofibroblast` = '#E15759',
    myofibroblast = '#E15759',
    

    
    `HAS1_PLIN2` = '#F1CE63',
    `CD34` = '#499894',
    `CCL19` = '#FABFD2',
    `Myofibroblasts` = '#E15759',
    `WNT5B+ 1` = '#FF9D9A',
    `WNT2B+ Fos-lo 2` = '#79706E',
    `Inflammatory Fibroblasts` = '#BAB0AC',
    `WNT2B+ Fos-hi` = '#D37295',
    `RSPO3+` = '#86BCB6',
    `WNT5B+ 2` = '#B07AA1',
    
    ## Joint fibroblast clusters
    `5` = '#4E79A7',
    `2` = '#A0CBE8',
    `4` = '#F28E2B',
    `1` = '#FFBE7D',
    `9` = '#59A14F',
    `11` = '#8CD17D',
    `7` = '#B6992D',
    `6` = '#F1CE63',
    `8` = '#499894',
    `3` = '#86BCB6',
    `12` = '#E15759',
    `10` = '#FF9D9A',
    `0` = '#79706E',
    `13` = '#BAB0AC',
    `FBLN1+ C5` = '#4E79A7',
    `C5` = '#4E79A7',
    `C2` = '#A0CBE8',
    `C4` = '#F28E2B',
    `SPARC+COL3A1+ C4` = '#F28E2B',
    `C1` = '#FFBE7D',
    `C9` = '#59A14F',
    `CD34+MFAP5+ C9` = '#59A14F',
    `C11` = '#8CD17D',
    `CXCL10+CCL19+ C11` = '#8CD17D',
    `C7` = '#B6992D',
    `C6` = '#F1CE63',
    `C8` = '#499894',
    `PTGS2+SEMA4A+ C8` = '#499894',
    `C3` = '#86BCB6',
    `C12` = '#E15759',
    `C10` = '#FF9D9A',
    `C0` = '#79706E',
    `C13` = '#BAB0AC',
    `MYH11+ C13` = '#BAB0AC',

    ## CellDive Subtypes
    `Other` = 'lightgrey', 
    `CCL19+` = '#8CD17D',
    `CCL19+ Fibroblast` = '#8CD17D',
    `SPARC+` = '#F28E2B',
    `SPARC+ Fibroblast` = '#F28E2B',
    `Other Fibroblast` = 'black', 
    `Other Cell` = 'lightgrey', 
    

    
    
    ## Inflam Subtype pathways
    `C4_GO` = '#F28E2B',
    `Gene Ontology (C4)` = '#F28E2B',
    `C11_GO` = '#8CD17D',
    `Gene Ontology (C11)` = '#8CD17D',
    Cytokines = '#9D7660',
    `Response to Cytokines` = '#9D7660',
    Morphogens = '#D7B5A6',
    `Response to Morphogens` = '#D7B5A6',
    
    ## Mouse time data
    Early = '#fec44f',
    early = '#fec44f',
    Acute = '#d95f0e',
    acute = '#d95f0e',
    Recovery = '#1c9099',

    ## Donors
    `BRI124` = '#F6FBF7',
    `BRI186` = '#EDF7F0',
    `BRI188` = '#E4F3E8',
    `BRI189` = '#DBEFE1',
    `BRI190` = '#D2EBD9',
    `BRI191` = '#C9E7D2',
    `BRI192` = '#C0E3CA',
    `BRI193` = '#B7DFC3',
    `BRI194` = '#AEDBBC',
    `BRI195` = '#A5D7B4',
    `BRI201` = '#9CD3AD',
    `BRI202` = '#93CFA5',
    `BRI204R` = '#8ACB9E',
    `BRI205` = '#81C696',
    `BRI207` = '#78C38F',
    `BRI208` = '#6FBF88',
    `BRI209` = '#66BB80',
    `BRI210` = '#5DB779',
    `BRI212` = '#54B371',
    `BRI213` = '#4BAF6A',
    `BRI214` = '#42AB62',
    `BRI215` = '#39A65B',
    `BRI216` = '#31A354',
#     `BRI013` = '#F6F6FA',
#     `BRI015` = '#EEEDF5',
#     `BRI019` = '#E6E4F1',
#     `BRI021` = '#DEDCEC',
#     `BRI111` = '#D6D3E8',
#     `BRI112` = '#CECAE3',
#     `BRI114` = '#C6C2DE',
#     `BRI133` = '#BEB9DA',
#     `BRI134` = '#B5B0D5',
#     `BRI135` = '#ADA7D1',
#     `BRI136` = '#A59FCC',
#     `BRI137` = '#9D96C7',
#     `BRI138` = '#958DC3',
#     `BRI139` = '#8D85BE',
#     `BRI140` = '#857CBA',
#     `BRI152` = '#7D73B5',
#     `BRI153` = '#756BB1',
    `BRI013` = "#FFFCF4",
    `BRI015` = "#FFF8EA",
    `BRI019` = "#FFF5DF",
    `BRI021` = "#FFF2D5",
    `BRI111` = "#FFEFCA",
    `BRI112` = "#FEEBC0",
    `BRI114` = "#FDE8B5",
    `BRI133` = "#FCE5AB",
    `BRI134` = "#FBE2A1",
    `BRI135` = "#FADF96",
    `BRI136` = "#F8DB8C",
    `BRI137` = "#F7D881",
    `BRI138` = "#F5D576",
    `BRI139` = "#F3D26B",
    `BRI140` = "#F1CF60",
    `BRI152` = "#EFCC54",
    `BRI153` = "#EDC948",
    `180116A` = '#FDF5F4',
    `180116B` = '#FBEBEA',
    `180123A` = '#FAE1E0',
    `BRI001` = '#F8D7D5',
    `BRI003` = '#F7CDCB',
    `BRI004` = '#F5C3C1',
    `BRI006` = '#F4B9B6',
    `BRI007` = '#F2AFAC',
    `BRI009` = '#F0A5A2',
    `BRI011` = '#EF9B97',
    `BRI039` = '#ED918D',
    `BRI041` = '#EC8783',
    `BRI072` = '#EA7D78',
    `BRI073` = '#E9736E',
    `BRI075` = '#E76964',
    `BRI077` = '#E55F59',
    `BRI083` = '#E4554F',
    `BRI106` = '#E24B45',
    `BWH075` = '#E1413A',
    `BWH076CD45n` = '#DF3730',
    `BWH078` = '#DE2D26',
    `GX44` = '#EFF5F9',
    `GX46` = '#DFEBF4',
    `GX47` = '#CFE2EF',
    `GX48` = '#BFD8EA',
    `GX50` = '#AFCEE5',
    `GX09` = '#9FC5E0',
    `GX21` = '#90BBDB',
    `GX33` = '#80B2D6',
    `GX45` = '#70A8D1',
    `GX57` = '#609ECC',
    `GX69` = '#5095C7',
    `GX81` = '#408BC2',
    `GX93` = '#3182BD',    
    
    ## CellDive Donors
    GI6645 = '#E69F00',
    GI6717 = '#56B4E9',
    GI6846 = '#009E73',
    JPR118 = '#F0E442',
    JPR125 = '#0072B2',
    S445250 = '#D55E00',    
    
    ## Public cross-tissue datasets
	Endothelial = '#4E79A7',
	Stromal = '#F28E2B',
	Epithelial = '#E15759',
	Immune = '#76B7B2',
	Mural_Muscle = 'black',
	Vasculature = '#4E79A7',
	Muscle = '#A0CBE8',
	Bladder = '#F28E2B',
	Lung = '#FFBE7D',
	Thymus = '#59A14F',
	Trachea = '#8CD17D',
	Large_Intestine = '#B6992D',
	Kidney = '#F1CE63',
	Ovary = '#499894',
	Pancreas = '#86BCB6',
	Uterus = '#E15759',
	Esophagus = '#FF9D9A',
	Heart = '#79706E',
	Rectum = '#BAB0AC',
	Skin = '#D37295',
	Mammary_Gland = '#FABFD2',
	Limb_Muscle = '#B07AA1',
    
    ## Stim conditions
    ECs = '#fff7bc', 
    Tcells = '#2b8cbe', 
    
    ## Niches
    `Other Niche` = 'lightgrey', 
    `Lymphoid Niche` = '#4E79A7', 
    `Vascular Niche` = '#F28E2B', 
    `Mural Niche` = '#76B7B2',
    
    ## Generic
    Yes = '#de2d26', ## Red
    No = '#bdbdbd', ## Light grey 

    `NULL` = '#000000'
    
)


## Data structures for ComplexHeatmap
palette_heatmap <- list(
    Tissue = palette_global[c(
        ## Our tissues
        'Lung', 'Gut', 'Synovium', 'SalivaryGland', 'Lip', 
        'lung', 'gut', 'synovium', 'salivarygland', 'lip',
        
        ## AHCA, TM, TS, and MCA tissues 
        'Vasculature','Muscle','Bladder','Lung','Thymus','Trachea',
        'Large_Intestine','Kidney','Ovary','Pancreas','Uterus','Esophagus',
        'Heart','Rectum','Skin','Mammary_Gland','Limb_Muscle'
        
    )],
    Cluster = palette_global[c(
        '5', '2', '4', '1', '9', '11', '7', '6', '8', '3', '12', '10', '0', '13', 
        'FBLN1+ C5', 'C2', 'SPARC+COL3A1+ C4', 'C1', 'CD34+MFAP5+ C9', 
        'CXCL10+CCL19+ C11', 'C7', 'C6', 'PTGS2+SEMA4A+ C8', 'C3', 'C12', 'C10', 'C0', 'MYH11+ C13',
        'Lining (SC-F4)',
        'Sublining',
        'CD34+ sublining (SC-F1)',
        'DKK3+ sublining (SC-F3)',
        'HLA-DRAhi sublining (SC-F2)',
        'Fibroblast',
        'Myofibroblast',
        'HAS1_PLIN2',
        'CD34',
        'CCL19',
        'Myofibroblasts',
        'WNT5B+ 1',
        'WNT2B+ Fos-lo 2',
        'Inflammatory Fibroblasts',
        'WNT2B+ Fos-hi',
        'RSPO3+',
        'WNT5B+ 2'
        
    )],
    Cluster_name = palette_global[c(
        '5', '2', '4', '1', '9', '11', '7', '6', '8', '3', '12', '10', '0', '13', 
        'FBLN1+ C5', 'C2', 'SPARC+COL3A1+ C4', 'C1', 'CD34+MFAP5+ C9', 
        'CXCL10+CCL19+ C11', 'C7', 'C6', 'PTGS2+SEMA4A+ C8', 'C3', 'C12', 'C10', 'C0', 'MYH11+ C13'
    )],
    Label = palette_global[c(
        'C4_GO', 'C11_GO', 'Cytokines', 'Morphogens',
        'Gene Ontology (C11)', 'Response to Cytokines', 'Gene Ontology (C4)', 'Response to Morphogens'
    )],
    CellType = palette_global[c(  
        'Endothelial','Stromal','Epithelial','Immune','Mural_Muscle'
    )],
    Condition = palette_global[c(
        'ECs', 'Tcells'
    )],
    Niche = palette_global[c(
        'Other Niche', 'Lymphoid Niche', 'Vascular Niche', 'Mural Niche'
    )], 
    LibraryID = palette_global[c(
        'GI6645', 'GI6717', 'GI6846', 'JPR118', 'JPR125', 'S445250'    
    )],
    Subtype = palette_global[c(
        'CCL19+', 'CCL19+ Fibroblast', 'SPARC+ Fibroblast', 'SPARC+', 'Other Fibroblast', 'Other', 'Other Cell'
    )]
)


scale_colour_discrete <- function(...) {
#   scale_colour_manual(..., values = viridis::viridis(2))
    scale_colour_manual(..., values = palette_global)
}

scale_fill_discrete <- function(...) {
#   scale_colour_manual(..., values = viridis::viridis(2))
    scale_fill_manual(..., values = palette_global)
}



# x <- unique(dplyr::select(obj$meta_data, Cluster, Cluster_name))
# colors <- tableau_color_pal('Tableau 20')(14)

# writeLines(glue("`{x$Cluster}` = \'{colors}\',"))
# writeLines(glue("`{x$Cluster_name}` = \'{colors}\',"))


# writeLines(Reduce(paste0, glue("\'{x$Cluster}\', ")))

# writeLines(Reduce(paste0, glue("\'{x$Cluster_name}\', ")))

