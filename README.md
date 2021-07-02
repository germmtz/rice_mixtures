# rice_mixtures

Data & Code for the study "Crop mixtures: does niche complementarity hold for belowground resources? an experimental test using rice genotypic pairs"

Data:
"Rice_traits.csv": this file contains trait and productivity data measured at the individual plant level. It has one row per plant and one column per trait.

Column headers:
"IDplant": unique plant identifier (1 to 200)
"IDpot": pot identifier with two plants per pot (1 to 100)
"Bloc": bloc identifier, with 20 pots per bloc (A, B, C, D, E)
"Treatment": P0 vs P+ = no P supply vs P supply
"Asso": pot type, either monoculture (M) or mixture (P)
"IDcouple": concatenation of the identifiers of the two genotypes in a pot (I64 = IR64, I64+=IR64 introgressed with QTL9, Pdi=Padi, Ktn=Ketan)
"IDgeno": focal genotype identifier (I64 = IR64, I64+=IR64 introgressed with QTL9, Pdi=Padi, Ktn=Ketan)
"IDnei": neighbour genotype identifier (I64 = IR64, I64+=IR64 introgressed with QTL9, Pdi=Padi, Ktn=Ketan)
"BIOM_above": aboveground biomass (g)
"Tillers": number of tillers
"PH": Plant height (cm)
"Biovolume": biovolume (m3)
"SLA": Specific Leaf Area (m2/kg)
"RB_top": Root biomass between 0 and 20 cm below the soil surface(g)
"RB_deep: Root biomass between 20 and 60 cm below the soil surface(g) (!!! Only measured at the pot-level)
"D_ad"/"D_bas": Mean root diameter (mm) of adventitious/basal roots, respectively
"SRL_ad"/"SRL_bas": Specific Root Length (m/g) of adventitious/basal roots, respectively
"RTD_ad"/"RTD_bas": Root Tissue Density (mg/cm3) of adventitious/basal roots, respectively
"RBI_ad"/"RBI_bas": Root Branching Intensity (nb tips/cm) of adventitious/basal roots, respectively
"PfR_ad"/"PfR_bas": Proportion of fine roots (diameter < 0.1 mm) (%) in adventitious/basal roots, respectively

Code:
"Rice_mixtures_analysis.R": this file contains the main statisticl analysis presented in the study. It uses "Rice_traits.csv" as an input. 
