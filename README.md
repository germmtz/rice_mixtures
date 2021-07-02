# rice_mixtures

Data & Code for the study "Crop mixtures: does niche complementarity hold for belowground resources? an experimental test using rice genotypic pairs"

Data:<br>
"Rice_traits.csv": this file contains trait and productivity data measured at the individual plant level. It has one row per plant and one column per trait.

Column headers:<br>
"IDplant": unique plant identifier (1 to 200)<br>
"IDpot": pot identifier with two plants per pot (1 to 100)<br>
"Bloc": bloc identifier, with 20 pots per bloc (A, B, C, D, E)<br>
"Treatment": P0 vs P+ = no P supply vs P supply<br>
"Asso": pot type, either monoculture (M) or mixture (P)<br>
"IDcouple": concatenation of the identifiers of the two genotypes in a pot (I64 = IR64, I64+=IR64 introgressed with QTL9, Pdi=Padi, Ktn=Ketan)<br>
"IDgeno": focal genotype identifier (I64 = IR64, I64+=IR64 introgressed with QTL9, Pdi=Padi, Ktn=Ketan)<br>
"IDnei": neighbour genotype identifier (I64 = IR64, I64+=IR64 introgressed with QTL9, Pdi=Padi, Ktn=Ketan)<br>
"BIOM_above": aboveground biomass (g)<br>
"Tillers": number of tillers<br>
"PH": Plant height (cm)<br>
"Biovolume": biovolume (m3)<br>
"SLA": Specific Leaf Area (m2/kg)<br>
"RB_top": Root biomass between 0 and 20 cm below the soil surface(g)<br>
"RB_deep: Root biomass between 20 and 60 cm below the soil surface(g) (!!! Only measured at the pot-level)<br>
"D_ad"/"D_bas": Mean root diameter (mm) of adventitious/basal roots, respectively<br>
"SRL_ad"/"SRL_bas": Specific Root Length (m/g) of adventitious/basal roots, respectively<br>
"RTD_ad"/"RTD_bas": Root Tissue Density (mg/cm3) of adventitious/basal roots, respectively<br>
"RBI_ad"/"RBI_bas": Root Branching Intensity (nb tips/cm) of adventitious/basal roots, respectively<br>
"PfR_ad"/"PfR_bas": Proportion of fine roots (diameter < 0.1 mm) (%) in adventitious/basal roots, respectively<br>

Code:<br>
"Rice_mixtures_analysis.R": this file contains the main statisticl analysis presented in the study. It uses "Rice_traits.csv" as an input. 
