This README file was generated on 2026-02-10 by Kelsey McGuire and updated on 2026-08-14 by Kelsey McGuire. 

**STUDY**: Aquatic vegetation and CH4 emission and production responses <br>

**DESCRIPTION/SUMMARY** <br>
Our study was conducted in the summer of 2025 in Kaska Ancestral Territory (Northern British Columbia, Canada) and our main objectives were to **(1)** understand how emergent vegetation influences overall methane flux and **(2)** how benthic vegetation influences potential methane production. 

Our field study was mainly to address objective 1, where we collected data across 18 sampling points, 9 vegetated and 9 open (i.e., non-vegetated). We sampled for both bubble and diffusive/plant-mediated flux using bubble traps and floating chambers respectively. Alongside flux measurements, we measured environmental conditions during each measurement, this included water temperature, % DO, and pH. The data from this sampling can be found in [DATA FILE]. 

Objective 2 was addressed through a laboratory incubation, where we collected sediment, benthic vegetation, and lake water samples during our field study and proceeded to run an incubation 6 months later. This incubation was run over the course of two-weeks where we collected headspace concentrations, and then nutrient and organic matter information post-hoc. The data from this experiment can be found in [DATA FILE].

This repository contains the data and code to reproduce the results from chapter 2 of Kelsey McGuire's thesis (2026): A little sedge goes a long way: emergent vegetation as a key driver of littoral methane flux.

**RESEARCH METHODOLOGY** <br>
Our study lake was sampled 2-3 times a week between July 2025 and August 2025. At each sampling point within the three zones, we captured both diffusive and plant-mediated flux alongside ebullitive flux using floating chambers and bubble traps respectively. Alongside measurements of flux, we also captured environmental variables such as depth, pH, % DO, temperature, and dissolved CH4 following the neonDissGas R package.

For chamber-based fluxes, we collected headspace samples from each sampling point once to twice a week across the sampling period. Samples were taken at 0, 15, and 30 minutes post chamber sealing using 30mL syringes with air-tight stopcocks, and for the final time step a 20mL vial that had been flushed with N2. We calculated slopes of CH4 accumulation and calculated an eventual flux following the equations from DelSontro et al., (2016). Plant-mediated flux was estimated following Desrosiers et al., (2022), where we subtracted open chamber flux values from vegetated chambers to estimate a value of plant-only flux. 

For ebullitive fluxes, bubble traps were left at all study locations throughout the sampling period and periodically sampled every one to nine days. Traps were checked for accumulated volume based on visual inspection, and when there was sufficient accumulated air, concentration samples were pulled from the traps and stored in 20mL N2 flushed vials. 

All concentration samples were run on a LI-7810 Trace Gas Analyzer (LICOR, Nebraska, USA), with chamber-based samples analyzed within 48 hours and ebullitive samples within 3 months due to the need for dilution to process samples properly. 

A subset of chamber and bubble trap samples were sent to Stanford University for d13C-CH4 composition on the Picarro G2201-I Isotopic Analyzer. 

To understand potential CH4 production in the presence of aquatic vegetation (benthic mosses), we established three treatments, sediment only, sediment + benthic, and benthic only, which was exposed to 20C and 30C. Headspace concentrations were sampled periodically across 10 days, and post-incubation two replicates from each treatment and temperature were analyzed for organic matter content and C/N ratios.

**SOFTWARE/TOOLS USED** <br>
- R version 4.6.0
- R packages used: here(), tidyverse(), janitor(), ggplot2(), lme4(), lmerTest(), glmmTMB(), emmeans(), ggsignif(), patchwork(), scales()

**DOI/PERSISTENT IDENTIFIERS** <br>
Hyperlinked DOI for this dataset: TBD

**LICENSE INFORMATION** <br>
The terms under which this dataset may be used follows the CC BY-SA 4.0 license. 

**AUTHOR/CREATOR INFORMATION** <br>
Position: MSc Student, Author <br>
Name: Kelsey McGuire <br>
ORCID: 0009-0003-0016-3838 <br>
Institution: University of British Columbia <br>
Address: https://kelsey-mcguire.github.io <br>
Email: kmcguire.9@outlook.com <br>

Position: Principal investigator, corresponding author <br>
Name: McKenzie Kuhn <br>
ORCID: 0000-0003-3871-1548 <br>
Institution: University of British Columbia <br>
Address: https://sites.chem.utoronto.ca/murphygroup/pi <br>
Email: mckenzie.kuhn@ubc.ca <br>

All experimental results can be reproduced using the code and data in this repository. Feel free to contact McKenzie Kuhn (corresponding author) by email at mckenzie.kuhn@ubc.ca if you have any questions about our work. 

**RELATED PUBLICATIONS** <br>
Journal article citation for this dataset: In preparation. <br>

Citation for Kelsey McGuire's thesis: TBD 

**FILE LIST/DIRECTORY** <br>
- `code_DataCleaning/` contains source code for the data cleaning from raw csv's which have field collected data. Multiple scripts are found in this directory including:
  - `EVCH4_CleanData.R` which contains all functions that are used in independently processing relevant csv's such as bubble trap, chamber-based, and vegetation specific data.
  - `EVCH4_Data.R` contains all subsets and data used during the processing of data for stats and figures, alongside the cleaning for incubation and isotope data which were treated separately to the field collected data.
  - `EVCH4_SlopeCalc.R` contains the function that calculated slopes from chamber-based LICOR data to find eventual flux measurements. 
- `code_PaperFiguresStats/` contains all the data used in this project. This subdirectory contains a data dictionary. Multiple scripts are found in this directory including:
  - `EVCH4_Figures.R` contains the code for the production of figures found within KM's thesis and the in-prep paper.
  - `EVCH4_ModelRunning.R` contains the functions for running repetitive model structures and filtering outputs based on significance values.
  - `EVCH4_Stats.R` contains the final model outputs and summaries for the statistics used in KM's thesis and in-prep paper.
- `data_csvs/` contains all the data used in this project. This subdirectory contains a Metadata file `Metadata_EmergentVegCH4_L&O_2026.docx` which holds all relevant information pertaining to data collection methodologies and explanations. Alongside the metadata file, you can find the csv's with relevant data from the project, named as follows: <br>
  - `Data1_ChamberCH4.csv`: Chamber-based flux data from the 2025 field campaign.
  - `Data2_BubbleTrapCH4.csv`: Bubble trap flux data from the 2025 field campaign.
  - `Data3_VegetationData.csv`: Vegetation data (i.e., biomass, greenness) from the 2025 field campaign.
  - `Data4_EnvironmentalData.csv`: Environmental  data (i.e., depth, pH, temperature) from the 2025 field campaign.
  - `Data5_CH4Isotopes.csv`: Isotope data from the 2025 field campaign.
  - `Data6_WaterChem.csv`: Lab analysed water chemistry data from the 2025 field campaign.
  - `Data7_Incubation.csv`: Incubation data from the potential CH4 production incubation in February 2026.
- `data_xlsx/` contains excel spreadsheets for the data used in this project. Files align with the above csv's, but have self-contained ReadME's as the first sheet within the .xlsx file. 

**ACKNOWLEDGEMENTS** <br>
We would like to acknowledge the help of Robbie Porter, Lisa Shephard, and Tanya Ball from the Dena Kayeh Institute, alongside Kaitlyn Chen, Graeme Morey, Samantha Solvis de Ovando and other members of the BAB Lab and UBC’s Department of Geography during field and laboratory work. The authors acknowledge the use of facilities and technical assistance of the Stanford Stable Isotope Biogeochemistry Laboratory (SIBL) Core Facility (RRID: SCR_023253). The research utilized the Picarro G2201-i Isotopic Analyzer, which was purchased with funding from Stanford University Community of Shared Advanced Research Platforms (C-ShaRP) (RRID:SCR_022986). Funding for this study was supported by the Natural Sciences and Engineering Research Council, University of British Columbia initiatives, and external partners. K.M. received funding through NSERC CGS-M, the Weston Family and Wildlife Conservation Society of Canada Boreal Research Fellowship, the Northern Scientist Training program, and UBC Graduate Student Research Awards. Further funding was provided through M.A.K.’s NSERC Discovery Grant and Northern Supplements.
