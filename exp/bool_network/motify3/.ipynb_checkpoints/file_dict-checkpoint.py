from motify_functions import *
import sys
sys.path.append("../..")
import json
file_dict = {}
directory = Path('cellcollective')
counter = 0
for file_path in directory.rglob('*'):
    if file_path.is_file():
        file_dict[counter] = str(file_path)
        counter += 1
file_dicts = {0: 'cellcollective/Apoptosis Network_19422837.txt', 
              1: 'cellcollective/Arabidopsis thaliana Cell Cycle_26340681.txt', 
              2: 'cellcollective/Aurora Kinase A in Neuroblastoma_26616283.txt', 
              3: 'cellcollective/Body Segmentation in Drosophila 2013_23520449.txt', 
              4: 'cellcollective/B cell differentiation_26751566.txt', 
              5: 'cellcollective/B bronchiseptica and T retortaeformis coinfection_22253585.txt', 
              6: 'cellcollective/Bordetella bronchiseptica_22253585.txt', 
              7: 'cellcollective/Bortezomib Responses in U266 Human Myeloma Cells_26163548.txt', 
              8: 'cellcollective/BT474 Breast Cell Line Long-term ErbB Network_24970389.txt', 
              9: 'cellcollective/BT474 Breast Cell Line Short-term ErbB Network_24970389.txt', 
              10: 'cellcollective/Budding Yeast Cell Cycle_23049686.txt', 
              11: 'cellcollective/Budding Yeast Cell Cycle 2009_19185585.txt', 
              12: 'cellcollective/Cardiac development_23056457.txt', 
              13: 'cellcollective/CD4 T cell signaling_25538703.txt', 
              14: 'cellcollective/CD4+ T Cell Differentiation and Plasticity_26090929.txt', 
              15: 'cellcollective/Cell Cycle Transcription by Coupled CDK and Network Oscillators_18463633.txt', 
              16: 'cellcollective/CD4+ T cell Differentiation_22871178.txt', 
              17: 'cellcollective/Cholesterol Regulatory Pathway_19025648.txt', 
              18: 'cellcollective/Cortical Area Development_20862356.txt', 
              19: 'cellcollective/Death Receptor Signaling_20221256.txt', 
              20: 'cellcollective/Colitis-associated colon cancer_26446703.txt', 
              21: 'cellcollective/Differentiation of T lymphocytes_23743337.txt', 
              22: 'cellcollective/EGFR & ErbB Signaling_19662154.txt', 
              23: 'cellcollective/ErbB (1-4) Receptor Signaling_23637902.txt', 
              24: 'cellcollective/FA BRCA pathway_22267503.txt', 
              25: 'cellcollective/Fanconi anemia and checkpoint recovery_26385365.txt', 
              26: 'cellcollective/FGF pathway of Drosophila Signaling Pathways_23868318.txt', 
              27: 'cellcollective/fMRI - Regulation of the Lac Operon_25790483.txt', 
              28: 'cellcollective/Glucose Repression Signaling 2009_19144179.txt', 
              29: 'cellcollective/Guard Cell Abscisic Acid Signaling_16968132.txt', 
              30: 'cellcollective/HCC1954 Breast Cell Line Long-term ErbB Network_24970389.txt', 
              31: 'cellcollective/HCC1954 Breast Cell Line Short-term ErbB Network_24970389.txt', 
              32: 'cellcollective/HGF Signaling in Keratinocytes_22962472.txt', 
              33: 'cellcollective/HH Pathway of Drosophila Signaling Pathways_23868318.txt', 
              34: 'cellcollective/HIV-1 interactions with T-Cell Signaling Pathway_25431332.txt', 
              35: 'cellcollective/Human Gonadal Sex Determination_26573569.txt', 
              36: 'cellcollective/IGVH mutations in chronic lymphocytic leukemia_26088082.txt', 
              37: 'cellcollective/IL-1 Signaling_21968890.txt', 
              38: 'cellcollective/IL-6 Signaling_21968890.txt', 
              39: 'cellcollective/Inflammatory Bowel Disease (IBD) Model_29513758.txt', 
              40: 'cellcollective/Influenza A Virus Replication Cycle_23081726.txt', 
              41: 'cellcollective/Iron acquisition and oxidative stress response in aspergillus fumigatus_25908096.txt', 
              42: 'cellcollective/Lac Operon_21563979.txt', 
              43: 'cellcollective/Lymphopoiesis Regulatory Network_26408858.txt', 
              44: 'cellcollective/Lymphoid and myeloid cell specification and transdifferentiation_28584084.txt', 
              45: 'cellcollective/Mammalian Cell Cycle 2006_16873462.txt', 
              46: 'cellcollective/Mammalian Cell Cycle_19118495.txt', 
              47: 'cellcollective/manual_mod_Treatment of Castration-Resistant Prostate Cancer_28361666.txt', 
              48: 'cellcollective/MAPK Cancer Cell Fate Network_24250280.txt', 
              49: 'cellcollective/Neurotransmitter Signaling Pathway_17010384.txt', 
              50: 'cellcollective/Metabolic Interactions in the Gut Microbiome_26102287.txt', 
              51: 'cellcollective/Oxidative Stress Pathway_23134720.txt', 
              52: 'cellcollective/Predicting Variabilities in Cardiac Gene_26207376.txt', 
              53: 'cellcollective/PC12 Cell Differentiation_27148350.txt', 
              54: 'cellcollective/Processing of Spz Network from the Drosophila Signaling Pathway_23868318.txt', 
              55: 'cellcollective/Regulation of the L-arabinose operon of Escherichia coli_28639170.txt', 
              56: 'cellcollective/Pro-inflammatory Tumor Microenvironment in Acute Lymphoblastic Leukemia_27594840.txt', 
              57: 'cellcollective/Senescence Associated Secretory Phenotype_29206223.txt', 
              58: 'cellcollective/Septation Initiation Network_26244885.txt', 
              59: 'cellcollective/Signal Transduction in Fibroblasts_18250321.txt', 
              60: 'cellcollective/SKBR3 Breast Cell Line Short-term ErbB Network_24970389.txt', 
              61: 'cellcollective/Signaling in Macrophage Activation_18433497.txt', 
              62: 'cellcollective/SKBR3 Breast Cell Line Long-term ErbB Network_24970389.txt', 
              63: 'cellcollective/T cell differentiation_16542429.txt', 
              64: 'cellcollective/T Cell Receptor Signaling_17722974.txt', 
              65: 'cellcollective/Stomatal Opening Model_27542373.txt', 
              66: 'cellcollective/T-LGL Survival Network 2011 Reduced Network_22102804.txt', 
              67: 'cellcollective/T-Cell Signaling 2006_16464248.txt', 
              68: 'cellcollective/T-LGL Survival Network 2008_18852469.txt', 
              69: 'cellcollective/Toll Pathway of Drosophila Signaling Pathway_23868318.txt', 
              70: 'cellcollective/TOL Regulatory Network_23171249.txt', 
              71: 'cellcollective/T-LGL Survival Network 2011_22102804.txt', 
              72: 'cellcollective/VEGF Pathway of Drosophila Signaling Pathway_23868318.txt', 
              73: 'cellcollective/Trichostrongylus retortaeformis_22253585.txt', 
              74: 'cellcollective/Tumour Cell Invasion and Migration_26528548.txt', 
              75: 'cellcollective/Wg Pathway of Drosophila Signalling Pathways_23868318.txt', 
              76: 'cellcollective/Yeast Apoptosis_23233838.txt'}
if file_dict == file_dicts:
    print("两个字典相同")
else:
    print("两个字典不相同")
#print(file_dict)