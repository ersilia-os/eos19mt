import numpy as np

def aminocoumarins(chebis):
    if "CHEBI:85128" in chebis:
        return 1
    else:
        return 0
    
def aminoglycosides_aminocyclitols(chebis):
    if ("CHEBI:47779" in chebis or "CHEBI:22507" in chebis 
        or "CHEBI:36974" in chebis or "CHEBI:38012" in chebis 
        or "CHEBI:61689" in chebis or "CHEBI:22479" in chebis):
        return 1
    else:
        return 0
    
def ansamycins_rifamycins_macrolides(chebis):
    if ("CHEBI:22565" in chebis or "CHEBI:26580" in chebis 
        or "CHEBI:25106" in chebis or "CHEBI:25105" in chebis 
        or "CHEBI:145565" in chebis):
        return 1
    else:
        return 0
    
def b_lactams_carbapenems(chebis):
    if "CHEBI:46633" in chebis:
        return 1
    else:
        return 0
    
def b_lactams_cephalosporins_cephems(chebis):
    if ("CHEBI:23066" in chebis or "CHEBI:38311" in chebis or
        "CHEBI:55506" in chebis or "CHEBI:55504" in chebis):
        return 1
    else:
        return 0
    
def b_lactams_penicillins(chebis):
    if ("CHEBI:17334" in chebis or "CHEBI:25865" in chebis or
        "CHEBI:88187" in chebis or "CHEBI:51212" in chebis):
        return 1
    else:
        return 0
    
def b_lactams_monobactams(chebis):
    if "CHEBI:50695" in chebis:
        return 1
    else:
        return 0
    
def b_lactamase_inhibitors(chebis):
    if "CHEBI:35625" in chebis:
        return 1
    else:
        return 0
    
def b_lactams_all(chebis):
    if "CHEBI:35627" in chebis or "CHEBI:27933" in chebis or "CHEBI:88225" in chebis:
        return 1
    else:
        return 0
    
def aminopyrimidines_trimethoprim_der(chebis):
    if "CHEBI:38338" in chebis or "CHEBI:45924" in chebis:
        return 1
    else:
        return 0
    
def glycopeptides(chebis):
    if "CHEBI:24396" in chebis:
        return 1
    else:
        return 0
    
def lipopeptides(chebis):
    if "CHEBI:46895" in chebis:
        return 1
    else:
        return 0
    
def nitrofurans(chebis):
    if "CHEBI:87230" in chebis:
        return 1
    else:
        return 0
    
def nucleosides(chebis):
    if "CHEBI:33838" in chebis:
        return 1
    else:
        return 0
    
def oxazolidinones(chebis):
    if "CHEBI:55374" in chebis:
        return 1
    else:
        return 0
    
def phenols(chebis):
    if "CHEBI:33853" in chebis or "CHEBI:26195" in chebis:
        return 1
    else:
        return 0


def polyketides_type_I(chebis):
    if ("CHEBI:26188" in chebis and
        ("CHEBI:25106" in chebis or "CHEBI:25105" in chebis)):
        return 1
    else:
        return 0
    
def polyketides_type_II(chebis):
    if "CHEBI:26188" in chebis and "CHEBI:26895" in chebis:
        return 1
    else:
        return 0
    
def polypeptides(chebis):
    if "CHEBI:15841" in chebis:
        return 1
    else:
        return 0
    
def fluoroquinolones(chebis):
    if "CHEBI:87211" in chebis:
        return 1
    else:
        return 0
    
def quinolines(chebis):
    if ("CHEBI:26513" in chebis or "CHEBI:26512" in chebis or "CHEBI:38921" in chebis or
        "CHEBI:53665" in chebis or "CHEBI:36709" in chebis or "CHEBI:38774" in chebis or
        "CHEBI:38775" in chebis):
        return 1
    else:
        return 0
    
def quat_ammonium_cpds(chebis):
    if "CHEBI:35267" in chebis or "CHEBI:35273" in chebis:
        return 1
    else:
        return 0
    
def sulfonamides(chebis):
    if "CHEBI:35358" in chebis or "CHEBI:87228" in chebis:
        return 1
    else:
        return 0
    
def arsenic_cpds(chebis):
    if ("CHEBI:33406" in chebis or "CHEBI:22632" in chebis or
        "CHEBI:33407" in chebis or "CHEBI:27563" in chebis):
        return 1
    else:
        return 0
    
def amino_acid_der(chebis):
    if "CHEBI:83821" in chebis:
        return 1
    else:
        return 0
    
def phosphonic_acid_der(chebis):
    if "CHEBI:37588" in chebis or "CHEBI:26069" in chebis:
        return 1
    else:
        return 0
    
def hydrazides(chebis):
    if "CHEBI:35362" in chebis:
        return 1
    else:
        return 0

def mc_fatty_acids(chebis):
    if "CHEBI:59554" in chebis or "CHEBI:59558" in chebis:
        return 1
    else:
        return 0
    
def lc_fatty_acids(chebis):
    if "CHEBI:15904" in chebis or "CHEBI:57560" in chebis:
        return 1
    else:
        return 0
    
def pyrazines(chebis):
    if "CHEBI:38314" in chebis:
        return 1
    else:
        return 0
    
def biguanides(chebis):
    if "CHEBI:53662" in chebis:
        return 1
    else:
        return 0
    
def depsipeptides(chebis):
    if "CHEBI:23643" in chebis or "CHEBI:35213" in chebis:
        return 1
    else:
        return 0
    
def benzenesulfonyls(chebis):
    if "CHEBI:52916" in chebis:
        return 1
    else:
        return 0
    
def peroxides(chebis):
    if "CHEBI:25940" in chebis:
        return 1
    else:
        return 0
    
def pyridinium(chebis):
    if "CHEBI:38188" in chebis or "CHEBI:50334" in chebis:
        return 1
    else:
        return 0
    
def antifungal(chebis):
    if "CHEBI:86478" in chebis:
        return 1
    else:
        return 0
    
def heterocyclic_antibiotic(chebis):
    if "CHEBI:24531" in chebis:
        return 1
    else:
        return 0
    
def quinolones(chebis):
    if "CHEBI:86324" in chebis or "CHEBI:23765" in chebis:
        return 1
    else:
        return 0
    
def penams(chebis):
    if "CHEBI:35992" in chebis:
        return 1
    else:
        return 0
    
def anthracyclines(chebis):
    if "CHEBI:48120" in chebis:
        return 1
    else:
        return 0