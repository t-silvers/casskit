# import omnipath as op
# op.options.cache = (CACHE_DIR / "omnipathdb").as_posix()
 
import pypath

# Set cache directory
print(pypath.share.settings.settings.cachedir)
from ...config import CACHE_DIR # TEMP

from pypath import omnipath


# See https://github.com/saezlab/pypath#integrated-databases

# cor.complex.Annotation
protein_sources_default = {
    'Dgidb',
    'Membranome',
    'Exocarta',
    'Vesiclepedia',
    'Matrisome',
    'Surfaceome',
    'CellSurfaceProteinAtlas',
    'CellSurfaceProteinAtlasCellType',
    'HumanPlasmaMembraneReceptome',
    'Matrixdb',
    'Locate',
    'GOIntercell',
    'CellPhoneDB',
    'Ramilowski2015',
    'Ramilowski2015Location',
    'Kirouac2010',
    'GuideToPharmacology',
    'Adhesome',
    'Integrins',
    'Opm',
    'Topdb',
    'Hgnc',
    'Zhong2015',
    'HumanProteinAtlas',
    'HumanProteinAtlasSubcellular',
    'HumanProteinAtlasSecretome',
    'Comppi',
    'SignorPathways',
    'SignalinkPathways',
    'SignalinkFunctions',
    'KeggPathways',
    'KeggPathwaysPC',
    'NetpathPathways',
    'Cpad',
    'Disgenet',
    'Kinasedotcom',
    'Phosphatome',
    'Tfcensus',
    'Intogen',
    'CancerGeneCensus',
    'Cancersea',
    'Msigdb',
    'Lrdb',
    'Baccin2019',
    'Almen2009',
    'Phobius',
    'Icellnet',
    'Cellcellinteractions',
    'Italk',
    'Embrace',
    'UniprotLocations',
    'UniprotFamilies',
    'UniprotTopologies',
    'UniprotTissues',
    'UniprotKeywords',
    'Tcdb',
    'Mcam',
    'Gpcrdb',
    'Celltalkdb',
    'Cellchatdb',
    'Connectomedb',
    'Talklr',
    'Humancellmap',
    'Cellcall',
    #'Biogps',
    'Cellinker',
    'Scconnect',
    'Cancerdrugsdb',
    'Progeny',
    'Celltypist',
    'Cytosig',
    'Wang',
}

#TODO this should be part of json files
complex_sources_default = {
    'CellPhoneDBComplex',
    'CorumFuncat',
    'CorumGO',
    'IcellnetComplex',
    'CellchatdbComplex',
    'CellinkerComplex',
    'ScconnectComplex',
}

#TODO this should be part of json files
default_fields = {
    'Matrisome': ('mainclass', 'subclass'),
    'Locate': ('location',),
    'Vesiclepedia': ('vesicle',),
    'Exocarta': ('vesicle',),
    'Ramilowski_location': ('location',),
    'HPA': ('tissue', 'level'),
    'CellPhoneDB': (
        'receptor',
        'adhesion',
        'cytoplasm',
        'peripheral',
        'secretion',
        'secreted',
        'transporter',
        'transmembrane',
        'extracellular',
    ),
    'CellPhoneDB_Complex': (
        'receptor',
        'adhesion',
        'cytoplasm',
        'peripheral',
        'secretion',
        'secreted',
        'transporter',
        'transmembrane',
        'extracellular',
    ),
    'Cpad': (
        'cancer',
        'effect_on_cancer',
    ),
    'Disgenet': (
        'disease',
    ),
}

# cor.complex.Complexes
complex_resources = (
    'Signor',
    'Corum',
    'CellPhoneDB',
    'Havugimana',
    'Compleat',
    'ComplexPortal',
    'Pdb',
    'GuideToPharmacology',
    'Humap',
    'Humap2',
    'Icellnet',
    'Kegg',
    'Cellchatdb',
    'Cellinker',
    'Spike',
)

pypath.core.complex.Complexes.get_db('Corum')