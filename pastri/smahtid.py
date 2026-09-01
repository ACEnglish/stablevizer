import re
import sys
import argparse
import dataclasses
import pandas as pd

pattern = r'SMHT(?:[^-]+-){5}[^-_.]+'
SMHTIDRE = re.compile(pattern)
############
# Metadata #
############

PROTOCOLS = {
    "3A": {"tissue": "Whole blood", "preservation": "Snap Frozen", "layer": "Clinically_accessible",
           "tissue_abv": "BLOO", "DedupKey": "D01", "color": "#FF382E"},
    "3B": {"tissue": "Buccal Swab", "preservation": "Fresh", "layer": "Clinically_accessible",
           "tissue_abv": "BUCC", "DedupKey": "D02", "color": "#75001F"},
    "3C": {"tissue": "Esophagus", "preservation": "Snap Frozen", "layer": "Endoderm",
           "tissue_abv": "ESOP", "DedupKey": "D03", "color": "#FFB161"},
    "3D": {"tissue": "Esophagus", "preservation": "Fixed", "layer": "Endoderm",
           "tissue_abv": "ESOP", "DedupKey": "D03", "color": "#FFB161"},
    "3E": {"tissue": "Colon-Ascending", "preservation": "Snap Frozen", "layer": "Endoderm",
           "tissue_abv": "COAS", "DedupKey": "D04", "color": "#C0B5E7"},
    "3F": {"tissue": "Colon-Ascending", "preservation": "Fixed", "layer": "Endoderm",
           "tissue_abv": "COAS", "DedupKey": "D04", "color": "#C0B5E7"},
    "3G": {"tissue": "Colon-Descending", "preservation": "Snap Frozen", "layer": "Endoderm",
           "tissue_abv": "COAD", "DedupKey": "D04", "color": "#9D4BBF"},
    "3H": {"tissue": "Colon-Descending", "preservation": "Fixed", "layer": "Endoderm",
           "tissue_abv": "COAD", "DedupKey": "D04", "color": "#9D4BBF"},
    "3I": {"tissue": "Liver", "preservation": "Snap Frozen", "layer": "Endoderm",
           "tissue_abv": "LIVR", "DedupKey": "D05", "color": "#FF4AA9"},
    "3J": {"tissue": "Liver", "preservation": "Fixed", "layer": "Endoderm",
           "tissue_abv": "LIVR", "DedupKey": "D05", "color": "#FF4AA9"},
    "3K": {"tissue": "Adrenal Gland, L", "preservation": "Snap Frozen", "layer": "Mesoderm",
           "tissue_abv": "ADGL", "DedupKey": "D06", "color": "#EA7A00"},
    "3L": {"tissue": "Adrenal Gland, L", "preservation": "Fixed", "layer": "Mesoderm",
           "tissue_abv": "ADGL", "DedupKey": "D06", "color": "#EA7A00"},
    "3M": {"tissue": "Adrenal Gland, R", "preservation": "Snap Frozen", "layer": "Mesoderm",
           "tissue_abv": "ADGR", "DedupKey": "D06", "color": "#EA7A00"},
    "3N": {"tissue": "Adrenal Gland, R", "preservation": "Fixed", "layer": "Mesoderm",
           "tissue_abv": "ADGR", "DedupKey": "D06", "color": "#EA7A00"},
    "3O": {"tissue": "Aorta", "preservation": "Snap Frozen", "layer": "Mesoderm",
           "tissue_abv": "AORT", "DedupKey": "D07", "color": "#F99B9B"},
    "3P": {"tissue": "Aorta", "preservation": "Fixed", "layer": "Mesoderm",
           "tissue_abv": "AORT", "DedupKey": "D07", "color": "#F99B9B"},
    "3Q": {"tissue": "Lung", "preservation": "Snap Frozen", "layer": "Endoderm",
           "tissue_abv": "LUNG", "DedupKey": "D08", "color": "#1DDD30"},
    "3R": {"tissue": "Lung", "preservation": "Fixed", "layer": "Endoderm",
           "tissue_abv": "LUNG", "DedupKey": "D08", "color": "#1DDD30"},
    "3S": {"tissue": "Heart, LV", "preservation": "Snap Frozen", "layer": "Mesoderm",
           "tissue_abv": "HART", "DedupKey": "D09", "color": "#DDD824"},
    "3T": {"tissue": "Heart, LV", "preservation": "Fixed", "layer": "Mesoderm",
           "tissue_abv": "HART", "DedupKey": "D09", "color": "#DDD824"},
    "3U": {"tissue": "Testis, L", "preservation": "Snap Frozen", "layer": "Germinal",
           "tissue_abv": "TESL", "DedupKey": "D10", "color": "#AAF8B4"},
    "3V": {"tissue": "Testis, L", "preservation": "Fixed", "layer": "Germinal",
           "tissue_abv": "TESL", "DedupKey": "D10", "color": "#AAF8B4"},
    "3W": {"tissue": "Testis, R", "preservation": "Snap Frozen", "layer": "Germinal",
           "tissue_abv": "TESR", "DedupKey": "D10", "color": "#AAF8B4"},
    "3X": {"tissue": "Testis, R", "preservation": "Fixed", "layer": "Germinal",
           "tissue_abv": "TESR", "DedupKey": "D10", "color": "#AAF8B4"},
    "3Y": {"tissue": "Ovary, L", "preservation": "Snap Frozen", "layer": "Germinal",
           "tissue_abv": "OVAL", "DedupKey": "D11", "color": "#FFAFD4"},
    "3Z": {"tissue": "Ovary, L", "preservation": "Fixed", "layer": "Germinal",
           "tissue_abv": "OVAL", "DedupKey": "D11", "color": "#FFAFD4"},
    "3AA": {"tissue": "Ovary, R", "preservation": "Snap Frozen", "layer": "Germinal",
            "tissue_abv": "OVAR", "DedupKey": "D11", "color": "#FFAFD4"},
    "3AB": {"tissue": "Ovary, R", "preservation": "Fixed", "layer": "Germinal",
            "tissue_abv": "OVAR", "DedupKey": "D11", "color": "#FFAFD4"},
    "3AC": {"tissue": "Dermal Fibroblast", "preservation": "Cultured Cells",  "layer": "Cultured_Cells",
            "tissue_abv": "FBRO", "DedupKey": "D12", "Notes": "Isolated from fresh calf skin", "color": "#FFE98D"},
    "3AD": {"tissue": "Skin-Calf", "preservation": "Snap Frozen", "layer": "Ectoderm",
            "tissue_abv": "SKSE", "DedupKey": "D13", "color": "#0EC1B8"},
    "3AE": {"tissue": "Skin-Calf", "preservation": "Fixed", "layer": "Ectoderm",
            "tissue_abv": "SKSE", "DedupKey": "D13", "color": "#0EC1B8"},
    "3AF": {"tissue": "Skin-Abdomen", "preservation": "Snap Frozen", "layer": "Ectoderm",
            "tissue_abv": "SKNE", "DedupKey": "D13", "color": "#91E5DB"},
    "3AG": {"tissue": "Skin-Abdomen", "preservation": "Fixed", "layer": "Ectoderm",
            "tissue_abv": "SKNE", "DedupKey": "D13", "color": "#91E5DB"},
    "3AH": {"tissue": "Skeletal Muscle", "preservation": "Snap Frozen", "layer": "Mesoderm",
            "tissue_abv": "MUSC", "DedupKey": "D14", "color": "#AA9A41"},
    "3AI": {"tissue": "Skeletal Muscle", "preservation": "Fixed", "layer": "Mesoderm",
            "tissue_abv": "MUSC", "DedupKey": "D14", "color": "#AA9A41"},
    "3AJ": {"tissue": "Brain", "preservation": "Fresh", "layer": "Ectoderm",
            "tissue_abv": "BRAI", "DedupKey": "D15", "color": "#BBE7FF"},
    "3AK": {"tissue": "Brain-Frontal lobe", "preservation": "Snap Frozen", "layer": "Ectoderm",
            "tissue_abv": "BRFL", "DedupKey": "D15", "color": "#BBE7FF"},
    "3AL": {"tissue": "Brain-Temporal lobe", "preservation": "Snap Frozen", "layer": "Ectoderm",
            "tissue_abv": "BRTL", "DedupKey": "D15", "color": "#76AEFF"},
    "3AM": {"tissue": "Brain-Cerebellum", "preservation": "Snap Frozen", "layer": "Ectoderm",
            "tissue_abv": "BRCE", "DedupKey": "D15", "color": "#1655C4"},
    "3AN": {"tissue": "Brain-Hippocampus, L", "preservation": "Snap Frozen", "layer": "Ectoderm",
            "tissue_abv": "BRHL", "DedupKey": "D15", "color": "#002A66"},
    "3AO": {"tissue": "Brain-Hippocampus, R", "preservation": "Snap Frozen", "layer": "Ectoderm",
            "tissue_abv": "BRHR", "DedupKey": "D15", "color": "#002A66"},
}

DEDUP_LOOKUP = {
    'D01': {'tissue': 'Blood', 'tissue_abv': 'BLO'},
    'D02': {'tissue': 'Buccal Swab', 'tissue_abv': 'BUC'},
    'D03': {'tissue': 'Esophagus', 'tissue_abv': 'ESO'},
    'D04': {'tissue': 'Colon', 'tissue_abv': 'COA'},
    'D05': {'tissue': 'Liver', 'tissue_abv': 'LIV'},
    'D06': {'tissue': 'Adrenal Gland', 'tissue_abv': 'ADG'},
    'D07': {'tissue': 'Aorta', 'tissue_abv': 'AOR'},
    'D08': {'tissue': 'Lung', 'tissue_abv': 'LUN'},
    'D09': {'tissue': 'Heart', 'tissue_abv': 'HAR'},
    'D10': {'tissue': 'Testis', 'tissue_abv': 'TES'},
    'D11': {'tissue': 'Ovary', 'tissue_abv': 'OVA'},
    'D12': {'tissue': 'Fibroblast', 'tissue_abv': 'FBR'},
    'D13': {'tissue': 'Skin', 'tissue_abv': 'SKI'},
    'D14': {'tissue': 'Muscle', 'tissue_abv': 'MUS'},
    'D15': {'tissue': 'Brain', 'tissue_abv': 'BRA'},
}

LAYER_SHORT = {
    'Clinically_accessible': 'Blood',
    'Cultured_Cells': 'Fibro',
    'Ectoderm': 'Ecto',
    'Endoderm': 'Endo',
    'Germinal': 'Germ',
    'Mesoderm': 'Meso',
}

PLATCODES = {
    "A": {"full": "Illumina NovaSeq X, Illumina NovaSeq X Plus", "name": "Ill"},
    "B": {"full": "PacBio Revio HiFi", "name": "HiFi"},
    "C": {"full": "Illumina NovaSeq 6000", "name": "Ill"},
    "D": {"full": "ONT PromethION 24", "name": "ONT"},
    "E": {"full": "ONT PromethION 2 Solo", "name": "ONT"},
    "F": {"full": "ONT MinION Mk1B", "name": "ONT"},
    "G": {"full": "Illumina HiSeq X", "name": "Ill"},
    "I": {"full": "BGI DNBSEQ-G400", "name": "BGI"},
    "J": {"full": "Element AVITI", "name": "Element"},
    "K": {"full": "Illumina NextSeq 2000", "name": "Ill"},
    "L": {"full": "PacBio Sequel IIe", "name": "PacBio"},
    "M": {"full": "Ultima Genomics UG 100", "name": "UG"},
    "X": {"full": "Multiplatform", "name": "Multi"},
}

ASSAYCODES = {
    "000": "Null",
    "001": "WGS",
    "002": "PCR WGS",
    "003": "Ultra-Long WGS",
    "004": "Fiber-seq",
    "005": "Hi-C",
    "007": "CODEC",
    "008": "Bot-seq",
    "009": "NanoSeq",
    "010": "scNanoSeq",
    "011": "DLP+",
    "012": "Microbulk MALBAC WGS",
    "013": "Single-cell MALBAC WGS",
    "014": "Microbulk PTA WGS",
    "015": "Single-cell PTA WGS",
    "016": "scDip-C",
    "017": "CompDuplex-seq",
    "018": "scCompDuplex-seq",
    "019": "Strand-seq",
    "020": "scStrand-seq",
    "021": "HiDEF-seq",
    "022": "HAT-seq",
    "023": "Microbulk HAT-seq",
    "024": "scHAT-seq",
    "025": "VISTA-seq",
    "026": "Microbulk VISTA-seq",
    "027": "scVISTA-seq",
    "028": "TEnCATS",
    "029": "L1-ONT",
    "030": "ppmSeq",
    "101": "RNA-seq",
    "102": "Kinnex",
    "103": "snRNA-seq",
    "104": "STORM-Seq",
    "105": "Tranquil-Seq",
    "201": "ATAC-seq",
    "202": "CUT&Tag",
    "203": "varCUT&Tag",
    "204": "sc-varCUT&Tag",
    "X": "Multi",
}


##################
# Helper Objects #
##################

@dataclasses.dataclass
class SMaHTid():
    """
    """
    full_name: str
    donor: str
    protocol: str
    aliquot: str
    core: str
    sex: str
    age: str
    platform: str
    assay_code: str
    center: str

    preservation: str
    tissue: str
    tissue_abv: str
    layer: str
    layer_abv: str
    dedup_tissue: str
    dedup_tissue_abv: str

    platform_full: str
    platform_name: str
    assay_name: str

    proto_core: str
    color: str

    def __init__(self, full_name):
        self.full_name = full_name

        pieces = full_name.split('-')

        assert len(
            pieces) == 6, f"Invalid SMaHTid '{full_name}'. Expected Donor-Protocol-AliquotCore-SexAge-PlatAssay-Center"
        self.donor, self.protocol, ali_core, sex_age, plat_assay, self.center = pieces

        assert self.protocol in PROTOCOLS, f"Invalid Protocol {self.protocol}"

        self.aliquot = ali_core[:3]
        self.core = ali_core[3:]
        self.sex = sex_age[0]
        self.age = sex_age[1:]

        self.platform = plat_assay[0]
        assert self.platform in PLATCODES, f"Invalid Platform {self.platform}"
        self.platform_full = PLATCODES[self.platform]['full']
        self.platform_name = PLATCODES[self.platform]['name']

        self.assay_code = plat_assay[1:]
        assert self.assay_code in ASSAYCODES, f"Invalid Assay Code {self.assay_code}"
        self.assay_name = ASSAYCODES[self.assay_code]

        proto = PROTOCOLS[self.protocol]
        self.tissue = proto['tissue']
        self.preservation = proto['preservation']
        self.layer = proto['layer']
        self.layer_abv = LAYER_SHORT[self.layer]
        self.tissue_abv = proto['tissue_abv']
        dedup = DEDUP_LOOKUP[proto['DedupKey']]
        self.dedup_tissue = dedup['tissue']
        self.dedup_tissue_abv = dedup['tissue_abv']

        self.proto_core = self.protocol
        if self.core:
            self.proto_core += '_' + self.core

        self.color = proto['color']
    
    @staticmethod
    def from_re(name):
        match = SMHTIDRE.findall(name)
        assert len(match) == 1,  f"Unable to parse string {name}"
        return SMaHTid(match[0])

    def to_dict(self):
        """
        Transform to dict
        """
        return dataclasses.asdict(self)

    def __str__(self):
        """
        """
        return self.full_name

    def __repr__(self):
        """
        """
        return f"SMaHTid(full_name={self.full_name!r})"

    def __eq__(self, other):
        """
        """
        if not isinstance(other, SMaHTid):
            return NotImplemented
        return repr(self) == repr(other)

    def __hash__(self):
        """
        """
        return hash(self.full_name)


def quick_smahtid_main(args):
    parser = argparse.ArgumentParser(prog="qid", description="SMaHTid Lister",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    args = parser.parse_args(args)
    for k,v in PROTOCOLS.items():
        print(k, v['tissue_abv'], v['tissue'])

def smahtid_main(args):
    """
    Command line SMaHTid parser
    """
    parser = argparse.ArgumentParser(prog="smhtid", description="SMaHTid Deconvolution",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("name", metavar="NAME", type=str, nargs="?",
                        help="SMaHTid (use `-` to read multiple ids from stdin)")
    parser.add_argument("keys", metavar="KEYS", nargs=argparse.REMAINDER,
                        help="Keys to report")
    parser.add_argument("-l", "--list", action="store_true",
                        help="List available keys and exit")
    parser.add_argument("-n", "--no-keys", action="store_true",
                        help="Only print the values")
    parser.add_argument("-t", "--transpose", action="store_true",
                        help="Output deconvolution in single row instead of columns")
    parser.add_argument("--re", action="store_true",
                        help="Use regex to parse out name")

    args = parser.parse_args(args)

    if not args.list and not args.name:
        parser.error("NAME is required when not using --list")

    if args.list:
        # name = SMaHTid("SMHT001-3AD-001B1-M42-B001-broad")
        for field in dataclasses.fields(SMaHTid):
            print(field.name)
        sys.exit(0)

    if args.name == '-':
        names = [line.strip() for line in sys.stdin if line.strip()]
    else:
        names = [args.name]

    if not args.keys:
        args.keys = [_.name for _ in dataclasses.fields(SMaHTid)]
    rows = []
    method = SMaHTid.from_re if args.re else SMaHTid
    for name in names:
        smaht_name = method(name)
        rows.append({
            field.name: getattr(smaht_name, field.name)
            for field in dataclasses.fields(SMaHTid)
            if field.name in args.keys
        })

    d = pd.DataFrame(rows)
    if args.transpose:
        print(d.T.to_csv(sep='\t', index=not args.no_keys, header=False), end='')
    else:
        
        print(d.to_csv(sep='\t', index = False, header=not args.no_keys), end='')

if __name__ == '__main__':
    smahtid_main(sys.argv[1:])
