import os
from dnachisel import EnforceTranslation, CodonOptimize, DnaOptimizationProblem
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

class ReverseTranslate:
    """
    Description:
        Converts an amino acid sequence to codon-optimized DNA and saves as a GenBank file.
    Input:
        aa_sequence (str): Protein sequence.
        host (str): Target organism.
        filename (str): Name for the output GenBank file (defaults to 'synthetic_pks.gb').
    Output:
        dict: Confirmation message and the path to the generated GenBank file.
    """

    def initiate(self) -> None:
        self.data_dir = "modules/pks/data"
        if not os.path.exists(self.data_dir):
            os.makedirs(self.data_dir)
            
        self.host_map = {
            "e_coli": "e_coli",
            "ecoli": "e_coli",
            "s_cerevisiae": "s_cerevisiae",
            "yeast": "s_cerevisiae",
            "p_putida": "160488",         
            "pseudomonas": "160488",
            "s_coelicolor": "100226",     
            "streptomyces": "100226", 
            "s_albus": "318586",            
            "s_venezuelae": "54571"      
        }
        
        # A basic lookup table just to generate a valid starting DNA string
        self.naive_codons = {
            'A': 'GCG', 'C': 'TGC', 'D': 'GAT', 'E': 'GAA',
            'F': 'TTT', 'G': 'GGC', 'H': 'CAT', 'I': 'ATT',
            'K': 'AAA', 'L': 'CTG', 'M': 'ATG', 'N': 'AAC',
            'P': 'CCG', 'Q': 'CAG', 'R': 'CGC', 'S': 'AGC',
            'T': 'ACC', 'V': 'GTG', 'W': 'TGG', 'Y': 'TAT',
            '*': 'TAA'
        }

    def run(self, aa_sequence: str, host: str = "e_coli", filename: str = "synthetic_pks.gb") -> dict:
        if not aa_sequence:
            raise ValueError("Amino acid sequence cannot be empty.")
            
        clean_host = host.lower().replace(" ", "_").replace(".", "")
        mapped_species = self.host_map.get(clean_host)
        
        if not mapped_species:
            mapped_species = "e_coli" 
            actual_host_used = "E. coli (Fallback: Requested host not found)"
        else:
            actual_host_used = clean_host

        # 1. Generate a naive DNA sequence so DnaChisel doesn't crash on AA letters
        try:
            initial_dna = "".join([self.naive_codons[aa] for aa in aa_sequence.upper()])
        except KeyError as e:
            raise ValueError(f"Invalid amino acid character found: {e}")

        # 2. Define Optimization Problem using the naive DNA
        problem = DnaOptimizationProblem(
            sequence=initial_dna,
            constraints=[EnforceTranslation()], # This locks the amino acids in place
            objectives=[CodonOptimize(species=mapped_species)] # This fixes the codons
        )
        
        # Perform optimization
        problem.optimize()
        optimized_dna = problem.sequence
        
        # 3. Create GenBank Record with BioPython
        protein_seq = str(Seq(optimized_dna).translate(to_stop=True))
        record = SeqRecord(
            Seq(optimized_dna),
            id="SYNTH_PKS",
            name="Synthetic_PKS",
            description=f"Codon optimized for {actual_host_used} using DnaChisel",
            annotations={"molecule_type": "DNA"}
        )

        cds_feature = SeqFeature(
            FeatureLocation(0, len(optimized_dna)),
            type="CDS",
            qualifiers={
                "label": "Optimized PKS Sequence",
                "product": "synthetic PKS module",
                "transl_table": "11",
                "codon_start": "1",
                "translation": protein_seq,
            }
        )
        record.features.append(cds_feature)
        
        # 4. Save the file
        file_path = os.path.join(self.data_dir, filename)
        SeqIO.write(record, file_path, "genbank")
        
        return {
            "status": "success",
            "dna_sequence": optimized_dna,
            "file_saved_at": file_path,
            "message": f"Sequence optimized for {actual_host_used} and saved to {filename}."
        }

_instance = ReverseTranslate()
_instance.initiate()
reverse_translate = _instance.run