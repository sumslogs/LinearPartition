from models.vienna import ViennaLinearPartition, kT
from models.contrafold import ContrafoldLinearPartition
from utils.dump_forest import dump_forest
from utils.mea_structure import PairProb_MEA
from utils.threshknot import threshknot_pairs, output_bpseq, dotbracket_to_pairs
from utils.dump_bpp import output_bpp_to_file
from utils.shape_data import load_shape_data

class BeamCKYParser:
    def __init__(self, beam_size=100, nosharpturn=True, is_verbose=False, bppfile="", bppfileindex="",
                 pf_only=False, bpp_cutoff=0.0, forestfile="", mea_=False, gamma=3.0, mea_file_index="",
                 bpseq=False, threshknot_=False, threshknot_threshold=0.3, threshknot_file_index="",
                 shape_file_path="", is_fasta=False, model_to_use='Vienna'):
        self.beam_size = beam_size
        self.no_sharp_turn = nosharpturn
        self.is_verbose = is_verbose
        self.bpp_file = bppfile
        self.bpp_file_index = bppfileindex
        self.pf_only = pf_only
        self.bpp_cutoff = bpp_cutoff
        self.forest_file = forestfile
        self.mea_ = mea_
        self.gamma = gamma
        self.mea_file_index = mea_file_index
        self.bpseq = bpseq
        self.threshknot_ = threshknot_
        self.threshknot_threshold = threshknot_threshold
        self.threshknot_file_index = threshknot_file_index
        self.is_fasta = is_fasta

        self.shape_file_path = shape_file_path
        self.using_model = model_to_use

    def parse(self, seq: str):
        if self.using_model == 'Vienna':
            self.model = ViennaLinearPartition(
                seq=seq,
                beam_size=self.beam_size,
                bpp_cutoff = self.bpp_cutoff,
                no_sharp_turn=self.no_sharp_turn,
                SHAPE_data=load_shape_data(self.shape_file_path) if len(self.shape_file_path) > 0 else []
            )
        elif self.using_model == 'Contrafold':
            self.model = ContrafoldLinearPartition(
                seq=seq,
                beam_size=self.beam_size,
                bpp_cutoff = self.bpp_cutoff,
                no_sharp_turn=self.no_sharp_turn)
        else:
            raise Exception(f'Model {self.using_model} not implemented.')

        self.model.inside()
        self.model.outside()
        self.model.cal_PairProb()

if __name__ == "__main__":
    with open('../smallseq', 'r') as f:
        rna_seq_list = f.read().split("\n")

    # parser = BeamCKYParser(model_to_use='Contrafold')
    # parser = BeamCKYParser(model_to_use='Vienna',
    #                        shape_file_path='../example.shape')
    parser = BeamCKYParser(model_to_use='Vienna')
    for i,seq in enumerate(rna_seq_list):
        parser.parse(seq)

        # Interact with the computed data structures:
        output_bpp_to_file(parser.model, 'scratch/test_bpp_out.txt')

        if parser.using_model == 'Vienna':
            print(f"Free Energy of Ensemble: { -kT * parser.model.log_partition_coeff() / 100.0 } kcal/mol")
        else:
            print(f"Log partition coeff: { parser.model.log_partition_coeff() } kcal/mol")

        # Forest dump
        dump_forest(parser.model, 'scratch/test_dump_forest.txt', False)

        # MEA structure:
        structure = PairProb_MEA(parser.model, parser.gamma)
        print(seq)
        print(structure)

        # Output of that into bpseq format"
        pairs_from_mea = dotbracket_to_pairs(structure)
        output_bpseq('scratch/mea_bpseq.txt', pairs_from_mea, seq)

        # Threshknot, and output as bpseq:
        tk_pairs = threshknot_pairs(parser.model, threshknot_threshold=parser.threshknot_threshold)
        output_bpseq('scratch/threshknot_bpseq.txt', tk_pairs, seq)

        break