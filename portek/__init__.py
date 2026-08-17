from .portek_utils import encode_seq
from .portek_utils import decode_kmer
from .portek_utils import assign_kmer_group_ava, assign_kmer_group_ovr
from .portek_utils import check_exclusivity

# Strand-aware helpers (new in this fork)
from .portek_utils import reverse_complement_kmer_id
from .portek_utils import canonical_kmer_id
from .portek_utils import reverse_complement_fasta

from .portek_findk import KmerFinder
from .portek_findk import FindOptimalKPipeline
from .portek_enriched import EnrichedKmersPipeline
from .portek_map import MappingPipeline
from .portek_tree import KmerPhyloTreeConstructor
