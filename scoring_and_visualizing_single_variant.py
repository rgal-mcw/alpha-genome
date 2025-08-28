'''
Ryan Gallagher

Scoring and Visualizing a Single Variant

Adapted from https://www.alphagenomedocs.com/colabs/variant_scoring_ui.html

Here, we go through the process of predicting the effect of a single variant on different modalities, such as gene expression and chromatin accessibility.

For my application, I'm going to substitute the tutorial variant for a real one that we see in MCW SVI 0162 UDD.

Utilizing Geneyx - we find that there is an intronic deletion on the DEAF1 gene - which seems to match our patient's pheno well.

'''

from secret.api_key import ALPHA_GENOME_API_KEY

from alphagenome.data import gene_annotation, genome, transcript
from alphagenome.models import dna_client, variant_scorers
from alphagenome.visualization import plot_components
import pandas as pd

# Load the model
dna_model = dna_client.create(ALPHA_GENOME_API_KEY())

HG38_GTF_FEATHER = (
    'https://storage.googleapis.com/alphagenome/reference/gencode/'
    'hg38/gencode.v46.annotation.gtf.gz.feather'
)

# Initialize an empty directory to serve as a variant effect prediction cache.
_prediction_cache = {}
_transcript_extractor_cache = {}

print('\n ---- Score Variant ---- \n')

organism = 'human'
organism_map = {'human': dna_client.Organism.HOMO_SAPIENS}
organism = organism_map[organism]

# A real deletion on the DEAF1 gene of SVI 0162 UDD
variant_chromosome = 'chr11'
variant_position = 647_807
variant_reference_bases = 'GTGCTGGGAGCAGGTGGGTGGACTTGACCCAGGCGGTGCTGGGAGCAGGTGGGTGGACTTGACCCAGGCGGTGCTGGGAGCAGGTGGGTGGACTTGACCCAGGCG'
variant_alternate_bases = ''

variant = genome.Variant(
        chromosome = variant_chromosome,
        position = variant_position,
        reference_bases = variant_reference_bases,
        alternate_bases = variant_alternate_bases,
        )

# Specify length of sequence around variant to predict (this is taking HG38 surrounding the event and predicting when there is vs. isnt a deletion there)

sequence_length = '1MB' # options = '2KB', '16KB', '100KB', '500KB', '1MB' as string
sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
        f'SEQUENCE_LENGTH_{sequence_length}'
        ]

# The input interval is derived from the variant (centered on it).
interval = variant.reference_interval.resize(sequence_length)

# Additonal settings

variant_scores = dna_model.score_variant(
        interval=interval,
        variant=variant,
        variant_scorers=list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values()),
        )

df_scores = variant_scorers.tidy_scores(variant_scores)

download_predictions = True
if download_predictions:
    df_scores.to_csv(f'chr11_DEAF1_DEL_scores.csv', index=False)

columns = [
        c for c in df_scores.columns if c not in ['variant_id', 'scored_intervals']
        ]
#print(df_scores[columns])


print('\n ---- Visualize Variant Effects ---- \n')

# Specify list of cell and tissue ontologies

# Since I don't know a single thing about what I would select here,
# I asked Gemini what cell types I should use given the sample phenotype / gene I'm looking at.
# It suggested the Cerebellum, Brain, and a bunch that weren't supported... I found a third in Caudate Nucleus

# It output Ontology CURIE, and I double checked them from the documentation: https://www.alphagenomedocs.com/colabs/tissue_ontology_mapping.html

ontology_terms = ['UBERON:0002037', 'UBERON:0000955', 'UBERON:0001873'] # Cerebellum, Brain, Caudate Nucleus

plot_gene_annotation = True
plot_longest_transcript_only = True

# Output types to plot:
# There are 13 different output options - each have their own scoring and interpretation
# https://www.alphagenomedocs.com/exploring_model_metadata.html

# Once again, I'm not exactly sure what I'm looking at, so I asked Gemini which would be the most useful.
# It said - CAGE (Cap Analysis of Gene Expression): This would give us an idea of the deletion increased or decreased expression of the DEAF1 gene.

# And ATAC (Assay for Transposase-Accessible Chromatin): This would tell us how "Accessible" the DNA is. This could help us interpret the change we see in the CAGE track. EXAMPLE: If the deletion removes an enhancer element, we would expect to see the ATAC peak disappear at that position. --- UPDATE: THESE DONT HAVE DATA FOR ATAC

# We can use ChIP-Histone instead: This tracks measure chemical marks on the DNA that reveal the function of the region. Could reveal a mechanistic link.o --- UPDATE: THIS DIDNT WORK EITHER. 

# I'll try with only CAGE

plot_rna_seq = True # Could be useful if CAGE shows something
plot_cage = True
plot_atac = False
plot_dnase = False
plot_chip_histone = False
plot_chip_tf = False
plot_splice_sites = False
plot_splice_site_usage = False
plot_contact_maps = False
plot_splice_junctions = False


# Options to filter tracks to only a specific DNA strand
## DEAF1 IS ON THE NEGATIVE STRAND
filter_to_positive_strand = False
filter_to_negative_strand = True

if filter_to_positive_strand and filter_to_negative_strand:
    raise ValueError(
            'Cannot specify both filter_to_positive_strand and '
            'filter_to_negative strand.'
            )

# Other visualization options:
ref_color = 'blue'
alt_color = 'red'
ref_alt_colors = {'REF': ref_color, 'ALT': alt_color}
plot_interval_width = 43008
plot_interval_shift = 0

# Load gene annotation
if organism in _transcript_extractor_cache:
    transcript_extractor, longest_transcript_extractor = (
            _transcript_extractor_cache[organism]
            )
else:
    match organism:
        case dna_client.Organism.HOMO_SAPIENS:
            gtf_path = HG38_GTF_FEATHER
        case _:
            raise ValueError(f'Unsupported organism: {organism}')

    gtf = pd.read_feather(gtf_path)

    # Filter to protein-coding genes and highly supported transcripts
    gtf_transcript = gene_annotation.filter_transcript_support_level(
            gene_annotation.filter_protein_coding(gtf), ['1']
            )

    # Extractor for identifying transcripts in the region
    transcript_extractor = transcript.TranscriptExtractor(gtf_transcript)
    
    # Also define an extractor that fetches only the longest transcript gene.
    gtf_longest_transcript = gene_annotation.filter_to_longest_transcript(
            gtf_transcript
            )

    longest_transcript_extractor = transcript.TranscriptExtractor(
            gtf_longest_transcript
            )
    _transcript_extractor_cache[organism] = (
            transcript_extractor,
            longest_transcript_extractor
            )

def _predict_variant_cached(
    interval, variant, organism, requested_outputs, ontology_terms
):
  """Cache wrapper of dna_model.predict_variant."""
  # Create a unique key from the function arguments.
  cache_key = (
      str(interval),
      str(variant),
      str(organism),
      tuple(requested_outputs),
      tuple(ontology_terms),
  )

  # Check if the result is already in the cache.
  if cache_key in _prediction_cache:
    return _prediction_cache[cache_key]

  # If not, compute the prediction and store it in the cache.
  result = dna_model.predict_variant(
      interval=interval,
      variant=variant,
      organism=organism,
      requested_outputs=requested_outputs,
      ontology_terms=ontology_terms,
  )
  _prediction_cache[cache_key] = result
  return result

output = _predict_variant_cached(
        interval = interval,
        variant = variant,
        organism = organism,
        requested_outputs = [*dna_client.OutputType],
        ontology_terms = ontology_terms,
        )

# Filter to DNA strand if requested
ref, alt = output.reference, output.alternate

if filter_to_positive_strand:
    ref = ref.filter_to_strand(strand = "+")
    alt = alt.filter_to_strand(strand = "+")
elif filter_to_negative_strand:
    ref = ref.filter_to_strand(strand='-')
    alt = alt.filter_to_strand(strand='-')


# Build plot.
components = []

# Gene and transcript annotation
if plot_gene_annotation:
    if plot_longest_transcript_only:
        transcripts = longest_transcript_extractor.extract(interval)
    else:
        transcripts = transcript_extractor.extract(interval)
    components.append(plot_components.TranscriptAnnotation(transcripts))


# Individual output type plots.
plot_map = {
        'plot_atac': (ref.atac, alt.atac, 'ATAC'),
        'plot_cage': (ref.cage, alt.cage, 'CAGE'),
        'plot_chip_histone': (ref.chip_histone, alt.chip_histone, 'CHIP_HISTONE'),
        'plot_chip_tf': (ref.chip_tf, alt.chip_tf, 'CHIP_TF'),
        'plot_contact_maps': (ref.contact_maps, alt.contact_maps, 'CONTACT_MAPS'),
        'plot_dnase': (ref.dnase, alt.dnase, 'DNASE'),
        'plot_rna_seq': (ref.rna_seq, alt.rna_seq, 'RNA_SEQ'),
        'plot_splice_junctions': (
            ref.splice_junctions,
            alt.splice_junctions,
            'SPLICE_JUNCTIONS',
        ),
        'plot_splice_sites': (ref.splice_sites, alt.splice_sites, 'SPLICE_SITES'),
        'plot_splice_site_usage': (
            ref.splice_site_usage,
            alt.splice_site_usage,
            'SPLICE_SITE_USAGE',
    ),
}

for key, (ref_data, alt_data, output_type) in plot_map.items():
  if eval(key) and ref_data is not None and ref_data.values.shape[-1] == 0:
    print(
        f'Requested plot for output {output_type} but no tracks exist in'
        ' output. This is likely because this output does not exist for your'
        ' ontologies or requested DNA strand.'
    )
  if eval(key) and ref_data and alt_data:
    match output_type:
      case 'CHIP_HISTONE':
        ylabel_template = (
            f'{output_type}: {{biosample_name}} ({{strand}})\n{{histone_mark}}'
        )
      case 'CHIP_TF':
        ylabel_template = (
            f'{output_type}: {{biosample_name}}'
            ' ({strand})\n{transcription_factor}'
        )
      case 'CONTACT_MAPS':
        ylabel_template = f'{output_type}: {{biosample_name}} ({{strand}})'
      case 'SPLICE_SITES':
        ylabel_template = f'{output_type}: {{name}} ({{strand}})'
      case _:
        ylabel_template = (
            f'{output_type}: {{biosample_name}} ({{strand}})\n{{name}}'
        )

    if output_type == 'CONTACT_MAPS':
      component = plot_components.ContactMapsDiff(
          tdata=alt_data - ref_data,
          ylabel_template=ylabel_template,
      )
      components.append(component)
    elif output_type == 'SPLICE_JUNCTIONS':
      ref_plot = plot_components.Sashimi(
          ref_data,
          ylabel_template='REF: ' + ylabel_template,
      )
      alt_plot = plot_components.Sashimi(
          alt_data,
          ylabel_template='ALT: ' + ylabel_template,
      )
      components.extend([ref_plot, alt_plot])
    else:
      component = plot_components.OverlaidTracks(
          tdata={'REF': ref_data, 'ALT': alt_data},
          colors=ref_alt_colors,
          ylabel_template=ylabel_template,
          alpha = 0.6 #ADDED FOR TRANSPARENCY
      )
      components.append(component)

if plot_interval_width > interval.width:
  raise ValueError(
      f'plot_interval_width ({plot_interval_width}) must be less than '
      f'interval.width ({interval.width}).'
  )

plot = plot_components.plot(
    components=components,
    interval=interval.shift(plot_interval_shift).resize(plot_interval_width),
    annotations=[
        plot_components.VariantAnnotation([variant]),
    ],
)

plot.savefig('deaf1_variant_plot.png')
