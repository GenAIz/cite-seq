docker-all: load_data/docker cell_type_annotation/docker differential_gene_expression_analysis/docker functional_state_profile/docker

load_data/docker:
	cd load_data && make docker

cell_type_annotation/docker:
	cd cell_type_annotation && make docker

differential_gene_expression_analysis/docker:
	cd differential_gene_expression_analysis && make docker

functional_state_profile/docker:
	cd functional_state_profile && make docker

genaiz-all: add_metadata/genaiz load_data/genaiz cell_type_annotation/genaiz differential_gene_expression_analysis/genaiz functional_state_profile/genaiz

add_metadata/genaiz:
	cd add_metadata && make genaiz

load_data/genaiz:
	cd load_data && make genaiz

cell_type_annotation/genaiz:
	cd cell_type_annotation && make genaiz

differential_gene_expression_analysis/genaiz:
	cd differential_gene_expression_analysis && make genaiz

functional_state_profile/genaiz:
	cd functional_state_profile && make genaiz

genaiz-publish:
	genaiz sn publish

drift:
	@python3 -c "$$DRIFT_PY"

define DRIFT_PY
import re, sys

PAIRS = [
    ('load_data/load_data.R', 'alternate_workflow/script/load_data_updated.R',
     ['.py_mapping_keys', '.read_h5ad_metadata', '.h5ad_x_dims',
      '.guard_square_h5ad', 'load_h5mu_to_seurat', '.extract_h5ad_adt']),
    ('cell_type_annotation/cell_type_annotation.R',
     'alternate_workflow/script/cell_type_annotation.R',
     ['cell_type_annotation', '.load_reference_dataset',
      '.run_singler_annotation', '.run_singler_fine_percell',
      '.vocab_rescue_signatures', '.rescue_missing_vocab',
      '.prediction_harmonization_map', '.harmonize_prediction_labels',
      '.drop_working_columns', '.is_tcell_label',
      '.compute_tcell_subset_embedding', '.ensure_python',
      '.run_stcat_annotation', '.print_annotation_summary']),
    ('add_metadata/add_metadata.R',
     'alternate_workflow/script/add_metadata.R', None),
    ('functional_state_profile/classify_functional_state.R',
     'alternate_workflow/script/process_seurat_pipeline.R',
     ['.adt_match_features']),
    ('load_data/error_codes.R',
     'functional_state_profile/error_codes.R', None),
    ('add_metadata/error_codes.R',
     'load_data/error_codes.R', None),
    ('cell_type_annotation/error_codes.R',
     'add_metadata/error_codes.R', None),
    ('differential_gene_expression_analysis/error_codes.R',
     'cell_type_annotation/error_codes.R', None),
    ('load_data/portability.R',
     'functional_state_profile/portability.R', None),
    ('add_metadata/portability.R',
     'load_data/portability.R', None),
    ('cell_type_annotation/portability.R',
     'add_metadata/portability.R', None),
    ('differential_gene_expression_analysis/portability.R',
     'cell_type_annotation/portability.R', None),
    ('load_data/disk_backed.R',
     'functional_state_profile/disk_backed.R', None),
    ('add_metadata/disk_backed.R',
     'load_data/disk_backed.R', None),
    ('cell_type_annotation/disk_backed.R',
     'add_metadata/disk_backed.R', None),
    ('differential_gene_expression_analysis/disk_backed.R',
     'cell_type_annotation/disk_backed.R', None),
]

def code_only(lines):
    return [l for l in lines if l and not l.startswith('#')]

def extract(path, name):
    """Function body from its definition to the brace that closes it."""
    src = open(path, encoding='utf-8', errors='replace').read().split('\n')
    pat = re.compile(r'^\s*' + re.escape(name) + r'\s*<-\s*function')
    for i, l in enumerate(src):
        if pat.match(l):
            depth, out = 0, []
            for l2 in src[i:]:
                out.append(l2.strip())
                depth += l2.count('{') - l2.count('}')
                if depth == 0 and len(out) > 1:
                    return out
            return out
    return None

bad = 0
for a, b, names in PAIRS:
    if names is None:
        ta, tb = open(a).read(), open(b).read()
        if ta != tb:
            print('DRIFT  %s' % a); print('       is no longer byte-identical to %s' % b); bad += 1
        else:
            print('ok     %s (byte-identical)' % a.split('/')[0])
        continue
    for n in names:
        fa, fb = extract(a, n), extract(b, n)
        if fa is None and fb is None:
            print('MISSING %s: absent from BOTH copies' % n); bad += 1
        elif fa is None or fb is None:
            print('DRIFT  %s: present in only one copy' % n); bad += 1
        elif code_only(fa) != code_only(fb):
            print('DRIFT  %s: CODE differs between' % n)
            print('       %s' % a); print('       %s' % b); bad += 1
        elif fa != fb:
            print('note   %s: comments differ, code identical' % n)
        else:
            print('ok     %s' % n)

print()
if bad:
    print('%d divergence(s). Port the change to both copies before committing.' % bad)
    sys.exit(1)
print('no drift')
endef
export DRIFT_PY
