"""S07.6 frozen functions and independent small list oracles; never production.

The gzip is the complete pre-change source. Load its functions with the current
unchanged classes so typed encoders exercise the real class registry; a private
class the live module has since changed or removed comes from the gzip. Every
function-to-function call stays in the frozen namespace.
"""
import ast
from contextlib import ExitStack, contextmanager
from dataclasses import replace
import gzip
import math
from pathlib import Path
import sys
from types import FunctionType, ModuleType
from unittest.mock import patch

from gbdraw.analysis import collinearity as cc, protein_colinearity as pc
from gbdraw.core.numeric import compensated_sum


def _moved_runtime_global(node):
    """Keep the frozen LOSAT runtime copies runnable.

    The runtime left protein_colinearity for gbdraw.comparisons.losat_runtime
    after S07.6, so the live module no longer provides its stdlib imports,
    constants, or runtime class.
    """
    if isinstance(node, ast.Import):
        return all(alias.name.split('.')[0] in sys.stdlib_module_names for alias in node.names)
    if isinstance(node, ast.ImportFrom):
        return node.module != '__future__' and (node.module or '').split('.')[0] in sys.stdlib_module_names
    if isinstance(node, ast.ClassDef):
        return node.name == 'ProteinBlastpRuntime'
    return isinstance(node, ast.Assign) and any(
        getattr(target, 'id', None) in {'_DEFAULT_LOSATP_BIN', '_BUNDLED_LOSATP_DIR'}
        for target in node.targets)


class _FitSumsToCompensated(ast.NodeTransformer):
    """Add the frozen fit model's floats with compensated_sum, as the runtime does since OV-1237.

    Built-in sum() compensates floats only from Python 3.12, so the frozen
    copy would otherwise differ from the runtime in the last bits on 3.10/3.11.
    """

    def visit_Call(self, node):
        self.generic_visit(node)
        if isinstance(node.func, ast.Name) and node.func.id == 'sum':
            node.func = ast.copy_location(ast.Name('compensated_sum', ast.Load()), node.func)
        return node


def _changed_private_class(node, live_classes):
    return (isinstance(node, ast.ClassDef) and node.name.startswith('_')
            and (node.name not in live_classes or ast.dump(node) != ast.dump(live_classes[node.name])))


def frozen_module():
    source = gzip.decompress(Path(__file__).with_name('protein_colinearity_s076.py.gz').read_bytes())
    module = ModuleType('protein_colinearity_s076')
    module.__dict__.update(vars(pc))
    tree = ast.parse(source)
    live_classes = {node.name: node for node in ast.parse(Path(pc.__file__).read_bytes()).body
                    if isinstance(node, ast.ClassDef)}
    tree.body = [_FitSumsToCompensated().visit(node) if getattr(node, 'name', None) == '_fit_expected_bitscore_model'
                 else node
                 for node in tree.body
                 if isinstance(node, ast.FunctionDef) or _moved_runtime_global(node)
                 or _changed_private_class(node, live_classes)]
    ast.fix_missing_locations(tree)
    exec(compile(tree, 'protein_colinearity_s076.py', 'exec'), module.__dict__)
    return module


frozen = frozen_module()


# Since OV-2370 the live hit selectors take the protein map that orders a
# score tie by CDS; the frozen copies order it by ID. Live callers pass the map,
# which a frozen copy drops; a differential gives IDs in CDS order.
_POSITIONAL_PROTEIN_MAP = {'select_best_hits_per_query': 1, 'select_reciprocal_best_hit_edges': 2}


def without_protein_map(name, function):
    if not isinstance(function, FunctionType):
        return function
    arity = _POSITIONAL_PROTEIN_MAP.get(name)

    def call(*args, protein_map=None, **kwargs):
        return function(*(args if arity is None else args[:arity]), **kwargs)
    return call


def ids_in_cds_order(pm, tables):
    """Rename the proteins so that ID order is CDS order.

    The frozen copies break a score tie on the protein ID, the live code on the
    CDS (OV-2370); with both orders equal, a differential compares the rest.
    """
    order = sorted(pm, key=lambda pid: pc._protein_sort_key(pm[pid]))
    names = {pid: f'c{rank:04d}_{pid}' for rank, pid in enumerate(order)}
    renamed = {names[pid]: replace(pm[pid], protein_id=names[pid]) for pid in order}
    return renamed, {pair: table.assign(query=table['query'].map(names), subject=table['subject'].map(names))
                     for pair, table in tables.items()}


@contextmanager
def frozen_collinear_callers():
    with ExitStack() as stack:
        for name, value in vars(cc).items():
            if name in vars(frozen) and value is vars(pc).get(name) and callable(value):
                stack.enter_context(patch.object(cc, name, without_protein_map(name, vars(frozen)[name])))
        yield


def id_key(frame, column, position):
    series = frame[column]
    if hasattr(series, 'cat'):
        return int(series.cat.codes.iloc[position])
    return series.iloc[position]


def fit_oracle(frame, fraction):
    """List ordinals, scalar logs and ordinary least squares in selected order."""
    order = sorted(range(len(frame)), key=lambda i: (
        frame['length_product'].iloc[i], id_key(frame, 'query', i), id_key(frame, 'subject', i)))
    width = max(4, math.ceil(len(order) / 8))
    selected = []
    for start in range(0, len(order), width):
        group = order[start:start + width]
        group.sort(key=lambda i: (-frame['bitscore'].iloc[i],
                                 id_key(frame, 'query', i), id_key(frame, 'subject', i)))
        selected.extend(group[:max(1, math.ceil(len(group) * fraction))])
    points = [(math.log10(float(frame['length_product'].iloc[i])),
               math.log10(float(frame['bitscore'].iloc[i]))) for i in selected
              if float(frame['length_product'].iloc[i]) > 0 and float(frame['bitscore'].iloc[i]) > 0]
    model = None
    if len(frame) >= 12 and len(points) >= 3:
        xs, ys = zip(*points)
        if len({round(x, 9) for x in xs}) >= 2:
            mx, my = compensated_sum(xs) / len(xs), compensated_sum(ys) / len(ys)
            denominator = compensated_sum((x - mx) ** 2 for x in xs)
            if denominator > 1e-12:
                slope = compensated_sum((x - mx) * (y - my) for x, y in points) / denominator
                intercept = my - slope * mx
                if math.isfinite(slope) and math.isfinite(intercept):
                    model = slope, intercept
    return selected, points, model


def hit_ordinals(frame, limit=1, *, subject_first=False):
    """Exhaustively rank original ordinals; count distinct pairs in each bucket."""
    first, last = ('subject', 'query') if subject_first else ('query', 'subject')
    order = sorted(range(len(frame)), key=lambda i: (
        id_key(frame, first, i), -float(frame['bitscore'].iloc[i]),
        float(frame['evalue'].iloc[i]), -float(frame['identity'].iloc[i]),
        -float(frame['alignment_length'].iloc[i]), id_key(frame, last, i)))
    seen = {}
    selected = []
    for i in order:
        q, s = frame[first].iloc[i], frame[last].iloc[i]
        targets = seen.setdefault(q, set())
        if s not in targets and len(targets) < limit:
            targets.add(s)
            selected.append(i)
    return selected
