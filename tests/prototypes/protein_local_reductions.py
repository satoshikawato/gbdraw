"""S07.6 frozen functions and independent small list oracles; never production.

The gzip is the complete pre-change source. Load its functions with the current
unchanged classes so typed encoders exercise the real class registry. Every
function-to-function call stays in the frozen namespace.
"""
import ast
from contextlib import ExitStack, contextmanager
import gzip
import math
from pathlib import Path
from types import ModuleType
from unittest.mock import patch

from gbdraw.analysis import collinearity as cc, protein_colinearity as pc


def frozen_module():
    source = gzip.decompress(Path(__file__).with_name('protein_colinearity_s076.py.gz').read_bytes())
    module = ModuleType('protein_colinearity_s076')
    module.__dict__.update(vars(pc))
    tree = ast.parse(source)
    tree.body = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
    exec(compile(tree, 'protein_colinearity_s076.py', 'exec'), module.__dict__)
    return module


frozen = frozen_module()


@contextmanager
def frozen_collinear_callers():
    with ExitStack() as stack:
        for name, value in vars(cc).items():
            if name in vars(frozen) and value is vars(pc).get(name) and callable(value):
                stack.enter_context(patch.object(cc, name, vars(frozen)[name]))
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
            mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
            denominator = sum((x - mx) ** 2 for x in xs)
            if denominator > 1e-12:
                slope = sum((x - mx) * (y - my) for x, y in points) / denominator
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
