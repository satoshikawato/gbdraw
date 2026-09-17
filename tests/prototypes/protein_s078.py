"""Exact 02ca8f95 owners, isolated functions/private classes; test use only."""
import ast
from contextlib import ExitStack, contextmanager
import gzip
from pathlib import Path
import sys
from types import ModuleType
from unittest.mock import patch

from gbdraw.analysis import collinearity as cc, collinearity_units as cu, protein_colinearity as pc


def load_frozen(current):
    name = current.__name__.rsplit('.', 1)[-1]
    source = gzip.decompress(Path(__file__).with_name(f'{name}_s077.py.gz').read_bytes())
    tree = ast.parse(source)
    live = {n.name: n for n in ast.parse(Path(current.__file__).read_bytes()).body
            if isinstance(n, ast.ClassDef)}
    module = ModuleType(f'{__name__}.{name}_s077')
    sys.modules[module.__name__] = module
    exec(compile(tree, f'{name}_s077.py', 'exec'), module.__dict__)
    # Typed encoders require the public class registry. Reuse only AST-identical
    # public classes; private accumulators and their methods remain frozen.
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and not node.name.startswith('_'):
            assert ast.dump(node) == ast.dump(live[node.name]), node.name
            setattr(module, node.name, getattr(current, node.name))
    return module


frozen = load_frozen(pc)
frozen_units = load_frozen(cu)


@contextmanager
def frozen_collinear_callers():
    with ExitStack() as stack:
        for name, value in vars(cc).items():
            if name in vars(frozen) and value is vars(pc).get(name) and callable(value):
                stack.enter_context(patch.object(cc, name, vars(frozen)[name]))
        stack.enter_context(patch.object(cc, 'build_collinearity_unit_index',
                                        frozen_units.build_collinearity_unit_index))
        yield
