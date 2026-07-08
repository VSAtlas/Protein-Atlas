from __future__ import annotations

import ast
import shutil
from pathlib import Path

SOURCE_DIRS = ["src", "analysis"]
OUTPUT_DIR = "_skeletons"
IGNORE_DIRS = ["__pycache__", "tests", "data", "egg-info"]


class Skeletonizer(ast.NodeTransformer):
    """Strip function and async-function bodies to ellipsis while keeping docstrings."""

    def visit_FunctionDef(self, node):  # type: ignore[override]
        return self._prune_body(node)

    def visit_AsyncFunctionDef(self, node):  # type: ignore[override]
        return self._prune_body(node)

    def visit_ClassDef(self, node):  # type: ignore[override]
        self.generic_visit(node)
        return node

    def _prune_body(
        self, node: ast.FunctionDef | ast.AsyncFunctionDef
    ) -> ast.FunctionDef | ast.AsyncFunctionDef:
        new_body: list[ast.stmt] = []

        if (
            node.body
            and isinstance(node.body[0], ast.Expr)
            and isinstance(node.body[0].value, (ast.Str, ast.Constant))
        ):
            new_body.append(node.body[0])
        new_body.append(ast.Expr(value=ast.Constant(value=Ellipsis)))

        node.body = new_body
        return node


def make_skeleton(src_path: Path, dest_path: Path) -> None:
    try:
        tree = ast.parse(src_path.read_text(encoding="utf-8"))
        transformer = Skeletonizer()
        clean_tree = transformer.visit(tree)
        ast.fix_missing_locations(clean_tree)
        skeleton_code = ast.unparse(clean_tree)

        header = f"# SKELETON OF: {src_path}\n# BODIES STRIPPED FOR CONTEXT\n\n"
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        dest_path.write_text(header + skeleton_code, encoding="utf-8")
    except Exception as exc:
        print(f"Skipping {src_path}: {exc}")


def main() -> int:
    root = Path(".")
    out_root = Path(OUTPUT_DIR)

    if out_root.exists():
        shutil.rmtree(out_root)
    out_root.mkdir()

    print(f"?? Generating skeletons in {OUTPUT_DIR}/ ...")
    for src_dir in SOURCE_DIRS:
        for path in root.glob(f"{src_dir}/**/*.py"):
            if any(token in IGNORE_DIRS for token in path.parts):
                continue
            rel_path = path.relative_to(root)
            dest_path = out_root / rel_path
            make_skeleton(path, dest_path)
            print(f"  - {rel_path}")

    print("\n? Done! Add '_skeletons/' to your .gitignore")
    return 0

