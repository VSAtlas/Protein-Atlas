import os
import ast
import shutil
from pathlib import Path

# CONFIGURATION
SOURCE_DIRS = ["src", "analysis"]  # Directories to scan
OUTPUT_DIR = "_skeletons"          # Where to put the skinny files
IGNORE_DIRS = ["__pycache__", "tests", "data", "egg-info"]

class Skeletonizer(ast.NodeTransformer):
    """
    Rewrites Python AST to strip function/class bodies,
    replacing them with '...' (Ellipsis) while keeping docstrings.
    """
    def visit_FunctionDef(self, node):
        return self._prune_body(node)

    def visit_AsyncFunctionDef(self, node):
        return self._prune_body(node)

    def visit_ClassDef(self, node):
        # We still want to traverse inside classes to prune methods
        self.generic_visit(node)
        return node

    def _prune_body(self, node):
        new_body = []
        
        # 1. Keep the Docstring (if it exists)
        if (node.body and isinstance(node.body[0], ast.Expr) and 
            isinstance(node.body[0].value, (ast.Str, ast.Constant))):
            new_body.append(node.body[0])
        
        # 2. Add '...' to make it valid syntax
        new_body.append(ast.Expr(value=ast.Constant(value=Ellipsis)))
        
        node.body = new_body
        return node

def make_skeleton(src_path, dest_path):
    try:
        with open(src_path, "r", encoding="utf-8") as f:
            tree = ast.parse(f.read())
        
        # Transform the tree (strip bodies)
        transformer = Skeletonizer()
        clean_tree = transformer.visit(tree)
        ast.fix_missing_locations(clean_tree)
        
        # Convert back to text
        skeleton_code = ast.unparse(clean_tree)
        
        # Add a header so the Agent knows this is a skeleton
        header = f"# SKELETON OF: {src_path}\n# BODIES STRIPPED FOR CONTEXT\n\n"
        
        dest_path.parent.mkdir(parents=True, exist_ok=True)
        with open(dest_path, "w", encoding="utf-8") as f:
            f.write(header + skeleton_code)
            
    except Exception as e:
        print(f"Skipping {src_path}: {e}")

def main():
    root = Path(".")
    out_root = Path(OUTPUT_DIR)
    
    # Clean previous run
    if out_root.exists():
        shutil.rmtree(out_root)
    out_root.mkdir()

    print(f"?? Generating skeletons in {OUTPUT_DIR}/ ...")

    for src_dir in SOURCE_DIRS:
        for path in root.glob(f"{src_dir}/**/*.py"):
            # Check ignore list
            parts = path.parts
            if any(x in IGNORE_DIRS for x in parts):
                continue
                
            # Calculate destination path
            rel_path = path.relative_to(root)
            dest_path = out_root / rel_path
            
            make_skeleton(path, dest_path)
            print(f"  - {rel_path}")

    print("\n? Done! Add '_skeletons/' to your .gitignore")

if __name__ == "__main__":
    main()