import subprocess
import activesite

def detect_active_site(cleaned_pdb):
    center, box_size = activesite.main(cleaned_pdb)
    if center is None:
        print("❌ Active site detection failed.")
    return center, box_size
def convert_to_pdbqt(cleaned_pdb, receptor_pdbqt, mgltools_python, prepare_script):
    try:
        subprocess.run([
            mgltools_python,
            prepare_script,
            "-r", cleaned_pdb,
            "-o", receptor_pdbqt
        ], check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to convert protein: {e}")
        return False