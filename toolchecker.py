import subprocess

def check_tool(tool_name, version_arg="--version"):
    try:
        result = subprocess.run([tool_name, version_arg], capture_output=True, text=True, check=True)
        print(f"{tool_name} found: {result.stdout.splitlines()[0]}")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        print(f"{tool_name} NOT found or not on PATH.")
        return False

def main():
    tools = ["phenix", "p2rank", "vina", "gnina", "obabel"]
    for tool in tools:
        check_tool(tool)

if __name__ == "__main__":
    main()