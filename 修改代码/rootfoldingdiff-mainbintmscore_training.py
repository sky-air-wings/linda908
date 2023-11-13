把sample移到tmscore-training相同目录下

def run_tmalign(query: str, reference: str, fast: bool = False) -> float:
    """
    Run TMalign on the two given input pdb files
    """
    assert os.path.isfile(query)
    assert os.path.isfile(reference)

    # Check if TMalign is installed
    exec = shutil.which("TMalign")
    if not exec:
        raise FileNotFoundError("TMalign not found in PATH")

    # Build the command
    cmd = f"{exec} {query} {reference}"
    if fast:
        cmd += " -fast"
    try:
        output = subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError:
        logging.warning(f"Tmalign failed on {query}|{reference}, returning NaN")
        return np.nan

    # Parse the outpu
    score_lines = []
    for line in output.decode().split("\n"):
        if line.startswith("TM-score"):
            score_lines.append(line)

    # Fetch the chain number
    key_getter = lambda s: re.findall(r"Chain_[12]{1}", s)[0]
    score_getter = lambda s: float(re.findall(r"=\s+([0-9.]+)", s)[0])
    results_dict = {key_getter(s): score_getter(s) for s in score_lines}
    return results_dict["Chain_2"]  # Normalize by reference length    这些代码改成下面代码

def run_tmalign(query: str, reference: str, fast: bool = False) -> float:
    """
    Run TMalign on the two given input pdb files
    """
    assert os.path.isfile(query)
    assert os.path.isfile(reference)

    # Check if TMalign is installed
    exec = "/root/foldingdiff-main/foldingdiff/tmalign.py"

    # Build the command
    cmd = f"{exec} {query} {reference}"
    if fast:
        cmd += " -fast"
    try:
        output = subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError:
        logging.warning(f"Tmalign failed on {query}|{reference}, returning NaN")
        return np.nan

    # Parse the outpu
    score_lines = []
    for line in output.decode().split("\n"):
        if line.startswith("TM-score"):
            score_lines.append(line)

    # Fetch the chain number
    key_getter = lambda s: re.findall(r"Chain_[12]{1}", s)[0]
    score_getter = lambda s: float(re.findall(r"=\s+([0-9.]+)", s)[0])
    results_dict = {key_getter(s): score_getter(s) for s in score_lines}
    return results_dict["Chain_2"]  # Normalize by reference length   这是改完后代码

