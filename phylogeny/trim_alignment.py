"""trimAl -automated1, via the binary if present, otherwise via pytrimal.

Called by 02_align_trim_tree.sh; can also be run directly:

    python phylogeny/trim_alignment.py <aligned.fasta> <trimmed.fasta>

WHY THIS EXISTS
    `trimal` is not packaged in Ubuntu 24.04, so `apt-get install trimal`
    fails there.  pytrimal is the same C++ code as a Python module and
    installs with pip on any platform, so the trimming step does not depend
    on the distribution.

    Either route is reported explicitly and recorded in output/phylo_versions.txt.
    What must never happen is the silent fallback the original notebook had -
    copying the untrimmed alignment forward on failure, leaving the Methods
    describing a step that did not run.
"""
import shutil, subprocess, sys

def ncols(path):
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                n += 1
                if n == 2:
                    break
            elif n == 1:
                yield len(line.strip())

def width(path):
    return sum(ncols(path))

def main(src, dst):
    binary = shutil.which("trimal")
    if binary:
        r = subprocess.run([binary, "-in", src, "-out", dst,
                            "-automated1", "-fasta"],
                           capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"trimal failed:\n{r.stderr.strip()}")
        ver = subprocess.run([binary, "--version"], capture_output=True,
                             text=True).stdout.strip().splitlines()
        route = f"trimal binary ({ver[0] if ver else 'version unknown'})"
    else:
        try:
            from pytrimal import Alignment, AutomaticTrimmer
        except ImportError:
            sys.exit("neither the trimal binary nor pytrimal is available.\n"
                     "  pip install pytrimal        (works everywhere)\n"
                     "  apt-get install trimal      (not packaged on Ubuntu 24.04)")
        import pytrimal
        aln = Alignment.load(src)
        AutomaticTrimmer(method="automated1").trim(aln).dump(dst, format="fasta")
        route = f"pytrimal {pytrimal.__version__} (automated1)"

    before, after = width(src), width(dst)
    if after == 0:
        sys.exit("trimming removed every column - check the alignment")
    print(f"  {route}")
    print(f"  alignment columns: {before} -> {after} "
          f"({100 * after / before:.1f}% retained)")
    with open("output/phylo_versions.txt", "a") as f:
        f.write(f"trimming   {route}\n")
        f.write(f"alignment  {before} columns before trimming, {after} after\n")
    return before, after


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
