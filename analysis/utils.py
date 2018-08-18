import subprocess


def pack_the_dataset(outdir="output"):
    subprocess.call("mkdir -p " + outdir, shell=True)
    subprocess.call("mv *.eps " + outdir, shell=True)
    subprocess.call("mv *.pdf " + outdir, shell=True)
    subprocess.call("mv *.csv " + outdir, shell=True)
    subprocess.call("mv *.tex " + outdir, shell=True)
    subprocess.call(
        "zip -r {out}.zip {out}".format(out=outdir), shell=True)