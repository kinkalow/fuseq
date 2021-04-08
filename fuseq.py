from fuseq.option import Option
from fuseq.genomon import Genomon
from fuseq.blat import Blat
from fuseq.checker import Checker


def main():
    Checker.has_tools()
    opts = Option().create()
    genomon = Genomon(opts)
    for mf_dir, mf_path in genomon.mf_dic.items():
        work_dir = f'{opts.out_base_dir}/{mf_dir}/_work'
        blat = Blat(mf_path, genomon.jun_dic, work_dir, opts)
        blat.run()


if __name__ == "__main__":
    main()
