from fuseq.option import Option
from fuseq.genomon import Genomon
from fuseq.preprocess import Preprocess


def main():
    opts = Option().create()
    genomon = Genomon(opts.genomon_dir)
    for mf_dir, mf_path in genomon.mf_dic.items():
        p = Preprocess(mf_dir, mf_path, genomon.jun_dic, opts)
        data = p.create_target_list()
        import time
        start = time.time()
        p.process(data)
        print(f"{time.time() - start:.3f}[s]")


if __name__ == "__main__":
    main()
