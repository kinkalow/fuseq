from fuseq.option import Option
from fuseq.genomon import Genomon
from fuseq.preprocess import Preprocess


def main():
    opts = Option().create()
    genomon = Genomon(opts)
    for mf_dir, mf_path in genomon.mf_dic.items():
        p = Preprocess(mf_dir, mf_path, genomon.jun_dic, opts)
        data = p.get_data()
        p.process(data)


if __name__ == "__main__":
    main()
