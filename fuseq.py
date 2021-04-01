from fuseq.option import Option
from fuseq.genomon import Genomon
from fuseq.preprocess import Preprocess


def main():
    opts = Option().opts
    genomon = Genomon(opts['genomon_dir'])
    Preprocess(genomon)


if __name__ == "__main__":
    main()
