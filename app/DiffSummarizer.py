
import glob

from YoosCufflinks import CuffDiff
from YoosCufflinks import CuffNorm
from YoosENS85 import EnsMysqlCore

class DiffSummarizer(CuffDiff, CuffNorm, EnsMysqlCore):
    def __init__(self):
        self.deg_dic = dict()

    def load_diff(self, diff_fn, deg_name):
        diff_dic = self.diff_to_dic(diff_fn)
        for test_id, info_dic in diff_dic.items():
            self.deg_dic.setdefault(test_id, {}).setdefault(deg_name, info_dic)

    def load_anno(self, fn):
        self.anno_dic = self.gene_to_dic(fn)

    def load_norm(self, fn):
        self.norm_dic = self.norm_to_dic(fn)
        self.sample_id_s = self.extract_sample_s(fn)

    def make_result(self, out_fn, deg_name_s):
        out_fh = open(out_fn, 'w')

        anno_headers = ['gene_id','locus','gene_name','biotype','gene_desc']
        norm_headers = ['EXP:{0}:RPKM'.format(sample_id) for sample_id in self.sample_id_s]
        deg_headers = []
        for deg_name in deg_name_s:
            deg_headers.append('DEG:{0}:log2fc'.format(deg_name))
            deg_headers.append('DEG:{0}:p'.format(deg_name))
            deg_headers.append('DEG:{0}:q'.format(deg_name))
            deg_headers.append('DEG:{0}:UPDOWN'.format(deg_name))
            deg_headers.append('DEG:{0}:YN'.format(deg_name))

        out_fh.write('{0}\t{1}\t{2}\n'.format('\t'.join(anno_headers),
                                              '\t'.join(norm_headers),
                                              '\t'.join(deg_headers)))

        for test_id, deg_name_dic in self.deg_dic.items():
            deg_items = list()
            for deg_name in deg_name_s:
                info_dic = deg_name_dic[deg_name]
                log2fc = info_dic['log2fc']
                p = info_dic['p']
                q = info_dic['q']
                #
                deg_items.append(log2fc)
                deg_items.append(p)
                deg_items.append(q)
                if float(log2fc) < 0:
                    deg_items.append('DOWN')
                elif float(log2fc) > 0:
                    deg_items.append('UP')
                else:
                    deg_items.append('FLAT')
                if abs(float(log2fc)) >= 0.584 and float(p) < 0.01:
                    deg_items.append('Y')
                else:
                    deg_items.append('N')
                #
            anno_items = [test_id, info_dic['locus'], info_dic['gene_name'],
                          self.anno_dic[test_id]['biotype'], self.anno_dic[test_id]['desc']]
            norm_items = [self.norm_dic[test_id][sample_name] for sample_name in self.sample_id_s]
            out_fh.write('{0}\t{1}\t{2}\n'.format('\t'.join(anno_items),
                                                  '\t'.join(norm_items),
                                                  '\t'.join(deg_items)))
        out_fh.close()



def main():
    ds = DiffSummarizer()

    ds.load_diff('./do_deg/DE_AE-AKTE/gene_exp.diff', 'AKTP +/M vs. AKTP-null')
    ds.load_diff('./do_deg/DE_AE-AKTP_LOH/gene_exp.diff', 'AKTP-LOH vs. AKTP-null')
    ds.load_diff('./do_deg/DE_AE-AKTP_M/gene_exp.diff', 'AKTP +/M vs. AKTP-LOH')
    ds.load_diff('./do_deg/DE_AE-AKTP_null/gene_exp.diff', 'AKT vs. AKTP +/M')
    ds.load_diff('./do_deg/DE_AKTE-AKTP_LOH/gene_exp.diff', 'AKT vs. AKTP-LOH')
    ds.load_diff('./do_deg/DE_AKTE-AKTP_M/gene_exp.diff', 'AKT vs. AKTP -null')
    ds.load_diff('./do_deg/DE_AKTE-AKTP_null/gene_exp.diff', 'A vs. AKT')
    ds.load_diff('./do_deg/DE_AKTP_LOH-AKTP_null/gene_exp.diff', 'A vs. AKTP +/M')
    ds.load_diff('./do_deg/DE_AKTP_M-AKTP_LOH/gene_exp.diff', 'A vs. AKTP-null')
    ds.load_diff('./do_deg/DE_AKTP_M-AKTP_null/gene_exp.diff', 'A vs. AKTP-LOH')

    ds.load_anno('/BiO/BioProjects/HCP-Nakayama-RNASeq-20190422/anno/gene.txt')

    ds.load_norm('./do_cuffnorm/genes.fpkm_table')

    deg_name_s = ['AKTP +/M vs. AKTP-null', 'AKTP-LOH vs. AKTP-null', 'AKTP +/M vs. AKTP-LOH',
                  'AKT vs. AKTP +/M', 'AKT vs. AKTP-LOH', 'AKT vs. AKTP -null',
                  'A vs. AKT', 'A vs. AKTP +/M', 'A vs. AKTP-null', 'A vs. AKTP-LOH']

    ds.make_result('./genes.xls', deg_name_s)


if __name__=='__main__':
    main()

