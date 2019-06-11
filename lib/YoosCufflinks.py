

import math

class CuffDiff:

    def diff_to_dic(self, diff_fn):
        diff_dic = dict()
        for line in open(diff_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['test_id']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item, idx)
                continue
            #
            test_id   = items[0]
            gene_id   = items[1]
            gene_name = items[2]
            locus     = items[3]
            sample_1  = items[4]
            sample_2  = items[5]
            value_1   = items[7]
            value_2   = items[8]
            _log2fc    = items[9]
            p         = items[11]
            q         = items[12]
            #
            log2fc = str(math.log(float(value_2)+1, 2) - math.log(float(value_1)+1, 2))
            #
            diff_dic.setdefault(test_id, {}).setdefault('gene_id', gene_id)
            diff_dic.setdefault(test_id, {}).setdefault('gene_name', gene_name)
            diff_dic.setdefault(test_id, {}).setdefault('locus', locus)
            diff_dic.setdefault(test_id, {}).setdefault('sample_1', sample_1)
            diff_dic.setdefault(test_id, {}).setdefault('sample_2', sample_2)
            diff_dic.setdefault(test_id, {}).setdefault('value_1', value_1)
            diff_dic.setdefault(test_id, {}).setdefault('value_2', value_2)
            diff_dic.setdefault(test_id, {}).setdefault('log2fc', log2fc)
            diff_dic.setdefault(test_id, {}).setdefault('p', p)
            diff_dic.setdefault(test_id, {}).setdefault('q', q)
        return diff_dic


class CuffNorm:
    def norm_to_dic(self, norm_fn, trim_tag='_0'):
        norm_dic = dict()
        for line in open(norm_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['tracking_id']:
                idx_dic = dict()
                for idx, item in enumerate(items):
                    idx_dic.setdefault(item.rstrip(trim_tag), idx)
                continue
            gene_name = items[0]
            for sample_id, idx in idx_dic.items():
                if sample_id in ['tracking_id']:
                    continue
                norm_dic.setdefault(gene_name, {}).setdefault(sample_id, items[idx])
        return norm_dic

    def extract_sample_s(self, norm_fn, trim_tag='_0'):
        sample_id_s = list()
        for line in open(norm_fn):
            items = line.rstrip('\n').split('\t')
            if items[0] in ['tracking_id']:
                for item in items[1:]:
                    sample_id_s.append(item.rstrip(trim_tag))
                break
        return sample_id_s


