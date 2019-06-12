
## Reference
# https://asia.ensembl.org/info/docs/api/core/core_schema.html
# ftp://ftp.ensembl.org/pub/release-85/mysql/bos_taurus_core_85_31/gene.txt.gz
# http://rest.ensembl.org/
# https://github.com/Ensembl/ensembl-rest

import gzip
import sys


class EnsMysqlCore:

    def gene_to_dic(self, fn):
        anno_dic = dict()

        if fn.endswith('.gz'):
            fh = gzip.open(fn)
        else:
            fh = open(fn)

        for line in fh:
            items = line.rstrip('\n').split('\t')

            gene_id = items[0]
            biotype = items[1]
            analysis_id = items[2]
            seq_region_id = items[3]
            seq_region_start = items[4]
            seq_region_end = items[5]
            seq_region_strand = items[6]
            display_xref_id = items[7]
            source = items[8]
            status = items[9]
            description = items[10]
            is_current = items[11]
            canonical_transcript_id = items[12]
            stable_id = items[13]
            version = items[14]
            created_date = items[15]
            modified_date = items[16]
            #
            if is_current in ['1']:
                anno_dic.setdefault(stable_id, {}).setdefault('biotype', biotype)
                anno_dic.setdefault(stable_id, {}).setdefault('desc', description)
            else:
                pass

        fh.close()
        return anno_dic



