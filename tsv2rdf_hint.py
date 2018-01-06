#coding: UTF-8
'''
Created on 2018/01/06
'''
import sys,os
from jinja2 import Template, Environment, FileSystemLoader
from optparse import OptionParser
import json
import logging.handlers
import logging

IN_DATA_DIR = 'data_h'
RDF_TEMPLATE='templ_hint.ttl'
MOUNT_PATH = '/mnt'
TEMPLATE_BODY = 'templ_hint.ttl'
TEMPLATE_EVI = 'templ_hint.ttl.evi'
TEMPLATE_PREFIX = 'templ_hint.ttl.prefix'

class Message(object):
    def __init__(self, subject_no, data ):
        '''
        @param subject_no: 主語に付与する連番
        @param data: HiNTのデータ（1レコード）
        @summary: HiNTのデータ（1レコード）をMessageオブジェクトに格納する
        '''

        self.SNo = subject_no
        self.UniprotA = data["Uniprot_A"]
        self.UniprotB = data["Uniprot_B"]
        self.GeneA = data["Gene_A"]
        self.GeneB = data["Gene_B"]
        self.ORFA = data["ORF_A"]
        self.ORFB = data["ORF_B"]
        self.AliasA = data["Alias_A"]
        self.AliasB = data["Alias_B"]

        #Publications
        self.Publications = [ Publications(comma, publication)
                             for comma, publication in enumerate(data["pmid:method:quality"].split('|')) ]

        self.hilo = data["hi/lo"]
        self.proteinForm = data["proteinForm"]

        return

class Publications(object):
    def __init__(self, comma, publication ):
        '''
        @param comma: 複数のPublicationを出力する際のカンマの出力制御用の変数
        @param publication: 文献のID
        @summary: 複数のPublicationをPublicationsオブジェクトに格納する
        '''

        #セミコロンを含むデータがあるためエスケープする
        publication.replace(";","\;")

        wk_publication = publication.split(':')

        self.Comma = comma
        self.DB = "pdb" if (wk_publication[0][0:4] == "PDB_") else "pubmed"
        self.ID = wk_publication[0][4:].lower() if (wk_publication[0][0:4] == "PDB_") else wk_publication[0]
        self.Method = wk_publication[1]
        self.Quality = "high-throughput" if (wk_publication[2] == "HT") else "literature-curated" if (wk_publication[2] == "LC") else "unknown"   #I列[2]

        return

def read_tsv(bi_all, bi_hq, cc_al, cc_hq):
    import pandas as pd

    hint_bi_allDF=pd.read_csv(bi_all, delimiter="\t")
    hint_cc_allDF=pd.read_csv(cc_al, delimiter="\t")

    if cc_hq is not None:
        hint_bi_hqDF=pd.read_csv(bi_hq, delimiter="\t")
    else:
        hint_bi_hqDF = pd.DataFrame(columns=hint_bi_allDF.columns)

    if cc_hq is not None:
        hint_cc_hqDF=pd.read_csv(cc_hq, delimiter="\t")
    else:
        hint_cc_hqDF = pd.DataFrame(columns=hint_cc_allDF.columns)

    qdiffDF=pd.merge(hint_bi_allDF, hint_bi_hqDF, on=["Uniprot_A", "Uniprot_B"], how="outer")
    hint_bi_onlyhqDF=qdiffDF[(qdiffDF["pmid:method:quality_y"].fillna("nan") == "nan")][[
        "Uniprot_A",
        "Uniprot_B",
        "Gene_A_x",
        "Gene_B_x",
        "ORF_A_x",
        "ORF_B_x",
        "Alias_A_x",
        "Alias_B_x",
        "pmid:method:quality_x"]]
    hint_bi_onlyhqDF.columns=[
        "Uniprot_A",
        "Uniprot_B",
        "Gene_A",
        "Gene_B",
        "ORF_A",
        "ORF_B",
        "Alias_A",
        "Alias_B",
        "pmid:method:quality"]
    hint_bi_onlyhqDF["hi/lo"]="high"
    hint_bi_onlyhqDF["proteinForm"]="binary"

    hint_bi_onlylqDF=qdiffDF[(
        qdiffDF["pmid:method:quality_y"].fillna("nan") != "nan")
        & (qdiffDF["pmid:method:quality_x"].fillna("nan") != "nan")
        & (qdiffDF["pmid:method:quality_x"].fillna("nan") == qdiffDF["pmid:method:quality_y"].fillna("nan"))][
            ["Uniprot_A",
             "Uniprot_B",
             "Gene_A_x",
             "Gene_B_x",
             "ORF_A_x",
             "ORF_B_x",
             "Alias_A_x",
             "Alias_B_x",
             "pmid:method:quality_x"]]
    hint_bi_onlylqDF.columns=[
        "Uniprot_A",
        "Uniprot_B",
        "Gene_A",
        "Gene_B",
        "ORF_A",
        "ORF_B",
        "Alias_A",
        "Alias_B",
        "pmid:method:quality"]
    hint_bi_onlylqDF["hi/lo"]="low"
    hint_bi_onlylqDF["proteinForm"]="binary"


    hint_bi_onlyhq2DF=qdiffDF[(qdiffDF["pmid:method:quality_x"].fillna("nan") == "nan")][[
        "Uniprot_A",
        "Uniprot_B",
        "Gene_A_y",
        "Gene_B_y",
        "ORF_A_y",
        "ORF_B_y",
        "Alias_A_y",
        "Alias_B_y",
        "pmid:method:quality_y"]]
    hint_bi_onlyhq2DF.columns=[
        "Uniprot_A",
        "Uniprot_B",
        "Gene_A",
        "Gene_B",
        "ORF_A",
        "ORF_B",
        "Alias_A",
        "Alias_B",
        "pmid:method:quality"]
    hint_bi_onlyhq2DF["hi/lo"]="high"
    hint_bi_onlyhq2DF["proteinForm"]="binary"

    hint_biDF=pd.concat([hint_bi_onlyhqDF, hint_bi_onlyhq2DF, hint_bi_onlylqDF]).sort_index()

    qdiff2DF=pd.merge(hint_cc_allDF, hint_cc_hqDF, on=["Uniprot_A", "Uniprot_B"], how="outer")
    hint_cc_onlyhqDF=qdiff2DF[(qdiff2DF["pmid:method:quality_y"].fillna("nan") == "nan")][[
        "Uniprot_A",
        "Uniprot_B",
        "Gene_A_x",
        "Gene_B_x",
        "ORF_A_x",
        "ORF_B_x",
        "Alias_A_x",
        "Alias_B_x",
        "pmid:method:quality_x"]]
    hint_cc_onlyhqDF.columns=[
        "Uniprot_A",
        "Uniprot_B",
        "Gene_A",
        "Gene_B",
        "ORF_A",
        "ORF_B",
        "Alias_A",
        "Alias_B",
        "pmid:method:quality"]
    hint_cc_onlyhqDF["hi/lo"]="low"
    hint_cc_onlyhqDF["proteinForm"]="co-complex"

    hint_cc_onlylqDF=qdiff2DF[(qdiff2DF["pmid:method:quality_y"].fillna("nan") != "nan") & (qdiff2DF["pmid:method:quality_x"].fillna("nan") != "nan")
       & (qdiff2DF["pmid:method:quality_x"].fillna("nan") == qdiff2DF["pmid:method:quality_y"].fillna("nan"))][
           ["Uniprot_A",
            "Uniprot_B",
            "Gene_A_y",
            "Gene_B_y",
            "ORF_A_y",
            "ORF_B_y",
            "Alias_A_y",
            "Alias_B_y",
            "pmid:method:quality_y"]]
    hint_cc_onlylqDF.columns=[
        "Uniprot_A",
        "Uniprot_B",
        "Gene_A",
        "Gene_B",
        "ORF_A",
        "ORF_B",
        "Alias_A",
        "Alias_B",
        "pmid:method:quality"]
    hint_cc_onlylqDF["hi/lo"]="high"
    hint_cc_onlylqDF["proteinForm"]="co-complex"

    hint_ccDF=pd.concat([hint_cc_onlyhqDF, hint_cc_onlylqDF]).sort_index()

    hint_DF=pd.concat([hint_biDF, hint_ccDF]).reset_index()
    del hint_DF["index"]

    return hint_DF

def main(fw, template, organism=None, bi_all = None, bi_hq = None, cc_al = None, cc_hq = None ):

    # テンプレートの読み込み
    env = Environment(loader = FileSystemLoader(".", encoding='utf8'), autoescape = False)
    imageTemplate = env.get_template(template['body'])

    evidenceList=[] #エビデンス

    hint_DF = read_tsv(bi_all, bi_hq, cc_al, cc_hq)

    for row in hint_DF.iterrows():

        #Message Body
        msgObject = Message(row[0], row[1])
        evidenceList.extend([ '\t'.join([pub.DB,pub.ID,pub.Method,pub.Quality]) for pub in msgObject.Publications])

        namespace = dict(message=msgObject)
        FeedContent = imageTemplate.render(namespace, organism=organism).encode('utf-8')

        fw.write(FeedContent)

    # テンプレートの読み込み(evidence)
    env = Environment(loader = FileSystemLoader(".", encoding='utf8'), autoescape = False)
    imageTemplate = env.get_template(template['evidence'])

    evidences_uq_wk = list(set(evidenceList))
    evidences_uq_wk2 = map(lambda x: x.split('\t'), evidences_uq_wk)
    evidences = [dict(DB=x[0],ID=x[1],Method=x[2],Quality=x[3]) for x in  evidences_uq_wk2]

    namespace = dict(Publications=evidences)
    FeedContent = imageTemplate.render(namespace).encode('utf-8')
    fw.write(FeedContent)

    return


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option( "-c", type="string", help="config_file", dest="config", default= 'tsv2rdf_hint.json')
    (options, args) = parser.parse_args()

    #ログ関連の設置
    argvs = sys.argv
    LOG_FILENAME = '%s.log' % argvs[0].split('.')[0]
    logging.basicConfig(format='%(asctime)s:%(levelname)s:%(message)s',level=logging.DEBUG)
    logger = logging.getLogger()

    #コンソールログ
    '''
    stream_log = logging.StreamHandler()
    stream_log.setLevel(logging.DEBUG)
    stream_log.setFormatter(logging.Formatter('%(asctime)s;%(levelname)s:%(message)s'))
    '''

    #ファイルログ
    file_log = logging.handlers.RotatingFileHandler(filename=LOG_FILENAME)
    file_log.setLevel(logging.DEBUG)
    file_log.setFormatter(logging.Formatter('%(asctime)s;%(levelname)s:%(message)s'))

    #ハンドラの設定
    #logging.getLogger().addHandler(stream_log)
    logging.getLogger().addHandler(file_log)

    #定義ファイルの読み込み
    f = open(options.config ,"r")
    config = json.load(f)
    logger.info("Organism: %s" % ",".join([ k for k,v in config['organism'].items() ]))

    output_file = config['output_file']
    template = config['template']
    data_path = config['data_path']

    #出力先ファイル
    if os.path.exists(output_file) != True:
        output_file = os.path.join(MOUNT_PATH, output_file)
        fw = open(output_file, 'w')
    else:
        fw = open(output_file, 'w')

    #テンプレートファイルのチェック
    if os.path.exists(template['body']) != True:
        template['body'] = os.path.join(MOUNT_PATH, template['body'])
        logger.info("Set new template file path: %s" % template['body'])
    if os.path.exists(template['evidence']) != True:
        template['evidence'] = os.path.join(MOUNT_PATH, template['evidence'])
        logger.info("Set new template file path: %s" % template['evidence'])
    if os.path.exists(template['prefix']) != True:
        template['prefix'] = os.path.join(MOUNT_PATH, template['prefix'])
        logger.info("Set new template file path: %s" % template['prefix'])

    if os.path.exists(data_path) != True:
        data_path = os.path.join(MOUNT_PATH, data_path)
        logger.info("Set new data_path: %s" % data_path)

    #prefixの出力
    fp = open(template['prefix'], 'r')
    prefixes = fp.readlines()
    for prefix in prefixes:
        fw.write(prefix)

    for k,v in config['organism'].items():

        if 'binary_hq' in v:
            bi_hq = os.path.join(data_path, v['binary_hq'])
        else:
            bi_hq = None

        if 'cocomp_hq' in v:
            cc_hq = os.path.join(data_path, v['cocomp_hq'])
        else:
            cc_hq = None

        main(
            fw,
            template,
            organism = k,
            bi_all = os.path.join(data_path, v['binary_all']),
            bi_hq = bi_hq,
            cc_al = os.path.join(data_path, v['cocomp_all']),
            cc_hq = cc_hq
            )
    pass
