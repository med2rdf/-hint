import copy
import json
import logging
import logging.handlers
import os
import sys
from optparse import OptionParser
from typing import Hashable

import pandas as pd

message_json_base = {
    "@context": "",  # tsv2jsonld_hint.jsonに定義したファイルパスを指定.
    "@id": "",  # f"hnt:{message.UniprotA}-{message.UniprotB}"
    "@type": "bp3:MolecularInteraction",
    "bp3:dataSource": "http://identifiers.org/HINT",
    "bp3:displayName": "",  # f"{message.UniprotA}-{message.UniprotB}"
    "bp3:evidence": [],  # f"pint:evidence-{Publication.DB}_{Publication.ID}-MI_{Publication.Method}-{Publication.Quality}". message.Publications分カンマ区切りで出力する。
    "bp3:name": "",  # f"{message.UniprotA}-{message.UniprotB}"
    "bp3:participant": {
        "obo:BFO_0000051": []
    },  # [{"@id": f"uni:{message.UniprotA}"}, {"@id": f"uni:{message.UniprotB}"}]
    "pint:networkVariety": "",  # message.proteinForm
    "pint:accuracy": "",  # message.hilo
}

evidence_json_base = {
    "@context": "",  # tsv2jsonld_hint.jsonに定義したファイルパスを指定.
    "@id": "",  # f"pint:evidence-{Publication.DB}_{Publication.ID}-MI_{Publication.Method}-{Publication.Quality}"
    "@type": "bp3:Evidence",
    "bp3:evidenceCode": "",  # f"obo:MI_{Publication.Method}"
    "dcterms:references": "",  # f"http://identifiers.org/{{Publication.DB}}/{{Publication.ID}}"
    "pint:dataSourceType": "",  # Publication.Quality
}


class Message(object):
    def __init__(self, subject_no: Hashable, data: pd.Series):
        """
        @param subject_no: Subject serial number
        @param data: HiNT data (1 row)
        @summary: Store HiNT data (1 row) in the Message object
        """

        self.SNo: Hashable = subject_no
        self.UniprotA: str = data["Uniprot_A"]
        self.UniprotB: str = data["Uniprot_B"]
        self.GeneA: str = data["Gene_A"]
        self.GeneB: str = data["Gene_B"]
        self.ORFA: str = data["ORF_A"]
        self.ORFB: str = data["ORF_B"]
        self.AliasA: str = data["Alias_A"]
        self.AliasB: str = data["Alias_B"]

        # Publications
        self.Publications = []
        comma = 0
        for publication in data["pmid:method:quality"].split("|"):

            pmids, method, quality = publication.split(":")
            if ";" in pmids or "," in pmids:
                # if publication has multiple pmid
                pmids = pmids.split(";") if ";" in pmids else pmids.split(",")
                for pmid in pmids:
                    new_publication = "%s:%s:%s" % (
                        pmid.strip(),
                        method.strip(),
                        quality.strip(),
                    )
                    self.Publications.append(
                        Publication(comma, new_publication)
                    )
                    comma += 1
            else:
                self.Publications.append(Publication(comma, publication))
                pass

            comma += 1

        """
        self.Publications = [ Publications(comma, publication)
                             for comma, publication in enumerate(data["pmid:method:quality"].split('|')) ]
        """

        self.hilo: str = data["hi/lo"]
        self.proteinForm: str = data["proteinForm"]

    def to_json(message, context_path: str) -> dict:
        message_json = copy.deepcopy(message_json_base)
        message_json["@context"] = context_path
        message_json["@id"] = f"hnt:{message.UniprotA}-{message.UniprotB}"
        message_json["bp3:displayName"] = (
            f"{message.UniprotA}-{message.UniprotB}"
        )
        evidences = [
            f"pint:evidence-{publication.DB}_{publication.ID}-MI_{publication.Method}-{publication.Quality}"
            for publication in message.Publications
        ]
        message_json["bp3:evidence"] = [
            {"@id": evidence} for evidence in evidences
        ]
        message_json["bp3:name"] = f"{message.UniprotA}-{message.UniprotB}"
        message_json["bp3:participant"]["obo:BFO_0000051"] = [
            {"@id": f"uni:{message.UniprotA}"},
            {"@id": f"uni:{message.UniprotB}"},
        ]
        message_json["pint:networkVariety"] = message.proteinForm
        message_json["pint:accuracy"] = message.hilo

        return message_json


class Publication(object):
    def __init__(self, comma: int, publication: list):
        """
        @param comma: Variables for controlling output of commas when outputting multiple publications
        @param publication: ID of document
        @summary: Store multiple Publications in the Publications object
        """

        wk_publication = publication.split(":")

        self.Comma: int = comma
        self.DB: str = (
            "pdb" if (wk_publication[0][0:4] == "PDB_") else "pubmed"
        )
        self.ID: str = (
            wk_publication[0][4:].lower()
            if (wk_publication[0][0:4] == "PDB_")
            else wk_publication[0]
        )
        self.Method: str = wk_publication[1]
        self.Quality: str = (
            "high-throughput"
            if (wk_publication[2] == "HT")
            else (
                "literature-curated"
                if (wk_publication[2] == "LC")
                else "unknown"
            )
        )  # I列[2]


def to_evidence_json(publication: dict, context_path: str) -> dict:
    evidence_json = copy.deepcopy(evidence_json_base)
    db = publication["DB"]
    id = publication["ID"]
    method = publication["Method"]
    quality = publication["Quality"]
    evidence_json["@context"] = context_path
    evidence_json["@id"] = f"pint:evidence-{db}_{id}-MI_{method}-{quality}"
    evidence_json["bp3:evidenceCode"] = f"obo:MI_{method}"
    evidence_json["dcterms:references"] = f"http://identifiers.org/{db}/{id}"
    evidence_json["pint:dataSourceType"] = quality

    return evidence_json


def read_tsv(bi_all: str, bi_hq: str, cc_al: str, cc_hq: str) -> pd.DataFrame:
    try:
        hint_bi_allDF = pd.read_csv(bi_all, delimiter="\t")
        hint_cc_allDF = pd.read_csv(cc_al, delimiter="\t")

        if cc_hq is not None:
            hint_bi_hqDF = pd.read_csv(bi_hq, delimiter="\t")
        else:
            hint_bi_hqDF = pd.DataFrame(columns=hint_bi_allDF.columns)

        if cc_hq is not None:
            hint_cc_hqDF = pd.read_csv(cc_hq, delimiter="\t")
        else:
            hint_cc_hqDF = pd.DataFrame(columns=hint_cc_allDF.columns)
    except Exception as e:
        logger.error("%s" % e.message)
        sys.exit()

    #
    # binary
    #
    qdiffDF = pd.merge(
        hint_bi_allDF, hint_bi_hqDF, on=["Uniprot_A", "Uniprot_B"], how="outer"
    )

    hint_bi_onlylqDF = qdiffDF[
        (qdiffDF["pmid:method:quality_y"].fillna("nan") == "nan")
    ][
        [
            "Uniprot_A",
            "Uniprot_B",
            "Gene_A_x",
            "Gene_B_x",
            "ORF_A_x",
            "ORF_B_x",
            "Alias_A_x",
            "Alias_B_x",
            "pmid:method:quality_x",
        ]
    ]
    hint_bi_onlylqDF.columns = [
        "Uniprot_A",
        "Uniprot_B",
        "Gene_A",
        "Gene_B",
        "ORF_A",
        "ORF_B",
        "Alias_A",
        "Alias_B",
        "pmid:method:quality",
    ]
    hint_bi_onlylqDF["hi/lo"] = "low"
    hint_bi_onlylqDF["proteinForm"] = "binary"

    hint_bi_onlyhqDF = qdiffDF[
        (qdiffDF["pmid:method:quality_y"].fillna("nan") != "nan")
        & (qdiffDF["pmid:method:quality_x"].fillna("nan") != "nan")
        & (
            qdiffDF["pmid:method:quality_x"].fillna("nan")
            == qdiffDF["pmid:method:quality_y"].fillna("nan")
        )
    ][
        [
            "Uniprot_A",
            "Uniprot_B",
            "Gene_A_x",
            "Gene_B_x",
            "ORF_A_x",
            "ORF_B_x",
            "Alias_A_x",
            "Alias_B_x",
            "pmid:method:quality_x",
        ]
    ]
    hint_bi_onlyhqDF.columns = [
        "Uniprot_A",
        "Uniprot_B",
        "Gene_A",
        "Gene_B",
        "ORF_A",
        "ORF_B",
        "Alias_A",
        "Alias_B",
        "pmid:method:quality",
    ]
    hint_bi_onlyhqDF["hi/lo"] = "high"
    hint_bi_onlyhqDF["proteinForm"] = "binary"

    hint_biDF = pd.concat([hint_bi_onlyhqDF, hint_bi_onlylqDF]).sort_index()

    #
    # co-complex
    #

    qdiff2DF = pd.merge(
        hint_cc_allDF, hint_cc_hqDF, on=["Uniprot_A", "Uniprot_B"], how="outer"
    )

    hint_cc_onlylqDF = qdiff2DF[
        (qdiff2DF["pmid:method:quality_y"].fillna("nan") == "nan")
    ][
        [
            "Uniprot_A",
            "Uniprot_B",
            "Gene_A_x",
            "Gene_B_x",
            "ORF_A_x",
            "ORF_B_x",
            "Alias_A_x",
            "Alias_B_x",
            "pmid:method:quality_x",
        ]
    ]
    hint_cc_onlylqDF.columns = [
        "Uniprot_A",
        "Uniprot_B",
        "Gene_A",
        "Gene_B",
        "ORF_A",
        "ORF_B",
        "Alias_A",
        "Alias_B",
        "pmid:method:quality",
    ]
    hint_cc_onlylqDF["hi/lo"] = "low"
    hint_cc_onlylqDF["proteinForm"] = "co-complex"

    hint_cc_onlyhqDF = qdiff2DF[
        (qdiff2DF["pmid:method:quality_y"].fillna("nan") != "nan")
        & (qdiff2DF["pmid:method:quality_x"].fillna("nan") != "nan")
        & (
            qdiff2DF["pmid:method:quality_x"].fillna("nan")
            == qdiff2DF["pmid:method:quality_y"].fillna("nan")
        )
    ][
        [
            "Uniprot_A",
            "Uniprot_B",
            "Gene_A_y",
            "Gene_B_y",
            "ORF_A_y",
            "ORF_B_y",
            "Alias_A_y",
            "Alias_B_y",
            "pmid:method:quality_y",
        ]
    ]
    hint_cc_onlyhqDF.columns = [
        "Uniprot_A",
        "Uniprot_B",
        "Gene_A",
        "Gene_B",
        "ORF_A",
        "ORF_B",
        "Alias_A",
        "Alias_B",
        "pmid:method:quality",
    ]
    hint_cc_onlyhqDF["hi/lo"] = "high"
    hint_cc_onlyhqDF["proteinForm"] = "co-complex"

    hint_ccDF = pd.concat([hint_cc_onlyhqDF, hint_cc_onlylqDF]).sort_index()

    hint_DF = pd.concat([hint_biDF, hint_ccDF]).reset_index()
    del hint_DF["index"]

    return hint_DF


def main(config: dict):
    organism: dict = config["organism"]
    input_folder: str = os.path.join(get_script_dir(), config["data_path"])
    output_file: str = os.path.join(get_script_dir(), config["output_file"])
    context_path: str = config["context_path"]

    jsonls = []
    for k, v in organism.items():
        if "binary_hq" in v:
            bi_hq = os.path.join(input_folder, v["binary_hq"])
        else:
            bi_hq = None

        if "cocomp_hq" in v:
            cc_hq = os.path.join(input_folder, v["cocomp_hq"])
        else:
            cc_hq = None

        bi_all = os.path.join(input_folder, v["binary_all"])
        cc_al = os.path.join(input_folder, v["cocomp_all"])

        logger.info("Processing organism: %s" % k)

        hint_DF = read_tsv(bi_all, bi_hq, cc_al, cc_hq)
        hint_DF.to_csv("output.csv", index=False)

        evidenceList = []
        for row in hint_DF.iterrows():
            msgObject = Message(row[0], row[1])
            evidenceList.extend(
                [
                    "\t".join([pub.DB, pub.ID, pub.Method, pub.Quality])
                    for pub in msgObject.Publications
                ]
            )
            jsonls.append(msgObject.to_json(context_path))

        evidences_uq_wk = list(set(evidenceList))
        evidences_uq_wk2 = map(lambda x: x.split("\t"), evidences_uq_wk)
        evidences = [
            dict(DB=x[0], ID=x[1], Method=x[2], Quality=x[3])
            for x in evidences_uq_wk2
        ]
        for evidence in evidences:
            jsonls.append(to_evidence_json(evidence, context_path))

        # with open('context.json', 'r') as file:
        #     context = json.load(file)

        # Write out in JSONL format to a file.
        with open(output_file, "w") as f:
            for jsonl in jsonls:
                # jsonl["@context"] = context
                line = json.dumps(jsonl)
                f.write(line + "\n")


def get_script_dir() -> str:
    script_path = os.path.abspath(__file__)
    return os.path.dirname(script_path)


if __name__ == "__main__":
    print("tsv2jsonld_hint.py start")
    parser = OptionParser()
    parser.add_option(
        "-c",
        type="string",
        help="config_file",
        dest="config",
        default="tsv2jsonld_hint.json",
    )
    (options, args) = parser.parse_args()

    # Log setting
    argvs = sys.argv
    LOG_FILENAME = "%s.log" % argvs[0].split(".")[0]
    logging.basicConfig(
        format="%(asctime)s:%(levelname)s:%(message)s", level=logging.DEBUG
    )
    logger = logging.getLogger()

    # File log
    file_log = logging.handlers.RotatingFileHandler(filename=LOG_FILENAME)
    file_log.setLevel(logging.DEBUG)
    file_log.setFormatter(
        logging.Formatter("%(asctime)s;%(levelname)s:%(message)s")
    )

    # Log handler setting
    logging.getLogger().addHandler(file_log)

    # Reading the configuration file

    try:
        with open(os.path.join(get_script_dir(), options.config), "r") as f:
            config = json.load(f)
    except Exception as e:
        logger.error("%s: %s" % (str(e), options.config))
        sys.exit()

    logger.info(
        "Organism: %s" % ",".join([k for k, v in config["organism"].items()])
    )

    main(config)

    print("tsv2jsonld_hint.py end")
