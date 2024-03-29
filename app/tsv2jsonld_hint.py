import copy
import json
import logging
import logging.handlers
import os
import sys
import time
from collections import Counter
from optparse import OptionParser
from typing import Hashable
from urllib.parse import urljoin, urlparse

import pandas as pd
import requests
from bs4 import BeautifulSoup

message_json_base = {
    "@context": "",  # tsv2jsonld_hint.jsonに定義したファイルパスを指定.
    "@id": "",  # f"hint:{message.UniprotA}-{message.UniprotB}"
    "label": "",  # f"{message.UniprotA}-{message.UniprotB}"
    "participant": [],  # [f"up:{message.UniprotA}", f"up:{message.UniprotB}"]
    "evidence": [],  # [{"reference": f"pmid:{Publication.ID}", "method": f"obo:{Publication.Method}", "quality": Publication.Quality}] message.Publications分のdictを出力する。
    "category": "",  # message.proteinForm
    "accuracy": "",  # message.hilo
    "taxonomy": "",  # f"taxid:{message.taxid}"
    "provenance": "",  # f"taxid:{message.dataSource}"
}

evidence_json_base = {
    "@context": "",  # tsv2jsonld_hint.jsonに定義したファイルパスを指定.
    "@id": "",  # f"pint:evidence-{Publication.DB}_{Publication.ID}-MI_{Publication.Method}-{Publication.Quality}"
    "@type": "bp3:Evidence",
    "evidenceCode": "",  # f"obo:MI_{Publication.Method}"
    "reference": "",  # f"http://identifiers.org/{{Publication.DB}}/{{Publication.ID}}"
    "dataSourceType": "",  # Publication.Quality
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
        self.dataSource: str = data["dataSource"]

    def to_json(message, context_path: str, taxid: str) -> dict:
        message_json = copy.deepcopy(message_json_base)
        message_json["@context"] = context_path
        message_json["@id"] = f"hint:{message.UniprotA}-{message.UniprotB}"
        message_json["label"] = f"{message.UniprotA}-{message.UniprotB}"
        message_json["participant"] = [
            f"up:{message.UniprotA}",
            f"up:{message.UniprotB}",
        ]
        evidences = [
            {
                "reference": f"pmid:{publication.ID}",
                "method": f"obo:MI_{publication.Method}",
                "quality": publication.Quality,
            }
            for publication in message.Publications
        ]
        message_json["evidence"] = evidences
        message_json["category"] = message.proteinForm
        message_json["accuracy"] = message.hilo
        message_json["taxonomy"] = f"taxid:{taxid}"
        message_json["provenance"] = message.dataSource

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
        # self.Quality: str = "high-throughput" if (wk_publication[2] == "HT") else "literature-curated" if (wk_publication[2] == "LC") else "unknown"   #I列[2]
        self.Quality: str = wk_publication[2]


def to_evidence_json(publication: dict, context_path: str) -> dict:
    evidence_json = copy.deepcopy(evidence_json_base)
    db = publication["DB"]
    id = publication["ID"]
    method = publication["Method"]
    quality = publication["Quality"]
    evidence_json["@context"] = context_path
    evidence_json["@id"] = f"pint:evidence-{db}_{id}-MI_{method}-{quality}"
    evidence_json["evidenceCode"] = f"obo:MI_{method}"
    evidence_json["reference"] = f"http://identifiers.org/{db}/{id}"
    evidence_json["dataSourceType"] = quality

    return evidence_json


def read_tsv(
    bi_all: str, bi_hq: str, cc_al: str, cc_hq: str, data_source: str
) -> pd.DataFrame:
    try:
        hint_bi_allDF = pd.read_csv(bi_all, delimiter="\t")
        hint_cc_allDF = pd.read_csv(cc_al, delimiter="\t")
        data_source_bi_all = (
            f"{data_source}/binaly/all/{os.path.basename(bi_all)}"
        )
        data_source_cc_all = (
            f"{data_source}/cocomp/all/{os.path.basename(cc_al)}"
        )

        if cc_hq is not None:
            hint_bi_hqDF = pd.read_csv(bi_hq, delimiter="\t")
            data_source_bi_hq = (
                f"{data_source}/binaly/hq/{os.path.basename(bi_hq)}"
            )
        else:
            hint_bi_hqDF = pd.DataFrame(columns=hint_bi_allDF.columns)
            data_source_bi_hq = ""

        if cc_hq is not None:
            hint_cc_hqDF = pd.read_csv(cc_hq, delimiter="\t")
            data_source_cc_hq = (
                f"{data_source}/cocomp/hq/{os.path.basename(cc_hq)}"
            )
        else:
            hint_cc_hqDF = pd.DataFrame(columns=hint_cc_allDF.columns)
            data_source_cc_hq = ""

    except Exception as e:
        logger.error("%s" % str(e))
        sys.exit()

    # set data source
    hint_bi_allDF["dataSource"] = data_source_bi_all
    hint_cc_allDF["dataSource"] = data_source_cc_all
    hint_bi_hqDF["dataSource"] = data_source_bi_hq
    hint_cc_hqDF["dataSource"] = data_source_cc_hq

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
            "dataSource_x",
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
        "dataSource",
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
            "dataSource_y",
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
        "dataSource",
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
            "dataSource_x",
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
        "dataSource",
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
            "dataSource_y",
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
        "dataSource",
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
    output_dir: str = os.path.join(get_script_dir(), config["output_dir"])
    output_file_prefix: str = config["output_file_prefix"]
    context_path: str = config["context_path"]

    for k, v in organism.items():
        jsonls = []

        data_source: str = v["data_source"]
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

        hint_DF = read_tsv(bi_all, bi_hq, cc_al, cc_hq, data_source)
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
            jsonls.append(msgObject.to_json(context_path, v["taxid"]))

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
        with open(
            os.path.join(output_dir, f"{output_file_prefix}{k}.jsonl"), "w"
        ) as f:
            for jsonl in jsonls:
                # jsonl["@context"] = context
                line = json.dumps(jsonl)
                f.write(line + "\n")


def get_script_dir() -> str:
    script_path = os.path.abspath(__file__)
    return os.path.dirname(script_path)


def extract_links_with_keyword_and_suffix(
    url: str, keyword: list[int | str], suffix: list[str]
) -> list[str]:
    """
    Extracts links from a specified URL that contain a given keyword and
    end with a specified suffix.

    This function sends an HTTP GET request to retrieve the content of
    the URL provided. It then parses the HTML content to find all hyperlinks,
    filters those that meet the criteria of containing the keyword and
    ending with the suffix, and returns a list of qualifying URLs.

    Args:
        url (str): The target URL from which to extract links.
        keyword (list[int | str]): The keyword for which to filter the links.
        suffix (list[str]): The suffix that the links should end with.

    Returns:
        list[str]: A list of strings, each representing a link that matches
        the specified conditions.

    Raises:
        requests.exceptions.HTTPError: For HTTP errors during URL retrieval.
        requests.exceptions.ConnectionError: For connection errors during
        URL retrieval.
    """

    filtered_links = []

    try:
        response = requests.get(url)
        response.raise_for_status()  # HTTPエラーチェック
    except requests.exceptions.HTTPError as e:
        logger.exception(f"HTTPエラーが発生しました: {e}")
        raise

    except requests.exceptions.ConnectionError as e:
        logger.exception(f"接続エラーが発生しました: {e}")
        raise

    soup = BeautifulSoup(response.text, "html.parser")
    links = soup.find_all("a")

    for link in links:
        href = link.get("href")
        if not href:
            continue

        # 相対URLを絶対URLに変換
        href = urljoin(url, href)

        # パスとクエリを取得
        parsed_url = urlparse(href)
        path = parsed_url.path

        # パスが指定したサフィックスで終わり、かつキーワードを含む場合
        if any(path.endswith(s) for s in suffix) and any(
            str(kw) in path for kw in keyword
        ):
            filtered_links.append(href)

    return filtered_links


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

    url: str = config["hint_url"]
    keyword: list[int | str] = config["organism"].keys()
    suffix: list[str] = [
        "/binary/all/",
        "/binary/hq/",
        "/cocomp/all/",
        "/cocomp/hq/",
    ]

    links = extract_links_with_keyword_and_suffix(url, keyword, suffix)

    link_aggr_list = copy.deepcopy(links)

    for i in range(len(link_aggr_list)):
        for s in suffix:
            if link_aggr_list[i].endswith(s):
                link_aggr_list[i] = link_aggr_list[i].removesuffix(s)

    link_counts = Counter(link_aggr_list)

    for item, count in link_counts.items():
        if count == 4:
            org = item.split("/")[-1]

            org_conf = {}

            org_conf["taxid"] = config["organism"][org]["taxid"]

            org_conf["data_source"] = item

            org_conf["binary_all"] = f"{org}_binary_all.txt"

            org_conf["binary_hq"] = f"{org}_binary_hq.txt"

            org_conf["cocomp_all"] = f"{org}_cocomp_all.txt"

            org_conf["cocomp_hq"] = f"{org}_cocomp_hq.txt"

            config["organism"][org] = org_conf

        else:
            raise Exception("")

    for link in links:
        splitted_link = link.split("/")

        org = splitted_link[-4]

        file_type = "_".join([splitted_link[-3], splitted_link[-2]])

        with requests.get(link, stream=True) as r:
            with open(
                os.path.join(
                    config["data_path"], config["organism"][org][file_type]
                ),
                "wb",
            ) as f:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    f.write(chunk)

        time.sleep(3)

    logger.info(
        "Organism: %s" % ",".join([k for k, v in config["organism"].items()])
    )

    main(config)

    print("tsv2jsonld_hint.py end")
