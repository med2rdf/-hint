# setting base image
FROM continuumio/anaconda:latest

ADD tsv2rdf_hint.py /tmp
ADD templ_hint.ttl /tmp
ADD templ_hint.ttl.prefix /tmp
ADD templ_hint.ttl.evi /tmp
WORKDIR "/mnt"
ENTRYPOINT ["python", "/tmp/tsv2rdf_hint.py"]
CMD ["-c","/mnt/tsv2rdf_hint.json"]
