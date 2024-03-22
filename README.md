# HiNT - tsv2jsonld

JSON-LD Converter for HiNT

## Downloading data

Download data from the following web site.

    http://hint.yulab.org/download/

Store the source file(tsv) to be converted in the following directory.

    data/input

![image](https://github.com/med2rdf/document/blob/master/hint_file_download.png)

## Cinfiguration

Edit the following definition file.

    app/tsvjsonld_hint.json

## Usage

### Build Docker Image

Run the following command to execute a Docker image.
```bash
docker build --no-cache -t tsv2jsonld_hint .
```

### Execute the conversion
Run the following command to execute the conversion process.
```bash
docker run -v "$(pwd)/data:/usr/src/app/data" tsv2jsonld_hint
```
