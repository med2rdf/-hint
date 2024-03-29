FROM python:3.11-slim-buster

# 作業ディレクトリを設定
WORKDIR /usr/src/app

# アプリケーションのコードをコピー
COPY ./app .
COPY ./requirements.txt .

# 必要なパッケージがあればここでインストール
RUN pip install --no-cache-dir -r requirements.txt

# コンテナ起動時に実行するコマンドexit
CMD ["python", "tsv2jsonld_hint.py"]