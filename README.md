# ConnecTF Backend

## Introduction

This is the API backend for ConnecTF.

## Install

```bash
pip install -r requirements.txt
```

A MySQL database instance is required. Details of the database should be configured in the [`config.yaml`](#configyaml) file in the 
`./connectf/` folder, alongside `settings.py`. You can edit `settings.py` directly if you are more comfortable in
configuring a Django project.

```bash
python manage.py migrate
```

## Import Data

Import data before starting the server.

```bash
python manage.py import_annotation annotation.csv  # import gene annotations
python manage.py import_data data.csv metadata.txt  # import data/metadata
python manage.py import_edges additional_edges.txt  # import additional edges
```

Sample files can be found at:

## Configuration

### config.yaml

A sample `config.yaml` file:

```yaml
SECRET_KEY: 'django_secret_key'  # see https://docs.djangoproject.com/en/2.2/ref/settings/#secret-key
DEBUG: True
DATABASE:
  NAME: 'db_name'
  USER: 'db_username'
  PASSWORD: 'db_password'
  HOST: 'localhost'
MOTIF_ANNOTATION: '/path/to/file'  # path to cluster motif file motifs.csv.gz
MOTIF_TF_ANNOTATION: '/path/to/file'  # path to tf motif file motifs_tf.csv.gz
MOTIF_CLUSTER_INFO: '/path/to/file'  # path to cluster_info.csv.gz
GENE_LISTS: '/path/to/folder'
TARGET_NETWORKS: '/path/to/folder'
```

## Deploying

### Development

```bash
python manage.py runserver 8000
```

### Production

Deploy using Gunicorn

```bash
pip install gunicorn
gunicorn --workers 5 --timeout 200 --bind unix:connectf_backend.sock -m 007 connectf.wsgi
# use "nohup gunicorn [OPTIONS] &" to run in background

```

This binds the server to a unix socket, which can then be connected to from a reverse proxy such as nginx.

### Sample Nginx Server Configuration

This listens to an HTTPS connection. Remember to include certificates and private keys in the configuration, 
or use an HTTP configuration instead.

```text
server {
        listen [::]:443 ssl http2;
        listen 443 ssl http2;
        server_name example.com;
        ssl_certificate /path/to/cert.cer;
        ssl_certificate_key /path/to/private_key.pem;
        ssl_protocols TLSv1.2;
        ssl_ciphers HIGH:!aNULL:!MD5;

        add_header Strict-Transport-Security "max-age=86400; includeSubDomains" always;

        client_max_body_size 20M;

        root /var/www/html; # path to html files

        index index.html;

        location / {
                # First attempt to serve request as file, then
                # as directory, then fall back to displaying a 404.
                try_files $uri /index.html;
        }

        location ~* ^/(api|queryapp)/ {
            # include proxy_params;
            proxy_set_header Host $http_host;
            proxy_set_header X-Real-IP $remote_addr;
            proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
            proxy_set_header X-Forwarded-Proto $scheme;
            proxy_connect_timeout       3600;
            proxy_send_timeout          3600;
            proxy_read_timeout          3600;
            send_timeout                3600;
            proxy_pass http://unix:/path/to/connectf_backend.sock;
        }
}

```
