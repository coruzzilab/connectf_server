# ConnecTF Backend

## Introduction

This is the API backend for ConnecTF. The instructions for the for the [react](#front-end-interface) server should be followed after these instructions.

## Pre-requisites
​
### Conda

Following instructions are using conda

```bash
# On linux
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
​
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
​
conda create -n connectf
conda activate connectf
conda install python
conda install nodejs
sudo apt-get install build-essential
​
```
​
### Mysql
​
A MySQL database instance is required. 
​
​
```bash
sudo apt install mariadb-server
sudo apt install mariadb-client
sudo apt-get install libmariadbclient-dev
​
#TODO: confirm the following are required
#sudo apt-get install gcc
#pip install mysql-connector-python
​
sudo system mysql start
sudo mysql_secure_installation
# set a password for user root on mysql
```
​
Now create database with username and password that will be accessed for connectf.
Details of the database should be configured in the [`config.yaml`](#configyaml) file in the `./connectf/` folder, alongside `settings.py`. You can edit `settings.py` directly if you are more comfortable in configuring a Django project.
​
​
- db_name=connectf
- db_username=connectfuser
- db_passwd=connectfpwd
​
```bash
sudo mysqladmin create connectf
sudo mysql -u root
grant all privileges on connectf.* to 'connectfuser'@'localhost' identified by 'connectfpwd';
flush privileges;
​
```


## Install

```bash
pip install -r requirements.txt
```

Create a [`config.yaml`](#configyaml) in the `./connectf` folder, alongside the `./connectf/settings.py` file.

To set up data files for the server to read from, it is recommended you create a `./data` folder at the top level of the project. Put required data files within the `./data` folder and edit [`config.yaml`](#configyaml) to reflect the changes.

A MySQL database instance is required. Details of the database should be configured in the [`config.yaml`](#configyaml) file. You can edit `settings.py` directly if you are more comfortable in
configuring a Django project.

```bash
python manage.py migrate
```

## Import Data

Import data before starting the server.

```bash
python manage.py import_annotation -i annotation.csv  # import gene annotations
python manage.py import_data data.csv metadata.txt  # import data/metadata
python manage.py import_edges additional_edges.txt  # import additional edges
```

Sample files can be found at:

## Configuration

### config.yaml

A sample `config.yaml` file:

*N.B.* If the file does not exist, create a new `config.yaml` in the same folder as `connectf/settings.py`, with the contents similar to the one seen in the sample.

```yaml
SECRET_KEY: 'django_secret_key'  # see https://docs.djangoproject.com/en/2.2/ref/settings/#secret-key
DEBUG: True
DATABASE:
  NAME: 'db_name'
  USER: 'db_username'
  PASSWORD: 'db_password'
  HOST: 'localhost'
MOTIF_ANNOTATION: '/path/to/file'  # path to cluster motif file motifs.csv.gz
MOTIF_TF_ANNOTATION: '/path/to/file'  # path to tf motif file motifs_indv.csv.gz
MOTIF_CLUSTER_INFO: '/path/to/file'  # path to cluster_info.csv.gz
GENE_LISTS: '/path/to/folder'  # optional gene list folder
TARGET_NETWORKS: '/path/to/folder' # optional target network folder
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

This listens to an HTTPS connection. Remember to include certificates and private keys in the configuration, or use an HTTP configuration instead.

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

        client_max_body_size 100M; # ensure file size is big enough for user upload

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

### Front End Interface

The front end of this project is built with ReactJS. You can find it at [connectf_react](https://github.com/coruzzilab/connectf_react).

### Containerization

To run all three -databse, backend and frontend- in a containerized environment, use the provided `Dockerfiles` and `docker-compose.yml` files.

- 1. The docker compose script expects all gene data to be located at `connectf_server/connectf_data_release_v1/` folder, following the structure:
```bash
connectf_server
├── connectf_data_release_v1
│   ├── arabidopsis
│   │   ├── additional_edges/
│   │   ├── data/
│   │   ├── gene_lists/
│   │   ├── metadata/
│   │   ├── motifs/
│   │   ├── networks/
│   │   ├── annotation.csv.gz
│   ├── common
│   │   ├── cluster_info.csv.gz
│   ├── maize
│   │   ├── same structure as arabidopsis
│   ├── rice
│   │   ├── same structure as arabidopsis
```
- 2. Make sure the folders for projects `connectf_react` and `connectf_server` are siblings.
- 2. cd into connectf_react and build a production version of the React app
```bash
cd ../connectf_react
npm run build
```
- 3. cd back to connectf_server and make sure the `./connectf/config.yaml` file is set up correctly.
```yaml
SECRET_KEY: 'django_secret_key' # see https://docs.djangoproject.com/en/2.2/ref/settings/#secret-key
DEBUG: True
DATABASE:
  NAME: 'connectf'
  USER: 'connectfuser'
  PASSWORD: 'connectfpwd'
  HOST: 'localhost'

MOTIF_ANNOTATION: './connectf_data_release_v1/$IMPORT/motifs/motifs.csv.gz'  # path to cluster motif file motifs.csv.gz
MOTIF_TF_ANNOTATION: './connectf_data_release_v1/$IMPORT/motifs/motifs_indv.csv.gz'  # path to tf motif file motifs_indv.csv.gz
MOTIF_CLUSTER_INFO: './connectf_data_release_v1/common/cluster_info.csv.gz'  # path to cluster_info.csv.gz
GENE_LISTS: './connectf_data_release_v1/$IMPORT/gene_lists'  # optional gene list folder
TARGET_NETWORKS: './connectf_data_release_v1/$IMPORT/networks' # optional target network folder
```

- 4. Make sure there are no running containers:
```bash
docker compose down
```

- 5. Finally, build and run images:
```bash
COMPOSE_PROJECT_NAME=connectf IMPORT=arabidopsis docker compose up --build
```
The env var `IMPORT` is always required and signals to install all your gene data into the mysql database if this is the first time running the containers. `IMPORT` is a string that specifies the folder that contains your gene data, for example, 'arabiopsis' or 'maize' or 'rice'.

- 5. The api should be running on `http://localhost:8001/api` and the frontend should be running on `http://localhost:80`

<img width="1155" alt="image" src="https://github.com/coruzzilab/connectf_server/assets/1679438/cab90ef1-a3d5-47b5-af47-de85cb6debf8">

### Data

All data for this project can be found at https://connectf.s3.amazonaws.com/connectf_data_release_v1.tar.gz

