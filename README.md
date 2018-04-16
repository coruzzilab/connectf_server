# TGDBbackend

## Introduction

This is the API backend for TargetDB.

## Install

```bash
pip install -r requirements.txt
python manage.py migrate
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
gunicorn --workers 5 --timeout 200 --bind unix:tgdbbackend.sock -m 007 tgdbbackend.wsgi
```

This binds the server to a unix socket, which can then be connected to from a reverse proxy such as nginx.
