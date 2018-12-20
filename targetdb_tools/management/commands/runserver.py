import atexit
import datetime
import os
import shutil
from itertools import filterfalse
from operator import methodcaller

from django.conf import settings
from django.contrib.staticfiles.management.commands.runserver import Command as RunServerCommand
from django.core.files.storage import FileSystemStorage

storage = FileSystemStorage(settings.QUERY_CACHE)


def clear_query_cache():
    dirs, files = storage.listdir('.')
    now = datetime.datetime.now()

    for d in filterfalse(methodcaller('startswith', '.'), dirs):
        if now - datetime.datetime.fromtimestamp(os.path.getmtime(storage.path(d))) > datetime.timedelta(minutes=360):
            shutil.rmtree(storage.path(d), ignore_errors=True)

    for f in filterfalse(methodcaller('startswith', '.'), files):
        if now - datetime.datetime.fromtimestamp(os.path.getmtime(storage.path(f))) > datetime.timedelta(minutes=360):
            os.remove(storage.path(f))


class Command(RunServerCommand):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        atexit.register(clear_query_cache)
