import atexit
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

    for d in filterfalse(methodcaller('startswith', '.'), dirs):
        shutil.rmtree(storage.path(d))

    for f in filterfalse(methodcaller('startswith', '.'), files):
        os.remove(storage.path(f))


class Command(RunServerCommand):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        atexit.register(clear_query_cache)
